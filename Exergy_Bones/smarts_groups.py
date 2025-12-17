# exergy/smarts_groups.py

from rdkit import Chem
from typing import Callable, Optional, List
import pandas as pd

class SMARTSGroup:
    def __init__(
        self,
        name: str,
        smarts: str,
        value_h: float = 0.0,
        value_s: float = 0.0,
        filter_func: Optional[Callable] = None,
        extra_count: Optional[Callable] = None,
        binary: bool = False,
        padelcol: str = '',
        padelfp: str = '',        
    ):
        self.name = name
        self.smarts = smarts
        self.value_h = value_h
        self.value_s = value_s
        self.pattern = Chem.MolFromSmarts(smarts)
        self.filter_func = filter_func
        self.extra_count = extra_count
        self.binary = binary
        self.padelcol = padelcol
        self.padelfp = padelfp

    def count_in_molecule(self, mol: Chem.Mol,padelcount: pd.DataFrame,padelhit: pd.DataFrame) -> int:
        """
        Generic counter:
        - finds base SMARTS matches
        - applies an optional per-match filter
        - applies optional extra_count logic
        - if binary=True, caps count at 1
        """
        if self.pattern is None:
            raise ValueError(f"Invalid SMARTS pattern: {self.smarts}")
        # 0. check if calculated by padel
        if not pd.isnull(self.padelcol): #padelcount,padelhit
            
            if self.padelfp == 'AtomPairs2DCount':
                base_count = int(padelcount[self.padelcol].values[0])
            elif self.padelfp == 'AtomPairs2D':
                base_count = int(padelhit[self.padelcol].values[0])
            else:
                raise ValueError(f"Invalid fingerprint name")
        else:
            # 1. base SMARTS matches
            matches = mol.GetSubstructMatches(self.pattern)

            # 2. optional per-match filter
            if self.filter_func:

                matches = [m for m in matches if self.filter_func(mol, m)]

            base_count = len(matches)

            # 3. optional extra group-level logic (can add to the count)
            if self.extra_count:
                extra = self.extra_count(mol)
                return base_count + int(extra)

            # 4. optional binary behavior
            if self.binary and base_count > 1:
                base_count = 1

        return base_count

    # ========== specialized helper logics ==========

    @staticmethod
    def group_11_author_logic(mol: Chem.Mol) -> int:
        nitril_logic = "[C;D2;!a]#[N;D1;!a]"
        pattern = Chem.MolFromSmarts(nitril_logic)
        new_matches = mol.GetSubstructMatches(pattern)
        return len(new_matches)

    @staticmethod
    def filter_Y_hetero_terminal(mol: Chem.Mol, match: tuple) -> bool:
        """
        Group 2: hetero atom attached to one carbon only (terminal).
        match[0] = hetero, match[1] = carbon.
        """
        hetero_idx = match[0]
        connected_c_idx = match[1]

        hetero_atom = mol.GetAtomWithIdx(hetero_idx)
        for neighbor in hetero_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != connected_c_idx:
                return False  # Found extra C attached — reject
        return True

    @staticmethod
    def filter_Y_non_terminal(mol: Chem.Mol, match: tuple) -> bool:
        """
        Group 3: Y non-terminal logic.
        """
        y_idx = match[0]
        c1_idx = match[1]

        y_atom = mol.GetAtomWithIdx(y_idx)

        if y_atom.GetAtomicNum() == 6:
            return True  # Carbon always allowed

        heavy_neighbors = [
            nbr
            for nbr in y_atom.GetNeighbors()
            if nbr.GetAtomicNum() != 1 and nbr.GetIdx() != c1_idx
        ]
        return len(heavy_neighbors) > 0

    @staticmethod
    def filter_Y2N_C_Al(mol: Chem.Mol, match: tuple) -> bool:
        """
        Group 8: Y2N_C_Al logic.
        match = (n_idx, c_idx, o_idx, al_idx)
        """
        n_idx, c_idx, o_idx, al_idx = match
        n_atom = mol.GetAtomWithIdx(n_idx)

        # Get neighbors excluding the carbonyl carbon
        y_neighbors = [
            nbr
            for nbr in n_atom.GetNeighbors()
            if nbr.GetIdx() != c_idx and nbr.GetAtomicNum() != 1
        ]

        if len(y_neighbors) != 2:
            return False

        for y in y_neighbors:
            if y.GetAtomicNum() == 1:
                return False  # No hydrogens
            if y.GetSymbol() == "C":
                for nbr in y.GetNeighbors():
                    if nbr.GetIdx() == n_idx:
                        continue
                    if nbr.GetAtomicNum() == 8 and nbr.GetDegree() == 1:  # C=O group
                        return False

        return True

    @staticmethod
    def filter_group_24(mol: Chem.Mol, match: tuple) -> bool:
        """
        Group 24: overcount correction when SMARTS can't easily
        capture second carbon behavior.
        """
        X_set = {"O", "N", "S", "F", "Cl", "Br", "I", "Si"}

        for idx in match:
            atom = mol.GetAtomWithIdx(idx)

            if atom.GetSymbol() != "C":
                continue

            # Count X neighbors of this carbon
            x_neighbors = [
                nbr for nbr in atom.GetNeighbors() if nbr.GetSymbol() in X_set
            ]

            if len(x_neighbors) == 0:
                # Likely main carbon: check neighboring carbons for X groups
                carbon_to_x_counter = 0
                for n in atom.GetNeighbors():
                    if n.GetSymbol() == "C":
                        new_atom_x = [
                            nbr
                            for nbr in n.GetNeighbors()
                            if nbr.GetSymbol() in X_set
                        ]
                        if len(new_atom_x) > 1:
                            return False
                        if len(new_atom_x) == 1:
                            carbon_to_x_counter += 1
                if carbon_to_x_counter > 1:
                    return False

            if len(x_neighbors) > 2:
                return False

        return True
    @staticmethod
    def filter_no_end_wildcard_connection(mol: Chem.Mol, match: tuple) -> bool:
        """
        Enforce that wildcard atoms (* in the SMARTS) form a strictly linear chain:
        - Pattern is assumed: Anchor ~ * ~ * ~ ... ~ * ~ Anchor
          i.e. match = (anchor_left, A1, A2, ..., Ak, anchor_right)
        - Any wildcard Ai may only be bonded to adjacent wildcards in this sequence:
            A1 ↔ A2, A2 ↔ A1/A3, A3 ↔ A2/A4, ..., Ak ↔ A(k-1)
        - Bonds from Ai to atoms outside this wildcard set are allowed (side chains).
        - If any Ai is bonded to a non-adjacent Aj (e.g. A1–A3, A2–A5), reject.
        """

        # There must be at least one anchor on each side; wildcards are in between
        if len(match) < 3:
            # Not the pattern we expect, but don't block it here
            return True

        # Indices of wildcard atoms in the molecule
        star_indices = list(match)#[1:-1]      # A1..Ak
        star_set = set(star_indices)

        # Map each wildcard atom index -> its position in the chain (0..k-1)
        position = {idx: i for i, idx in enumerate(star_indices)}

        # For each wildcard atom, ensure any wildcard neighbor is only adjacent in sequence
        for idx in star_indices:
            atom = mol.GetAtomWithIdx(idx)
            i_pos = position[idx]

            for nbr in atom.GetNeighbors():
                j = nbr.GetIdx()
                if j in star_set:
                    j_pos = position[j]
                    # Only allow neighbors that are at positions i±1
                    if abs(i_pos - j_pos) != 1:
                        # Found a bond like A1–A3, A2–A4, etc. → reject
                        return False

        # If we didn't violate the rule, this match is OK
        return True
    

    @staticmethod
    def filtergroup_26(mol: Chem.Mol, match: tuple) -> bool:
        # Get atom indices
        cl_idx = match[0]  # Chlorine atom
        connected_idx = match[1]  # Atom connected to chlorine
        
        cl_atom = mol.GetAtomWithIdx(cl_idx)
        connected_atom = mol.GetAtomWithIdx(connected_idx)
        if connected_atom.GetSymbol() == 'C':
            # Get hybridization
            hybridization = connected_atom.GetHybridization()
            hybridization_str = str(hybridization)
            if hybridization_str != 'SP3':
                neibhorcount = 1
                for neighbor in connected_atom.GetNeighbors():
                    if neighbor.GetIdx() == cl_idx:
                        continue  # Skip the chlorine
                    elif neighbor.GetSymbol() == 'C':
                        neibhorcount += 0
                    elif neighbor.GetSymbol() == 'H':  
                        neibhorcount += -1 
                    elif (neighbor.GetSymbol() == 'O')|(neighbor.GetSymbol() == 'S'):
                        neibhorcount += 2
                    elif (neighbor.GetSymbol() in ['F','Cl','Br','I']):
                        neibhorcount += 1
                if (neibhorcount > 1) & (hybridization_str == 'SP2'):
                    return True
                elif (neibhorcount == 1) & (hybridization_str == 'SP'):
                    return True
                elif (neibhorcount == 4) & (hybridization_str == 'SP'):
                    return True
                else:
                    return False   
            else:
                return False
        else:
            return True    


class SMARTSCollection:
    """
    Store definitions of SMARTSGroup (e.g. from your table).
    """

    def __init__(self):
        self.groups: List[SMARTSGroup] = []

    def add_group(
        self,
        name: str,
        smarts: str,
        value_h: float = 0.0,
        value_s: float = 0.0,
        filter_func: Optional[Callable] = None,
        extra_count: Optional[Callable] = None,
        binary: bool = False,
        padelcol: str = '',
        padelfp: str = '',
    ):
        """
        Create a SMARTSGroup and add it to the collection.
        """
        group = SMARTSGroup(
            name,
            smarts,
            value_h,
            value_s,
            filter_func,
            extra_count,
            binary,
            padelcol,
            padelfp,
        )
        self.groups.append(group)

    def analyze_smiles(self, smiles: str, padelcount: pd.DataFrame, padelhit: pd.DataFrame):
        """
        Analyze a SMILES string against all stored SMARTS groups.
        Returns a list of dicts: name, smarts, count, value_h, value_s
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string")
        mol = Chem.AddHs(mol)

        results = []
        for group in self.groups:
            count = group.count_in_molecule(mol,padelcount,padelhit)
            results.append(
                {
                    "name": group.name,
                    "smarts": group.smarts,
                    "count": count,
                    "value_h": group.value_h,
                    "value_s": group.value_s,
                }
            )

        return results
