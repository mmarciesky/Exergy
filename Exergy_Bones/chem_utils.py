# exergy/chem_utils.py

from rdkit import Chem

# constants for equation 4
Epsilon = {
    "C": 410.26,
    "H": 118.05,
    "N": 0.36,
    "O": 1.985,
    "S": 609.6,
    "F": 233.15,
    "Cl": 61.8,
    "Br": 50.6,
    "I": 87.35,
    "Si": 854.6,
}

def element_counts(smiles: str) -> dict:
    """
    Return a dict of element counts for a SMILES string.
    Includes H after adding explicit hydrogens.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:  # stop the script if the SMILES cant be found
        raise ValueError(f"Problem with SMILES: {smiles}")
    mol = Chem.AddHs(mol)  # add hydrogens to count them

    Z_to_symbol = {
        6: "C",
        1: "H",
        7: "N",
        8: "O",
        16: "S",
        9: "F",
        17: "Cl",
        35: "Br",
        53: "I",
        14: "Si",
    }

    counts = {sym: 0 for sym in Z_to_symbol.values()}
    for atom in mol.GetAtoms():
        Z = atom.GetAtomicNum()
        if Z in Z_to_symbol:
            counts[Z_to_symbol[Z]] += 1
    return counts


def equation_4_part_2(counts: dict, epsilon: dict = Epsilon) -> float:
    """
    Given a dict of element counts, compute the sum over counts * epsilon.
    """
    return sum(counts[atom] * epsilon[atom] for atom in counts)
#wrap
def eq4_from_smiles(smiles: str) -> float:
    counts = element_counts(smiles)
    return equation_4_part_2(counts)

if __name__ == "__main__":
    counts = element_counts("C=C=C")
    print("Test eq4:", equation_4_part_2(counts))
