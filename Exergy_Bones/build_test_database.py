# scripts/build_test_database.py

import os
import time
import numpy as np
import pandas as pd
import pubchempy as pcp
import re
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from pathlib import Path

# --------- PubChem helper ---------

def cas2smiles_robust(cas, name=None, max_retries=5, base_sleep=1):
    """
    CAS + name → SMILES with retries on PUGREST.ServerBusy.
    Tries CAS first, then name.
    Returns SMILES or None.
    """
    cas  = str(cas).strip() if cas  else ""
    name = str(name).strip() if name else ""

    def _try_query(query, kind):
        sleep = base_sleep
        for attempt in range(1, max_retries + 1):
            try:
                comps = pcp.get_compounds(query, "name")
                if not comps:
                    return None
                props = pcp.get_properties("smiles", [comps[0].cid], "cid")
                if not props:
                    return None
                return props[0]["SMILES"]
            except Exception as e:
                msg = str(e)
                if "PUGREST.ServerBusy" in msg:
                    print(
                        f"[{kind}] ServerBusy for {query!r}, attempt {attempt}/{max_retries}. "
                        f"Sleeping {sleep:.1f} s..."
                    )
                    time.sleep(sleep)
                    sleep *= 2
                    continue
                else:
                    print(f"[{kind}] Error for {query!r}: {e}")
                    return None
        print(f"[{kind}] Giving up on {query!r} after {max_retries} retries.")
        return None

    if cas:
        smi = _try_query(cas, "CAS")
        if smi:
            return smi

    if name:
        smi = _try_query(name, "NAME")
        if smi:
            return smi

    return None

# --------- formula parsing helper ---------

element_pattern = re.compile(r"([A-Z][a-z]?)(\d*)")

def formula_to_counts(formula):
    if not isinstance(formula, str):
        return None
    formula = formula.strip()
    if not formula:
        return None

    counts = {}
    for elem, num in element_pattern.findall(formula):
        n = int(num) if num else 1
        counts[elem] = counts.get(elem, 0) + n
    return counts

# --------- main build function ---------

def build_test_database(raw_path: Path, partial_path: Path, clean_out_path: Path):
    """
    Build the cleaned exergy test database:
    - add SMILES from CAS/name via PubChem
    - check formula consistency via RDKit
    - keep only matching rows
    - save final Excel file
    """
    # ==== Load data (resume if partial file exists) ====
    if partial_path.exists():
        print(f"Resuming from existing file: {partial_path}")
        df_ref = pd.read_excel(partial_path)  # header already flat
    else:
        print(f"Loading original file: {raw_path}")
        df_ref = pd.read_excel(raw_path, header=[0, 1])

        flat_cols = []
        for top, bottom in df_ref.columns:
            top = "" if (isinstance(top, float) and np.isnan(top)) else str(top).strip()
            bottom = "" if (isinstance(bottom, float) and np.isnan(bottom)) else str(bottom).strip()
            if top and not bottom:
                flat_cols.append(top)
            elif bottom and not top:
                flat_cols.append(bottom)
            else:
                flat_cols.append(bottom)
        df_ref.columns = flat_cols

        df_ref.to_excel(partial_path, index=False)
        print(f"Initialized partial file at {partial_path}")

    cas_list  = df_ref["CASN"].tolist()
    name_list = df_ref["Chemical Name"].tolist()
    n_rows    = len(df_ref)

    if "SMILES" in df_ref.columns:
        smiles_from_cas = df_ref["SMILES"].tolist()
    else:
        smiles_from_cas = [None] * n_rows
        df_ref["SMILES"] = smiles_from_cas

    GLOBAL_SLEEP = 0.2

    # ==== PubChem loop ====
    try:
        for i, (cas, name) in enumerate(zip(cas_list, name_list)):
            existing = smiles_from_cas[i]

            if existing is not None and not (isinstance(existing, float) and np.isnan(existing)):
                continue

            smi = cas2smiles_robust(cas, name)
            smiles_from_cas[i] = smi
            df_ref.at[i, "SMILES"] = smi

            if (i + 1) % 50 == 0 or i == n_rows - 1:
                print(f"Processed {i + 1}/{n_rows}")

            if (i + 1) % 250 == 0 or i == n_rows - 1:
                df_ref.to_excel(partial_path, index=False)
                print(f"Saved checkpoint at row {i + 1} → {partial_path}")

            time.sleep(GLOBAL_SLEEP)

    except KeyboardInterrupt:
        print("Interrupted by user, saving partial results...")
        df_ref.to_excel(partial_path, index=False)
        print(f"Partial results saved to {partial_path}")

    df_ref.to_excel(partial_path, index=False)
    print(f"Done. Final results saved to {partial_path}")

    # Stats on SMILES found
    found = sum(
        1
        for x in df_ref["SMILES"]
        if x is not None and not (isinstance(x, float) and np.isnan(x)) and str(x).strip() != ""
    )
    print("SMILES found so far:", found)
    print("Missing:", n_rows - found)

    # ==== formula consistency check ====
    rdkit_formula_list = []
    formula_match_list = []
    formula_note_list  = []

    for smi, formula in zip(df_ref.get("SMILES", []), df_ref.get("Formula", [])):
        rdkit_f = None
        match   = False
        note    = ""

        if not isinstance(smi, str) or not smi.strip():
            note = "no_smiles"
        elif not isinstance(formula, str) or not formula.strip():
            note = "no_formula"
        else:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                note = "bad_smiles_parse"
            else:
                rdkit_f = rdMolDescriptors.CalcMolFormula(mol)

                ref_counts   = formula_to_counts(formula)
                rdkit_counts = formula_to_counts(rdkit_f)

                if ref_counts is None:
                    note = "bad_ref_formula"
                elif ref_counts == rdkit_counts:
                    match = True
                else:
                    note = f"counts_mismatch: ref={ref_counts}, rdkit={rdkit_counts}"

        rdkit_formula_list.append(rdkit_f)
        formula_match_list.append(match)
        formula_note_list.append(note)

    df_ref["RDKit_Formula"] = rdkit_formula_list
    df_ref["Formula_Match"] = formula_match_list
    df_ref["Formula_Note"]  = formula_note_list

    partial_path.parent.mkdir(parents=True, exist_ok=True)
    df_ref.to_excel(partial_path, index=False)
    print("Formula check complete and saved.")

    total      = len(df_ref)
    matches    = df_ref["Formula_Match"].sum()
    mismatches = total - matches

    print(f"Total rows: {total}")
    print(f"Formula matches: {matches}")
    print(f"Formula mismatches: {mismatches}")
    print(f"Percent correct: {matches/total*100:.2f}%")

    df_matches = df_ref[df_ref["Formula_Match"] == True].reset_index(drop=True))

    clean_out_path.parent.mkdir(parents=True, exist_ok=True)
    df_matches.to_excel(clean_out_path, index=False)
    print(f"Clean test set saved to {clean_out_path}")

    return {
        "total": int(total),
        "matches": int(matches),
        "mismatches": int(mismatches),
        "percent_correct": float(matches/total*100.0),
    }

if __name__ == "__main__":
    PROJECT_ROOT   = Path(__file__).resolve().parents[1]
    raw_path       = PROJECT_ROOT / "Data" / "raw" / "True_set.xlsx"
    partial_path   = PROJECT_ROOT / "Data" / "processed" / "True_set_with_smiles_partial.xlsx"
    clean_out_path = PROJECT_ROOT / "Data" / "processed" / "Exergy_Clean_Test.xlsx"

    stats = build_test_database(raw_path, partial_path, clean_out_path)
    print("Summary stats:", stats)
