# exergy_bones/group_tester.py

from typing import List, Tuple
import pandas as pd

from .smarts_groups import SMARTSCollection


def summarize_groups(smarts_db: SMARTSCollection) -> None:
    """
    Print a quick summary of all groups loaded from Excel.
    """
    for group in smarts_db.groups:
        print(f"Name:    {group.name}")
        print(f"SMARTS:  {group.smarts}")
        print(f"ΔH:      {group.value_h}, ΔS: {group.value_s}")
        print(f"Filter:  {group.filter_func.__name__ if group.filter_func else None}")
        print(f"Extra:   {group.extra_count.__name__ if group.extra_count else None}")
        print(f"Binary:  {group.binary}")
        print("RDKit SMARTS valid:", group.pattern is not None)
        print("-" * 40)


def compare_groups_on_dataframe(
    df_test: pd.DataFrame,
    smarts_db: SMARTSCollection,
    group_name_prefix: str = "Group "
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Compare original group counts in df_test (e.g. 'Group 1', 'Group 2', ...)
    against counts computed by the SMARTSCollection.

    Returns:
        df_work     : copy of df_test with added '(calc)' and '(diff)' columns
        df_mismatch : long-form table listing each mismatched group per molecule
    """

    # 1. Get list of implemented group names (e.g. "Group 1", "Group 2", ...)
    implemented_groups: List[str] = [g.name for g in smarts_db.groups]
    print("Implemented groups:", implemented_groups)

    # 2. Work on a copy
    df_work = df_test.copy()

    # 3. Add "(calc)" columns initialized to 0
    for g in implemented_groups:
        df_work[g + " (calc)"] = 0

    # 4. For each molecule, compute group counts from SMARTS
    for idx, row in df_work.iterrows():
        smiles = row["SMILES"]
        results = smarts_db.analyze_smiles(smiles)
        # results is list of dicts: {"name", "smarts", "count", "value_h", "value_s"}
        count_map = {r["name"]: r["count"] for r in results}

        for g in implemented_groups:
            df_work.at[idx, g + " (calc)"] = count_map.get(g, 0)

    # 5. Compute "(diff)" columns = calc - original
    for g in implemented_groups:
        orig_col = g              # e.g. "Group 1"
        calc_col = g + " (calc)"  # e.g. "Group 1 (calc)"
        diff_col = g + " (diff)"  # e.g. "Group 1 (diff)"

        df_work[diff_col] = df_work[calc_col] - df_work[orig_col].fillna(0)

    # 6. Find rows where ANY group has a non-zero diff
    diff_cols = [g + " (diff)" for g in implemented_groups]
    mask_any_diff = (df_work[diff_cols] != 0).any(axis=1)
    df_diff_rows = df_work[mask_any_diff].copy()

    print(f"Rows with at least one group mismatch: {len(df_diff_rows)}")

    # 7. Build long-form table listing each mismatched group per molecule
    records = []

    for idx, row in df_diff_rows.iterrows():
        smiles = row.get("SMILES", None)
        name   = row.get("Chemical Name", row.get("Name", None))

        for g in implemented_groups:
            orig_val = row.get(g, 0)
            calc_val = row.get(g + " (calc)", 0)
            diff_val = row.get(g + " (diff)", 0)

            if diff_val != 0:
                records.append({
                    "index": idx,
                    "Name": name,
                    "SMILES": smiles,
                    "Group": g,
                    "Orig_Count": orig_val,
                    "Calc_Count": calc_val,
                    "Diff": diff_val,
                })

    df_mismatch = pd.DataFrame(records)
    if not df_mismatch.empty:
        df_mismatch["AbsDiff"] = df_mismatch["Diff"].abs()
        df_mismatch = df_mismatch.sort_values(
            ["AbsDiff", "Group", "Name"], ascending=[False, True, True]
        )

    return df_work, df_mismatch
