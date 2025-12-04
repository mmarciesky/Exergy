# Exergy_Bones/exergy_calc.py

import pandas as pd
from . import chem_utils  # assumes chem_utils.py is in the same package


def attach_exergy_columns(
    df_work: pd.DataFrame,
    smarts_db,
    T: float = 298.15,
    H_const: float = -23.9527,  # your coefficient constant
    S_const: float = 0.0205,    # your coefficient constant
) -> pd.DataFrame:
    """
    From df_work (which already has Group i (calc) columns and SMILES),
    build df_calc with:
      - ID columns
      - Eq4_term
      - Sum_H_groups, Sum_S_groups
      - Exergy (kJ/mol)
    """

    # 1. Identify columns
    id_cols = ["No.", "Chemical Name", "SMILES", "CASN", "Formula"]
    id_cols_present = [c for c in id_cols if c in df_work.columns]

    # only group count columns that end with ' (calc)'
    calc_cols = [c for c in df_work.columns if c.endswith(" (calc)")]

    # 2. Start df_calc
    df_calc = df_work[id_cols_present + calc_cols].copy()

    # 3. Eq. 4 term from SMILES
    df_calc["Eq4_term"] = df_calc["SMILES"].apply(chem_utils.eq4_from_smiles)

    # 4. Build group → (ΔH, ΔS) lookup
    group_params = {
        g.name: (float(g.value_h), float(g.value_s))
        for g in smarts_db.groups
    }

    # 5. Sum Σ n_i ΔH_i and Σ n_i ΔS_i (+ your constant coefficients)
    def _sum_group_H_S(row):
        H_sum = H_const
        S_sum = S_const

        for col in calc_cols:
            group_name = col[:-7]  # strip " (calc)"

            n_i = row[col]
            if pd.isna(n_i):
                n_i = 0.0

            H_i, S_i = group_params.get(group_name, (0.0, 0.0))
            H_sum += n_i * H_i
            S_sum += n_i * S_i

        return pd.Series({"Sum_H_groups": H_sum, "Sum_S_groups": S_sum})

    df_calc[["Sum_H_groups", "Sum_S_groups"]] = df_calc.apply(_sum_group_H_S, axis=1)

    # 6. Final exergy: H - T * S + Eq4
    df_calc["Exergy (kJ/mol)"] = (
        df_calc["Sum_H_groups"] - T * df_calc["Sum_S_groups"] + df_calc["Eq4_term"]
    )

    return df_calc
