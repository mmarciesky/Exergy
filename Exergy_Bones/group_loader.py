# exergy_bones/group_loader.py

from pathlib import Path
from typing import Optional

import pandas as pd

from .smarts_groups import SMARTSGroup, SMARTSCollection


def load_smarts_groups_from_excel(path: Optional[Path] = None) -> SMARTSCollection:
    """
    Load SMARTS groups from an Excel file and return a SMARTSCollection.

    Expected columns in the Excel file:
      - 'SMARTS'      : SMARTS pattern string
      - 'delta_H'     : group enthalpy contribution
      - 'delta_S'     : group entropy contribution
      - 'filter func' : (optional) name of a staticmethod on SMARTSGroup
      - 'extra func'  : (optional) name of a staticmethod on SMARTSGroup
      - 'binary'      : (optional) truthy value (1/True) if group is binary-counted

    Rows with missing SMARTS are skipped.
    """
    # Default path if none provided
    if path is None:
        project_root = Path(__file__).resolve().parents[1]
        path = project_root / "Data" / "groups" / "Exergy_SMARTS_Groups.xlsx"

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Group Excel file not found at: {path}")

    df = pd.read_excel(path, engine="openpyxl")

    smarts_db = SMARTSCollection()

    for idx, row in df.iterrows():
        smarts = row.get("SMARTS", None)

        # Skip completely empty or NaN SMARTS
        if not isinstance(smarts, str) or not smarts.strip():
            continue

        name = f"Group {idx + 1}"
        value_h = row.get("delta_H", 0.0)
        value_s = row.get("delta_S", 0.0)

        # filter function name from Excel
        filter_name = row.get("filter func", None)
        if isinstance(filter_name, str) and filter_name.strip() != "":
            # Look up a staticmethod on SMARTSGroup with that name
            filter_func = getattr(SMARTSGroup, filter_name)
        else:
            filter_func = None

        # extra counting function name from Excel
        extra_name = row.get("extra func", None)
        if isinstance(extra_name, str) and extra_name.strip() != "":
            extra_count = getattr(SMARTSGroup, extra_name)
        else:
            extra_count = None

        # binary flag
        binary_val = row.get("binary", None)
        if pd.notna(binary_val) and bool(binary_val):
            binary_count = True
        else:
            binary_count = False  # better than None: explicit False

        smarts_db.add_group(
            name=name,
            smarts=smarts,
            value_h=value_h,
            value_s=value_s,
            filter_func=filter_func,
            extra_count=extra_count,
            binary=binary_count,
        )

    return smarts_db
