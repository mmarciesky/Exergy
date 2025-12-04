# exergy_bones/loaders.py

from pathlib import Path
import pandas as pd

def load_exergy_test_db(project_root: Path = None) -> pd.DataFrame:
    if project_root is None:
        project_root = Path(__file__).resolve().parents[1]
    path = project_root / "data" / "processed" / "Exergy_Clean_Test.xlsx"
    return pd.read_excel(path)