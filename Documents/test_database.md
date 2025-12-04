# Exergy Test Database

## Source

- Original data: `data/raw/True_set.xlsx` taken form 10.1016/j.energy.2018.05.186
- SMILES generated from CAS or chemical name using PubChem (via PubChemPy).
- Formula consistency checked with RDKit (`CalcMolFormula`).

The test set was built by running:

```bash
python scripts/build_test_database.py
Total rows: 3148
Formula matches: 2609
Formula mismatches: 539
Percent correct: 82.88%
Final Test database was created with these matching molecules. 
```python
import json

stats_path = PROJECT_ROOT / "data" / "processed" / "Exergy_Clean_Test_stats.json"
with open(stats_path, "w") as f:
    json.dump(stats, f, indent=2)