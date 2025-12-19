# Exergy-PFAS: Python Implementation of the Gharagheizi Chemical Exergy Model

This repository contains a full Python implementation of the Standard Molar Chemical Exergy model developed by:

Gharagheizi, F.; Ilani-Kashkouli, P.; Hedden, R. C.  
"Standard Molar Chemical Exergy: A New Accurate Model."  
Energy 2017, 121, 384–394.  
https://doi.org/10.1016/j.energy.2017.01.010

The model computes exergy at 298.15 K from molecular structure alone (SMILES).  
It includes:

- Eq. 4 element-based term  
- Group-based enthalpy contributions (ΔHf)  
- Group-based entropy contributions (ΔSf)  
- Final exergy prediction

Users supply a SMILES string or an Excel file containing a SMILES column, and the model outputs:

- Eq.4 term  
- Group contributions  
- Total ΔH  
- Total ΔS  
- Final exergy (kJ/mol)

---

# Repository Structure
```
Exergy_Project/
│
├── Exergy_Bones/
│   ├── build_test_database.py        # Retrieve/test SMILES from CASN-style databases
│   ├── chem_utils.py                 # Eq.4 utilities: element counting & epsilon sum
│   ├── exergy_calc.py                # Compute Eq4, Sum(H), Sum(S), final exergy
│   ├── smarts_groups.py              # SMARTSGroup + SMARTSCollection classes
│   ├── group_loader.py               # Load SMARTS groups & values from Excel
│   ├── group_tester.py               # Tools for validating/tuning group definitions
│   ├── run_exergy.py                 # High-level function to compute exergy from Excel
│   ├── test_molecule_loader.py       # Load and parse test datasets
│   └── __init__.py
│
├── data/
│   ├── groups/
│   │   └── Exergy_Smarts_Groups.xlsx      # Editable Group definitions & values
│   │
│   ├── processed/
│   │   ├── Exergy_clean_test.xlsx         # Curated test set with resolved SMILES
│   │   └── Test_output.xlsx               # Example output from run_exergy()
│   │
│   └── raw/
│       ├── True_set.xlsx                  # Original test database (as in publication)
│       └── Test_sample_excel.xlsx         # Example input for run_exergy()
│   └── various fingerprint files
│
├── notebooks/
│   └── group_tuning.ipynb                 # For debugging and tuning SMARTS groups
│   └──group_with_fingerprints.ipynb # Testing the original database
│
└── README.md

```

---

# Installation
Requirements:
- Python 3.9+
- RDKit
- pandas
- NumPy

# Citation

If you use this code in academic work, please cite:

Gharagheizi, F.; Ilani-Kashkouli, P.; Hedden, R. C.
Standard Molar Chemical Exergy: A New Accurate Model.
Energy, 2017, 121, 384–394.
https://doi.org/10.1016/j.energy.2017.01.010

This project was developed for PFAS-related computational chemistry work in the Ng lab, University of Pittsburgh Pittsburgh, PA, USA by Dr. Joyce Li and Mel Marciesky
