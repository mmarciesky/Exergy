Python Implementation for Automated Exergy Prediction from SMILES based of the work by Gharagheizi, F.; Ilani-Kashkouli, P.; Hedden, R. C. Standard Molar Chemical Exergy: A New Accurate Model. Energy 2017, 121, 384–394. https://doi.org/10.1016/j.energy.2017.01.010.

This repository contains a full Python implementation of the Exergy Group Contribution Model, allowing you to compute:

Eq. 4 (element-based term)

Group-based enthalpy & entropy contributions

Exergy at 298.15 K

from the SMILES

The model reads SMARTS patterns and group values from a customizable Excel file, applies RDKit substructure matching, performs optional filtering, computes group contributions, and returns the full exergy value. Excel file must have a "SMILES" column.

Exergy_Project/
│
├── Exergy_Bones/
│   ├── build_test_database.py # for searching through CASN database and retrieving SMILES for making Exergy_clean_test
│   ├── chem_utils.py          # Eq.4 utilities: element counting, epsilon sum
│   ├── exergy_calc.py         # Functions to compute H, S, Eq4, and Exergy
│   ├── smarts_groups.py       # SMARTSGroup + SMARTSCollection classes
│   ├── group_loader.py        # Load SMARTS groups from Excel
│   ├── group_tester.py        # for group tuning notebooks
│   ├── run_exergy.py          # High-level helper to process SMILES files
│   ├── test_molecule_loader.py # load in the test data set
│   └── __init__.py
│
├── data/
|   ├── groups
|       └──Exergy_Smarts_Groups.xslx # The Group parameters this can be edited if new Groups are needed
│   ├── processed 
|       └──Exergy_clean_test.xslx # The processed True_set that contains all the SMILES available from the CASN
|       └──Test_output.xslx # The final output of the Testing file for Run_Exergy
│   └── raw     
|       └──True_set.xslx  # contains the Test database with the groups found in the study
|       └──Test_sample_excel.xslx  # Testing file for Run_Exergy
│
├── notebooks/
│   └── group_tuning.ipynb   #for testing Group parameters against Exergy_clean_test
│
└── README.md

Installing requirements
Python 3.9+
RDKit
pandas
NumPy

