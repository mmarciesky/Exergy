# Exergy SMARTS Group Definitions

## Source

The group contribution scheme is defined in:

- `data/groups/Exergy_SMARTS_Groups.xlsx`

Columns:

- `SMARTS`      : RDKit SMARTS pattern for the group
- `delta_H`     : Enthalpy contribution (units: ...)
- `delta_S`     : Entropy contribution (units: ...)
- `filter func` : (optional) name of a static method on `SMARTSGroup` used as a per-match filter
- `extra func`  : (optional) name of a static method on `SMARTSGroup` used to adjust the count
- `binary`      : (optional) 0/1 or False/True, whether the group count is capped at 1

## How code uses this table

The loader in `exergy_bones/group_loader.py` reads the Excel file and:

1. Creates a `SMARTSCollection`.
2. For each row, builds a `SMARTSGroup` with:
   - name (currently "Group i")
   - SMARTS pattern from `SMARTS`
   - `delta_H`, `delta_S`
   - optional filter function: looked up from `SMARTSGroup.<filter func>`
   - optional extra counting function: `SMARTSGroup.<extra func>`
   - binary flag controlling whether the group count is capped at 1
3. Returns the `SMARTSCollection`, which is used to analyze molecules by SMILES.
