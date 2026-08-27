# Table S1 Triples Correction Data

This repository contains computational chemistry data and input files related to the triple corrections for META, ORTHO, and PARA configurations. The calculations were performed using CCSD(T) and DLPNO-CCSD(T) methods at the 6-311G* basis set level of theory.

## Contents

The data is organized into two main directories based on the computational method used:

### `CCSDpT/`
Contains results from standard **CCSD(T)** calculations.
- **`META/`**: Output and input files for the META configuration.
- **`ORTHO/`**: Output and input files for the ORTHO configuration.
- **`PARA/`**: Output and input files for the PARA configuration.

### `DLPNO_CCSDpT/`
Contains results from **DLPNO-CCSD(T)** calculations (Domain Based Local Pair Natural Orbital).
- **`META/`**: Output and input files for the META configuration.
- **`ORTHO/`**: Output and input files for the ORTHO configuration.
- **`PARA/`**: Output and input files for the PARA configuration.

## File Descriptions

Within each configuration directory, you will find the following types of files:

- `*.out`: The primary output file containing the results of the calculation.
- `input.inp`: The ORCA input file used to perform the calculation (defines method, basis set, and molecular geometry).
- `input.bibtex`: BibTeX entries for the references cited in the calculations.
- `input.densities` / `input.densitiesinfo`: Files related to electron density information.
- `*.gbw`: Gaussian Basis Set (Orbital) files used by ORCA.
- `input.property.txt`: Extracted properties from the calculation output.
- `input.loc`: (In DLPNO directories) Localized orbital information.

## Computational Details

All calculations were performed using the **ORCA** quantum chemistry software package. The specific level of theory used for these entries is:
- **Method**: CCSD(T) or DLPNO-CCSD(T)
- **Basis Set**: 6-311G*
- **Additional Features**: RIJCOSX, AutoAux, TightSCF

## Usage for Reviewers and Collaborators

This dataset is intended to support the findings presented in the associated publication (Table S1). Researchers can use these files to verify the reported values or to replicate the calculations using the provided input files.
