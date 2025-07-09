[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15118996.svg)](https://doi.org/10.5281/zenodo.15118996)

This repository contains the datasets and scripts needed to reproduce the analysis of the manuscript: [Pintado-Grima, C., Bárcenas, O., Arribas-Ruiz, E. et al. Comprehensive protein datasets and benchmarking for liquid–liquid phase separation studies. Genome Biol 26, 198 (2025)](https://doi.org/10.1186/s13059-025-03668-6). Datasets are also publicly available at: [https://llpsdatasets.ppmclab.com/](https://llpsdatasets.ppmclab.com/).

# LLPS Datasets

This file has been created on 05-02-2024 by Carlos Pintado-Grima & and last modified on 09-07-2025

**Access link to the online dataset**: [LLPS Datasets](https://llpsdatasets.ppmclab.com/)

## GENERAL INFORMATION

### 1. Dataset title

Confident protein datasets for liquid-liquid phase separation studies

### 2. Authorship

- Carlos Pintado-Grima
  - Institution: Institut de Biotecnologia i de Biomedicina, Universitat Autònoma de Barcelona
  - Email: <Carlos.Pintado@uab.cat>
  - ORCID: [0000-0002-8544-959X](https://orcid.org/0000-0002-8544-959X)

- Oriol Bárcenas
  - Institution: Institut de Biotecnologia i de Biomedicina, Universitat Autònoma de Barcelona
  - Email: <Oriol.Barcenas@uab.cat>
  - ORCID: [0000-0002-8439-4005](https://orcid.org/0000-0002-8439-4005)

- Eva Arribas-Ruiz
  - Institution: Institut de Biotecnologia i de Biomedicina, Universitat Autònoma de Barcelona
  - Email: <eva.arribas00@gmail.com>
  - ORCID: [0009-0006-3110-0057](https://orcid.org/0009-0006-3110-0057)

- Valentín Iglesias
  - Institution: Institut de Biotecnologia i de Biomedicina, Universitat Autònoma de Barcelona
  - Email: <Valentin.Iglesias@uab.cat>
  - ORCID: [0000-0002-6133-0869](https://orcid.org/0000-0002-6133-0869)


- Michał Burdukiewicz
  - Institution: Institut de Biotecnologia i de Biomedicina, Universitat Autònoma de Barcelona and Medical University of Białystok (Poland)
  - Email: <michaljan.burdukiewicz@uab.cat>
  - ORCID: [0000-0001-8926-582X](https://orcid.org/0000-0001-8926-582X)

- Salvador Ventura
  - Institution: Institut de Biotecnologia i de Biomedicina, Universitat Autònoma de Barcelona
  - Email: <Salvador.Ventura@uab.cat>
  - ORCID: [0000-0002-9652-6351](https://orcid.org/0000-0002-9652-6351)

## DESCRIPTION

### 1. Dataset language

English

### 2. Abstract

Confident dataset of client, driver and negative proteins in liquid-liquid phase separation. The datasets collects unambiguous drivers, unambiguous clients, both clients and drivers and potential non-LLPS proteins. Careful data curation is applied to obtain integrated data confidence.

### 3. Keywords

Liquid-liquid phase separation; client; driver; negative; proteins; prediction; benchmark

### 4. Date of data collection (single date or date range)

29-01-2024

### 5. Publication Date

July 8th, 2025 (Genome Biology)

### 6. Grant information


- Grant Agency: Spanish Ministry of Science and Innovation
  - Grant Number: PID2019 – 105017RB-I00 to SV
  - Grant Number: FPU22/03656 to OB

- Grant Agency: AGAUR (Generalitat de Catalunya)
  - Grant Number: 2021-SGR-00635 to SV
  - Grant Number: 2023 FI_3 00018 to CPG 

- Grant Agency: Polish National Agency for Academic Exchange
  - Grant Number: BPN/ULM/2023/1/00189/U/00001 to VI

- Grant Agency: EU-Maria Zambrano to MB



### 7. Geographical location/s of data collection

Barcelona, Spain

## ACCESS INFORMATION

1. Open source License of the dataset:
MIT.

2. Dataset DOI:
[10.5281/zenodo.15118996](https://doi.org/10.5281/zenodo.15118996)

3. Related publication:
[10.1186/s13059-025-03668-6](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-025-03668-6)

4. Link to related datasets:
See source databases in the next section.

## VERSIONING AND PROVENANCE

### 1. Last modification date

05-06-2025

### 2. Were data derived from another source?

Yes, curated data derived from source databases of liquid-liquid phase separation.

- [PhaSePro](https://doi.org/10.1093/nar/gkz848)
- [PhaSepDB](https://doi.org/10.1093/nar/gkac783)
- [LLPSDB](https://doi.org/10.1093/bioinformatics/btac026)
- [CD-CODE](https://doi.org/10.1038/s41592-023-01831-0)
- [DrLLPS](https://doi.org/10.1093/nar/gkz1027)
- [DisProt](https://doi.org/10.1093/nar/gkad928)
- [PDB](https://doi.org/10.1093/nar/28.1.235)

### 3. Additional related data not included in this dataset

No

## METHODOLOGICAL INFORMATION

### 1. Description of the methods used to collect and generate the data

Computational script to automatically collect and generate protein datasets.

### 2. Data processing methods

Integrative rational filtering applied to independent souce databases.

### 3. Software or instruments needed to interpret the data

Spreadsheet software like Excel or LibreOffice. Programming languages software for further analyses (e.g. Python, R).

### 4. Information about instruments, calibration and standards

Not applicable

### 5. Environmental or experimental conditions

Not applicable

### 6. Quality-assurance procedures performed on the data

Rational data quality assessment and technical validation on data properties.

## FILE OVERVIEW

```{bash}
.
├── datasets.tsv
├── sequential_elements.json
├── README.MD
├── LICENSE
├── figures/
│   ├── 'Figure-0 Graphical abstract.png'
│   ├── 'Figure-1 Datasets.png'
│   ├── 'Figure-2 Integrative filters.png'
│   ├── 'Figure-3 Crosstalk.png'
│   ├── 'Figure-4 Order distribution.pdf'
│   └── 'Figure-5 Properties.png'
├── properties/
│   └── calculations_properties_datasets.csv
├── biogrid_interactors/
│   └── proteins.tsv
└── scripts/
    ├── generate_datasets.py
    ├── json_create.py
    ├── physicochemical_properties_analyses.py
    ├── properties-figure.R
    └── tsv_generation.py
```

### 1. Explain the file naming conversion, if applicable

Not applicable

### 2. File list

- File name: `datasets.tsv`
  - Description: main dataset file with all proteins classified into driver, client or negative.

- File name: `sequential_elements.json`
  - Description: json dictionary with all IDRs and PrLDs for each protein.

- Directory name: `./scripts`
  - Description: Python scripts to automatically generate protein datsets, generate tsv/json files and physicochemical analysis.

- Directory name: `./figures`
  - Description: Figures generated for the manuscript.

- Directory name: `./biogrid_interactors`
  - Description: Potential interactors of positive LLPS proteins obtained with BioGRID.

- File name: `./properties/calculations_properties_datasets.csv`
  - Description: property calculation for all sequences.  

### 3. Relationship between files

The json file collects subsequences of the annotated proteins from the main tsv file. The python scripts are code to generate independent files.

### 4. File format

`.tsv` (tab separated values), `.json` (JavaScript object notation) and `.py` (Python script)

### 5. If the dataset includes multiple files, specify the directory structure and relationships between the files

Main files in main directory. Python scripts in `./scripts` folder. Figures in `./figures`. Property calculations in `./properties`. BioGRID data in `./biogrid_interactors`

## SPECIFIC INFORMATION FOR TABULAR DATA

### 1. Name file

`datasets.tsv
`

### 2. Number of rows and columns

2876 rows, 6 columns

### 3. Variables list

- `UniProt Acc`.
  - Description: UniProt accession number of each protein
  - Units of measure or value labels: text string

- `Datasets`
 Description: list of datasets where each protein is located
 Units of measure or value labels: text string, each dataset separated by semicolon

- `GO`
  - Description: GO cellular location term from UniProt
  - Units of measure or value labels: N (nuclear), C (cytoplasmatic), both (CN), none(_)

- `Frac.Order`
  - Description: fraction of PDB-derived ordered residues from MobiDB
  - Units of measure or value labels: numerical float

- `Frac.Disorder`
  - Description: fraction of curated disordered residues from MobiDB
  - Units of measure or value labels: numberical float

- `Full.seq`
  - Description: residue sequence of the full length protein
  - Units of measure or value labels: aminoacid string

### 4. Codes or symbols for missing data

`NA`

## SPECIFIC INFORMATION FOR JSON DATA

### 1. Name file

`sequential_elements.json`

### 2. Number of keys

2876

### 3. Variables list

- `IDRs`
  - Description: list of intrinsically disordered regions' sequences obtained from MobiDB-AlphaFold
  - Units of measure or value labels: list of strings

- `PrLDs`
  - Description: list of prion-like domains' sequences obtained with the PLAAC algorithm
  - Units of measure or value labels: list of strings

- `Full_seq`
  - Description: full protein sequence
  - Units of measure or value labels: string

- `Dataset`
- Description: protein's dataset
- Units of measure or value labels: string

### 4. Codes or symbols for missing data

Empty lists

### 5. Special formats or abbreviations used

`NA`
