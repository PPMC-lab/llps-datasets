This file has been created on 05-02-2024 by Carlos Pintado-Grima and last modified on 09-04-2024
 

GENERAL INFORMATION
------------------

1. Dataset title:
Confident datasets of client, driver and negative proteins in liquid-liquid phase separation


2. Authorship:

	Name: Carlos Pintado-Grima
	Institution: Institut de Biotecnologia i de Biomedicina, Universitat Autònoma de Barcelona
	Email: Carlos.Pintado@uab.cat
	ORCID: 0000-0002-8544-959X

	Name: Oriol Bárcenas
	Institution: Institut de Biotecnologia i de Biomedicina, Universitat Autònoma de Barcelona
	Email: Oriol.Barcenas@uab.cat
	ORCID: 0000-0002-8439-4005

	Name: Valentín Iglesias
	Institution: Institut de Biotecnologia i de Biomedicina, Universitat Autònoma de Barcelona
	Email: Valentin.Iglesias@uab.cat
	ORCID: 0000-0002-6133-0869

	Name: Eva Arribas-Ruiz
	Institution: Institut de Biotecnologia i de Biomedicina, Universitat Autònoma de Barcelona
	Email: Eva.Arribas@uab.cat
	ORCID: XYZ

	Name: Michał Burdukiewicz
	Institution: Institut de Biotecnologia i de Biomedicina, Universitat Autònoma de Barcelona and Medical University of Białystok (Poland)
	Email: michaljan.burdukiewicz@uab.cat
	ORCID: 0000-0001-8926-582X

	Name: Salvador Ventura
	Institution: Institut de Biotecnologia i de Biomedicina, Universitat Autònoma de Barcelona
	Email: Salvador.Ventura@uab.cat
	ORCID: 0000-0002-9652-6351




DESCRIPTION
----------

1. Dataset language:
English


2. Abstract:
Confident dataset of client, driver and negative proteins in liquid-liquid phase separation. The datasets collects unambiguous drivers, unambiguous clients and different dataset labels according to the independent source databases used to filter data. Careful data curation is applied to obtain integrated data confidence.


3. Keywords: 
Liquid-liquid phase separation; client; driver; negative; proteins; prediction; benchmark


4. Date of data collection (single date or date range):
29-01-2024

5. Publication Date:
Not applicable


6. Grant information:

	Grant Agency: Spanish Ministry of Science and Innovation 
	Grant Number: PID2019– 105017RB-I00

	Grant Agency: European Union’s Horizon 2020 research and innovation programme
	Grant Number: GA 952334 (PhasAGE)

	Grant Agency: Secretariat of Universities and Research of the Catalan Government and the European Social Fund
	Grant Number: 2023 FI_3 00018


7. Geographical location/s of data collection:
Barcelona, Spain



ACCESS INFORMATION
------------------------

1. Creative Commons License of the dataset:
CC-BY SA. 


2. Dataset DOI:
NA


3. Related publication:
Submitted

4. Link to related datasets:
Not applicable



VERSIONING AND PROVENANCE
---------------

1. Last modification date:
09-04-2024


2. Were data derived from another source?:
Yes, curated data derived from source databases of liquid-liquid phase separation.
- PhaSePro: 10.1093/nar/gkz848
- PhaSepDB:10.1093/nar/gkac783 
- LLPSDB: 10.1093/bioinformatics/btac026 
- CD-CODE:10.1038/s41592-023-01831-0
- DrLLPS: 10.1093/nar/gkz1027 
- DisProt: 10.1093/nar/gkad928 
- PDB: 10.1093/nar/28.1.235


3. Additional related data not included in this dataset:
No


METHODOLOGICAL INFORMATION
-----------------------

1. Description of the methods used to collect and generate the data:
Computational script to automatically collect and generate protein datasets.


2. Data processing methods:
Integrative rational filtering applied to independent souce databases.


3. Software or instruments needed to interpret the data:
Spreadsheet software like Excel or LibreOffice. Programming languages software for further analyses (e.g. Python, R).

4. Information about instruments, calibration and standards:
Not applicable


5. Environmental or experimental conditions:
Not applicable


6. Quality-assurance procedures performed on the data:
Rational data quality assessment and technical validation on data properties.



FILE OVERVIEW
----------------------
(main)
dataset_llps_29_01_24.tsv
sequential_elements.json
README.txt

(subdirectories)
./scripts/python_scripts
./figures/pngs
./properties/.csv

1. Explain the file naming conversion, if applicable:
dataset_llps_generation-date, for the main dataset file (tsv)


2. File list:

	File name: dataset_llps_29_01_24.tsv
	Description: main dataset file with all proteins classified into driver, client or negative.

	File name: sequential_elements.json
	Description: json dictionary with all IDRs, CARs and PrLDs for each protein.   

	File name: ./scripts/python_scripts
	Description: Python scripts to automatically generate protein datsets, generate tsv/json files and physicochemical analysis.

	File name: ./figures/pngs
	Description: Figures generated for the manuscript.   

	File name: ./properties/calculations_properties_datasets.csv
	Description: propoerty calculation for all sequences.  


3. Relationship between files:
The json file collects subsequences of the annotated proteins from the main tsv file. The python scripts are code to generate independent files.


4. File format:
.tsv (tab separated values), .json (JavaScript object notation) and .py (Python script)


5. If the dataset includes multiple files, specify the directory structure and relationships between the files:
Main files in main directory. Python scripts in ./scripts folder. Figures in ./figures. Property calculations in ./properties.



SPECIFIC INFORMATION FOR TABULAR DATA
-------------------------------------------

1. Name file:
dataset_llps_29_01_24.tsv


2. Number of rows and columns: 
4527 rows, 6 columns

3. Variables list:

	Variable name: UniProt Acc.
	Description: UniProt accession number of each protein
	Units of measure or value labels: text string

	Variable name: Datasets
	Description: list of datasets where each protein is located
	Units of measure or value labels: text string, each dataset separated by semicolon

	Variable name: GO
	Description: GO cellular location term from UniProt
	Units of measure or value labels: N (nuclear), C (cytoplasmatic), both (CN), none(_)

	Variable name: Derived.Order
	Description: fraction of PDB-derived ordered residues from MobiDB
	Units of measure or value labels: numerical float

	Variable name: Curated.Disorder
	Description: fraction of curated disordered residues from MobiDB
	Units of measure or value labels: numberical float

	Variable name: Full seq
	Description: residue sequence of the full length protein
	Units of measure or value labels: aminoacid string

4. Codes or symbols for missing data:
NA


SPECIFIC INFORMATION FOR JSON DATA
-------------------------------------------
1. Name file:
sequential_elements.json


2. Number of keys:
4527

3. Variables list:

	Variable name: IDRs.
	Description: list of intrinsically disordered regions' sequences obtained from MobiDB-AlphaFold
	Units of measure or value labels: list of strings

	Variable name: CARs
	Description: list of cryptic amyloidogenic regions' sequences obtained from the Waltz algorithm at threshold 85
	Units of measure or value labels: list of strings

	Variable name: PrLDs
	Description: list of prion-like domains' sequences obtained with the PLAAC algorithm
	Units of measure or value labels: list of strings


4. Codes or symbols for missing data:
Empty lists

        
5. Special formats or abbreviations used:
NA


MORE INFORMATION
--------------
The datasets here presented are still unbuplished but are aimed to be submitted to a scientific journal in the following weeks.