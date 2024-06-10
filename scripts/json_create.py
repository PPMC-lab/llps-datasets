import sys
import json

def write_idrs(d, uniprot):
	with open("consensus_idrs.json", "rt") as json_file: #From MobiDB
		d_idrs = json.load(json_file)
		for key, value in d_idrs.items():
			if key.strip() == uniprot:
				d[uniprot]["IDRs"] = value


def write_cars(d, uniprot):
	#generar archivo de nuevo CARs
	d[uniprot]["CARs"]=[]
	with open("./CARs_85/CARs_85_datasets.csv", "rt") as csv: #CARs outfile
		for line in csv:
			uniprot_csv=line.split(",")[0].split("_")[0].strip()
			waltz=float(line.split(",")[4].strip())
			car=line.split(",")[3].strip()
			if uniprot == uniprot_csv:
				if waltz >= 85 and len(car) >= 7:
					d[uniprot]["CARs"].append(car)

def write_prlds(d, uniprot):
	d[uniprot]["PrLDs"]=[]
	with open("./PrLDs/prlds_plaac.txt", "rt") as txt: #PLAAC outfile
		txt.readline()
		for line in txt:
			uniprot_txt=line.split("\t")[0].strip()
			corescore=float(line.split("\t")[11].strip())
			prdaa=line.split("\t")[25].strip()
			if uniprot == uniprot_txt:
				if corescore > 0 and len(prdaa)>=20:
					d[uniprot]["PrLDs"].append(prdaa)


def write_dataset(d, uniprot, line):
		d[uniprot]["Dataset"]=[]
		dataset_cell=line.split("\t")[1]
		if "CO" in dataset_cell:
			dataset="CE"
		elif "DO" in dataset_cell:
			dataset="DE"
		elif "C_D" in dataset_cell:
			dataset="C_D"
		elif "DisProt" in dataset_cell:
			dataset="ND"
		elif "PDB" in dataset_cell:
			dataset="NP"
		else:
			pass

		d[uniprot]["Dataset"]=dataset

def write_fullseq(d, uniprot, line):
	d[uniprot]["Full_seq"]=[]
	seq=line.split("\t")[3].strip()
	d[uniprot]["Full_seq"]=seq



#generate json
d={}
with open("datasets_info_3.tsv", "rt") as tsv: #Dataset tsv
	tsv.readline()
	for line in tsv:
		uniprot=line.split("\t")[0].strip()
		d[uniprot]={}
		write_fullseq(d, uniprot, line)
		write_idrs(d, uniprot)
		#write_cars(d, uniprot)
		write_prlds(d, uniprot)
		write_dataset(d, uniprot, line) #CO, DO, C_D, DisProt, PDB

#print(d)
with open("element_annotations_final.json", "w") as outfile: #final json
	json.dump(d, outfile)		