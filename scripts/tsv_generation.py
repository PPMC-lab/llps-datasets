# import required module
import os
import sys
import subprocess
import re
import json
import time

# assign directory
directory = './Datasets_26_01_24' #datasets directory


def check_datasets(uniprot, directory):
	datasets=""
	for filename in os.listdir(directory):
		f = os.path.join(directory, filename)
		with open(f, "rt") as txt:
			for line in txt:
				uniprot2=line.strip("\n")
				if uniprot == uniprot2:
					datasets+=filename[:-4]+";"
				else:
					continue

	return datasets[:-1]

def json_parse(d):
	# Convert and write JSON object to file
	with open("annotations_3.json", "w") as outfile: 
   		json.dump(d, outfile)

def uniprot_parse(uniprot, d):
	GO=""
	
	#Obtain UniProt's txt file for protein X
	os.system("curl -H 'Accept: text/plain; format=flatfile' 'https://rest.uniprot.org/uniprotkb/"+uniprot+"' > temporary_file.txt")

	
	try:
		output_GO= str(subprocess.check_output("grep 'DR   GO;' temporary_file.txt", shell=True)) #Obtain GOcc
	except:
		return d, "fail", "fail"
	output_fasta= str(subprocess.check_output("sed -n '/SQ   SEQUENCE/,$p' temporary_file.txt | tail -n +2 | head -n -1", shell=True)) #Obtain seq

	#GO parse
	GO_list=output_GO.split(';')
	for GO_term in GO_list:
		if "C:" in GO_term: #GO_CC term
			match_cyto = re.search(r'cyto',GO_term)
			match_nucl = re.search(r'nucl',GO_term)
			if match_cyto:
				if "C" in GO: #cytoplasmatic protein
					continue
				GO+="C"
			if match_nucl: #nuclear protein
				if "N" in GO:
					continue
				GO+="N"
	
	if GO == "":
		GO = "_" #Not nuclear nor cytoplasmatic
	
	d[uniprot]["GO"] = GO





	#protein_name=output_description.strip().split("=")[1][:-4]
	#d["Protein_name"] = protein_name

	#Seq parse
	sequence=output_fasta.strip("b").strip("'").split("\\n") #Format sequence string
	final_seq=""
	for subseq in sequence: #Remove blank spaces
		final_seq+=subseq.replace(" ", "")
	d[uniprot]["Sequence"] = final_seq


	return d, GO, final_seq


def write_tsv(uniprot, datasets, GO, sequence):
	with open("datasets_info_3.tsv", "at") as tab:
		if "GO" != "fail":
			tab.write(uniprot+"\t"+datasets+"\t"+GO+"\t"+sequence+"\n")



# iterate
uniprots_list=[]
d={}

for filename in os.listdir(directory):
	f = os.path.join(directory, filename)
	with open(f, "rt") as txt:
		for line in txt:
			uniprot=line.strip("\n")
			if uniprot not in uniprots_list: #Already computed
				d[uniprot]={}
				uniprots_list.append(uniprot)
				datasets=check_datasets(uniprot, directory)
				print(uniprot, datasets)

			else:
				continue

			time.sleep(0.5)
			d, GO, final_seq=uniprot_parse(uniprot, d)

			write_tsv(uniprot, datasets, GO, final_seq)
			#print(d)


json_parse(d)



