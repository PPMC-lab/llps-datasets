#physiochemical_properties_analysis_heatmatplot
import os
import json
import sys
import subprocess

def kappa(seq):
	from localcider.sequenceParameters import SequenceParameters
	# create a SequenceParameters object with your amino acid sequence
	try:
		score = SequenceParameters(seq)
		return score.get_kappa()
	except:
		return "NA"

def kappa_ss(seq):
	from michalcider.sequenceParameters import SequenceParameters
	# adapted version of localcider ~ michalcider for kappa ss
	
	try:
		score = SequenceParameters(seq)
		return score.get_kappa()
	except:
		return "NA"

def yr(seq):
	yr_count=seq.count("Y")+seq.count("R")
	score=(yr_count/len(seq))*100

	return score

def get_pKa(res):
	"""
	Return pKa values for titratable residues. """

	d_libro_nuevo={'H': 6,'E': 4.25,'D': 3.65,'K': 10.53,'R': 12.48, }

	return d_libro_nuevo[res]

def charge_equation(res, pH=7):
	'''Apply H-H equation. pKa+pKb=14
	Start with charge at pH 7'''

	acids=['E', 'D',] # 'Ct']
	bases=['H', 'K', 'R', ] #'Nt']

	if res in acids:
		pKa=get_pKa(res)
		carga=-(10**(pH-pKa))/(1+ 10**(pH-pKa))

	elif res in bases:
		pKa=get_pKa(res)
		carga=(10**(pKa-pH))/(1+ 10**(pKa-pH))
	else:
		carga=0

	return carga


def ncpr(seq, pH):
	prot_charge=0
	for res in seq:
		prot_charge=prot_charge+charge_equation(res, pH,)

	return (prot_charge)/len(seq)

def aggrescan(seq):
	with open("tmp.fasta", "wt") as fasta:
		fasta.write(">tmp\n%s"% (seq))
	try:
		score=float(subprocess.check_output('./aggrescan tmp.fasta 1', shell=True))
	except:
		score="NA"

	return score

def anchor(seq):
	os.system("python3 ./iupred3/iupred3.py tmp.fasta long -a > anchor_out.txt")
	anchor_mean=[]
	with open("anchor_out.txt", "rt") as txt:
		for line in txt:
			if line[0] == "#":
				continue
			else:
				anchor_res=float(line.split("\t")[3].strip())
				anchor_mean.append(anchor_res)
	try:
		return sum(anchor_mean)/len(anchor_mean)
	except:
		return "NA"

def write_csv(uniprot, dataset, datatype, kappa, kappa_ss, yr, ncpr, aggrescan, anchor):
	with open("./properties/properties_v3.csv", "at") as csv_final:
		csv_final.write("%s,%s,%s,%s,%s\n" % (uniprot, dataset, datatype, "kappa", kappa))
		csv_final.write("%s,%s,%s,%s,%s\n" % (uniprot, dataset, datatype, "kappa_ss", kappa_ss))
		csv_final.write("%s,%s,%s,%s,%s\n" % (uniprot, dataset, datatype, "yr", yr))
		csv_final.write("%s,%s,%s,%s,%s\n" % (uniprot, dataset, datatype, "ncpr", ncpr))
		csv_final.write("%s,%s,%s,%s,%s\n" % (uniprot, dataset, datatype, "aggrescan", aggrescan))
		csv_final.write("%s,%s,%s,%s,%s\n" % (uniprot, dataset, datatype, "anchor", anchor))


with open("element_annotations_final.json", "rt") as json_file:
	d = json.load(json_file)

for uniprot, element in d.items():
	for element, value in element.items():

		if element == "IDRs":
			if not value: #empty list, value is 'NA'
				write_csv(uniprot, d[uniprot]["Dataset"], "IDRs", "NA", "NA", "NA", "NA", "NA", "NA")

			else:
				#create empty lists for all properties:
				score_idrs_list_kappa=[]
				score_idrs_list_kappa_ss=[]
				score_idrs_list_yr=[]
				score_idrs_list_ncpr=[]
				score_idrs_list_aggrescan=[]
				score_idrs_list_anchor=[]

				for seq in value:
					seq=seq.strip()
					#kappa
					scorep_idrs_kappa=kappa(seq)
					score_idrs_list_kappa.append(scorep_idrs_kappa)

					#kappa_ss
					scorep_idrs_kappa_ss=kappa_ss(seq)
					score_idrs_list_kappa_ss.append(scorep_idrs_kappa_ss)

					#yr
					scorep_idrs_yr=yr(seq)
					score_idrs_list_yr.append(scorep_idrs_yr)

					#ncpr
					scorep_idrs_ncpr=ncpr(seq, 7)
					score_idrs_list_ncpr.append(scorep_idrs_ncpr)


					#aggrescan
					scorep_idrs_aggrescan=aggrescan(seq)
					score_idrs_list_aggrescan.append(scorep_idrs_aggrescan)


					#anchor
					scorep_idrs_anchor=anchor(seq)
					score_idrs_list_anchor.append(scorep_idrs_anchor)
				
				#final averaged values for each protein/property
				try:
					score_idrs_kappa=sum(score_idrs_list_kappa) / len(score_idrs_list_kappa)
				except:
					score_idrs_kappa="NA"
				try:
					score_idrs_kappa_ss=sum(score_idrs_list_kappa_ss) / len(score_idrs_list_kappa_ss)
				except:
					score_idrs_kappa_ss="NA"

				score_idrs_yr=sum(score_idrs_list_yr) / len(score_idrs_list_yr)

				score_idrs_ncpr=sum(score_idrs_list_ncpr) / len(score_idrs_list_ncpr)

				try:
					score_idrs_aggrescan=sum(score_idrs_list_aggrescan) / len(score_idrs_list_aggrescan)
				except:
					score_idrs_aggrescan="NA"

				try:
					if "NA" in score_idrs_list_anchor:
						score_idrs_list_anchor=score_idrs_list_anchor.remove("NA")

					score_idrs_anchor=sum(score_idrs_list_anchor) / len(score_idrs_list_anchor)

				except:
					score_idrs_anchor="NA"


				#write results
				write_csv(uniprot, d[uniprot]["Dataset"], "IDRs", score_idrs_kappa, score_idrs_kappa_ss, score_idrs_yr, score_idrs_ncpr, score_idrs_aggrescan, score_idrs_anchor)

		elif element == "CARs":
			if not value:
				write_csv(uniprot, d[uniprot]["Dataset"], "CARs", "NA", "NA", "NA", "NA", "NA", "NA")

			else:
				#create empty lists for all properties:
				score_cars_list_kappa=[]
				score_cars_list_kappa_ss=[]
				score_cars_list_yr=[]
				score_cars_list_ncpr=[]
				score_cars_list_aggrescan=[]
				score_cars_list_anchor=[]

				for seq in value:
					#kappa
					scorep_cars_kappa=kappa(seq)
					score_cars_list_kappa.append(scorep_cars_kappa)

					#kappa_ss
					scorep_cars_kappa_ss=kappa_ss(seq)
					score_cars_list_kappa_ss.append(scorep_cars_kappa_ss)

					#yr
					scorep_cars_yr=yr(seq)
					score_cars_list_yr.append(scorep_cars_yr)

					#ncpr
					scorep_cars_ncpr=ncpr(seq, 7)
					score_cars_list_ncpr.append(scorep_cars_ncpr)


					#aggrescan
					scorep_cars_aggrescan=aggrescan(seq)
					score_cars_list_aggrescan.append(scorep_cars_aggrescan)


					#anchor
					scorep_cars_anchor=anchor(seq)
					score_cars_list_anchor.append(scorep_cars_anchor)
				
				#final averaged values for each protein/property
				try:
					score_cars_kappa=sum(score_cars_list_kappa) / len(score_cars_list_kappa)
				except:
					score_cars_kappa="NA"
				try:
					score_cars_kappa_ss=sum(score_cars_list_kappa_ss) / len(score_cars_list_kappa_ss)
				except:
					score_cars_kappa_ss="NA"

				score_cars_yr=sum(score_cars_list_yr) / len(score_cars_list_yr)

				score_cars_ncpr=sum(score_cars_list_ncpr) / len(score_cars_list_ncpr)

				try:
					score_cars_aggrescan=sum(score_cars_list_aggrescan) / len(score_cars_list_aggrescan)
				except:
					score_cars_aggrescan="NA"

				try:
					if "NA" in score_cars_list_anchor:
						score_cars_list_anchor=score_cars_list_anchor.remove("NA")
					score_cars_anchor=sum(score_cars_list_anchor) / len(score_cars_list_anchor)
				except:
					score_cars_anchor="NA"

				#write results
				write_csv(uniprot, d[uniprot]["Dataset"], "CARs", score_cars_kappa, score_cars_kappa_ss, score_cars_yr, score_cars_ncpr, score_cars_aggrescan, score_cars_anchor)

		elif element == "PrLDs":
			if not value:
				write_csv(uniprot, d[uniprot]["Dataset"], "PrLDs", "NA", "NA", "NA", "NA", "NA", "NA")

			else:
				#create empty lists for all properties:
				score_prlds_list_kappa=[]
				score_prlds_list_kappa_ss=[]
				score_prlds_list_yr=[]
				score_prlds_list_ncpr=[]
				score_prlds_list_aggrescan=[]
				score_prlds_list_anchor=[]

				for seq in value:
					#kappa
					scorep_prlds_kappa=kappa(seq)
					score_prlds_list_kappa.append(scorep_prlds_kappa)

					#kappa_ss
					scorep_prlds_kappa_ss=kappa_ss(seq)
					score_prlds_list_kappa_ss.append(scorep_prlds_kappa_ss)

					#yr
					scorep_prlds_yr=yr(seq)
					score_prlds_list_yr.append(scorep_prlds_yr)

					#ncpr
					scorep_prlds_ncpr=ncpr(seq, 7)
					score_prlds_list_ncpr.append(scorep_prlds_ncpr)


					#aggrescan
					scorep_prlds_aggrescan=aggrescan(seq)
					score_prlds_list_aggrescan.append(scorep_prlds_aggrescan)


					#anchor
					scorep_prlds_anchor=anchor(seq)
					score_prlds_list_anchor.append(scorep_prlds_anchor)
				
				#final averaged values for each protein/property
				try:
					score_prlds_kappa=sum(score_prlds_list_kappa) / len(score_prlds_list_kappa)
				except:
					score_prlds_kappa="NA"
				try:
					score_prlds_kappa_ss=sum(score_prlds_list_kappa_ss) / len(score_prlds_list_kappa_ss)
				except:
					score_prlds_kappa_ss="NA"

				score_prlds_yr=sum(score_prlds_list_yr) / len(score_prlds_list_yr)

				score_prlds_ncpr=sum(score_prlds_list_ncpr) / len(score_prlds_list_ncpr)

				try:
					score_prlds_aggrescan=sum(score_prlds_list_aggrescan) / len(score_prlds_list_aggrescan)
				except:
					score_prlds_aggrescan="NA"

				try:
					if "NA" in score_prlds_list_anchor:
						score_prlds_list_anchor=score_prlds_list_anchor.remove("NA")
					score_prlds_anchor=sum(score_prlds_list_anchor) / len(score_prlds_list_anchor)
				except:
					score_prlds_anchor="NA"

				#write results
				write_csv(uniprot, d[uniprot]["Dataset"], "PrLDs", score_prlds_kappa, score_prlds_kappa_ss, score_prlds_yr, score_prlds_ncpr, score_prlds_aggrescan, score_prlds_anchor)
		
		elif element == "Full_seq":
			seq=value

			score_fullseq_kappa=kappa(seq)

			#kappa_ss
			score_fullseq_kappa_ss=kappa_ss(seq)

			#yr
			score_fullseq_yr=yr(seq)

			#ncpr
			score_fullseq_ncpr=ncpr(seq, 7)

			#aggrescan
			score_fullseq_aggrescan=aggrescan(seq)

			#anchor
			score_fullseq_anchor=anchor(seq)
			#try:
			#	score_fullseq_anchor=sum(score_idrs_list_anchor) / len(score_idrs_list_anchor)
			#except:
			#	score_fullseq_anchor="NA"

			write_csv(uniprot, d[uniprot]["Dataset"], "Full_seq", score_fullseq_kappa, score_fullseq_kappa_ss, score_fullseq_yr, score_fullseq_ncpr, score_fullseq_aggrescan, score_fullseq_anchor)


		elif element == "Dataset":
			pass
		
		else:
			print("Unexpected key dict: " + str(element))

