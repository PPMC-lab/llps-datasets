import urllib.request, json
import os
import sys

def phasepro_D1():
    def json_parse(data):
        for key, value in data.items():
            if value["partner_dep"] == "N" and value["ptm_dep"] == "N" and value["rna_dep"] == "N":
                with open("D1.txt", "at") as out:
                    out.write(key+"\n")
    
    #Try request in PhaSePro website
    try:   
        with urllib.request.urlopen("https://phasepro.elte.hu/download_full.json") as url:
            # Variable 'data' will contain the full database as a nested dictionary
            data = json.loads(url.read().decode())
            json_parse(data)
    
    #Request failed, use the local file    
    except: 
        with open("phasepro.json", "rt") as json_file:
            data = json.load(json_file)
            json_parse(data)

def phsasepdb_D2():
    uniprots=[]
    with open("phasepdb.csv", "rt") as csv:
        csv.readline() #Skip header
        for line in csv:
            elements=line.split(",")
            #psself proteins with no partners and no regulation (oligomerization allowed)
            if "psself" in elements[0] and elements[7] == "_" and (elements[8] == "_" or elements[8] == "oligomerization"):
                uniprots.append(elements[5])
    
    #Obtain unique uniprots             
    uniprots_filtered=list(set(uniprots))
    for protein in uniprots_filtered:
        with open("D2.txt", "at") as out:
            out.write(protein+"\n")

def llpsdb_D3():
    llpsdb_ids=[]
    with open("llpsdb.csv", "rt") as csv:
        csv.readline()
        for line in csv:
            elements=line.split(",")
            #Natural proteins, one protein component without repears, PTMs or mutations
            if elements[5] == "N" and elements[2] == "protein(1)" and ";" not in elements[3] and elements[8] == "-" and elements[9] == "-" and elements[10] == "-":
                llpsdb_ids.append(elements[3])
    
    #Unique proteins
    llpsdb_ids_filtered=list(set(llpsdb_ids))
    for identifier in llpsdb_ids_filtered:
        with open("llpsdb_proteins.csv", "rt") as prots:
            prots.readline()
            for line in prots:
                identifier2=line.split(",")[0]
                #Parse protein info for uniprot accs
                if identifier == identifier2:
                    uniprot = line.split(",")[5]
                    if uniprot == "": #No associated uniprot
                        continue
                    with open("D3.txt", "at") as out:
                        out.write(uniprot+"\n")

def cdcode_D4_C1():
    #Ask for a token to CD-CODE (https://wiki.cd-code.org/api-doc) and replace it below, after -H 'Authorization: ' 
    os.system("GET   'https://cd-code.org/api/condensates?is_experimental=False&size=500' -H 'Authorization: your_token_here' > biocondensates_all.json")
    uniprots_D4=[]
    uniprots_C1=[]
    try:
        with open("biocondensates_all.json", "rt") as json_file:
            data_dict = json.load(json_file)
            data_list=data_dict["data"]
            for biocondensate in data_list:
                for key, value in biocondensate.items():
                    if key == "proteins":
                        uniprots=value
                        for uniprot in uniprots:
                            try:
                                protein_type=biocondensate["protein_functional_type"][uniprot.strip()]
                                confidence=biocondensate["protein_confidence_score"][uniprot.strip()]
                                #At least in vitro evidence
                                if protein_type == "driver" and confidence >= 3: #driver
                                    if uniprot not in uniprots_D4:
                                        uniprots_D4.append(uniprot)
                                        with open("D4.txt", "at") as out:
                                            out.write(uniprot+"\n")
                                    else:
                                        continue
                                elif "d" in str(protein_type) and confidence >= 3: #driver-wrongly annotated
                                    if uniprot not in uniprots_D4:
                                        uniprots_D4.append(uniprot)
                                        with open("D4.txt", "at") as out:
                                            out.write(uniprot+"\n")
                                    else:
                                        continue
                                elif protein_type == "member" and confidence >= 3: #member
                                    if uniprot not in uniprots_C1:
                                        uniprots_C1.append(uniprot)
                                        with open("C1.txt", "at") as out:
                                            out.write(uniprot+"\n")
                                    else:
                                        continue
                                elif "m" in str(protein_type) and confidence >= 3: #member-wrongly annotated
                                    if uniprot not in uniprots_C1:
                                        uniprots_C1.append(uniprot)
                                        with open("C1.txt", "at") as out:
                                            out.write(uniprot+"\n")
                                    else:
                                        continue
                                else:
                                    continue #does not meet the criteria for confidency
                            except:
                                continue #Not annotated
    #Token missing                             
    except ValueError:
        print("Request to generate CDCODE-related datasets (D4 and C1) has failed.\n\nAccess to CD-CODE datasets require previous notification to authors. You need to ask first for a private authorization token to CD-CODE developers (https://wiki.cd-code.org/api-doc) and substitute it in the first line of the function 'cdcode_D4_C1'.\n\nRun it again with a valid token to obtain the complete list of datasets.")
        sys.exit()

def drllps_D5_C2():
    ###D5
    with open("full-scaffold.txt", "rt") as tab:
        tab.readline()
        uniprot_list=[]
        for line in tab:
            uniprot_id=line.split("\t")[3].strip('"').strip("\n").rstrip('"')
            uniprot_list.append(uniprot_id)

    uniprot_unique=set(uniprot_list) #Unique UniProt Accs.

    for protein in uniprot_unique:
        with open("D5.txt", "at") as drllps:
            drllps.write(protein+"\n")

    ###C2
    from collections import Counter

    proteins=[] #List with tupples for each protein

    with open("full-clients.txt", "rt") as tab:
        description_list=[]
        tab.readline()
        for line in tab:
            proteins.append((line.split("\t")[0].strip('"'), line.split("\t")[1].strip('"'), line.split("\t")[2].strip('"'), line.split("\t")[3].strip('"').strip('\\').strip('"'), line.split("\t")[4].strip('"'), line.split("\t")[5].strip('\n').strip('"')))
            description_list.append(line.split("\t")[3].strip('"').strip('\\').strip('"'))
    
    c = Counter(description_list) #Count descriptions

    repeated_descriptions=[]
    for element in c.most_common(9000):
        if element[1] > 10: #We consider it high throuput or tissue not defined, next
            repeated_descriptions.append(element[0])
        else:
            continue #good one

    uniprots_final=[]
    for element in proteins:
        if element[3] in repeated_descriptions or element[4] == "N/A": #Tissue/Cell not defined
            continue
        else:
            with open("C2.txt", "at") as drllps:
                drllps.write(element[1]+"\n")

def right_side(C1, C2, D1, D2, D3, D4, D5):
    #client confidence
    elements=[]
    for element in C1:
        elements.append(element)
        if element in C2:
            with open("C+.txt", "at") as Cplus:
                Cplus.write(element)
        else:
            with open("C-.txt", "at") as Cminus:
                Cminus.write(element) 
    for element in C2:
        if element not in elements:   
            with open("C-.txt", "at") as Cminus:
                Cminus.write(element)

    #driver confidence
    uniprots=[]
    integration=D1+D2+D3+D4+D5
    for element in integration:
        if element not in uniprots:
            uniprots.append(element)
        else:
            continue

        n=integration.count(element)
        if n>=3:
            with open("D+.txt", "at") as Dplus:
                Dplus.write(element)
        else:
            with open("D-.txt", "at") as Dminus:
                Dminus.write(element)              


def confidence_clients_drivers():
    #Read all files
    with open("D1.txt", "rt") as d1:
        D1=d1.readlines()
    with open("D2.txt", "rt") as d2:
        D2=d2.readlines()
    with open("D3.txt", "rt") as d3:
        D3=d3.readlines()
    with open("D4.txt", "rt") as d4:
        D4=d4.readlines()
    with open("D5.txt", "rt") as d5:
        D5=d5.readlines()
    with open("C1.txt", "rt") as c1:
        C1=c1.readlines()
    with open("C2.txt", "rt") as c1:
        C2=c1.readlines()
    with open("disprot_thematic.txt", "rt") as dconsensus:
        disprot_thematic=dconsensus.readlines()

    #check clients
    with open("C1.txt", "rt") as c1: #Check CD-CODE
        for line in c1:
            if line not in D1 and line not in D2 and line not in D3 and line not in D4 and line not in D5: #clients only
                with open("CO.txt", "at") as co:
                    co.write(line)
            else: #client but also driver
                with open("C_D.txt", "at") as c_d:
                    c_d.write(line)
   
    with open("CO.txt", "rt") as co:
        CO=co.readlines()

    with open("C_D.txt", "rt") as cd:
        C_D=cd.readlines()


    with open("C2.txt", "rt") as c2: #Check DrLLPS
        for line in c2:
            if line not in D1 and line not in D2 and line not in D3 and line not in D4 and line not in D5 and line not in CO: #clients only
                with open("CO.txt", "at") as co:
                    co.write(line)
            else:
                if line not in C_D and line not in CO:#client but also driver
                    with open("C_D.txt", "at") as c_d:
                        c_d.write(line)

    #only drivers
    with open("D1.txt", "rt") as d1: #Check D1
        for line in d1:
            if line not in C1 and line not in C2: #driver only
                with open("DO.txt", "at") as do:
                    do.write(line)
    with open("DO.txt", "rt") as do:
        DO=do.readlines()

    with open("D2.txt", "rt") as d2: #Check D2
        for line in d2:
            if line not in C1 and line not in C2 and line and line not in DO: #driver only
                with open("DO.txt", "at") as do:
                    do.write(line)

    with open("DO.txt", "rt") as do:
        DO=do.readlines()

    with open("D3.txt", "rt") as d3: #Check D3
        for line in d3:
            if line not in C1 and line not in C2 and line not in DO: #driver only
                with open("DO.txt", "at") as do:
                    do.write(line)
    
    with open("DO.txt", "rt") as do:
        DO=do.readlines()

    with open("D4.txt", "rt") as d4: #Check D4
        for line in d4:
            if line not in C1 and line not in C2 and line not in DO: #driver only
                with open("DO.txt", "at") as do:
                    do.write(line)
    
    with open("DO.txt", "rt") as do:
        DO=do.readlines()

    with open("D5.txt", "rt") as d5: #Check D5
        for line in d5:
            if line not in C1 and line not in C2 and line not in DO: #driver only
                with open("DO.txt", "at") as do:
                    do.write(line)

    right_side(C1, C2, D1, D2, D3, D4, D5)

    return D1, D2, D3, D4, D5, C1, C2, disprot_thematic

def disprot(D1, D2, D3, D4, D5, C1, C2, PDB, disprot_thematic):
    try: #try download from disprot api
        with urllib.request.urlopen("https://disprot.org/api/search?release=2023_06&show_ambiguous=true&show_obsolete=false&format=json") as url:
            # Variable 'data' will contain the full database as a nested dictionary
            data = json.loads(url.read().decode())
    except: #used json file
        with open("disprot_2023_06.json", "rt") as json_file:
            data = json.load(json_file)

    uniprots=[] #all uniprot from disprot
    scaffolds=[] #scaffold list

    for key, value in data.items():
        try:
            for element in value:
                uniprots.append(element["acc"])
                for key2, element2 in element.items():
                    if key2 == "regions":
                        for element3 in element[key2]: #list
                            if element3["term_name"] == "molecular condensate scaffold activity":
                                scaffolds.append(element["acc"])
        except:
            continue

    uniprots_unique = list(set(uniprots))
    scaffolds_unique = list(set(scaffolds))

    final_list=[] #uniprots of DisProt proteins not associated with scaffold activity
    for up in uniprots_unique:
        if up not in scaffolds_unique:
            final_list.append(up)

    #Write final DisProt dataset
    for up in final_list:
        upn=up+"\n"
        if upn not in D1 and upn not in D2 and upn not in D3 and upn not in D4 and upn not in D5 and upn not in C1 and upn not in C2 and upn not in PDB and upn not in disprot_thematic:
            with open("DisProt.txt", "at") as txt:
                txt.write(up+"\n")

def pdb(D1, D2, D3, D4, D5, C1, C2):
    with open("PDB.txt", "rt") as pdb_file:
        for upn in pdb_file:
            if upn not in D1 and upn not in D2 and upn not in D3 and upn not in D4 and upn not in D5 and upn not in C1 and upn not in C2:
                with open("PDB_2.txt", "at") as txt:
                    txt.write(upn)
    with open("PDB_2.txt", "rt") as pdb_file2:
        content=pdb_file2.readlines()
        return content

#Filtering drivers and clients (D1, D2, D3, D4, D5 and C1, C2 datasets)
phasepro_D1()
phsasepdb_D2()
llpsdb_D3()
cdcode_D4_C1() 
drllps_D5_C2()

#Condifent clients and cofident drivers (CO/DO/C_D/C+/C-/D+/D- datasets)
D1, D2, D3, D4, D5, C1, C2, disprot_thematic = confidence_clients_drivers()

#Negative datasets (DisProt and PDB datasets)
PDB=pdb(D1, D2, D3, D4, D5, C1, C2)
disprot(D1, D2, D3, D4, D5, C1, C2, PDB, disprot_thematic)


#Done
print("All datasets have been successfully generated.\n\nThey are saved in your current directory as: 'dataset_name.txt'")