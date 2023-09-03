#%%
from requests import get, Session
import gzip
import os
import gzip
import json
import shutil
import xmltodict
from subprocess import Popen, PIPE, run, check_output
import time
import argparse
import re
#%%
argparser = argparse.ArgumentParser(description='Update refdb')
argparser.add_argument('-db', '--database', help='Database to update', default="refseq")
arg, _ = argparser.parse_known_args()
#%%
def gbffgz_download(gbff_URI, des):
    #URI = "https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Fungi/fungi.ITS.gbff.gz"
    name = gbff_URI.split("/")[-1]
    basename, _ = os.path.splitext(name)
    print(f"Downloading {name} database from NCBI refseq ftp...")
    r = get(gbff_URI, allow_redirects=True)
    open(f"{des}/{name}", 'wb').write(r.content)
    if name.endswith(".gz"):
        #Extract gbff.gz
        print("Extracting gbff.gz file...")
        with gzip.open(f"{des}/{name}", 'rb') as f_in:
            with open(f"{des}/{basename}", 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        #Delete gbff.gz
        os.remove(f"{des}/{name}")
        return f"{des}/{basename}"
    else:
        return f"{des}/{name}"
#%%
def gbffgz_to_taxfas(gbff_path, des, gene=None):
    name = gbff_path.split("/")[-1]
    name, _ = os.path.splitext(name)
    #Get taxinfo for each record
    print("Getting taxinfo for each record...")
    taxinfos = {}
    taxid_list = set()
    for rec in gbff_reader(open(gbff_path, 'r'), gene=gene):
        taxid_list.add(rec["taxid"])
        if len(taxid_list) % 1000 == 0:
            print(f"{len(taxid_list)} taxid loaded...")
    #Show how many taxid to process
    print(f"{len(taxid_list)} taxid to process...")
    #Retrieve taxon info by taxid
    batch = 150
    for i in range(0, len(taxid_list), batch):
        taxinfo_batch = lineage_by_taxid(list(taxid_list)[i:i+batch])
        if taxinfo_batch != None:
            taxinfos.update(taxinfo_batch)
        else:
            #Try again
            taxinfo_batch = lineage_by_taxid(list(taxid_list)[i:i+batch])
            if taxinfo_batch != None:
                taxinfos.update(taxinfo_batch)
            else:
                print("Error when retrieving taxinfo, skip this batch")
        print(f"{len(taxinfos)}/{len(taxid_list)} taxid processed...", end="\r")
    print(f"{len(taxinfos)}/{len(taxid_list)} taxid processed...")
    #write fasta
    with open(f"{des}/{name}.fas", 'w') as f:
        for rec in gbff_reader(open(gbff_path, 'r'), gene=gene):
            try:
                lineage = ";".join([taxinfos[rec["taxid"]][i] for i in ["kingdom", "phylum", "class", "order", "family", "genus"]])
            except Exception as e:
                lineage = ";".join(["Unclassified"]*6)
            title = "{}||{}||{}||{}".format(rec["accession"], rec["organism"], lineage, rec["taxid"])
            title = title.replace(" ", "_")
            try:
                f.write(">{}\n{}\n".format(title, rec["seq"]))  
            except Exception as e:
                print(rec)
                print(e)
                pass
    return f"{des}/{name}.fas"   
#%%
def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 
                  'N': 'N', 'R': 'Y', 'Y': 'R', 'S': 'S', 
                  'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 
                  'D': 'H', 'H': 'D', 'V': 'B', 'X': 'X',
                  'a': 't', 'c': 'g', 'g': 'c', 't': 'a',
                        'n': 'n', 'r': 'y', 'y': 'r', 's': 's',
                        'w': 'w', 'k': 'm', 'm': 'k', 'b': 'v',
                        'd': 'h', 'h': 'd', 'v': 'b', 'x': 'x'}
    return "".join([complement[base] for base in reversed(seq)])
#%%
def gbff_reader(handle, gene = None):
    line = handle.readline()
    data = {}
    gene_loc_tmp = ""
    while line:
        if line.startswith("ACCESSION"):
            data["accession"] = line.split()[1]
        if line.startswith("  ORGANISM  "):
            data["organism"] = line.replace("  ORGANISM  ", "").strip()
        if line.startswith('                     /db_xref="taxon:'):
            data["taxid"] = line.replace('                     /db_xref="taxon:', "").replace('"', "").strip()
        if gene != None:
            #Save gene location 
            if line.startswith("     gene            "):
                gene_loc_tmp = line.replace("     gene            ", "").strip()
            #Check if this gene is the gene we want
            if line.startswith("                     /gene="):
                gene_name = line.replace('                     /gene="', "").replace('"', "").strip()
                if gene_name == gene and gene_loc_tmp != "":
                    #Save gene location
                    data['gene_loc'] = gene_loc_tmp
                    #print(f"gene found for {data['accession']}")
        if line.startswith("ORIGIN"):
            line = handle.readline()
            seq = ""
            while not line.startswith("//"):
                seq += line.replace(" ", "").replace("\n", "").replace("\r", "").replace("0", "").replace("1", "").replace("2", "").replace("3", "").replace("4", "").replace("5", "").replace("6", "").replace("7", "").replace("8", "").replace("9", "").upper()
                line = handle.readline()
            data["seq"] = seq
        
        if line.startswith("//"):
            if gene != None:
                if "gene_loc" in data:
                    #Extract the gene sequence
                    if "complement" in data["gene_loc"]:
                        gene_loc = data["gene_loc"].replace("complement(", "").replace(")", "").replace(">","").replace("<","").split("..")
                        data['seq'] = reverse_complement(data["seq"][int(gene_loc[0])-1:int(gene_loc[1])])
                    else:
                        gene_loc = data["gene_loc"].replace(">", "").replace("<","").split("..")
                        data['seq'] = data["seq"][int(gene_loc[0])-1:int(gene_loc[1])]
                else:
                    #If this record does not contain the gene we want, skip it
                    data = {}
                    gene_loc_tmp = ""
                    line = handle.readline()
                    continue
            yield data
            data = {}    
            gene_loc_tmp = ""
        line = handle.readline()
#%%
def lineage_by_taxid( taxid=['3016022', '2888342'], tax_id_cache = f"./taxid_cache/taxid_cache.json"):
    #try load tax_id_cache from cache file
    ranks = {} 
    try:
        taxid_json = json.load(open(tax_id_cache, 'r'))
    except:
        taxid_json = {}
    #Add requested taxid to ranks if it is in cache
    #And remove it from taxid list
    new_taxid = []
    for t in taxid:
        if t in taxid_json:
            ranks[t] = taxid_json[t]
        else:
            new_taxid.append(t)
    #If all taxid are in cache, return ranks
    if len(new_taxid) == 0:
        return ranks
    
    # Get taxon info by taxid
    taxid_info_URI = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&rettype=xml&id={}'.format(",".join(new_taxid))
    retry = 4
    while True:
        #Wait for 1 second before retry
        time.sleep(1)
        try:
            r = get(taxid_info_URI)
            r = xmltodict.parse(r.text)
            break
        except Exception as e:
            print(taxid_info_URI,e)
            #print(r.text)
            if retry == 0:
                return None
            retry -= 1
    rank_base = {"kingdom":"incertae sedis", "phylum":"incertae sedis", "class":"incertae sedis", "order":"incertae sedis", "family":"incertae sedis", "genus":"incertae sedis"}
    try:
        
        #If there is only one taxid, the result will be a dict
        if isinstance(r["TaxaSet"]["Taxon"], dict):
            r["TaxaSet"]["Taxon"] = [r["TaxaSet"]["Taxon"]]
        for query in r["TaxaSet"]["Taxon"]:
            
            rank = rank_base.copy()
            #Get taxon of this taxid. Because the taxon of this taxid is not included in the LineageEx
            tid_name = query['ScientificName']
            tid_rank = query['Rank']
            if tid_rank == 'superkingdom':
                tid_rank = 'kingdom'
            if tid_rank in rank.keys():
                rank[tid_rank] = tid_name
            #If there is only one lineage, the result will be a dict
            if isinstance(query['LineageEx']['Taxon'], dict):
                query['LineageEx']['Taxon'] = [query['LineageEx']['Taxon']]
            for i in query['LineageEx']['Taxon']:
                if i['Rank'] == 'superkingdom':
                    #For Bacteria and Archaea
                    i['Rank'] = 'kingdom'
                if i['Rank'] in rank.keys():
                    rank[i['Rank']] = i['ScientificName']
            ranks[query['TaxId']] = rank
            #Save taxid info to cache
            taxid_json[query['TaxId']] = rank
            
    except Exception as e:
        print(taxid_info_URI)
        print(r, e)
        pass
    #Save taxid info to cache
    
    json.dump(taxid_json, open(tax_id_cache, 'w'))
    return ranks
def _exec(cmd,suppress_output=True):
    if suppress_output:
        with open(os.devnull, 'w') as DEVNULL:
            out = run(cmd, stdout=DEVNULL, stderr=DEVNULL, shell=True)
        return None
    else:
        out = run(cmd, stdout=PIPE, shell=True)
        print (">>", cmd)
        print("Output:")
        print(out.stdout.decode('utf-8'))
        print("Exception:")
        print(out.stderr)
        return out.stdout.decode('utf-8'), out.stderr   
# %%
gbff_list = ["https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Fungi/fungi.ITS.gbff.gz",
             "https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Fungi/fungi.28SrRNA.gbff.gz",
             "https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Fungi/fungi.18SrRNA.gbff.gz",
             "https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.gbff.gz",
             "https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.23SrRNA.gbff.gz",
             "https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.5SrRNA.gbff.gz",
             "https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Archaea/archaea.16SrRNA.gbff.gz",
             "https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Archaea/archaea.23SrRNA.gbff.gz",
             "https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Archaea/archaea.5SrRNA.gbff.gz",
             ]
refdb_path = f"./refdb/"
refdb_path = os.path.abspath(refdb_path)
# %%
if arg.database == "refseq":
    for gbff_URI in gbff_list[:]:
        gbff_path = gbffgz_download(gbff_URI=gbff_URI, 
                            des=refdb_path
                            )
        gbffgz_to_taxfas(gbff_path=gbff_path,
                            des = refdb_path
                                )
        os.remove(gbff_path)

# %%
#Download query result from ncbi
def gbff_from_query(query="internal transcribed spacer AND Olpidium[Organism] AND biomol_genomic [PROP] AND is _nuccore [filter] 300:1000[SLEN]",
                    gbff_path = "./refdb/nr_Olpidium.ITS.gbff"):
    URI = f"https://www.ncbi.nlm.nih.gov/nuccore/?term={query}"
    print(URI)
    session = Session()
    #Query ncbi
    query = session.get(URI)
    #save cookie
    cookies = query.cookies
    #Get ncbi_phid
    ncbi_phid = re.findall(r"ncbi_phid=(\w+)\"", query.text)[0]
    print(f"ncbi_phid:{ncbi_phid}")
    #Get gbff
    if ncbi_phid == None:
        print("ncbi_phid not found")
    else:
        print("ncbi_phid found")
        response = session.get(f"https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=gbwithparts&query_key=1&filter=all&extrafeat=undefined&withparts=on&ncbi_phid={ncbi_phid}",
                               cookies=cookies,
                               stream=True
                               )
        with open(gbff_path, 'wb') as f:
            downloaded = 0
            for chunk in response.iter_content(chunk_size=1024*10):
                if chunk:
                    f.write(chunk)
                    downloaded += len(chunk)
                    print(f"Downloaded {downloaded/1024/1024} MB", end="\r")


    return None

#%%
##Build database from ncbi query
def refdb_from_query(query="rbcL[Gene%20Name]%20AND%20plants[filter]%20AND%20biomol_genomic[PROP]%20AND%20srcdb_refseq[PROP]",
                     des ="./refdb",
                     name="plant_rbcl",
                     gene="rbcL",
                     batch = 200,
                     split_per_seq = 50000
                     ):
    #Download plant rbcl from refseq
    URI = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?RetMax=900000&db=nucleotide&term={query}"
    r = get(URI)
    r = xmltodict.parse(r.text)
    idlist = r["eSearchResult"]["IdList"]["Id"]
    id_no = len(idlist)
    print(f"{id_no} records found")
    #Split idlist into files
    splits = 0
    split_count = 0
    while splits < len(idlist):
        print(f"Split {split_count}/{int(len(idlist)/split_per_seq)}")
        idlist_split = idlist[splits:splits+split_per_seq]
        gbff_path = f"{des}/{name}_{split_count}.gbff"
        #Download sequences in batch
        debug_count=1000
        with open(gbff_path, 'w') as f:
            f.write("")
            for i in range(0, len(idlist_split), batch):
                idlist_batch = idlist_split[i:i+batch]
                idlist_batch = ",".join(idlist_batch)
                URI = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=gb&id={idlist_batch}"
                r = get(URI)
                with open(gbff_path, 'a') as f:
                    f.write(r.text)
                print(f"{i+batch}/{len(idlist_split)} processed...")
                debug_count -= 1
                if debug_count == 0:
                    pass
                    #break
            gbffgz_to_taxfas(gbff_path, "./refdb", gene=gene)
            os.remove(gbff_path)
        splits += split_per_seq
        split_count += 1
# %%
#Download plant rbcL gene from refseq
if arg.database == "plant_rbcl": #79058 records
    refdb_from_query(query="rbcL AND plants [filter] AND biomol_genomic [PROP] AND is _nuccore [filter] 1000:3000[SLEN] ",
                        des ="./refdb",
                        name="plant_rbcl",
                        gene=None)

if arg.database == "plant_rbcl_new":
    #Download gbff
    gbff_from_query(query="rbcL[Gene%20Name]%20AND%20plants[filter]%20AND%20biomol_genomic[PROP]%20AND%20srcdb_refseq[PROP]",
                    gbff_path = "./refdb/plant_rbcl.gbff")
    #Convert gbff to fasta
    gbffgz_to_taxfas(gbff_path="./refdb/plant_rbcl.gbff",
                        des = "./refdb",
                        gene="rbcL")
# %%
#Download plant rbcL gene from refseq
if arg.database == "nr_basidiomycota_ITS":
    refdb_from_query(query="ITS1 AND Basidiomycota[Organism] AND biomol_genomic [PROP] AND is _nuccore [filter] 300:1000[SLEN]",
                        des ="./refdb",
                        name="nr_basidiomycota_ITS",
                        gene=None)   
if arg.database == "nr_ascomycota_ITS":
    refdb_from_query(query="internal transcribed spacer AND Ascomycota[Organism] AND biomol_genomic [PROP] AND is _nuccore [filter] 300:1000[SLEN]",
                        des ="./refdb",
                        name="nr_ascomycota_ITS",
                        gene=None)   
if arg.database == "nr_oomycota_ITS":
    refdb_from_query(query="ITS1 AND Oomycota[Organism] AND biomol_genomic [PROP] AND is _nuccore [filter] 300:1000[SLEN]",
                        des ="./refdb",
                        name="nr_oomycota_ITS",
                        gene=None)
# %%
#Compress each taxonomic fasta file
for tax in os.scandir(refdb_path):
    if tax.name.endswith(".fas"):
        print("Compressing "+tax.name)
        _exec(f"gzip {tax.path} -c > {tax.path}.gz")
        #with open(refdb_path+tax.name, "rb") as f_in:
        #    with gzip.open(refdb_path+tax.name+".gz", "wb") as f_out:
        #        f_out.writelines(f_in)
        #    f_out.close()
        #f_in.close()
        os.remove(tax.path)


# %%
