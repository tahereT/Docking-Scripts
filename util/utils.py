import urllib.parse
import urllib.request
import requests
import subprocess
import json
import os


def get_zincid(inchikey):
    ### search in zinc15 database by inchikey and return zincid
    ### input - inchikey: str
    ### return - zinc_id: first match inchikey
    ###########################################################
    url_part1 = 'http://zinc.docking.org/substances/?inchikey='
    url_part3 = ''

    zinc_ids = []

    try:
        response = urllib.request.urlopen('{}{}{}'.format(url_part1, inchikey, url_part3))
    except urllib.error.HTTPError:
        print('Invalid inchikey string {}'.format(smile_str))
        response = []

    for line in response:
        line = line.decode(encoding='UTF-8').strip()

        if line.startswith('<a href="/substances/ZINC'):
            line = line.split('/')[-2]
            zinc_id = urllib.parse.unquote(line)
            zinc_ids.append(str(zinc_id))
    return zinc_ids[0]


def download_with_zincid(zincid, fr):
    ### download sample from zinc15
    ### zincid: str
    ### fr: str - format of sample that you want to downalod from zinc15
    ####################################################################
    if os.path.isdir(f'ligands/{zincid}.{fr}') == False:
        # os.mkdir('ligands')
        urllib.request.urlretrieve(f'https://zinc15.docking.org/substances/{zincid}.{fr}', f'ligands/{zincid}.{fr}')
        #print(f'download {zincid}.{fr} ligand')

def get_rcsb_id(gene):
    ### create a query for search in rcsb database and return rcsb_id
    ### this query is based on gene name in homo sapiens
    ### this function returns first output id for query
    #################################################################
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entity_source_organism.rcsb_gene_name.value",
                        "operator": "exact_match",
                        "value": str(gene)

                    }},
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "operator": "exact_match",
                        "value": "Homo sapiens",
                        "attribute": "rcsb_entity_source_organism.taxonomy_lineage.name"
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "operator": "exact_match",
                        "value": "UniProt",
                        "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name"
                    }
                }
            ]
        },
        "return_type": "polymer_entity"
    }

    myurl = "https://search.rcsb.org/rcsbsearch/v2/query?json="

    req = urllib.request.Request(myurl)
    req.add_header('Content-Type', 'application/json; charset=utf-8')
    jsondata = json.dumps(query)
    jsondataasbytes = jsondata.encode('utf-8')  # needs to be bytes
    req.add_header('Content-Length', len(jsondataasbytes))
    response = urllib.request.urlopen(req, jsondataasbytes)
    with urllib.request.urlopen(req) as f:
        response = f.read()
    my_json = response.decode('utf8').replace("'", '"')
    data = json.loads(my_json)
    return data["result_set"][0]['identifier'][0:4]


def download_pdb(rcsbid):
    ### download pbd file of protein from rcsb database
    rcsbid = str(rcsbid)
    print('download pdb receptor', rcsbid)
    if os.path.isdir(f'receptor/{rcsbid}.pdb') == False:
        #os.mkdir('receptor')
        urllib.request.urlretrieve(f'https://files.rcsb.org/download/{rcsbid}.pdb', f'receptor/{rcsbid}.pdb')

def prepare_receptor(pdb_list):
    for pdb in pdb_list:
        print('processing ', pdb)
        command = ["prepare_receptor","-A", "-r" , f"receptor/{pdb}.pdb" , "-o", f"receptor/{pdb}.pdbqt"]
        subprocess.run(command)

def prepare_ligand(ligand_list):
    for ligand in ligand_list:
        print('processing ', ligand)
        command = ["mk_prepare_ligand.py", "-i", f"ligands/{ligand}.sdf", "-o", f"ligands/{ligand}.pdbqt"]
        subprocess.run(command)

def calculate_docking(receptor, ligand):
    r = "receptor/" + receptor + '.pdbqt'
    l = "ligands/" + ligand + '.pdbqt'
    command = ["vina", l , "--exhaustiveness" , "32" , "--out", r]
    subprocess.run(command)