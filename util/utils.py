import urllib.parse
import urllib.request
import requests
import subprocess
import json
import os
from vina import Vina


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
    if len(zinc_ids) == 0:
        print('no match found')
    return zinc_ids[0]


def get_rcsb_id(gene):
    ### create a query for search in rcsb database and return rcsb_id
    ### this query is based on gene name in homo sapiens
    ### this function returns first output rcsb_id for query
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
    try:
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
    except:
        print('no match found')
        return []

def download_with_zincid(zincid, fr='sdf'):
    ### download sample from zinc15
    ### zincid: str
    ### fr: str - format of sample that you want to downalod from zinc15
    ###           default is "sdf"
    ####################################################################
    if os.path.isdir(f'Data/ligands/{zincid}.{fr}') == False:
        print('download ligand' , zincid)
        urllib.request.urlretrieve(f'https://zinc15.docking.org/substances/{zincid}.{fr}', f'Data/ligands/{zincid}.{fr}')

def download_pdb(rcsbid):
    ### download pbd file of protein from rcsb database with rcsb_id
    rcsbid = str(rcsbid)
    if os.path.isdir(f'Data/receptor/{rcsbid}.pdb') == False:
        print('download pdb receptor', rcsbid)
        urllib.request.urlretrieve(f'https://files.rcsb.org/download/{rcsbid}.pdb', f'Data/receptor/{rcsbid}.pdb')

def prepare_receptor(pdb_list):
    ### convert pdb to pdbqt
    ### input: list of pdb
    ### save outputs at Data/receptor/
    for pdb in pdb_list:
        if os.path.isdir(f"Data/receptor/{pdb}.pdbqt"):
            continue
        print('preparing ', pdb)
        command = ["prepare_receptor", "-r" , f"Data/receptor/{pdb}.pdb", "-o", f"Data/receptor/{pdb}.pdbqt"]
        subprocess.run(command)

def prepare_ligand(ligand_list):
    ### convert sdf to pdbqt
    ### input: list of liagands
    ### save outputs at Data/ligands/
    for ligand in ligand_list:
        if os.path.isdir(f"Data/ligands/{ligand}.pdbqt"):
            continue
        print('preparing ', ligand)
        command = ["mk_prepare_ligand.py", "-i", f"Data/ligands/{ligand}.sdf", "-o", f"Data/ligands/{ligand}.pdbqt"]
        subprocess.run(command)

def download_all_protein(genes):
    ### input list of genes, and download pdb file for each gene
    protiens = {}
    for g in genes:
        try:
            i = get_rcsb_id(g)
            download_pdb(i)
            protiens[g] = i
        except:
            pass
    print(list(protiens.values()))
    prepare_receptor(list(protiens.values()))
    return protiens

def download_all_ligands(inchikeys):
    ### input list of inchikeys, and download sdf file for each inchikey
    zincids = {}
    for c in inchikeys:
        try:
            zinc = get_zincid(c)
            download_with_zincid(zinc, 'sdf')
            zincids[c] = zinc
        except:
            pass
    prepare_ligand(list(zincids.values()))
    return zincids

def docking(receptor, ligand, exhaustiveness=32, n_poses=20, maps_center=[0,0,0], box_size=[20,20,20]):
    v = Vina(sf_name='vina')
    v.set_receptor(f'Data/receptor/{receptor}.pdbqt')

    v.set_ligand_from_file(f'Data/ligands/{ligand}.pdbqt')
    v.compute_vina_maps(center=maps_center, box_size=[20,20,20])

    # Score the current pose
    #energy = v.score()
    # print('Score before minimization: %.3f (kcal/mol)' % energy[0])

    # Minimized locally the current pose
    #energy_minimized = v.optimize()
    # print('Score after minimization : %.3f (kcal/mol)' % energy_minimized[0])
    # v.write_pose('1iep_ligand_minimized.pdbqt', overwrite=True)

    # Dock the ligand
    v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
    v.write_poses('Data/Docks/1iep_ligand_vina_out.pdbqt', n_poses=n_poses, overwrite=True)