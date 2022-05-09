import pandas as pd
from util import utils
from vina import Vina


# batch = "2020_11_04_CPJUMP1"
# experiment_df = (
#     pd.read_csv('output/experiment-metadata.tsv', sep='\t')
#     .query('Batch==@batch')
#     .query('Density==100')
#     .query('Antibiotics=="absent"')
# )

#metadata = pd.read_csv('../2021_Chandrasekaran_submitted/metadata/external_metadata/JUMP-Target-1_compound_metadata.tsv', sep='\t')

orf = pd.read_csv('Data/orf_pairs.csv')
crispr = pd.read_csv('Data/crispr_pairs.csv')

sorted_orf = orf.sort_values(by=['correlation'], ascending=False)
sorted_crispr = crispr.sort_values(by=['correlation'], ascending=False)
# print('-----------------------------------------------------')
def download_all_protein(genes):
    protiens = {}
    for g in genes:
        print(g)
        try:
            i = utils.get_rcsb_id(g)
            utils.download_pdb(i)
            protiens[g] = i
        except:
            i = []
        print(i)
        print('-------------------------')
    print(list(protiens.values()))
    utils.prepare_receptor(list(protiens.values()))
    return protiens

def download_all_ligands(inchikeys):
    zincids = {}
    for c in inchikeys:
        print(c)
        try:
            zinc = utils.get_zincid(c)
            utils.download_with_zincid(zinc, 'sdf')
            zincids[c] = zinc
        except:
            zinc = []
        print(zinc)
    utils.prepare_ligand(list(zincids.values()))
    return zincids

def generate_zincids_csv(metadata):
    metadata["ZINCid"] = None
    for i in range(len(metadata)):
        inchikey = metadata.iloc[i]["InChIKey"]
        try:
            zinc = utils.get_zincid(inchikey)
            print(inchikey, zinc)
            metadata.loc[i, "ZINCid"] = zinc
        except:
            pass
    metadata.to_csv('Data/metadata/zinc_metadata_compound.csv', index=False)

def generate_protien_csv(genes):
    protein_df = {'gene': [], 'rcsb_id': []}
    protein_df = pd.DataFrame(data=protein_df)
    for gene in set(genes):
        try:
            rcsb = utils.get_rcsb_id(gene)
            protein_df = protein_df.append({'gene': gene, 'rcsb_id': rcsb}, ignore_index=True)

        except:
            protein_df = protein_df.append({'gene': gene, 'rcsb_id': None}, ignore_index=True)
    protein_df.to_csv('Data/metadata/receptor_metadata.csv', index=False)

def docking(receptor, ligand):
    v = Vina(sf_name='vina')
    v.set_receptor(f'Data/receptor/{receptor}.pdbqt')

    v.set_ligand_from_file(f'Data/ligands/{ligand}.pdbqt')
    v.compute_vina_maps(center=[0, 0, 0], box_size=[20, 20, 20])

    # Score the current pose
    energy = v.score()
    # print('Score before minimization: %.3f (kcal/mol)' % energy[0])

    # Minimized locally the current pose
    energy_minimized = v.optimize()
    # print('Score after minimization : %.3f (kcal/mol)' % energy_minimized[0])
    # v.write_pose('1iep_ligand_minimized.pdbqt', overwrite=True)

    # Dock the ligand
    v.dock(exhaustiveness=32, n_poses=5)
    # v.write_poses('1iep_ligand_vina_out.pdbqt', n_poses=5, overwrite=True)

if __name__ == '__main__':
    ### download all proteins as receptors
    # print(list(set(genes))
    # download_all_protein(set(genes))
    # print('all prteins prepared')
    ### download all braod samples as ligans
    # inchikeys = [metadata[metadata["broad_sample"] == broad_sample]["InChIKey"].item() for broad_sample in sorted_crispr["broad_sample_modality_1"]]
    # download_all_ligands(set(inchikeys))
    # print('all ligands prepared')
    ### process for docking
    # metadata["ZINCid"] = None

    metadata = pd.read_csv('Data/metadata/zinc_metadata_compound.csv')
    receptor_meta = pd.read_csv('Data/metadata/receptor_metadata.csv')

    zincids = list(metadata["ZINCid"])
    proteins = list(receptor_meta["rcsb_id"])

    # for z in set(zincids):
    #     z = str(z)
    #     if z=='nan':
    #         continue
    #     utils.download_with_zincid(z , 'sdf')
    # utils.prepare_ligand(zincids)
    # for p in set(proteins):
    #     if p=='nan':
    #         continue
    #     utils.download_pdb(p)
    # utils.prepare_receptor(proteins)
    #
    # print('all is dowloaded........')

    for i in range(len(sorted_crispr)):
        broad_sample = sorted_crispr.iloc[i]["broad_sample_modality_1"]
        gene_1 = sorted_crispr.iloc[i]["target_gene_sample_1"].split("|")[0]
        gene_2 = sorted_crispr.iloc[i]["target_gene_sample_2"].split("|")[0]
        zincid = metadata[metadata["broad_sample"] == broad_sample]["ZINCid"].item()

        print(f'### {metadata[metadata["broad_sample"]==broad_sample]["pert_iname"].item()} targhet gene is {gene_1}\n### our result show high correlation with gene {gene_2}.')
    #
        receptor_true = receptor_meta[receptor_meta["gene"]==gene_1]["rcsb_id"].item()
        receptor_predicted = receptor_meta[receptor_meta["gene"]==gene_2]["rcsb_id"].item()
        ligand = zincid
        check = True

        if receptor_true==None:
            print('gene_1 is None. ')
            check=False
        if receptor_predicted==None:
            print('gene_2 is None. ')
            check = False
        if ligand==None:
            print('compound is None.')
            check = False


        if check:
            print('**** docking true pair ************************************')
            print('### ' ,receptor_true , ligand)
            print('***********************************************************')
            try:
                docking(receptor_true, ligand)
            except:
                pass

            print('***********************************************************')
            print('**** docking predicted pair *******************************')
            print('### ' ,receptor_predicted, ligand)
            print('***********************************************************')
            try:
                docking(receptor_predicted, ligand)
            except:
                pass

        print('##################################################################')
        print('##################################################################')

# metadata.to_csv('metadata_with_zinc.csv', index=False)
