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

metadata = pd.read_csv('../2021_Chandrasekaran_submitted/metadata/external_metadata/JUMP-Target-1_compound_metadata.tsv', sep='\t')

orf = pd.read_csv('Data/orf_pairs.csv')
crispr = pd.read_csv('Data/crispr_pairs.csv')

sorted_orf = orf.sort_values(by=['correlation'], ascending=False)
sorted_crispr = crispr.sort_values(by=['correlation'], ascending=False)
# print('-----------------------------------------------------')
def download_all_protein(genes):
    for g in genes:
        print(g)
        try:
            i = utils.get_rcsb_id(g)
        except:
            i = []
        print(i)
        print('-------------------------')
    # protiens = [utils.get_rcsb_id(g) for g in genes]
    # for p in protiens:
    #     utils.download_pdb(p)
    # utils.prepare_receptor(protiens)

def download_all_ligands(inchikeys):
    for c in inchikeys:
        print(c)
        try:
            zinc = utils.get_zincid(c)
        except:
            zinc = []
        print(zinc)
    # for zincid in zincids:
    #     utils.download_with_zincid(zincid, 'sdf')
    # utils.prepare_ligand(zincids)

def docking(receptor, ligand):
    v = Vina(sf_name='vina')
    v.set_receptor(f'receptor/{receptor}.pdbqt')

    v.set_ligand_from_file(f'ligands/{ligand}.pdbqt')
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
    genes = list(sorted_crispr["target_gene_sample_1"])
    genes += list(sorted_crispr["target_gene_sample_2"])
    # print(list(set(genes))
    download_all_protein(set(genes))
    print('all prteins prepared')
    ### download all braod samples as ligans
    inchikeys = [metadata[metadata["broad_sample"] == broad_sample]["InChIKey"].item() for broad_sample in sorted_crispr["broad_sample_modality_1"]]
    download_all_ligands(set(inchikeys))
    print('all ligands prepared')
    ### process for docking
    for i in range(len(sorted_crispr)):
        broad_sample = sorted_crispr.iloc[i]["broad_sample_modality_1"]
        gene_1 = sorted_crispr.iloc[i]["target_gene_sample_1"].split("|")[0]
        gene_2 = sorted_crispr.iloc[i]["target_gene_sample_2"].split("|")[0]
        inchikey = metadata[metadata["broad_sample"] == broad_sample]["InChIKey"].item()

        print(f'{metadata[metadata["broad_sample"]==broad_sample]["pert_iname"]} target gene is {gene_1}\nour result show high correlation with gene {gene_2}.')

        receptor_true = utils.get_rcsb_id(gene_1)
        receptor_predicted = utils.get_rcsb_id(gene_2)
        ligand = utils.get_zincid(inchikey)

        print('docking true pair ******************************')
        docking(receptor_true, ligand)

        print('************************************************')
        print('docking predicted pair *************************')
        docking(receptor_predicted, ligand)

