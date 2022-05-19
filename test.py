import pandas as pd
from util import utils
from vina import Vina


orf = pd.read_csv('Data/orf_pairs.csv')
crispr = pd.read_csv('Data/crispr_pairs.csv')

sorted_orf = orf.sort_values(by=['correlation'], ascending=False)
sorted_crispr = crispr.sort_values(by=['correlation'], ascending=False)
# print('-----------------------------------------------------')
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

if __name__ == '__main__':

    metadata = pd.read_csv('Data/metadata/zinc_metadata_compound.csv')
    receptor_meta = pd.read_csv('Data/metadata/receptor_metadata.csv')

    zincids = list(metadata["ZINCid"])
    proteins = list(receptor_meta["rcsb_id"])

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
                utils.docking(receptor_true, ligand)
            except:
                pass

            print('***********************************************************')
            print('**** docking predicted pair *******************************')
            print('### ' ,receptor_predicted, ligand)
            print('***********************************************************')
            try:
                utils.docking(receptor_predicted, ligand)
            except:
                pass

        print('##################################################################')
        print('##################################################################')

