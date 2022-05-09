# Docking-Scripts
scripts for download data and use autodock vina in python. <br />
in this repository I prepared codes to download and Dock compounds to protein, automaticly.

## download data
for Docking we need a receptor(likd proteins) and a ligand(like compound).<br />
in this repository I use [rcsb protein databank](https://www.rcsb.org/) for download receptros and [ZINC](https://zinc.docking.org/) database for dowloand ligands.<br />
I search in ZINC by inchikey, and in rcsb by gene_names. any query can be wirte and send to these database to search and find desired match.

## Docking
I use [Autodock Vina](https://vina.scripps.edu/) to dock protein-ligand pairs.
Autodock vina has [python implementation](https://autodock-vina.readthedocs.io/en/latest/introduction.html).

### requiremnets
- [vina](https://autodock-vina.readthedocs.io/en/latest/installation.html)
- [meeko](https://autodock-vina.readthedocs.io/en/latest/docking_requirements.html)
- [ADFR](https://ccsb.scripps.edu/adfr/downloads/)
- urllib
- requests
- subprocess
