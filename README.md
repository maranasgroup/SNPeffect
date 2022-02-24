# SNPeffect-

This code corresponds to doi: [10.1111/tpj.14746](https://onlinelibrary.wiley.com/doi/full/10.1111/tpj.14746)

### SNPeffect: identifying functional roles of SNPs using metabolic networks

Input data is stored under data/. main.py makes input files for the GAMS file findSNPs.gms, which (currently) uses CPLEX to identify SNP activities. The following inputs are needed for main.py - 

1) under data/
- AllGenotypesInStudy.txt - list of genotypes in the study
- Biomass_Cons_AllGenotypes.txt - relative growth rate of each genotype wrt the reference genotype
- Biomass_ConsList_AllGenotypes.txt - list of constraints for the above file


2) under data/gsm/
- smbl_metabolites.txt - list of metabolites in the GSM model
- smbl_reactions.txt - list of reactions in the GSM model
- smbl_reaction_bounds.txt - reaction bounds for reactions in the GSM model
- smbl_sij.txt - stoichiometric matrix from the GSM model
- FVA.txt - FVA results using nutrients from medium.txt
- maxbiomass_pFBA_FVA.txt - FVA results at maximal biomass value, followed by total flux limited to pFBA value at maximal biomass
- Poplar_GPRrxns_090419.txt - tab delimited file, with reactions (column 1) and associated list of genes (column 2, comma delimited)

3) under data/metabolomics/  - stores genotype-specific metabolite levels

4) under data/genomics/ - stores genotype-wise sequence data
