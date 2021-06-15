# ps2-npjBM

This repo contains all the files and code needed to reproduce analyses and figures from ""Microbial predictors of healing and short-term effect of debridement on the microbiome of chronic wounds" in npj B&M


Contents are as follows:
1. data
	‘ps2mapping.txt’ - metadata mapping file
	‘rep_set.tre’ - phylogenetic tree file
	‘silva.sorted.raw.biom’ - OTU table w/ Silva taxonomy in BIOM format
	‘staph.taxtable.forps.csv’ - Species-level taxonomic assignment from Blast for Staph OTUs 
2. sv_code
	‘PS2.DESeq2.O2req.CG.Rmd’ - R Notebook for association of O2 requirements
	‘PS2.DESeq2.Rmd’ - R Notebook for taxonomic association analysis using DESeq2
	‘PS2.Revised.Publication.Figures.Dryad.Rmd’ R Notebook of nearly all OTU table analyses 		and figure generation.
3. ys_code
	‘data’ - contains data from bglmm, blast, deseq2, otu tables, and sequence preparation
	‘scripts’ - contains a README for ‘ys_code’ and the figure generation Python notebook
	‘src’ - scripts
