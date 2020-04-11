# 

fish:
	Rscript analysis/01-simulate-data.R

sir: fish
	Rscript analysis/02-sir-run.R

figs: sir
	Rscript analysis/03-draw-graphs.R

supp:
	Rscript analysis/99-suppmat-declining-popn.R

