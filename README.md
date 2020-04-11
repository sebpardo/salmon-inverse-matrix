Estimating marine survival of Atlantic salmon using an inverse matrix approach
=====

This repository contains all the R code for the analyses and figures for the
manuscript titled "Estimating marine survival of Atlantic salmon using an
inverse matrix approach" by Sebastián A. Pardo and Jeffrey A. Hutchings,
currently under review in PLOS ONE.

To recreate the simulation, SIR algorirthm, and figures, run the files in the `analysis/` folder 
sequentially, or use the Makefile. This will populate the `figures/` and `data/` folders with 
the simulations' outputs.
These two folders are not under version control, which is why they're empty in the repository.

### Using the Makefile

To simulate the data:

```sh
make fish
```

To run the sample-importance-resampling algorithm:

```sh
make sir
```

To recreate the figures and LaTeX tables:

```sh
make figs
```

To run the Supporting Information analysis for a declining population and recreate the figures:

```sh
make supp
```
