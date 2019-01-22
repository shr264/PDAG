# PDAG
This contains all the simulation code for PDAG. The main function used in generating the tables is located in simulations.R. This requires installation of the partitionDAG package, which in turn requires the pdagDFS package.

## Files

A few of the key files in this repo include:
* amat.Rdata
* data_generating_functions.R
* metric_functions.R
* simulations.R
* pDag_dairy_cattle_data.pdf
* pDag_dairy_cattle_data.Rmd

## Details on key Files

### amat.Rdata

This file contains the adjacency matrices for Yeast1, Yeast2, Yeast3, Ecoli1, Ecoli2 from the DREAM3 challenge.

### data_generating_functions.R

This file contains the files to generate the true covariance matrix according the various adjacency matrices incluiding those from the DREAM3 challenge as well as random DAGs.

### metric_functions.R

This file contains functions that calculate various metrics of interest, especially the macro-averaged AUC for DAGs.

### simulations.R

This file contains the code to generate the simulations reported in the paper.

### pDag_dairy_cattle_data.Rmd

This file contains the code to conduct the analysis of the dairy cattle data.

### pDag_dairy_cattle_data.pdf

This file contains the output for the analysis of the dairy cattle data.

## Example

To generate all the data for the tables, open up bash and type:
```
Rscript simulations.R
```
In particular, the lines that generates Table 1 are:

```
# Table 1
generate_tables(Methods = c('pcalg_custom','ccdr_paper_t','partial2'),
              Btypes = c('genB_mult_Yeast1','genB_mult_Yeast2','genB_mult_Yeast3', 'genB_mult_Ecoli1','genB_mult_Ecoli2'),
              Ns = c(40,50,100,200),
              Seeds = 1:20,
              m = 2,
              m1 = 25)
```

# Authors
Syed Rahman
