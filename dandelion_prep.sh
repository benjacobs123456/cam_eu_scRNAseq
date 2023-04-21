#!/bin/bash
module load R/4.1.0-icelake
Rscript dandelion_preparation.R TCR &
Rscript dandelion_preparation.R BCR &
Rscript dandelion_preparation2.R TCR &
Rscript dandelion_preparation2.R BCR &

