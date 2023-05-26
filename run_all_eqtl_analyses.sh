#!/bin/bash

sbatch eqtl_icelake_ms_csf.sh
sbatch eqtl_icelake_ms_pbmc.sh
sbatch eqtl_icelake_oind_pbmc.sh
sbatch eqtl_icelake_oind_csf.sh
sbatch eqtl_icelake_nind_csf.sh
sbatch eqtl_icelake_nind_pbmc.sh

