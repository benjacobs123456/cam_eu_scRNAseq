#!/bin/bash
#!
#! Example SLURM job script for Peta4-Skylake (Skylake CPUs, OPA)
#! Last updated: Mon 13 Nov 12:25:17 GMT 2017
#!

#!#############################################################
#!#### Modify the options in this section as appropriate ######
#!#############################################################

#! sbatch directives begin here ###############################
#! Name of the job:
#SBATCH -J viral_alignment
#! Which project should be charged:
#SBATCH -A SAWCER-SL3-CPU
#! How many whole nodes should be allocated?
#SBATCH --nodes=1
#! How many (MPI) tasks will there be in total? (<= nodes*32)
#! The skylake/skylake-himem nodes have 32 CPUs (cores) each.
#SBATCH --ntasks=1
#! How much wallclock time will be required?
#SBATCH --time=8:00:00
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=NONE
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
#SBATCH --no-requeue
#SBATCH --export=ALL
#! For 6GB per CPU, set "-p skylake"; for 12GB per CPU, set "-p skylake-himem":
#SBATCH -p cclake-himem
#SBATCH --cpus-per-task=32
#SBATCH --array=1-70

#! sbatch directives end here (put any additional directives above this line)



#! Number of nodes and tasks per node allocated by SLURM (do not change):
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS
mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/')
#! ############################################################
#! Modify the settings below to specify the application's environment, location
#! and launch method:

#! Optionally modify the environment seen by the application
#! (note that SLURM reproduces the environment at submission irrespective of ~/.bashrc):
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment

# git: https://github.com/PierreBSC/Viral-Track
source /home/hpcjaco1/.bashrc
source activate /home/hpcjaco1/.conda/envs/umitools
module load R/4.0.3
module load samtools
module load star

# navigate to the wd
cd /rds/project/sjs1016/rds-sjs1016-msgen/10X_RawData/5prime_V1/SLX-21497

# get list of unique samples
sample=$(ls *.fastq.gz | awk -F "_" '{print $1}' | uniq | awk -v x=$SLURM_ARRAY_TASK_ID 'NR==x{print $1}')

cd /rds/project/sjs1016/rds-sjs1016-msgen/bj_scrna/Cambridge_EU_combined/viral_alignment/cam_v2
mkdir $sample
cd $sample

echo "Running alignment for:"
echo $sample

# combine fastQ lanes
cat /rds/project/sjs1016/rds-sjs1016-msgen/10X_RawData/5prime_V1/SLX-21497/$sample\_S1_L00?_R1_001.fastq.gz > read_1.gz
cat /rds/project/sjs1016/rds-sjs1016-msgen/10X_RawData/5prime_V1/SLX-21497/$sample\_S1_L00?_R2_001.fastq.gz > read_2.gz

# barcode whitelist
umi_tools whitelist --stdin read_1.gz \
                    --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
                    --log2stderr > whitelist.txt


# extract cells
umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
                  --stdin read_1.gz \
                  --stdout read_1_extracted.fastq.gz \
                  --read2-in read_2.gz \
                  --read2-out=read_2_extracted.fastq.gz \
                  --whitelist=whitelist.txt

# make target file
echo read_2_extracted.fastq.gz > target_file.txt

# parameter file
echo N_thread=20 > parameters.txt
echo Output_directory="temp_umi_tools" >> parameters.txt
echo Index_genome="/rds/project/sjs1016/rds-sjs1016-msgen/bj_scrna/references/index_for_star/" >> parameters.txt
echo Viral_annotation_file="/home/hpcjaco1/Viral-Track/Virusite_annotation_file.txt" >> parameters.txt
echo Name_run="EBV_test" >> parameters.txt
echo Load_STAR_module=FALSE >> parameters.txt
echo Load_samtools_module=FALSE >> parameters.txt
echo Load_stringtie_module=FALSE >> parameters.txt
echo Minimal_read_mapped=50 >> parameters.txt

# run
Rscript ~/Viral-Track/Viral_Track_scanning.R parameters.txt target_file.txt

# Rscript ~/Viral-Track/Viral_Track_cell_demultiplexing.R  parameters.txt target_file.txt




#! Work directory (i.e. where the job will run):
workdir="$SLURM_SUBMIT_DIR"  # The value of SLURM_SUBMIT_DIR sets workdir to the directory
                         # in which sbatch is run.

#! Are you using OpenMP (NB this is unrelated to OpenMPI)? If so increase this
#! safe value to no more than 32:
export OMP_NUM_THREADS=1

#! Number of MPI tasks to be started by the application per node and in total (do not change):
np=$[${numnodes}*${mpi_tasks_per_node}]

#! The following variables define a sensible pinning strategy for Intel MPI tasks -
#! this should be suitable for both pure MPI and hybrid MPI/OpenMP jobs:
export I_MPI_PIN_DOMAIN=omp:compact # Domains are $OMP_NUM_THREADS cores in size
export I_MPI_PIN_ORDER=scatter # Adjacent domains have minimal sharing of caches/sockets
#! Notes:
#! 1. These variables influence Intel MPI only.
#! 2. Domains are non-overlapping sets of cores which map 1-1 to MPI tasks.
#! 3. I_MPI_PIN_PROCESSOR_LIST is ignored if I_MPI_PIN_DOMAIN is set.
#! 4. If MPI tasks perform better when sharing caches/sockets, try I_MPI_PIN_ORDER=compact.


#! Uncomment one choice for CMD below (add mpirun/mpiexec options if necessary):

#! Choose this for a MPI code (possibly using OpenMP) using Intel MPI.
CMD="mpirun -ppn $mpi_tasks_per_node -np $np $application $options"

#! Choose this for a pure shared-memory OpenMP parallel program on a single node:
#! (OMP_NUM_THREADS threads will be created):
#CMD="$application $options"

#! Choose this for a MPI code (possibly using OpenMP) using OpenMPI:
#CMD="mpirun -npernode $mpi_tasks_per_node -np $np $application $options"


###############################################################
### You should not have to change anything below this line ####
###############################################################

cd $workdir
echo -e "Changed directory to `pwd`.\n"

JOBID=$SLURM_JOB_ID

echo -e "JobID: $JOBID\n======"
echo "Time: `date`"
echo "Running on master node: `hostname`"
echo "Current directory: `pwd`"

if [ "$SLURM_JOB_NODELIST" ]; then
        #! Create a machine file:
        export NODEFILE=`generate_pbs_nodefile`
        cat $NODEFILE | uniq > machine.file.$JOBID
        echo -e "\nNodes allocated:\n================"
        echo `cat machine.file.$JOBID | sed -e 's/\..*$//g'`
fi

echo -e "\nnumtasks=$numtasks, numnodes=$numnodes, mpi_tasks_per_node=$mpi_tasks_per_node (OMP_NUM_THREADS=$OMP_NUM_THREADS)"

echo -e "\nExecuting command:\n==================\n$CMD\n"

eval $CMD
