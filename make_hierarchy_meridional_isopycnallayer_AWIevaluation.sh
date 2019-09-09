#!/bin/bash
#SBATCH --job-name=hierarchy    # Specify job name
#SBATCH --partition=prepost     # Specify partition name
#SBATCH --ntasks=24           # Specify max. number of tasks to be invoked
#SBATCH --cpus-per-task=2     # Specify number of CPUs per task
#SBATCH --time=10:00:00        # Set a limit on the total run time
#SBATCH --mail-type=FAIL       # Notify user by email in case of job failure
#SBATCH --account=ab0995       # Charge resources on this project account
#SBATCH --output /mnt/lustre01/pf/a/a270046/hierarchy/myvisualization/out2.txt
#SBATCH --error /mnt/lustre01/pf/a/a270046/hierarchy/myvisualization/error2.txt

#module load python/2.7-ve0

#python -V
#which python
nproc

window=10

### --------------------------------- CORE --------------------------------- ##
### ------------------------------------------------------------------------ ## 
/work/ab0995/a270046/miniconda2-install/bin/python make_hierarchy_meridional_isopycnallayer_ALL_AWIevaluation.py LR 2008 2107 $window

### --------------------------------- BOLD --------------------------------- ##
### ------------------------------------------------------------------------ ##
/work/ab0995/a270046/miniconda2-install/bin/python make_hierarchy_meridional_isopycnallayer_ALL_AWIevaluation.py HR 2008 2107 $window

