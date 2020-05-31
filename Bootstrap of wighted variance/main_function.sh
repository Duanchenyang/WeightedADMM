#!/bin/sh

# see where the job is being run
hostname

# print date and time
date

module load R

Rscript main_function.R ${SGE_TASK_ID}
