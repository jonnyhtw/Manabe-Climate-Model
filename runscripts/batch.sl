#!/bin/bash
#SBATCH --job-name=run_model_simple
#SBATCH --partition=standard
#SBATCH --qos=short
#SBATCH --account=ncas_cms

./runscripts/run_model_simple

