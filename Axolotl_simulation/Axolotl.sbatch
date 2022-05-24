#!/bin/bash 
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name  
#SBATCH --job-name=ElectronHyrdoSims
################# 
#a file for job output, you can check job progress, append the job ID with %j to make it unique
#SBATCH --output=teslaValve.%j.out
#################
# a file for errors from the job
#SBATCH --error=teslaValve.%j.err
#################
#time you think you need; default is 2 hours
#format could be dd-hh:mm:ss, hh:mm:ss, mm:ss, or mm
#SBATCH --time=10-00:00:00
#################
#Quality of Service (QOS); think of it as job priority. --qos=MaxJobs<N> limits the number of jobs allowed for one user to N:[1,10]
#SBATCH --qos=normal

# To use nodes not owned by our group, specify the account and partition as --account=barnardlab-ckpt and --partition=ckpt

#SBATCH --account=barnardlab-ckpt
#SBATCH --partition=ckpt

#################
#number of nodes you are requesting, the more you ask for the longer you wait
#SBATCH --nodes=1
#SBATCH --ntasks=1                   # Run a single task	
#SBATCH --cpus-per-task=40         # Number of CPU cores per tasks

# module load system
# module load parallel

#################
# --mem is memory per node; default is 4000 MB per CPU, remember to ask for enough mem to match your CPU request, since 
# sherlock automatically allocates 8 Gigs of RAM/CPU, if you ask for 8 CPUs you will need 32 Gigs of RAM, so either 
# leave --mem commented out or request >= to the RAM needed for your CPU request.
 
#SBATCH --mem=30000
#################

# Have SLURM send you an email when the job ends or fails, careful, the email could end up in your clutter folder
# Also, if you submit hundreds of jobs at once you will get hundreds of emails.
 
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
# Remember to change this to your email 
#SBATCH --mail-user=hanhiller2@gmail.com


#now run normal batch commands
python ../runTeslaValveScript.py ./teslaValve_T2_DE_generateProbs.txt
