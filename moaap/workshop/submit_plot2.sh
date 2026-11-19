#!/bin/bash

#SBATCH --account CMPNS_ictpclim
#SBATCH -p dcgp_usr_prod
#SBATCH -N 1
#SBATCH --ntasks-per-node=112
#SBATCH -t 12:00:00
#SBATCH -J Plot
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=mda_silv@ictp.it
#SBATCH --qos=qos_lowprio

#__author__      = 'Leidinice Silva'
#__email__       = 'leidinicesilva@gmail.com'
#__date__        = 'Mar 02, 2026'
#__description__ = 'Submit plot'

{
set -eo pipefail

echo
echo "--------------- INIT PLOT ----------------"

python3 plot_maps_moaap_mcs_peak.py

echo
echo "--------------- THE END PLOT ----------------"

}
