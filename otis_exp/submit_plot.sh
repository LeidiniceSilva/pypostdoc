#!/bin/bash

#SBATCH -A ICT26_ESP
#SBATCH -p dcgp_usr_prod
#SBATCH -N 1
#SBATCH --ntasks-per-node=112
#SBATCH -t 1-00:00:00
#SBATCH -J Plot
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=mda_silv@ictp.it

#__author__      = 'Leidinice Silva'
#__email__       = 'leidinicesilva@gmail.com'
#__date__        = 'Jun 21, 2026'
#__description__ = 'Call plot script'

{
set -eo pipefail

echo
echo "--------------- INIT PLOT ----------------"

python3 plot_maps_hurricane_otis_p3_prec.py

echo
echo "--------------- THE END PLOT ----------------"

}
