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
#__date__        = 'JUl 28, 2026'
#__description__ = 'Submit plot'

{
set -eo pipefail

echo
echo "--------------- INIT PLOT ----------------"

python3 plot_graph_taylor_diagram.py

echo
echo "--------------- THE END PLOT ----------------"

}
