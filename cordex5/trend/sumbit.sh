#!/bin/bash

#SBATCH -A CMPNS_ictpclim
#SBATCH -p dcgp_usr_prod
#SBATCH -N 1
#SBATCH --ntasks-per-node=112
#SBATCH -t 1-00:00:00
#SBATCH -J Plot
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=mda_silv@ictp.it

#__author__      = 'Leidinice Silva'
#__email__       = 'leidinicesilva@gmail.com'
#__date__        = 'Jan 27, 2026'
#__description__ = 'Call python code to plot'

{
set -eo pipefail

echo
echo "--------------- INIT PLOT ----------------"

#python3 plot_maps_trend_v3.py --var pr --domain AUS-12 --idt 1970 --fdt 2023

python3 plot_maps_trend_v3.py --var pr --domain CSAM-3 --idt 2000 --fdt 2009
python3 plot_maps_trend_v3.py --var tas --domain CSAM-3 --idt 2000 --fdt 2009
python3 plot_maps_trend_v3.py --var tasmax --domain CSAM-3 --idt 2000 --fdt 2009
python3 plot_maps_trend_v3.py --var tasmin --domain CSAM-3 --idt 2000 --fdt 2009

echo
echo "--------------- THE END PLOT ----------------"

}
