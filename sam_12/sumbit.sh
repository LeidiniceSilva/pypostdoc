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

VAR_LIST=("pr" "tas")

RCM_LIST=("RegCM5-ECEarth_ICTP" "RegCM5-MPI_ICTP" "RegCM5-Nor_USP")
GCM_LIST=("EC-Earth3-Veg" "MPI-ESM1-2-HR" "NorESM2-MM")

for i in "${!RCM_LIST[@]}"; do
    RCM=${RCM_LIST[$i]}
    GCM=${GCM_LIST[$i]}

    for VAR in "${VAR_LIST[@]}"; do
        echo "Running VAR=${VAR} RCM=${RCM} GCM=${GCM}"

        python3 plot_maps_clim_srf.py --var "${VAR}" --rcm_ii "${RCM}" --gcm_iii "${GCM}"
        python3 plot_maps_bias_srf.py --var "${VAR}" --rcm_ii "${RCM}" --gcm_iii "${GCM}"

    done
done

echo
echo "--------------- THE END PLOT ----------------"

}
