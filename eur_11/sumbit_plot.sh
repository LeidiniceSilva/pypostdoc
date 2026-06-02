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
#__date__        = 'Jan 27, 2026'
#__description__ = 'Call python code to plot'

{
set -eo pipefail

echo
echo "--------------- INIT PLOT ----------------"

years=$1

# SRF

for var in pr tas clt; do
    python3 plot_maps_clim_srf.py --var "$var" --dt "$years"
    python3 plot_maps_bias_srf.py --var "$var" --dt "$years"
done

python3 plot_maps_bias_srf_int_p99.py  --dt "$years"
python3 plot_maps_bias_srf_int_freq.py --stats int  --dt "$years"
python3 plot_maps_bias_srf_int_freq.py --stats freq  --dt "$years"

python3 plot_graph_annual_cycle.py --dt "$years"
python3 plot_graph_diurnal_cycle.py --dt "$years"
python3 plot_graph_pdf_daily.py --dt "$years"
python3 plot_graph_pdf_hourly.py --dt "$years"

# ATM

for var in uv q; do
    for level in 850hPa 500hPa 200hPa; do
        python3 plot_maps_clim_atm.py --var "$var" --level "$level" --dt "$years"
        python3 plot_maps_bias_atm.py --var "$var" --level "$level" --dt "$years"
    done
done

for var in cl clw cli; do
    python3 plot_graph_vertical_profile.py --var "$var" --dt "$years"
done

python3 plot_graph_vertical_profile_wdm7.py --dt "$years"
python3 plot_graph_vertical_profile_wdm7-wsm7.py --dt "$years"

echo
echo "--------------- THE END PLOT ----------------"

}
