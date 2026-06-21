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

python3 plot_maps_hurricane_otis_p1_track.py --domain large --experiment exps_v2
python3 plot_maps_hurricane_otis_p2_vts.py --domain large --experiment exps_v2
python3 plot_maps_hurricane_otis_p3_prec.py --domain large --experiment exps_v2 --date 2023101900

python3 plot_maps_hurricane_otis_p1_track.py --domain small --experiment exps_v2
python3 plot_maps_hurricane_otis_p2_vts.py --domain small --experiment exps_v2
python3 plot_maps_hurricane_otis_p3_prec.py --domain small --experiment exps_v2 --date 2023101900

python3 plot_maps_hurricane_otis_p1_track.py --domain large --experiment exps_v3
python3 plot_maps_hurricane_otis_p2_vts.py --domain large --experiment exps_v3
python3 plot_maps_hurricane_otis_p3_prec.py --domain large --experiment exps_v3 --date 2023101900

python3 plot_maps_hurricane_otis_p1_track.py --domain small --experiment exps_v3
python3 plot_maps_hurricane_otis_p2_vts.py --domain small --experiment exps_v3
python3 plot_maps_hurricane_otis_p3_prec.py --domain small --experiment exps_v3 --date 2023101900

python3 plot_maps_hurricane_otis_p1_track.py --domain large --experiment exps_v4
python3 plot_maps_hurricane_otis_p2_vts.py --domain large --experiment exps_v4
python3 plot_maps_hurricane_otis_p3_prec.py --domain large --experiment exps_v4 --date 2023101200

python3 plot_maps_hurricane_otis_p1_track.py --domain large --experiment exps_v5
python3 plot_maps_hurricane_otis_p2_vts.py --domain large --experiment exps_v5
python3 plot_maps_hurricane_otis_p3_prec.py --domain large --experiment exps_v5 --date 2023102000

echo
echo "--------------- THE END PLOT ----------------"

}
