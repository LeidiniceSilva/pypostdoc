# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jul 28, 2026"
__description__ = "This script compute stats for Taylor diagram"

import pathlib
import glob
import numpy as np
import netCDF4


def import_regcm5_srf(path, param, domain, model, exp, season, dt, res):

    arq = (f"{path}/{param}_{domain}_{model}_{exp}_{season}_{dt}_{res}_box.nc")

    files = glob.glob(arq)

    if len(files) == 0:
        raise FileNotFoundError(f"No file found: {arq}")

    data = netCDF4.Dataset(files[0])

    var = data.variables[param][:]

    # Convert masked values to NaN
    if np.ma.isMaskedArray(var):
        var = var.filled(np.nan)

    lat = data.variables["lat"][:]
    lon = data.variables["lon"][:]

    data.close()

    return var, lat, lon



def import_cpc_srf(path, param, domain, obs_name, season, dt, res):

    arq = (
        f"{path}/{param}_{domain}_{obs_name}_{season}_{dt}_{res}_box.nc")

    files = glob.glob(arq)

    if len(files) == 0:
        raise FileNotFoundError(f"No file found: {arq}")

    data = netCDF4.Dataset(files[0])

    var = data.variables[param][:]

    # Convert masked values to NaN
    if np.ma.isMaskedArray(var):
        var = var.filled(np.nan)

    lat = data.variables["lat"][:]
    lon = data.variables["lon"][:]

    data.close()

    return var, lat, lon


def compute_taylor_stats(model, obs, lats):
    """
    Computes Taylor statistics:
    correlation coefficient,
    standard deviation ratio,
    centered RMSE.
    """

    # Remove extra dimensions
    while model.ndim > 2:
        model = model[0]

    while obs.ndim > 2:
        obs = obs[0]

    # Latitude weights
    if lats.ndim == 1:
        weights = np.cos(np.radians(lats))
        weights_2d = np.broadcast_to(weights[:, None], model.shape)
    else:
        weights_2d = np.cos(np.radians(lats))

    # Valid points
    valid_mask = (
        ~np.isnan(model)
        & ~np.isnan(obs)
        & ~np.isnan(weights_2d)
    )

    m = model[valid_mask]
    o = obs[valid_mask]
    w = weights_2d[valid_mask]

    # Weighted mean
    w_sum = np.sum(w)

    m_mean = np.sum(w * m) / w_sum
    o_mean = np.sum(w * o) / w_sum

    # Anomalies
    m_prime = m - m_mean
    o_prime = o - o_mean

    # Standard deviations
    sigma_m = np.sqrt(np.sum(w * m_prime**2) / w_sum)

    sigma_o = np.sqrt(np.sum(w * o_prime**2) / w_sum)

    # Correlation coefficient
    cc = (np.sum(w * m_prime * o_prime) / (w_sum * sigma_m * sigma_o))

    # Standard deviation ratio
    if sigma_o != 0:
        ratio = sigma_m / sigma_o
    else:
        ratio = np.nan

    # Centered RMSE
    rmse = np.sqrt(np.sum(w * (m_prime - o_prime)**2) / w_sum)

    return cc, ratio, rmse


def main():

    path = ("/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5/postproc/urban/paper")
    output_dir = pathlib.Path(path) / "txt_files"
    output_dir.mkdir(parents=True, exist_ok=True)

    seasons = ["ANN", "DJF", "MAM", "JJA", "SON"]

    variables = ["pr", "tasmax", "tasmin"]

    dt = "2000-2009"
    domain = "CSAM-3"
    res = "0.25"

    model = "RegCM5-ERA5"
    exp = "URBAN"
    obs_name = "CPC"

    for param in variables:

        print(f"\nProcessing variable: {param}")

        if param == "pr":
            obs_param = "precip"

        elif param == "tasmax":
            obs_param = "tmax"

        elif param == "tasmin":
            obs_param = "tmin"

        cc_vec = np.zeros(len(seasons))
        ratio_vec = np.zeros(len(seasons))
        rmse_vec = np.zeros(len(seasons))

        for idx, season in enumerate(seasons):
            print(f"Processing: {param} - {season}")

            # Model & obs
            v_mod, lat_mod, lon_mod = import_regcm5_srf(path, param, domain, model, exp, season, dt, res)
            v_obs, lat_obs, lon_obs = import_cpc_srf(path, obs_param, domain, obs_name, season, dt, res)

            # Check grid compatibility
            if v_mod.shape != v_obs.shape:
                raise ValueError(f"Different shapes: model {v_mod.shape}, obs {v_obs.shape}")

            cc, ratio, rmse = compute_taylor_stats(v_mod, v_obs, lat_obs)

            cc_vec[idx] = cc
            ratio_vec[idx] = ratio
            rmse_vec[idx] = rmse

        # Save output
        prefix = (f"{param}_{model}_{exp}_{obs_name}_{domain}_{dt}")

        np.savetxt(output_dir / f"{prefix}_cc.txt", cc_vec.reshape(1,-1), fmt="%12.4f")
        np.savetxt(output_dir / f"{prefix}_ratio.txt", ratio_vec.reshape(1,-1), fmt="%12.4f")
        np.savetxt(output_dir / f"{prefix}_rmse.txt", rmse_vec.reshape(1,-1), fmt="%12.4f")

        print("\nFinished:", param)
        print("CC:", cc_vec)
        print("Ratio:", ratio_vec)
        print("RMSE:", rmse_vec)

if __name__ == "__main__":
    main()

