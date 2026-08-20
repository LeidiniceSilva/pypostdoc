# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jul 28, 2026"
__description__ = "This script computes statistics for Taylor diagrams"

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


def import_obs_srf(path, param, domain, obs_name, season, dt, res):

    arq = (f"{path}/{param}_{domain}_{obs_name}_{season}_{dt}_{res}_box.nc")

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

    - correlation coefficient
    - standard deviation ratio
    - centered RMSE
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

    # Check if there are valid points
    if len(m) == 0:
        return np.nan, np.nan, np.nan

    # Sum of weights
    w_sum = np.sum(w)

    # Weighted means
    m_mean = np.sum(w * m) / w_sum
    o_mean = np.sum(w * o) / w_sum

    # Anomalies
    m_prime = m - m_mean
    o_prime = o - o_mean

    # Standard deviations
    sigma_m = np.sqrt(
        np.sum(w * m_prime**2) / w_sum
    )

    sigma_o = np.sqrt(
        np.sum(w * o_prime**2) / w_sum
    )

    # Correlation coefficient
    if sigma_m == 0 or sigma_o == 0:
        cc = np.nan
    else:
        cc = (
            np.sum(w * m_prime * o_prime)
            / (w_sum * sigma_m * sigma_o)
        )

    # Standard deviation ratio
    if sigma_o != 0:
        ratio = sigma_m / sigma_o
    else:
        ratio = np.nan

    # Centered RMSE
    rmse = np.sqrt(
        np.sum(w * (m_prime - o_prime)**2) / w_sum
    )

    return cc, ratio, rmse


def main():

    path = ("/leonardo/home/userexternal/mdasilva/leonardo_work/CORDEX5/postproc/urban/paper")
    output_dir = pathlib.Path(path) / "txt_files"
    output_dir.mkdir(parents=True, exist_ok=True)

    seasons = ["DJF", "MAM", "JJA", "SON", "ANN"]

    # All variables
    variables = [
        "pr",
        "tasmax",
        "tasmin",
        "huss",
        "hfss",
        "hfls",
        "sfcWind"
    ]

    dt = "2000-2009"
    domain = "CSAM-3"
    res = "0.25"

    model = "RegCM5-ERA5"
    exp = "CTRL"

    for param in variables:

        print("\n" + "=" * 70)
        print(f"Processing variable: {param}")
        print("=" * 70)

        # Define obs dataset and vars
        if param == "pr":
            obs_param = "precip"
            obs_name = "CPC"
        elif param == "tasmax":
            obs_param = "tmax"
            obs_name = "CPC"
        elif param == "tasmin":
            obs_param = "tmin"
            obs_name = "CPC"
        elif param == "huss":
            obs_param = "huss"
            obs_name = "ERA5"
        elif param == "hfss":
            obs_param = "msshf"
            obs_name = "ERA5"
        elif param == "hfls":
            obs_param = "mslhf"
            obs_name = "ERA5"
        elif param == "sfcWind":
            obs_param = "si10"
            obs_name = "ERA5"
        else:
            raise ValueError(f"Unknown variable: {param}")

        print(f"Model variable : {param}")
        print(f"Obs variable   : {obs_param}")
        print(f"Obs dataset    : {obs_name}")

        # Arrays for the five seasons
        cc_vec = np.full(len(seasons), np.nan)
        ratio_vec = np.full(len(seasons), np.nan)
        rmse_vec = np.full(len(seasons), np.nan)

        # Process each season
        for idx, season in enumerate(seasons):

            print(f"\nProcessing: {param} - {season}")

            # Model
            v_mod, lat_mod, lon_mod = import_regcm5_srf(
                path,
                param,
                domain,
                model,
                exp,
                season,
                dt,
                res
            )

            # Obsservation
            v_obs, lat_obs, lon_obs = import_obs_srf(
                path,
                obs_param,
                domain,
                obs_name,
                season,
                dt,
                res
            )

            # Check grid compatibility
            if v_mod.shape != v_obs.shape:

                raise ValueError(
                    f"Different shapes for {param} - {season}: "
                    f"model {v_mod.shape}, obs {v_obs.shape}"
                )

            # Check latitude compatibility
            if lat_mod.shape != lat_obs.shape:

                raise ValueError(
                    f"Different latitude shapes for "
                    f"{param} - {season}: "
                    f"model {lat_mod.shape}, obs {lat_obs.shape}"
                )

            # Compute Taylor statistics
            cc, ratio, rmse = compute_taylor_stats(
                v_mod,
                v_obs,
                lat_obs
            )

            cc_vec[idx] = cc
            ratio_vec[idx] = ratio
            rmse_vec[idx] = rmse

            print(f"CC   = {cc:.4f}")
            print(f"Ratio = {ratio:.4f}")
            print(f"RMSE = {rmse:.4f}")

        # Output prefix
        prefix = (f"{param}_{model}_{exp}_{obs_name}_{domain}_{dt}")

        # Save Taylor statistics
        np.savetxt(output_dir / f"{prefix}_cc.txt", cc_vec.reshape(1, -1), fmt="%12.4f")
        np.savetxt(output_dir / f"{prefix}_ratio.txt", ratio_vec.reshape(1, -1), fmt="%12.4f")
        np.savetxt(output_dir / f"{prefix}_rmse.txt", rmse_vec.reshape(1, -1), fmt="%12.4f")

        # Print results
        print("\nFinished:", param)
        print("Observation:", obs_name)

        print("CC:")
        print(cc_vec)

        print("Ratio:")
        print(ratio_vec)

        print("RMSE:")
        print(rmse_vec)


if __name__ == "__main__":
    main()
