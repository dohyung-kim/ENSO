"""Global OND ENSO precipitation teleconnection masks (Lenssen et al. 2020, IRI Data Library).

For each ENSO phase (El Nino, La Nina) builds:
  - wet mask : significant, elevated probability of an above-normal (upper tercile) OND season
  - dry mask : significant, elevated probability of a below-normal (lower tercile) OND season
  - combined : int8 raster, +1 = wetter, -1 = drier, 0 = no robust signal

Outputs (GeoTIFF, EPSG:4326, 0.5 deg) go to ./outputs, plus a NetCDF with everything and a PNG overview.
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import rioxarray  # noqa: F401  (registers .rio)
import xarray as xr
from matplotlib.colors import ListedColormap

URL = "http://iridl.ldeo.columbia.edu/home/.lenssen/.ensoTeleconnections/dods"
OND = 10.5  # T = centre month of the 3-month season, months since 1960-01 (0.5 = DJF)
PHASES = {"ElNino": "El Niño", "LaNina": "La Niña"}
CONDITIONS = {"wet": "Above_Normal", "dry": "Below_Normal"}


def load_ond(significance):
    ds = xr.open_dataset(URL, decode_times=False)
    # IRI stores the string coordinates as bytes (b'ElNino'); decode so .sel works with str labels
    ds = ds.assign_coords(
        ENSO=ds.ENSO.values.astype(str), anomaly=ds.anomaly.values.astype(str)
    )
    if significance not in ds.significance.values:
        raise SystemExit(f"significance must be one of {ds.significance.values.tolist()}")
    # Cells not significant at this level are already NaN
    return ds["prob"].sel(T=OND, significance=significance).load()


def to_geotiff(da, path, nodata=None):
    da = da.sortby("Y", ascending=False).rio.set_spatial_dims(x_dim="X", y_dim="Y")
    da = da.rio.write_crs("EPSG:4326")
    if nodata is not None:
        da = da.rio.write_nodata(nodata)
    da.rio.to_raster(path, compress="deflate")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--significance", type=int, default=90, help="85, 90, 95 or 99 (default 90)")
    parser.add_argument("--min-prob", type=float, default=0.5,
                        help="minimum tercile probability; climatology is 1/3 (default 0.5)")
    parser.add_argument("--out", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    print(f"Downloading OND probabilities (significance {args.significance}%)...")
    prob = load_ond(args.significance)
    tag = f"OND_sig{args.significance}_p{int(round(args.min_prob * 100))}"

    results = {}
    for phase, phase_label in PHASES.items():
        masks = {}
        for cond, anomaly in CONDITIONS.items():
            p = prob.sel(ENSO=phase, anomaly=anomaly, drop=True)
            masks[cond] = p.notnull() & (p >= args.min_prob)
            results[f"prob_{cond}_{phase}"] = p
            results[f"mask_{cond}_{phase}"] = masks[cond].astype("uint8")
            to_geotiff(p, args.out / f"{phase}_{cond}_prob_{tag}.tif", nodata=np.nan)
            to_geotiff(masks[cond].astype("uint8"), args.out / f"{phase}_{cond}_mask_{tag}.tif")

        overlap = masks["wet"] & masks["dry"]
        if overlap.any():
            print(f"  warning: {int(overlap.sum())} {phase} cells flagged both wet and dry; set to 0")
        combined = (masks["wet"].astype("int8") - masks["dry"].astype("int8")).where(~overlap, 0)
        results[f"combined_{phase}"] = combined.astype("int8")
        to_geotiff(combined.astype("int8"), args.out / f"{phase}_combined_{tag}.tif")

        land = prob.sel(ENSO=phase).notnull().any("anomaly").sum().item()  # cells with any signal
        print(f"  {phase_label}: wet {int(masks['wet'].sum())} cells, dry {int(masks['dry'].sum())} cells "
              f"(of {land} significant cells)")

    ds_out = xr.Dataset(results).assign_attrs(
        source="Lenssen et al. (2020) ENSO teleconnections, IRI Data Library",
        season="OND", significance_pct=args.significance, min_prob=args.min_prob,
        combined_codes="1 = wetter (above-normal tercile), -1 = drier (below-normal tercile), 0 = no signal",
    )
    ds_out.to_netcdf(args.out / f"enso_masks_{tag}.nc")

    # Overview figure
    cmap = ListedColormap(["#b35806", "#f0f0f0", "#2166ac"])
    fig, axes = plt.subplots(2, 1, figsize=(12, 11), constrained_layout=True)
    for ax, (phase, phase_label) in zip(axes, PHASES.items()):
        results[f"combined_{phase}"].plot(ax=ax, cmap=cmap, vmin=-1.5, vmax=1.5, add_colorbar=False)
        ax.set_title(f"{phase_label} OND precipitation teleconnection "
                     f"(sig ≥ {args.significance}%, P ≥ {args.min_prob:.2f})")
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_aspect("equal")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in ("#2166ac", "#b35806")]
    fig.legend(handles, ["Wetter (above-normal tercile)", "Drier (below-normal tercile)"],
               loc="outside lower center", ncol=2, frameon=False)
    fig.savefig(args.out / f"enso_masks_{tag}.png", dpi=150)
    print(f"Wrote outputs to {args.out.resolve()}")


if __name__ == "__main__":
    main()
