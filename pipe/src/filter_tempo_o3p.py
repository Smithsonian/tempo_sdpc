#!/usr/bin/env python3
import argparse
import sys
import shutil
from pathlib import Path
import numpy as np
from netCDF4 import Dataset

# use:
# python filter_tempo_o3prof.py "/data/tempo4/TEMPO/o3prof_test/L2/RAD/D16232/S001/G01/O3PROF/TEMPO_O3PROF_L2_V03_20240615T103025Z_S001G01.nc"

def _get_attr(var, name, default=None):
    try:
        return var.getncattr(name)
    except Exception:
        return default

def _safe_read(var, fillv):
    """
    Read variable data; if it's a masked array, fill masked values with `fillv`.
    """
    data = var[...]
    if isinstance(data, np.ma.MaskedArray):
        data = data.filled(fillv if fillv is not None else 0)
    return data

def _get_fill_value(var):
    """
    Determine a safe FillValue for `var`.
      1) Prefer `_FillValue` or `missing_value` if present and not NaN.
      2) Otherwise choose by dtype:
         - float  -> np.nan
         - int    -> dtype min (signed) / dtype max (unsigned)
         - other  -> None (skip masking)
    """
    fv = _get_attr(var, "_FillValue", None)
    if fv is None:
        fv = _get_attr(var, "missing_value", None)

    dt = var.dtype
    if isinstance(fv, float) and np.isnan(fv):
        fv = None

    if fv is not None:
        return fv

    if np.issubdtype(dt, np.floating):
        return np.nan
    if np.issubdtype(dt, np.signedinteger):
        return np.iinfo(dt).min
    if np.issubdtype(dt, np.unsignedinteger):
        return np.iinfo(dt).max

    # non-numeric types (e.g., char) — skip
    return None

def _ensure_2d_mask(mask, varname):
    if mask.ndim != 2:
        raise ValueError(f"Bad mask for {varname}: expected 2D (mirror_step,xtrack), got {mask.shape}")

def _expand_mask_to_var(bad_mask, var):
    """
    Broadcast a (mirror_step,xtrack) mask to `var`'s shape.
    """
    dims = var.dimensions
    if "mirror_step" not in dims or "xtrack" not in dims:
        raise ValueError(f"Target variable dims {dims} do not contain ('mirror_step','xtrack')")
    ms_axis = dims.index("mirror_step")
    xt_axis = dims.index("xtrack")
    shape = [1] * len(dims)
    shape[ms_axis] = bad_mask.shape[0]
    shape[xt_axis] = bad_mask.shape[1]
    return np.broadcast_to(bad_mask.reshape(shape), var.shape)

def compute_bad_pixel_mask(ds, avg_residuals_thresh=0.3, fit_rms_thresh=3.0,
                           exit_good_min=1, exit_good_max=99):
    """
    Build bad pixel mask (True=bad) from qa_statistics:
      - any(avg_residuals) > 0.3
      - any(fit_RMS) > 3
      - exit_status < 1 or > 99
    avg_residuals, fit_RMS: (mirror_step,xtrack,fitting_windows) or (mirror_step,xtrack).
    """
    qa = ds.groups["qa_statistics"]

    # exit_status: (mirror_step, xtrack)
    exit_status_var = qa.variables["exit_status"]
    exit_status = _safe_read(exit_status_var, _get_fill_value(exit_status_var))
    bad_exit = (exit_status < exit_good_min) | (exit_status > exit_good_max)

    # avg_residuals
    avg_res_var = qa.variables["avg_residuals"]
    avg_res = _safe_read(avg_res_var, _get_fill_value(avg_res_var))
    if avg_res.ndim == 3:
        bad_avg = np.any(avg_res > avg_residuals_thresh, axis=-1)
    elif avg_res.ndim == 2:
        bad_avg = (avg_res > avg_residuals_thresh)
    else:
        raise ValueError(f"Unexpected avg_residuals dims: {avg_res.shape}")

    # fit_RMS
    fit_rms_var = qa.variables["fit_RMS"]
    fit_rms = _safe_read(fit_rms_var, _get_fill_value(fit_rms_var))
    if fit_rms.ndim == 3:
        bad_rms = np.any(fit_rms > fit_rms_thresh, axis=-1)
    elif fit_rms.ndim == 2:
        bad_rms = (fit_rms > fit_rms_thresh)
    else:
        raise ValueError(f"Unexpected fit_RMS dims: {fit_rms.shape}")

    _ensure_2d_mask(bad_exit, "exit_status")
    _ensure_2d_mask(bad_avg, "avg_residuals")
    _ensure_2d_mask(bad_rms, "fit_RMS")

    bad = bad_exit | bad_avg | bad_rms
    return bad

def _mask_one_group(ds, group_name, bad_mask, skip_names=None):
    """
    Apply the QA mask to ALL numeric variables in `group_name`
    that include ('mirror_step','xtrack'). Each variable is filled
    with its own FillValue (safe fallback if missing).
    """
    if skip_names is None:
        skip_names = set()
    grp = ds.groups[group_name]
    for name, var in grp.variables.items():
        if name in skip_names:
            continue
        dims = var.dimensions
        if "mirror_step" not in dims or "xtrack" not in dims:
            continue

        fillv = _get_fill_value(var)
        if fillv is None:
            # non-numeric or no safe fill value → skip
            continue

        data = _safe_read(var, fillv)
        bad_for_var = _expand_mask_to_var(bad_mask, var)

        out = np.array(data, copy=True)

        # avoid writing NaN to integer arrays
        if np.issubdtype(out.dtype, np.integer) and (isinstance(fillv, float) and np.isnan(fillv)):
            fillv = np.iinfo(out.dtype).min if np.issubdtype(out.dtype, np.signedinteger) else np.iinfo(out.dtype).max

        out[bad_for_var] = fillv
        var[:] = out

def process_one_file(in_path: Path, overwrite=False):
    """Copy input file to filtered/<name>, apply mask, update history."""
    if not in_path.is_file():
        print(f"Not a file: {in_path}", file=sys.stderr)
        return

    out_path = in_path.parent.joinpath('filtered', in_path.name)
    if out_path.exists() and not overwrite:
        print(f"Exists, skip: {out_path.name}")
        return
    out_path.parent.mkdir (parents=True, exist_ok=True)

    shutil.copy2(in_path, out_path)

    with Dataset(str(out_path), mode="r+") as ds:
        # 1) build bad mask from qa_statistics
        bad = compute_bad_pixel_mask(
            ds,
            avg_residuals_thresh=0.3,
            fit_rms_thresh=3.0,
            exit_good_min=1,
            exit_good_max=99,
        )

        # 2) apply to ALL variables in product and support_data that are on the grid
        _mask_one_group(ds, "product", bad_mask=bad)
        _mask_one_group(ds, "support_data", bad_mask=bad)

        # 3) history
        try:
            hist = getattr(ds, "history", "")
            msg = ("filtered product/* and support_data/* with QA: "
                   "exit_status in [1,99], avg_residuals <= 0.3, fit_RMS <= 3 "
                   "→ masked grid variables with FillValue")
            ds.history = (hist + " " + msg).strip() if hist else msg
        except Exception:
            pass

        kept = (~bad).sum()
        total = bad.size
        print(f"{in_path.name} -> {out_path.name} | kept {kept}/{total} pixels ({kept/total:.1%}) after QA filters")

def expand_inputs(patterns):
    files = []
    for p in patterns:
        if any(ch in p for ch in "*?[]"):
            files.extend(sorted(Path().glob(p)))
        else:
            files.append(Path(p))
    # Restrict to O3PROF L2 NetCDFs
    return [f for f in files if f.is_file() and "TEMPO_O3PROF_L2_" in f.name and f.suffix == ".nc"]

def main():
    ap = argparse.ArgumentParser(
        description="Filter ALL grid variables in TEMPO O3PROF L2 (product & support_data) using QA thresholds."
    )
    ap.add_argument("inputs", nargs="+", help="Input files or glob patterns")
    ap.add_argument("--overwrite", action="store_true", help="Overwrite if filtered file exists")
    args = ap.parse_args()

    files = expand_inputs(args.inputs)
    if not files:
        print("No input files found.", file=sys.stderr)
        sys.exit(2)

    for f in files:
        process_one_file(f, overwrite=args.overwrite)

if __name__ == "__main__":
    main()
