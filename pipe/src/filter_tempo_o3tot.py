#!/usr/bin/env python3
import argparse
import sys
import shutil
from pathlib import Path
import numpy as np
from netCDF4 import Dataset

# use:
# python filter_tempo_o3tot.py "/path/to/TEMPO_O3TOT_L2_*.nc"

def _get_attr(var, name, default=None):
    """Safely retrieve an attribute from a NetCDF variable."""
    try:
        return var.getncattr(name)
    except Exception:
        return default

def _get_fill_value(var):
    """
    Return a safe FillValue for this variable:
      1) use _FillValue or missing_value if present and valid
      2) otherwise pick a dtype-appropriate sentinel
         - float: np.nan
         - signed int: dtype min
         - unsigned int: dtype max
         - other (e.g., char/strings): return None (caller should skip)
    """
    fv = _get_attr(var, "_FillValue", None)
    if fv is None:
        fv = _get_attr(var, "missing_value", None)

    dt = var.dtype

    # If fv is NaN-like, treat as missing
    if isinstance(fv, float) and np.isnan(fv):
        fv = None

    if fv is not None:
        return fv

    # Derive sensible default
    if np.issubdtype(dt, np.floating):
        return np.nan
    if np.issubdtype(dt, np.signedinteger):
        return np.iinfo(dt).min
    if np.issubdtype(dt, np.unsignedinteger):
        return np.iinfo(dt).max

    # For non-numeric types (e.g., S1/char), skip masking
    return None

def compute_bad_pixel_mask(ds):
    """
    Build bad-pixel mask (True=bad) from product/quality_flag.
    Good pixels are those with quality_flag in {0,1,2,5}.
    """
    prod = ds.groups["product"]
    qf = prod.variables["quality_flag"][...]
    allowed = (qf == 0) | (qf == 1) | (qf == 2) | (qf == 5)
    return ~allowed

def _safe_read(var, fillv):
    """Read variable as ndarray and fill masked values with fillv."""
    data = var[...]
    if isinstance(data, np.ma.MaskedArray):
        # If fillv is None (non-numeric), just convert via .filled(0) to allow reading;
        # but caller will skip masking anyway.
        data = data.filled(0 if fillv is None else fillv)
    return data

def _broadcast_bad_mask(bad_mask, var):
    """
    Broadcast a (mirror_step,xtrack) mask to the shape of `var`.
    Requires that var.dimensions contain both 'mirror_step' and 'xtrack'.
    """
    dims = var.dimensions
    if "mirror_step" not in dims or "xtrack" not in dims:
        raise ValueError(f"{getattr(var, 'name', 'var')} dims {dims} do not contain ('mirror_step','xtrack')")
    ms_axis = dims.index("mirror_step")
    xt_axis = dims.index("xtrack")
    shape = [1] * len(dims)
    shape[ms_axis] = bad_mask.shape[0]
    shape[xt_axis] = bad_mask.shape[1]
    return np.broadcast_to(bad_mask.reshape(shape), var.shape)

def _mask_one_group(ds, group_name, bad_mask, skip_names=None):
    """
    Apply mask to ALL numeric variables in group `group_name` whose dims include
    ('mirror_step','xtrack'). Each var gets masked with its own FillValue.
    """
    if skip_names is None:
        skip_names = set()
    grp = ds.groups[group_name]
    for name, var in grp.variables.items():
        if name in skip_names:
            continue
        dims = var.dimensions
        if "mirror_step" not in dims or "xtrack" not in dims:
            continue  # only grid variables

        fillv = _get_fill_value(var)
        if fillv is None:
            # non-numeric (e.g., char) or no safe fill → skip to avoid dtype issues
            continue

        data = _safe_read(var, fillv)
        bad_for_var = _broadcast_bad_mask(bad_mask, var)

        out = np.array(data, copy=True)

        # Ensure we don't try to write NaN into integer arrays
        if np.issubdtype(out.dtype, np.integer) and (isinstance(fillv, float) and np.isnan(fillv)):
            # pick a safe integer sentinel
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
        bad = compute_bad_pixel_mask(ds)

        # product: skip quality_flag itself
        _mask_one_group(ds, "product", bad_mask=bad, skip_names={"quality_flag"})
        # support_data: mask everything with (mirror_step,xtrack)
        _mask_one_group(ds, "support_data", bad_mask=bad)

        # update history
        try:
            hist = getattr(ds, "history", "")
            msg = "filtered product/* and support_data/* with quality_flag in {0,1,2,5}"
            ds.history = (hist + " " + msg).strip() if hist else msg
        except Exception:
            pass

    kept = (~bad).sum()
    total = bad.size
    print(f"{in_path.name} -> {out_path.name} | kept {kept}/{total} pixels ({kept/total:.1%}) across product/* and support_data/*")

def expand_inputs(patterns):
    """Expand patterns into a list of matching TEMPO_O3TOT_L2_ NetCDF files."""
    files = []
    for p in patterns:
        if any(ch in p for ch in "*?[]"):
            files.extend(sorted(Path().glob(p)))
        else:
            files.append(Path(p))
    return [f for f in files if f.is_file() and "TEMPO_O3TOT_L2_" in f.name and f.suffix == ".nc"]

def main():
    ap = argparse.ArgumentParser(description="Mask product and support_data groups in TEMPO O3TOT L2 using quality_flag ∈ {0,1,2,5}.")
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
