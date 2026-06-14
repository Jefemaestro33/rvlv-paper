"""Central path resolution for the reproducible pipeline.

Override the data/mask/output roots via environment variables. By default the
repo looks for local `raw/` and `masks/` directories next to this checkout:

    PLI_RAW    dir containing the Zenodo records  (default: <repo>/raw)
               -> expects  <PLI_RAW>/zenodo_11243763/...  and  <PLI_RAW>/zenodo_12636019/...
    PLI_MASKS  dir with the TotalSegmentator outputs (default: <repo>/masks)
               -> expects  <PLI_MASKS>/totalseg_11243763/  and  <PLI_MASKS>/totalseg_chambers_11243763/
    PLI_OUT    output dir (default: <repo>/results)

All scripts import this module and build their paths from RAW / MASKS / OUT.
Figures are written to OUT/plots and copied to paper/figures by run_all.sh.
"""
import os

_repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # code/ -> repo root
RAW   = os.environ.get("PLI_RAW", os.path.join(_repo, "raw"))
MASKS = os.environ.get("PLI_MASKS", os.path.join(_repo, "masks"))
OUT   = os.environ.get("PLI_OUT", os.path.join(_repo, "results"))

os.makedirs(OUT, exist_ok=True)
os.makedirs(os.path.join(OUT, "plots"), exist_ok=True)
