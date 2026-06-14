#!/usr/bin/env python3
import os as _o, sys as _s
_s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), '..'))
import config as cfg
"""Independent input verification for the blind RV/LV audit.
Writes outputs/input_manifest.json (mine) and prints a summary, including the
listmode t-distribution so I can design the lifetime model.
NO repo docs/results/derived are read.
"""
import os, json, hashlib, numpy as np, nibabel as nib

AUD = cfg.OUT
OUT = cfg.OUT
os.makedirs(OUT, exist_ok=True)
os.makedirs(os.path.join(OUT, "plots"), exist_ok=True)

LM   = f"{cfg.RAW}/zenodo_11243763/20230606_positronium_patient_evaluated_data_Histo_Out.l"
CT   = f"{cfg.RAW}/zenodo_11243763/20230606_lm_2_AC_CT_WB_1.65mm.nii.gz"
HIST = f"{cfg.RAW}/zenodo_11243763/20230606_histoimage3d_lm.nii.gz"
RV   = f"{cfg.MASKS}/totalseg_chambers_11243763/heart_ventricle_right.nii.gz"
LV   = f"{cfg.MASKS}/totalseg_chambers_11243763/heart_ventricle_left.nii.gz"
HRT  = f"{cfg.MASKS}/totalseg_11243763/heart.nii.gz"

man = {"_about": "independent input manifest, blind RV/LV audit"}

# ---------- listmode ----------
sz = os.path.getsize(LM)
n  = sz // 32
mm = np.memmap(LM, dtype=np.float64, mode="r", shape=(n, 4))
print("== LISTMODE ==")
print("size_bytes", sz, "size_mod_32", sz % 32, "event_count", n)
print("first 5 rows (x,y,z,t):")
print(np.array(mm[:5]))

# sample ~5M events as 20 contiguous blocks evenly spaced across the file
nblk, blk = 20, 250_000
offs = np.linspace(0, n - blk, nblk).astype(np.int64)
samp = np.concatenate([np.array(mm[o:o+blk]) for o in offs], axis=0)
cols = ["x", "y", "z", "t"]
stats = {}
for i, c in enumerate(cols):
    v = samp[:, i]
    stats[c] = {
        "min": float(v.min()), "max": float(v.max()),
        "mean": float(v.mean()), "std": float(v.std()),
        "p0.01": float(np.percentile(v, 0.01)), "p1": float(np.percentile(v, 1)),
        "p50": float(np.percentile(v, 50)), "p99": float(np.percentile(v, 99)),
        "p99.99": float(np.percentile(v, 99.99)),
    }
print("sample stats (n=%d):" % samp.shape[0])
for c in cols:
    s = stats[c]
    print(f"  {c}: min {s['min']:.3f} max {s['max']:.3f} mean {s['mean']:.3f} "
          f"std {s['std']:.3f} p1 {s['p1']:.3f} p50 {s['p50']:.3f} p99 {s['p99']:.3f}")

# t-histogram across a wide range to see the lifetime axis
tl, th = np.percentile(samp[:, 3], [0.005, 99.995])
edges = np.linspace(tl, th, 401)
hist, _ = np.histogram(samp[:, 3], bins=edges)
np.save(os.path.join(OUT, "tsample_hist.npy"), np.vstack([edges[:-1], hist]))
# print a coarse text view: 40 bins
cb = np.linspace(tl, th, 41)
ch, _ = np.histogram(samp[:, 3], bins=cb)
print("t-distribution (coarse 40 bins over [%.3f, %.3f]):" % (tl, th))
mxc = ch.max()
for i in range(40):
    bar = "#" * int(60 * ch[i] / mxc)
    print(f"  [{cb[i]:7.3f},{cb[i+1]:7.3f}) {ch[i]:9d} {bar}")

man["listmode"] = {
    "path": LM, "size_bytes": sz, "bytes_per_event": 32,
    "event_count": int(n), "size_mod_32": int(sz % 32),
    "sample_n": int(samp.shape[0]), "sample_blocks": nblk,
    "column_stats": stats,
    "first5": np.array(mm[:5]).tolist(),
}

# ---------- niftis ----------
def affine_summary(aff):
    aff = np.asarray(aff)
    vox = np.sqrt((aff[:3, :3] ** 2).sum(axis=0))
    return aff.tolist(), vox.tolist()

def mask_info(path, ct_affine=None):
    img = nib.load(path)
    data = np.asanyarray(img.dataobj)
    aff, vox = affine_summary(img.affine)
    nz = int((data > 0).sum())
    info = {"path": path, "shape": list(img.shape), "dtype": str(data.dtype),
            "affine": aff, "voxel_sizes": vox, "mask_voxels": nz}
    if nz > 0:
        ijk = np.argwhere(data > 0)
        info["voxel_bbox_min"] = ijk.min(0).tolist()
        info["voxel_bbox_max"] = ijk.max(0).tolist()
        cen_vox = ijk.mean(0)
        # world centroid via this image's own affine
        wc = img.affine @ np.append(cen_vox, 1.0)
        info["world_centroid_self"] = wc[:3].tolist()
        if ct_affine is not None:
            wc2 = np.asarray(ct_affine) @ np.append(cen_vox, 1.0)
            info["world_centroid_ct"] = wc2[:3].tolist()
    return info, img

print("\n== NIFTIs ==")
ct_img = nib.load(CT); hist_img = nib.load(HIST)
ct_aff, ct_vox = affine_summary(ct_img.affine)
hist_aff, hist_vox = affine_summary(hist_img.affine)
man["ct"]   = {"path": CT,   "shape": list(ct_img.shape),   "affine": ct_aff,   "voxel_sizes": ct_vox,   "dtype": str(ct_img.get_data_dtype())}
man["histoimage"] = {"path": HIST, "shape": list(hist_img.shape), "affine": hist_aff, "voxel_sizes": hist_vox, "dtype": str(hist_img.get_data_dtype())}
for nm, p in [("rv_mask", RV), ("lv_mask", LV), ("heart_mask", HRT)]:
    info, _ = mask_info(p, ct_affine=ct_img.affine)
    man[nm] = info
    print(f"{nm}: voxels {info['mask_voxels']} vox_bbox {info.get('voxel_bbox_min')}..{info.get('voxel_bbox_max')} "
          f"world_centroid_self {[round(x,1) for x in info['world_centroid_self']]}")

print("CT   affine:\n", np.array(ct_aff))
print("HIST affine:\n", np.array(hist_aff))
print("CT   voxel_sizes", [round(x,4) for x in ct_vox])
print("HIST voxel_sizes", [round(x,4) for x in hist_vox])

# world bbox of CT and histo (8 corners)
def world_bbox(img):
    sh = np.array(img.shape) - 1
    corners = np.array([[i*sh[0], j*sh[1], k*sh[2], 1]
                        for i in (0,1) for j in (0,1) for k in (0,1)], float)
    w = (img.affine @ corners.T).T[:, :3]
    return w.min(0).tolist(), w.max(0).tolist()
ct_bb = world_bbox(ct_img); hist_bb = world_bbox(hist_img)
man["ct"]["world_bbox"] = ct_bb
man["histoimage"]["world_bbox"] = hist_bb
print("CT   world bbox:", [round(x,1) for x in ct_bb[0]], "..", [round(x,1) for x in ct_bb[1]])
print("HIST world bbox:", [round(x,1) for x in hist_bb[0]], "..", [round(x,1) for x in hist_bb[1]])
print("LM xyz sample range: x[%.1f,%.1f] y[%.1f,%.1f] z[%.1f,%.1f]" % (
    stats['x']['min'], stats['x']['max'], stats['y']['min'], stats['y']['max'],
    stats['z']['min'], stats['z']['max']))

with open(os.path.join(OUT, "input_manifest.json"), "w") as f:
    json.dump(man, f, indent=2)
print("\nWROTE", os.path.join(OUT, "input_manifest.json"))
