#!/usr/bin/env python3
import os as _o, sys as _s
_s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), '..'))
import config as cfg
"""Step 2: derive the listmode->histo-voxel mapping by searching axis
permutations & signs, validated against the provided histoimage (which was
reconstructed FROM this listmode). Then test histo<->CT index correspondence
using the heart mask activity enrichment. No repo results are read.
"""
import os, json, itertools, numpy as np, nibabel as nib

AUD = cfg.OUT
OUT = cfg.OUT
LM   = f"{cfg.RAW}/zenodo_11243763/20230606_positronium_patient_evaluated_data_Histo_Out.l"
HIST = f"{cfg.RAW}/zenodo_11243763/20230606_histoimage3d_lm.nii.gz"
RV   = f"{cfg.MASKS}/totalseg_chambers_11243763/heart_ventricle_right.nii.gz"
LV   = f"{cfg.MASKS}/totalseg_chambers_11243763/heart_ventricle_left.nii.gz"
HRT  = f"{cfg.MASKS}/totalseg_11243763/heart.nii.gz"

# ---- histoimage (the reconstruction target) ----
himg = nib.load(HIST)
H = np.asanyarray(himg.dataobj).astype(np.float32)
SH = H.shape  # (512,512,644)
Nx, Ny, Nz = SH
A = himg.affine
# diagonal affine -> world = diag*vox + origin
sxx, syy, szz = A[0,0], A[1,1], A[2,2]
ox, oy, oz = A[0,3], A[1,3], A[2,3]
Hflat = H.reshape(-1)
Hsupport = Hflat > 0
print("histo support voxels:", int(Hsupport.sum()), "/", Hflat.size,
      "  H sum:", float(Hflat.sum()))

# ---- listmode sample ----
sz = os.path.getsize(LM); n = sz // 32
mm = np.memmap(LM, dtype=np.float64, mode="r", shape=(n, 4))
nblk, blk = 40, 250_000
offs = np.linspace(0, n - blk, nblk).astype(np.int64)
samp = np.concatenate([np.array(mm[o:o+blk]) for o in offs], axis=0)
XYZ = samp[:, :3]
T = samp[:, 3]
# drop sentinel-t events for spatial test (keep all spatially valid)
print("sample events:", XYZ.shape[0])

def world_to_idx(wx, wy, wz):
    i = np.round((wx - ox) / sxx).astype(np.int64)
    j = np.round((wy - oy) / syy).astype(np.int64)
    k = np.round((wz - oz) / szz).astype(np.int64)
    return i, j, k

flatlen = Nx * Ny * Nz
results = []
for perm in itertools.permutations(range(3)):          # which LM axis -> (world x,y,z)
    for sgn in itertools.product([1, -1], repeat=3):
        wx = sgn[0] * XYZ[:, perm[0]]
        wy = sgn[1] * XYZ[:, perm[1]]
        wz = sgn[2] * XYZ[:, perm[2]]
        i, j, k = world_to_idx(wx, wy, wz)
        inb = (i >= 0) & (i < Nx) & (j >= 0) & (j < Ny) & (k >= 0) & (k < Nz)
        frac_in = inb.mean()
        if frac_in < 0.5:
            results.append((perm, sgn, frac_in, 0.0, 0.0)); continue
        fl = (i[inb] * Ny * Nz + j[inb] * Nz + k[inb])
        in_support = Hsupport[fl].mean()
        meanH = Hflat[fl].mean()
        results.append((perm, sgn, float(frac_in), float(in_support), float(meanH)))

results.sort(key=lambda r: (r[3], r[4]), reverse=True)
print("\n top mappings by (frac events in histo-support, mean H at event voxel):")
print(" perm(LM->wxyz)  signs       frac_in  in_support  meanH")
for perm, sgn, fin, sup, mh in results[:6]:
    print(f"  {perm}  {sgn}   {fin:6.3f}   {sup:7.4f}   {mh:9.2f}")
best = results[0]
perm, sgn = best[0], best[1]
print("\nCHOSEN mapping: LM axes", perm, "-> world (x,y,z), signs", sgn)

# ---- apply chosen mapping to full sample, get histo voxel idx ----
wx = sgn[0]*XYZ[:,perm[0]]; wy = sgn[1]*XYZ[:,perm[1]]; wz = sgn[2]*XYZ[:,perm[2]]
i, j, k = world_to_idx(wx, wy, wz)
inb = (i>=0)&(i<Nx)&(j>=0)&(j<Ny)&(k>=0)&(k<Nz)
print("frac in-bounds (chosen):", round(float(inb.mean()),4))

# ---- histo<->CT index correspondence test via heart mask enrichment ----
def load_mask(p):
    return np.asanyarray(nib.load(p).dataobj) > 0
heart = load_mask(HRT); rv = load_mask(RV); lv = load_mask(LV)
ii, jj, kk = i[inb], j[inb], k[inb]
in_heart = heart[ii, jj, kk]
in_rv = rv[ii, jj, kk]; in_lv = lv[ii, jj, kk]
frac_heart = in_heart.mean()
heart_vol_frac = heart.sum() / heart.size
enrich = frac_heart / heart_vol_frac
print("\n== histo<->CT index-correspondence test (heart mask) ==")
print(f"heart mask vol fraction: {heart_vol_frac:.5f}")
print(f"events-in-heart fraction: {frac_heart:.5f}  ->  ENRICHMENT {enrich:.2f}x")
print(f"events-in-RV: {in_rv.mean():.6f} ({in_rv.sum()} of {ii.size})  enrich {in_rv.mean()/(rv.sum()/rv.size):.2f}x")
print(f"events-in-LV: {in_lv.mean():.6f} ({in_lv.sum()} of {ii.size})  enrich {in_lv.mean()/(lv.sum()/lv.size):.2f}x")
# extrapolated full-dataset chamber counts (sample fraction * N)
samp_n = XYZ.shape[0]
print(f"\nExtrapolated full-dataset counts (x{n/samp_n:.1f}):")
print(f"  RV ~ {int(in_rv.sum()*n/samp_n):,}   LV ~ {int(in_lv.sum()*n/samp_n):,}")

mapping = {
    "method": "axis permutation+sign search vs histoimage reconstruction support",
    "perm_LM_to_world_xyz": list(perm), "signs": list(sgn),
    "histo_affine_diag": [float(sxx),float(syy),float(szz)],
    "histo_affine_origin": [float(ox),float(oy),float(oz)],
    "frac_in_bounds": float(inb.mean()),
    "best_in_support": best[3], "best_meanH": best[4],
    "runner_up": {"perm": list(results[1][0]), "signs": list(results[1][1]),
                  "in_support": results[1][3]},
    "histo_ct_correspondence": "voxel-index (i,j,k) identity (masks live on CT grid, same shape)",
    "heart_enrichment": float(enrich), "heart_vol_frac": float(heart_vol_frac),
    "rv_enrichment": float(in_rv.mean()/(rv.sum()/rv.size)),
    "lv_enrichment": float(in_lv.mean()/(lv.sum()/lv.size)),
    "sample_n": int(samp_n),
}
with open(os.path.join(OUT, "coordinate_mapping.json"), "w") as f:
    json.dump(mapping, f, indent=2)
print("\nWROTE", os.path.join(OUT, "coordinate_mapping.json"))
