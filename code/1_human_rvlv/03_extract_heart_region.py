#!/usr/bin/env python3
import os as _o, sys as _s
_s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), '..'))
import config as cfg
"""Step 3 (Pass 1, the single 17GB streaming read): apply the derived mapping,
compute continuous histo-voxel coords, and cache every event whose rounded
baseline voxel lands in the heart bounding box expanded by MARGIN voxels.
Cache columns: (fi, fj, fk, t) float32 -> work/heart_region_events.npy
This lets the registration sweep run cheaply on a few-million-event subset.
"""
import os, json, time, numpy as np, nibabel as nib

AUD = cfg.OUT
OUT = cfg.OUT; WRK = cfg.OUT
LM   = f"{cfg.RAW}/zenodo_11243763/20230606_positronium_patient_evaluated_data_Histo_Out.l"
HIST = f"{cfg.RAW}/zenodo_11243763/20230606_histoimage3d_lm.nii.gz"
HRT  = f"{cfg.MASKS}/totalseg_11243763/heart.nii.gz"
RV   = f"{cfg.MASKS}/totalseg_chambers_11243763/heart_ventricle_right.nii.gz"
LV   = f"{cfg.MASKS}/totalseg_chambers_11243763/heart_ventricle_left.nii.gz"
MARGIN = 15

mapping = json.load(open(os.path.join(OUT, "coordinate_mapping.json")))
perm = mapping["perm_LM_to_world_xyz"]; sgn = mapping["signs"]
sxx, syy, szz = mapping["histo_affine_diag"]
ox, oy, oz = mapping["histo_affine_origin"]
print("mapping perm", perm, "signs", sgn)

himg = nib.load(HIST); Nx, Ny, Nz = himg.shape
heart = np.asanyarray(nib.load(HRT).dataobj) > 0
ijk = np.argwhere(heart)
bmin = np.maximum(ijk.min(0) - MARGIN, 0)
bmax = np.minimum(ijk.max(0) + MARGIN, [Nx-1, Ny-1, Nz-1])
print("heart bbox (vox):", ijk.min(0).tolist(), "..", ijk.max(0).tolist())
print("expanded bbox (+%d):" % MARGIN, bmin.tolist(), "..", bmax.tolist())

sz = os.path.getsize(LM); n = sz // 32
mm = np.memmap(LM, dtype=np.float64, mode="r", shape=(n, 4))
CH = 25_000_000
parts = []
t0 = time.time(); kept = 0; seen = 0
for s in range(0, n, CH):
    e = min(s + CH, n)
    a = np.asarray(mm[s:e])               # (m,4) float64
    wx = sgn[0]*a[:, perm[0]]; wy = sgn[1]*a[:, perm[1]]; wz = sgn[2]*a[:, perm[2]]
    fi = (wx - ox)/sxx; fj = (wy - oy)/syy; fk = (wz - oz)/szz
    ri = np.round(fi); rj = np.round(fj); rk = np.round(fk)
    m = ((ri >= bmin[0]) & (ri <= bmax[0]) &
         (rj >= bmin[1]) & (rj <= bmax[1]) &
         (rk >= bmin[2]) & (rk <= bmax[2]))
    if m.any():
        out = np.empty((int(m.sum()), 4), np.float32)
        out[:, 0] = fi[m]; out[:, 1] = fj[m]; out[:, 2] = fk[m]; out[:, 3] = a[m, 3]
        parts.append(out)
    seen += (e - s); kept += int(m.sum())
    if (s // CH) % 4 == 0:
        print(f"  chunk {s//CH:2d}: seen {seen:,} kept {kept:,}  ({time.time()-t0:.0f}s)")
cache = np.concatenate(parts, axis=0)
print(f"DONE read in {time.time()-t0:.0f}s  total kept {cache.shape[0]:,} of {n:,}")
np.save(os.path.join(WRK, "heart_region_events.npy"), cache)
print("WROTE", os.path.join(WRK, "heart_region_events.npy"), cache.shape, cache.dtype)

# quick baseline membership at rounded voxel (sanity)
rv = np.asanyarray(nib.load(RV).dataobj) > 0
lv = np.asanyarray(nib.load(LV).dataobj) > 0
ri = np.round(cache[:,0]).astype(int); rj = np.round(cache[:,1]).astype(int); rk = np.round(cache[:,2]).astype(int)
ib = (ri>=0)&(ri<Nx)&(rj>=0)&(rj<Ny)&(rk>=0)&(rk<Nz)
ri,rj,rk = ri[ib],rj[ib],rk[ib]
in_rv = rv[ri,rj,rk]; in_lv = lv[ri,rj,rk]
print(f"baseline cached-region RV events: {int(in_rv.sum()):,}   LV events: {int(in_lv.sum()):,}")
# t window stats for cached events
t = cache[:,3]
phys = (t>-50)&(t<50)
print(f"cached t in [-50,50]ns: {int(phys.sum()):,} ({phys.mean()*100:.1f}%)  median {np.median(t[phys]):.3f}")
json.dump({"margin":MARGIN,"bbox_min":bmin.tolist(),"bbox_max":bmax.tolist(),
           "kept":int(cache.shape[0]),"total":int(n),
           "baseline_rv":int(in_rv.sum()),"baseline_lv":int(in_lv.sum())},
          open(os.path.join(WRK,"extract_meta.json"),"w"), indent=2)
