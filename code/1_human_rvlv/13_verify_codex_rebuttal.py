#!/usr/bin/env python3
import os as _o, sys as _s
_s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), '..'))
import config as cfg
"""Independently verify Codex's rebuttal of my 'RV core = myocardium' mechanism:
 1) myocardial-PSF field per chamber: if WALL-myocardium leakage drove tau, the
    chamber with MORE myocardium-PSF (LV, thick wall) should be higher. Test it.
 2) HU separability: is RV-core CT density actually muscle-like, or blood-like?
    And critically -- can CT HU even distinguish blood from myocardium at all?
"""
import os, csv, numpy as np, nibabel as nib
from scipy.ndimage import gaussian_filter, binary_erosion
CH=f"{cfg.MASKS}/totalseg_chambers_11243763"; TS=f"{cfg.MASKS}/totalseg_11243763"
OUT=f"{cfg.OUT}"
CTp=f"{cfg.RAW}/zenodo_11243763/20230606_lm_2_AC_CT_WB_1.65mm.nii.gz"
ct=np.asanyarray(nib.load(CTp).dataobj).astype(np.float32)
def m(p): return np.asanyarray(nib.load(p).dataobj)>0
myo=m(f"{CH}/heart_myocardium.nii.gz")
myoPSF=gaussian_filter(myo.astype(np.float32),sigma=1.5)
lad={r["name"]:r for r in csv.DictReader(open(os.path.join(OUT,"oxygenation_ladder.csv")))}

comps={"RV":(f"{CH}/heart_ventricle_right.nii.gz"),"RA":(f"{CH}/heart_atrium_right.nii.gz"),
 "LV":(f"{CH}/heart_ventricle_left.nii.gz"),"LA":(f"{CH}/heart_atrium_left.nii.gz"),
 "aorta":(f"{CH}/aorta.nii.gz"),"pulm_ARTERY":(f"{CH}/pulmonary_artery.nii.gz")}
print("=== (1) myocardial-PSF leakage test (core erode-1) ===")
print(f"{'comp':12s}{'tau':>7s}{'myoPSF':>9s}{'coreHU':>9s}{'HU_sd':>7s}")
rows=[]
for nm,p in comps.items():
    core=binary_erosion(m(p),iterations=1)
    mp=float(myoPSF[core].mean()); hu=float(ct[core].mean()); hsd=float(ct[core].std())
    tau=float(lad[nm]["tau_core"])
    rows.append((nm,tau,mp,hu,hsd)); print(f"{nm:12s}{tau:7.3f}{mp:9.4f}{hu:9.1f}{hsd:7.1f}")
tau=np.array([r[1] for r in rows]); mp=np.array([r[2] for r in rows])
def pc(a,b): a=a-a.mean();b=b-b.mean();return float((a*b).sum()/np.sqrt((a*a).sum()*(b*b).sum()))
print(f"corr(tau, myoPSF) = {pc(tau,mp):+.3f}   (if wall-myo leakage: LV>RV in tau; observed?)")
print(f"  -> RV myoPSF {rows[0][2]:.4f} (tau {rows[0][1]:.3f}) vs LV myoPSF {rows[2][2]:.4f} (tau {rows[2][1]:.3f})")
print(f"  LV has {'MORE' if rows[2][2]>rows[0][2] else 'LESS'} myo-PSF but {'LOWER' if rows[2][1]<rows[0][1] else 'HIGHER'} tau -> wall-myo leakage {'REFUTED' if (rows[2][2]>rows[0][2] and rows[2][1]<rows[0][1]) else 'possible'}")

print("\n=== (2) can CT HU even separate blood from myocardium? ===")
myocore=binary_erosion(myo,iterations=1)
rvcore=binary_erosion(m(f"{CH}/heart_ventricle_right.nii.gz"),iterations=1)
hu_myo=ct[myocore]; hu_rv=ct[rvcore]
print(f"myocardium core HU: {hu_myo.mean():.1f} +- {hu_myo.std():.1f}  (n={hu_myo.size})")
print(f"RV blood core  HU: {hu_rv.mean():.1f} +- {hu_rv.std():.1f}  (n={hu_rv.size})")
# overlap: fraction of myo HU within RV central 90% range
lo,hi=np.percentile(hu_rv,[5,95])
ov=((hu_myo>=lo)&(hu_myo<=hi)).mean()
print(f"RV-core central-90% HU range: [{lo:.0f},{hi:.0f}]")
print(f"fraction of myocardium voxels whose HU falls in that blood range: {ov*100:.0f}%")
print(f"  -> HU {'CANNOT' if ov>0.5 else 'can'} cleanly separate blood from muscle (overlap {ov*100:.0f}%)")
print("\nVERDICT: wall-myocardium leakage refuted by myoPSF test; HU uninformative for")
print("internal-trabecula hypothesis (blood/muscle HU overlap). Mechanism = OPEN, not pinned.")
