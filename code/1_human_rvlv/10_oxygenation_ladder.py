#!/usr/bin/env python3
import os as _o, sys as _s
_s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), '..'))
import config as cfg
"""Definitive falsification test: does tau_oPs order along the whole oxygenation
ladder, not just RV>LV? Key built-in controls:
 - pulmonary ARTERY is DEoxygenated -> should sit with venous (long tau)
 - pulmonary VEIN is OXygenated     -> should sit with arterial (short tau)
 (an 'artery vs vein' anatomical confound predicts the opposite for these.)
Also records position (depth r, axial z) and CT HU per compartment as confound
covariates. One streaming pass over the full 529.8M-event listmode.
"""
import os, csv, numpy as np, nibabel as nib, time
from scipy.special import erfc
from scipy.optimize import curve_fit
from scipy.ndimage import binary_erosion

AUD=cfg.OUT
OUT=cfg.OUT; WRK=cfg.OUT
CH=f"{cfg.MASKS}/totalseg_chambers_11243763"
TS=f"{cfg.MASKS}/totalseg_11243763"
LM=f"{cfg.RAW}/zenodo_11243763/20230606_positronium_patient_evaluated_data_Histo_Out.l"
HIST=f"{cfg.RAW}/zenodo_11243763/20230606_histoimage3d_lm.nii.gz"
CT=f"{cfg.RAW}/zenodo_11243763/20230606_lm_2_AC_CT_WB_1.65mm.nii.gz"

# name, path, class, nominal SO2
COMPS=[
 ("SVC",        f"{TS}/superior_vena_cava.nii.gz","venous",0.72),
 ("IVC",        f"{TS}/inferior_vena_cava.nii.gz","venous",0.72),
 ("brachioceph_vein_L",f"{TS}/brachiocephalic_vein_left.nii.gz","venous",0.72),
 ("brachioceph_vein_R",f"{TS}/brachiocephalic_vein_right.nii.gz","venous",0.72),
 ("RA",         f"{CH}/heart_atrium_right.nii.gz","venous",0.72),
 ("RV",         f"{CH}/heart_ventricle_right.nii.gz","venous",0.72),
 ("pulm_ARTERY",f"{CH}/pulmonary_artery.nii.gz","venous",0.75),   # deoxygenated!
 ("pulm_VEIN",  f"{TS}/pulmonary_vein.nii.gz","arterial",0.98),    # oxygenated!
 ("LA",         f"{CH}/heart_atrium_left.nii.gz","arterial",0.98),
 ("LV",         f"{CH}/heart_ventricle_left.nii.gz","arterial",0.98),
 ("aorta",      f"{CH}/aorta.nii.gz","arterial",0.98),
 ("brachioceph_trunk",f"{TS}/brachiocephalic_trunk.nii.gz","arterial",0.98),
 ("carotid_L",  f"{TS}/common_carotid_artery_left.nii.gz","arterial",0.98),
 ("carotid_R",  f"{TS}/common_carotid_artery_right.nii.gz","arterial",0.98),
 ("subclavian_L",f"{TS}/subclavian_artery_left.nii.gz","arterial",0.98),
 ("subclavian_R",f"{TS}/subclavian_artery_right.nii.gz","arterial",0.98),
 ("myocardium", f"{CH}/heart_myocardium.nii.gz","muscle",np.nan), # structural control
]
mapping={"perm":[0,2,1],"sgn":[-1,1,-1]}
himg=nib.load(HIST); Nx,Ny,Nz=himg.shape; A=himg.affine
sxx,syy,szz=A[0,0],A[1,1],A[2,2]; ox,oy,oz=A[0,3],A[1,3],A[2,3]
ctd=np.asanyarray(nib.load(CT).dataobj)

# build label volumes (full + erode1), overlap-excluded
labels_full=np.zeros((Nx,Ny,Nz),np.int16)
labels_core=np.zeros((Nx,Ny,Nz),np.int16)
cnt=np.zeros((Nx,Ny,Nz),np.uint8)
meta=[]
for idx,(nm,p,cls,so2) in enumerate(COMPS,start=1):
    m=np.asanyarray(nib.load(p).dataobj)>0
    cnt+=m.astype(np.uint8)
    labels_full[m & (labels_full==0)]=idx
    core=binary_erosion(m,iterations=1)
    labels_core[core & (labels_core==0)]=idx
    # position (histo/scanner world) + HU
    ijk=np.argwhere(m); c=ijk.mean(0)
    wx=A[0,0]*c[0]+ox; wy=A[1,1]*c[1]+oy; wz=A[2,2]*c[2]+oz
    hu=float(ctd[m].mean())
    meta.append(dict(idx=idx,name=nm,cls=cls,so2=so2,nvox=int(m.sum()),
                     r_mm=float(np.hypot(wx,wy)),z_mm=float(wz),HU=hu))
    print(f"  {nm:20s} vox={int(m.sum()):7d} core_vox={int(core.sum()):7d} r={np.hypot(wx,wy):5.0f} z={wz:6.0f} HU={hu:6.0f}")
# exclude voxels claimed by >1 full mask
overlap=cnt>1
labels_full[overlap]=0; labels_core[overlap]=0
print("overlap voxels excluded:",int(overlap.sum()))
lf=labels_full.reshape(-1); lc=labels_core.reshape(-1)

# streaming pass -> per-label dt histograms
BW=0.1; edges=np.arange(-15,15+BW,BW); ctr=0.5*(edges[:-1]+edges[1:]); nb=len(ctr)
NL=len(COMPS)+1
Hf=np.zeros(NL*nb,np.int64); Hc=np.zeros(NL*nb,np.int64)
sz=os.path.getsize(LM); n=sz//32
mm=np.memmap(LM,dtype=np.float64,mode="r",shape=(n,4))
CHK=25_000_000; t0=time.time()
for s in range(0,n,CHK):
    e=min(s+CHK,n); a=np.asarray(mm[s:e])
    wx=mapping["sgn"][0]*a[:,mapping["perm"][0]]; wy=mapping["sgn"][1]*a[:,mapping["perm"][1]]; wz=mapping["sgn"][2]*a[:,mapping["perm"][2]]
    i=np.round((wx-ox)/sxx); j=np.round((wy-oy)/syy); k=np.round((wz-oz)/szz)
    ib=(i>=0)&(i<Nx)&(j>=0)&(j<Ny)&(k>=0)&(k<Nz)
    t=a[:,3]; tb=np.floor((t+15)/BW).astype(np.int64); tin=(tb>=0)&(tb<nb)
    good=ib&tin
    vflat=(i[good]*Ny*Nz+j[good]*Nz+k[good]).astype(np.int64); tbb=tb[good]
    Lf=lf[vflat]; Lc=lc[vflat]
    mf=Lf>0; Hf+=np.bincount((Lf[mf]*nb+tbb[mf]),minlength=NL*nb)
    mc=Lc>0; Hc+=np.bincount((Lc[mc]*nb+tbb[mc]),minlength=NL*nb)
print(f"pass done {time.time()-t0:.0f}s")
Hf=Hf.reshape(NL,nb); Hc=Hc.reshape(NL,nb)

# fit
neg=(ctr>-12)&(ctr<-6); FW=(ctr>=-3)&(ctr<=12); xw=ctr[FW]
def emg(tt,p0,sg,ta):
    x=tt-p0; z=(sg/ta-x/sg)/np.sqrt(2); return (1/(2*ta))*np.exp((sg**2)/(2*ta**2)-x/ta)*erfc(z)
def model(tt,B,Amp,f1,p0,sg,t1,t2): return B+Amp*(f1*emg(tt,p0,sg,t1)+(1-f1)*emg(tt,p0,sg,t2))
def fit(h):
    h=h.astype(float); N=h.sum()
    if N<6000: return np.nan,np.nan,int(N)
    y=h[FW]; w=1/np.sqrt(np.maximum(y,1)); B0=np.median(h[neg]); A0=max((y.sum()-B0*len(y))*BW,1)
    try:
        p,c=curve_fit(model,xw,y,p0=[B0,A0,0.55,0.066,0.108,0.39,1.8],sigma=1/w,absolute_sigma=True,
                      bounds=([0,0,0.05,-1,0.05,0.1,1],[np.inf,np.inf,0.98,1,0.6,0.9,5]),maxfev=120000)
        return float(p[6]),float(np.sqrt(c[6,6])),int(N)
    except Exception: return np.nan,np.nan,int(N)

rows=[]
for md in meta:
    tc,ec,nc=fit(Hc[md["idx"]]); tf,ef,nf=fit(Hf[md["idx"]])
    rows.append({**md,"tau_core":tc,"tau_core_err":ec,"n_core":nc,"tau_full":tf,"tau_full_err":ef,"n_full":nf})

rows_sorted=sorted([r for r in rows if np.isfinite(r["tau_core"])],key=lambda r:-r["tau_core"])
print("\n=== OXYGENATION LADDER (sorted by core tau, long->short) ===")
print(f"{'compartment':20s}{'class':9s}{'SO2':>5s}{'tau_core':>14s}{'n_core':>9s}{'r_mm':>6s}{'z_mm':>7s}{'HU':>6s}")
for r in rows_sorted:
    print(f"{r['name']:20s}{r['cls']:9s}{r['so2']:5.2f}{r['tau_core']:8.3f}+-{r['tau_core_err']:.3f}{r['n_core']:9d}{r['r_mm']:6.0f}{r['z_mm']:7.0f}{r['HU']:6.0f}")

def grp(cls):
    v=[r for r in rows if r["cls"]==cls and np.isfinite(r["tau_core"]) and r["n_core"]>20000]
    return v
ven=grp("venous"); art=grp("arterial")
def wmean(v):
    w=np.array([1/r["tau_core_err"]**2 for r in v]); t=np.array([r["tau_core"] for r in v])
    return (w*t).sum()/w.sum()
print("\n=== GROUP TEST (core, n>20k) ===")
print(f"venous   (deoxy) compartments n={len(ven)}: weighted mean tau = {wmean(ven):.3f} ns  [{', '.join(r['name'] for r in ven)}]")
print(f"arterial (oxy)   compartments n={len(art)}: weighted mean tau = {wmean(art):.3f} ns  [{', '.join(r['name'] for r in art)}]")
print(f"venous - arterial separation = {wmean(ven)-wmean(art):+.3f} ns")
# crossover controls
for nm in ["pulm_ARTERY","pulm_VEIN"]:
    r=next((x for x in rows if x["name"]==nm),None)
    if r and np.isfinite(r["tau_core"]):
        print(f"  CROSSOVER {nm}: tau={r['tau_core']:.3f} (class={r['cls']}, predict { 'long/venous' if r['cls']=='venous' else 'short/arterial'})")
myo=next((x for x in rows if x["name"]=="myocardium"),None)
if myo: print(f"  CONTROL myocardium (muscle): tau={myo['tau_core']:.3f}")
# correlations
good=[r for r in rows if np.isfinite(r["tau_core"]) and r["n_core"]>20000 and r["cls"]!="muscle"]
ta=np.array([r["tau_core"] for r in good]); so=np.array([r["so2"] for r in good])
rr=np.array([r["r_mm"] for r in good]); zz=np.array([r["z_mm"] for r in good]); hu=np.array([r["HU"] for r in good])
def pear(a,b):
    a=a-a.mean(); b=b-b.mean(); return float((a*b).sum()/np.sqrt((a*a).sum()*(b*b).sum()))
print(f"\n=== CONFOUND CORRELATIONS (blood compartments, n={len(good)}) ===")
print(f"  tau vs SO2  (oxygenation): r = {pear(ta,so):+.3f}   (expect strong NEGATIVE if oxygenation-driven)")
print(f"  tau vs r    (depth)      : r = {pear(ta,rr):+.3f}")
print(f"  tau vs z    (axial pos)  : r = {pear(ta,zz):+.3f}")
print(f"  tau vs HU   (CT density) : r = {pear(ta,hu):+.3f}")

with open(os.path.join(OUT,"oxygenation_ladder.csv"),"w",newline="") as f:
    k=["name","cls","so2","nvox","n_core","tau_core","tau_core_err","n_full","tau_full","tau_full_err","r_mm","z_mm","HU"]
    w=csv.DictWriter(f,fieldnames=k); w.writeheader()
    for r in rows: w.writerow({kk:r.get(kk) for kk in k})
print("\nWROTE outputs/oxygenation_ladder.csv")
np.save(os.path.join(WRK,"ladder_hists.npy"),{"ctr":ctr,"Hc":Hc,"Hf":Hf,"meta":meta},allow_pickle=True)
