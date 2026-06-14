#!/usr/bin/env python3
import os as _o, sys as _s
_s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), '..'))
import config as cfg
"""Mechanism test: is the venous>arterial / RV>LV tau elevation driven by
partial-volume contamination from adjacent LUNG AIR (air = huge free volume ->
very long o-Ps tau)? Build a lung-air-fraction field (lung lobe masks convolved
with the PET PSF), then:
 (A) per-compartment: correlate tau vs lung-air-fraction vs oxygenation;
 (B) decisive matched test: stratify cardiac blood by lung-air-fraction and,
     within each stratum, compare venous(RV,RA) vs arterial(LV,LA). If the
     contrast vanishes at matched air-fraction -> lung partial volume, not O2.
"""
import os, csv, numpy as np, nibabel as nib
from scipy.special import erfc
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter, binary_erosion

AUD=cfg.OUT
OUT=cfg.OUT; WRK=cfg.OUT
CH=f"{cfg.MASKS}/totalseg_chambers_11243763"; TS=f"{cfg.MASKS}/totalseg_11243763"
HIST=f"{cfg.RAW}/zenodo_11243763/20230606_histoimage3d_lm.nii.gz"
CTp=f"{cfg.RAW}/zenodo_11243763/20230606_lm_2_AC_CT_WB_1.65mm.nii.gz"
himg=nib.load(HIST); Nx,Ny,Nz=himg.shape; A=himg.affine
ox,oy,oz=A[0,3],A[1,3],A[2,3]; sxx,syy,szz=A[0,0],A[1,1],A[2,2]
ctd=np.asanyarray(nib.load(CTp).dataobj).astype(np.float32)

# lung-air fraction field: union of 5 lung lobes, convolved with PET PSF
lung=np.zeros((Nx,Ny,Nz),bool)
for f in ["lung_lower_lobe_left","lung_lower_lobe_right","lung_middle_lobe_right",
          "lung_upper_lobe_left","lung_upper_lobe_right"]:
    lung|=np.asanyarray(nib.load(f"{TS}/{f}.nii.gz").dataobj)>0
SIGMA=1.5  # voxels (~5mm FWHM PET PSF)
airf=gaussian_filter(lung.astype(np.float32),sigma=SIGMA)
print(f"lung voxels={int(lung.sum())}  air-frac field range [{airf.min():.3f},{airf.max():.3f}]  PSF sigma={SIGMA}vox")

BW=0.1; edges=np.arange(-15,15+BW,BW); ctr=0.5*(edges[:-1]+edges[1:])
neg=(ctr>-12)&(ctr<-6); FW=(ctr>=-3)&(ctr<=12); xw=ctr[FW]
def emg(tt,p0,s,ta):
    x=tt-p0; z=(s/ta-x/s)/np.sqrt(2); return (1/(2*ta))*np.exp((s**2)/(2*ta**2)-x/ta)*erfc(z)
def model(tt,B,Am,f1,p0,s,t1,t2): return B+Am*(f1*emg(tt,p0,s,t1)+(1-f1)*emg(tt,p0,s,t2))
def fit(tarr):
    if tarr.size<5000: return np.nan,np.nan,int(tarr.size)
    h,_=np.histogram(tarr,bins=edges);h=h.astype(float);y=h[FW];w=1/np.sqrt(np.maximum(y,1))
    B0=np.median(h[neg]);A0=max((y.sum()-B0*len(y))*BW,1)
    try:
        p,c=curve_fit(model,xw,y,p0=[B0,A0,0.55,0.066,0.108,0.39,1.8],sigma=1/w,absolute_sigma=True,
                      bounds=([0,0,0.05,-1,0.05,0.1,1],[np.inf,np.inf,0.98,1,0.6,0.9,5]),maxfev=120000)
        return float(p[6]),float(np.sqrt(c[6,6])),int(tarr.size)
    except Exception: return np.nan,np.nan,int(tarr.size)

# ---------- (A) per-compartment: tau vs air-fraction ----------
comps={"SVC":(f"{TS}/superior_vena_cava.nii.gz","venous"),"IVC":(f"{TS}/inferior_vena_cava.nii.gz","venous"),
 "RA":(f"{CH}/heart_atrium_right.nii.gz","venous"),"RV":(f"{CH}/heart_ventricle_right.nii.gz","venous"),
 "pulm_ARTERY":(f"{CH}/pulmonary_artery.nii.gz","venous"),"LA":(f"{CH}/heart_atrium_left.nii.gz","arterial"),
 "LV":(f"{CH}/heart_ventricle_left.nii.gz","arterial"),"aorta":(f"{CH}/aorta.nii.gz","arterial")}
lad={r["name"]:r for r in csv.DictReader(open(os.path.join(OUT,"oxygenation_ladder.csv")))}
print("\n=== (A) per-compartment: tau_core vs lung-air-fraction ===")
print(f"{'comp':13s}{'class':9s}{'tau':>7s}{'air_frac':>10s}{'meanHU':>8s}")
recs=[]
for nm,(p,cls) in comps.items():
    if nm not in lad: continue
    core=binary_erosion(np.asanyarray(nib.load(p).dataobj)>0,iterations=1)
    af=float(airf[core].mean()); hu=float(ctd[core].mean()); tau=float(lad[nm]["tau_core"])
    so2=0.72 if cls=="venous" else 0.98
    recs.append(dict(name=nm,cls=cls,so2=so2,tau=tau,airf=af,HU=hu))
    print(f"{nm:13s}{cls:9s}{tau:7.3f}{af:10.4f}{hu:8.0f}")
ta=np.array([r["tau"] for r in recs]); af=np.array([r["airf"] for r in recs])
so=np.array([r["so2"] for r in recs])
def pc(a,b): a=a-a.mean();b=b-b.mean();return float((a*b).sum()/np.sqrt((a*a).sum()*(b*b).sum()))
def partial(x,y,c):
    C=np.column_stack([np.ones_like(x),c]); bx=np.linalg.lstsq(C,x,rcond=None)[0]; by=np.linalg.lstsq(C,y,rcond=None)[0]
    return pc(x-C@bx,y-C@by)
print(f"  corr(tau, lung_air_frac)          = {pc(ta,af):+.3f}")
print(f"  corr(tau, oxygenation SO2)        = {pc(ta,so):+.3f}")
print(f"  partial corr(tau,SO2 | air_frac)  = {partial(so,ta,af):+.3f}   <-- does O2 survive controlling lung air?")
print(f"  partial corr(tau,air_frac | SO2)  = {partial(af,ta,so):+.3f}   <-- does lung air survive controlling O2?")

# ---------- (B) decisive matched stratification (cardiac blood from cache) ----------
cache=np.load(os.path.join(WRK,"heart_region_events.npy"))
ri=np.round(cache[:,0]).astype(int);rj=np.round(cache[:,1]).astype(int);rk=np.round(cache[:,2]).astype(int)
ib=(ri>=0)&(ri<Nx)&(rj>=0)&(rj<Ny)&(rk>=0)&(rk<Nz); t=cache[:,3]
def mask(p): return np.asanyarray(nib.load(p).dataobj)>0
rv=mask(f"{CH}/heart_ventricle_right.nii.gz"); lv=mask(f"{CH}/heart_ventricle_left.nii.gz")
ra=mask(f"{CH}/heart_atrium_right.nii.gz"); la=mask(f"{CH}/heart_atrium_left.nii.gz")
ven=np.zeros(cache.shape[0],bool); art=np.zeros(cache.shape[0],bool)
ven[ib]=rv[ri[ib],rj[ib],rk[ib]]|ra[ri[ib],rj[ib],rk[ib]]
art[ib]=lv[ri[ib],rj[ib],rk[ib]]|la[ri[ib],rj[ib],rk[ib]]
ev_af=np.full(cache.shape[0],np.nan); ev_af[ib]=airf[ri[ib],rj[ib],rk[ib]]
blood=(ven|art)&np.isfinite(ev_af)
print(f"\n=== (B) matched stratification by lung-air-fraction (cardiac blood, n={int(blood.sum())}) ===")
# overall tau vs air-fraction (pooled venous+arterial)
qs=np.quantile(ev_af[blood],[0,1/3,2/3,1.0])
print(f"air-frac tertile edges: {[round(x,4) for x in qs]}")
print(f"{'air-frac stratum':22s}{'tau_VENOUS(RV+RA)':>20s}{'tau_ARTERIAL(LV+LA)':>22s}{'V-A':>8s}")
strat=[]
for i in range(3):
    lo,hi=qs[i],qs[i+1]
    inb=blood&(ev_af>=lo)&(ev_af<=hi if i==2 else ev_af<hi)
    tv,ev,nv=fit(t[inb&ven]); tau_a,ea,na=fit(t[inb&art])
    d=tv-tau_a
    print(f"[{lo:.3f},{hi:.3f}] air {('LOW' if i==0 else 'MID' if i==1 else 'HIGH'):5s}{tv:8.3f}+-{ev:.3f}({nv//1000:>4}k){tau_a:9.3f}+-{ea:.3f}({na//1000:>4}k){d:+8.3f}")
    strat.append((0.5*(lo+hi),tv,tau_a,d))
# pooled tau vs air-fraction (both classes together)
print("\npooled cardiac blood, tau vs air-fraction (5 bins):")
qs5=np.quantile(ev_af[blood],np.linspace(0,1,6))
for i in range(5):
    lo,hi=qs5[i],qs5[i+1]
    inb=blood&(ev_af>=lo)&(ev_af<=hi if i==4 else ev_af<hi)
    tt,ee,nn=fit(t[inb]); print(f"  air[{lo:.3f},{hi:.3f}] <af>={ev_af[inb].mean():.3f}  tau={tt:.3f}+-{ee:.3f}  n={nn//1000}k")

print("\nINTERPRETATION:")
print("  (A) if corr(tau,air)>>corr(tau,SO2) and partial(SO2|air)->0: lung air explains it.")
print("  (B) if venous-arterial gap shrinks toward 0 in matched air strata: confound, not O2.")
print("      if gap persists at matched air: residual real contrast (O2 or other blood property).")

# plot
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig,ax=plt.subplots(1,2,figsize=(13,5))
cols={"venous":"#c0392b","arterial":"#2471a3"}
for r in recs:
    ax[0].scatter(r["airf"],r["tau"],c=cols[r["cls"]],s=70)
    ax[0].annotate(r["name"],(r["airf"],r["tau"]),fontsize=8,xytext=(4,3),textcoords="offset points")
ax[0].set_xlabel("mean lung-air fraction (PET-PSF)"); ax[0].set_ylabel("tau_oPs core (ns)")
ax[0].set_title(f"(A) tau vs lung-air: r={pc(ta,af):+.2f}  (vs O2 r={pc(ta,so):+.2f})")
sc=np.array(strat)
ax[1].plot(sc[:,0],sc[:,1],"o-",color="#c0392b",label="venous (RV+RA)")
ax[1].plot(sc[:,0],sc[:,2],"o-",color="#2471a3",label="arterial (LV+LA)")
ax[1].set_xlabel("lung-air fraction stratum"); ax[1].set_ylabel("tau_oPs (ns)")
ax[1].set_title("(B) matched strata: does venous>arterial survive?"); ax[1].legend()
fig.tight_layout(); fig.savefig(os.path.join(OUT,"plots","fig6_lung_air_confound.png"),dpi=130)
print("\nWROTE plots/fig6_lung_air_confound.png")
