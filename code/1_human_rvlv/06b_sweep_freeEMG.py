#!/usr/bin/env python3
import os as _o, sys as _s
_s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), '..'))
import config as cfg
"""Canonical sweep, take 2: per-chamber FREE 2-component EMG (background + short
+ o-Ps), seeded near global values but all params free (the fixed-IRF variant
was degenerate). Tail single-exp (lo=t0+1.5) as model-independent cross-check.
Builds the full uncertainty budget for the claim decision.
"""
import os, json, csv, numpy as np, nibabel as nib
from scipy.special import erfc
from scipy.optimize import curve_fit
from scipy.ndimage import binary_dilation

AUD=cfg.OUT
OUT=cfg.OUT; WRK=cfg.OUT
RV=f"{cfg.MASKS}/totalseg_chambers_11243763/heart_ventricle_right.nii.gz"
LV=f"{cfg.MASKS}/totalseg_chambers_11243763/heart_ventricle_left.nii.gz"
HRT=f"{cfg.MASKS}/totalseg_11243763/heart.nii.gz"
HIST=f"{cfg.RAW}/zenodo_11243763/20230606_histoimage3d_lm.nii.gz"

cache=np.load(os.path.join(WRK,"heart_region_events.npy"))
Nx,Ny,Nz=nib.load(HIST).shape
rv=np.asanyarray(nib.load(RV).dataobj)>0; lv=np.asanyarray(nib.load(LV).dataobj)>0
heart=np.asanyarray(nib.load(HRT).dataobj)>0
vox_mm=np.array([1.63282919,1.63713944,1.64513969])
ri0=np.round(cache[:,0]).astype(np.int32); rj0=np.round(cache[:,1]).astype(np.int32)
rk0=np.round(cache[:,2]).astype(np.int32); tt=cache[:,3].astype(np.float32)
region=binary_dilation(rv|lv,iterations=11)
ib0=(ri0>=0)&(ri0<Nx)&(rj0>=0)&(rj0<Ny)&(rk0>=0)&(rk0<Nz)
keep=np.zeros(cache.shape[0],bool); keep[ib0]=region[ri0[ib0],rj0[ib0],rk0[ib0]]
ri0,rj0,rk0,tt=ri0[keep],rj0[keep],rk0[keep],tt[keep]

BW=0.10; edges=np.arange(-15,15+BW,BW); ctr=0.5*(edges[:-1]+edges[1:])
neg=(ctr>-12)&(ctr<-6); FITW=(ctr>=-3)&(ctr<=12); xw=ctr[FITW]
def emg(t,t0,s,tau):
    x=t-t0; z=(s/tau-x/s)/np.sqrt(2)
    return (1.0/(2*tau))*np.exp((s**2)/(2*tau**2)-x/tau)*erfc(z)
def model(t,B,A,f1,t0,s,tau1,tau2):
    return B+A*(f1*emg(t,t0,s,tau1)+(1-f1)*emg(t,t0,s,tau2))
LO=[0,0,0.05,-1,0.05,0.1,1.0]; HI=[np.inf,np.inf,0.98,1.0,0.6,0.9,5.0]
def fit_emg(t):
    if t.size<5000: return np.nan,np.nan,np.nan
    h,_=np.histogram(t,bins=edges); h=h.astype(float)
    y=h[FITW]; w=1.0/np.sqrt(np.maximum(y,1)); B0=np.median(h[neg])
    A0=max((y.sum()-B0*len(y))*BW,1.0)
    try:
        p,c=curve_fit(model,xw,y,p0=[B0,A0,0.55,0.066,0.108,0.39,1.8],
                      sigma=1/w,absolute_sigma=True,bounds=(LO,HI),maxfev=120000)
        chi=np.sum(((y-model(xw,*p))*w)**2)/(len(xw)-7)
        return float(p[6]),float(np.sqrt(c[6,6])),float(chi)
    except Exception:
        return np.nan,np.nan,np.nan
def tau_tail(t,t0=0.066):
    if t.size<5000: return np.nan
    h,_=np.histogram(t,bins=edges); h=h.astype(float); B=np.median(h[neg]); sig=h-B
    v=(ctr>=t0+1.5)&(ctr<=9.0)&(sig>3*np.sqrt(max(B,1)))
    if v.sum()<6: return np.nan
    xx=ctr[v]; yy=np.log(np.maximum(sig[v],1e-6)); W=(1.0/np.maximum(np.sqrt(h[v]),1))**2
    X=np.vstack([np.ones_like(xx),xx]).T
    beta=np.linalg.solve((X*W[:,None]).T@X,(X*W[:,None]).T@yy)
    return float(-1.0/beta[1])
def chamber_t(di,dj,dk,mask):
    i=ri0+di;j=rj0+dj;k=rk0+dk
    ib=(i>=0)&(i<Nx)&(j>=0)&(j<Ny)&(k>=0)&(k<Nz)
    sel=np.zeros(ri0.size,bool); sel[ib]=mask[i[ib],j[ib],k[ib]]; return tt[sel]
def measure(di,dj,dk):
    tr=chamber_t(di,dj,dk,rv); tl=chamber_t(di,dj,dk,lv)
    ar,er,cr=fit_emg(tr); al,el,cl=fit_emg(tl)
    return dict(di=di,dj=dj,dk=dk,r_mm=float(np.sqrt(((np.array([di,dj,dk])*vox_mm)**2).sum())),
                n_rv=int(tr.size),n_lv=int(tl.size),tau_rv=ar,tau_lv=al,d_rv_lv=ar-al,
                tau_rv_err=er,tau_lv_err=el,chi_rv=cr,chi_lv=cl,
                tail_rv=tau_tail(tr),tail_lv=tau_tail(tl),tail_d=tau_tail(tr)-tau_tail(tl))
b=measure(0,0,0)
d_err=float(np.hypot(b['tau_rv_err'],b['tau_lv_err']))
print(f"BASELINE free-EMG: RV={b['tau_rv']:.4f}+-{b['tau_rv_err']:.4f} LV={b['tau_lv']:.4f}+-{b['tau_lv_err']:.4f}")
print(f"  RV-LV(EMG)={b['d_rv_lv']:+.4f}ns  ({b['d_rv_lv']/d_err:.2f} sigma_stat, d_err={d_err:.4f})  chi2 {b['chi_rv']:.2f}/{b['chi_lv']:.2f}")
print(f"  RV-LV(tail)={b['tail_d']:+.4f}ns  [RV {b['tail_rv']:.3f} LV {b['tail_lv']:.3f}]")

rows=[{**b,"sweep":"baseline"}]
S=range(-3,4)
for di in S:
  for dj in S:
    for dk in S:
      if di==dj==dk==0: continue
      rows.append({**measure(di,dj,dk),"sweep":"grid3"})
for ax in range(3):
  for d in range(-10,11):
    if d==0: continue
    sh=[0,0,0]; sh[ax]=d; rows.append({**measure(*sh),"sweep":f"axis{ax}"})

keys=["sweep","di","dj","dk","r_mm","n_rv","n_lv","tau_rv","tau_lv","d_rv_lv",
      "tau_rv_err","tau_lv_err","chi_rv","chi_lv","tail_rv","tail_lv","tail_d"]
with open(os.path.join(OUT,"registration_stability.csv"),"w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=keys); w.writeheader()
    [w.writerow({k:r.get(k) for k in keys}) for r in rows]

grid=[r for r in rows if r["sweep"]=="grid3" and np.isfinite(r["d_rv_lv"])]
dv=np.array([r["d_rv_lv"] for r in grid]); base=b["d_rv_lv"]
reg_sd=float(dv.std()); model_sys=abs(b["d_rv_lv"]-b["tail_d"])
total=float(np.sqrt(d_err**2+reg_sd**2+model_sys**2))
print(f"\n=== SMALL-GRID +-3vox free-EMG n={len(grid)} ===")
print(f"  baseline {base:+.4f} | median {np.median(dv):+.4f} | mean {dv.mean():+.4f} +- {reg_sd:.4f}")
print(f"  min/max {dv.min():+.4f}/{dv.max():+.4f} | frac>0 {100*(dv>0).mean():.1f}% | frac>+0.1 {100*(dv>0.1).mean():.1f}%")
for lo,hi in [(0,2),(2,4),(4,6),(6,9)]:
    s=np.array([r["d_rv_lv"] for r in grid if lo<r["r_mm"]<=hi])
    if s.size: print(f"    {lo}-{hi}mm n={s.size:3d} median {np.median(s):+.3f} range[{s.min():+.3f},{s.max():+.3f}] frac>0 {100*(s>0).mean():.0f}%")
tl=np.array([r["tail_d"] for r in grid if np.isfinite(r["tail_d"])])
print(f"  tail cross-check: median {np.median(tl):+.3f} frac>0 {100*(tl>0).mean():.0f}%")
print("\nPER-AXIS delta free-EMG -10..10:")
for ax,nm in [(0,"x"),(1,"y"),(2,"z")]:
    key={0:'di',1:'dj',2:'dk'}[ax]; v=[]
    for d in range(-10,11):
        rr=b if d==0 else next((r for r in rows if r["sweep"]==f"axis{ax}" and r[key]==d),None)
        v.append(f"{rr['d_rv_lv']:+.2f}" if rr and np.isfinite(rr['d_rv_lv']) else " nan")
    print(f"  {nm}: "+" ".join(v))
print(f"\n=== UNCERTAINTY BUDGET (RV-LV) ===")
print(f"  effect (baseline EMG)     : {base:+.3f} ns")
print(f"  statistical (1 sigma)     : {d_err:.3f} ns")
print(f"  model systematic EMG-tail : {model_sys:.3f} ns")
print(f"  registration sd (<=8mm)   : {reg_sd:.3f} ns")
print(f"  TOTAL combined 1 sigma    : {total:.3f} ns  -> effect/total = {base/total:.2f}")

with open(os.path.join(OUT,"rv_lv_results.csv"),"w",newline="") as f:
    wr=csv.writer(f); wr.writerow(["chamber","n_events","tau_oPs_ns","tau_err_ns","chi2_ndf","tau_tail_ns","method"])
    wr.writerow(["RV",b["n_rv"],round(b["tau_rv"],4),round(b["tau_rv_err"],4),round(b["chi_rv"],2),round(b["tail_rv"],4),"free 2-EMG"])
    wr.writerow(["LV",b["n_lv"],round(b["tau_lv"],4),round(b["tau_lv_err"],4),round(b["chi_lv"],2),round(b["tail_lv"],4),"free 2-EMG"])
    wr.writerow(["RV-LV","",round(b["d_rv_lv"],4),round(d_err,4),"",round(b["tail_d"],4),"delta"])
json.dump({"baseline":b,"d_err":d_err,"reg_sd":reg_sd,"model_sys":model_sys,"total":total,
           "effect_over_total":base/total,"grid_n":len(grid),"frac_pos":float((dv>0).mean()),
           "grid_median":float(np.median(dv)),"grid_min":float(dv.min()),"grid_max":float(dv.max())},
          open(os.path.join(WRK,"budget.json"),"w"),indent=2,default=float)
print("\nWROTE outputs/registration_stability.csv + rv_lv_results.csv")
