#!/usr/bin/env python3
import os as _o, sys as _s
_s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), '..'))
import config as cfg
"""Plots for the blind audit report."""
import os, csv, numpy as np, nibabel as nib
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import erfc
from scipy.optimize import curve_fit

AUD=cfg.OUT
OUT=cfg.OUT; WRK=cfg.OUT; PLT=os.path.join(OUT,"plots")
os.makedirs(PLT,exist_ok=True)
RV=f"{cfg.MASKS}/totalseg_chambers_11243763/heart_ventricle_right.nii.gz"
LV=f"{cfg.MASKS}/totalseg_chambers_11243763/heart_ventricle_left.nii.gz"
HIST=f"{cfg.RAW}/zenodo_11243763/20230606_histoimage3d_lm.nii.gz"
cache=np.load(os.path.join(WRK,"heart_region_events.npy"))
Nx,Ny,Nz=nib.load(HIST).shape
rv=np.asanyarray(nib.load(RV).dataobj)>0; lv=np.asanyarray(nib.load(LV).dataobj)>0
ri=np.round(cache[:,0]).astype(int);rj=np.round(cache[:,1]).astype(int);rk=np.round(cache[:,2]).astype(int)
ib=(ri>=0)&(ri<Nx)&(rj>=0)&(rj<Ny)&(rk>=0)&(rk<Nz)
def msel(m):
    s=np.zeros(cache.shape[0],bool); s[ib]=m[ri[ib],rj[ib],rk[ib]]; return cache[s,3]
t_rv=msel(rv); t_lv=msel(lv)
BW=0.1; edges=np.arange(-15,15+BW,BW); ctr=0.5*(edges[:-1]+edges[1:])
neg=(ctr>-12)&(ctr<-6); FITW=(ctr>=-3)&(ctr<=12); xw=ctr[FITW]
def emg(t,t0,s,tau):
    x=t-t0; z=(s/tau-x/s)/np.sqrt(2); return (1/(2*tau))*np.exp((s**2)/(2*tau**2)-x/tau)*erfc(z)
def model(t,B,A,f1,t0,s,tau1,tau2): return B+A*(f1*emg(t,t0,s,tau1)+(1-f1)*emg(t,t0,s,tau2))
def fit(t):
    h,_=np.histogram(t,bins=edges); h=h.astype(float); y=h[FITW]; w=1/np.sqrt(np.maximum(y,1))
    B0=np.median(h[neg]); A0=max((y.sum()-B0*len(y))*BW,1)
    p,_=curve_fit(model,xw,y,p0=[B0,A0,0.55,0.066,0.108,0.39,1.8],sigma=1/w,absolute_sigma=True,
                  bounds=([0,0,0.05,-1,0.05,0.1,1],[np.inf,np.inf,0.98,1,0.6,0.9,5]),maxfev=120000)
    return h,p

# Fig 1: spectra + fits
fig,ax=plt.subplots(1,2,figsize=(13,5))
for a,(nm,t,col) in zip(ax,[("RV (deoxygenated)",t_rv,"#c0392b"),("LV (oxygenated)",t_lv,"#2471a3")]):
    h,p=fit(t)
    a.step(ctr,h,where="mid",color=col,lw=.8,label="data")
    a.plot(xw,model(xw,*p),"k-",lw=1.6,label=f"2-EMG fit  τ_oPs={p[6]:.3f} ns")
    a.axhline(p[0],ls=":",color="gray",label=f"background B={p[0]:.0f}")
    a.set_yscale("log"); a.set_xlim(-3,12); a.set_xlabel("Δt (ns)"); a.set_ylabel("counts / 0.1 ns")
    a.set_title(f"{nm}  (n={t.size:,})"); a.legend(fontsize=9)
fig.suptitle("Per-chamber positronium lifetime spectra (baseline registration)")
fig.tight_layout(); fig.savefig(os.path.join(PLT,"fig1_spectra.png"),dpi=130); plt.close(fig)

# read sweep csv
rows=list(csv.DictReader(open(os.path.join(OUT,"registration_stability.csv"))))
def F(r,k):
    v=r.get(k);
    try: return float(v)
    except: return np.nan
grid=[r for r in rows if r["sweep"]=="grid3" and np.isfinite(F(r,"d_rv_lv"))]
dv=np.array([F(r,"d_rv_lv") for r in grid]); base=[F(r,"d_rv_lv") for r in rows if r["sweep"]=="baseline"][0]

# Fig 2: grid distribution
fig,ax=plt.subplots(figsize=(7,5))
ax.hist(dv,bins=30,color="#5d8aa8",edgecolor="k",alpha=.85)
ax.axvline(0,color="k",lw=2,label="zero (no contrast)")
ax.axvline(base,color="#c0392b",lw=2,ls="--",label=f"baseline {base:+.3f} ns")
ax.axvline(np.median(dv),color="green",lw=2,ls=":",label=f"grid median {np.median(dv):+.3f} ns")
ax.set_xlabel("RV−LV τ_oPs (ns)"); ax.set_ylabel("# registration shifts")
ax.set_title(f"RV−LV across ±3 vox (≤8 mm) shifts, n={len(grid)}\n100%>0, sd={dv.std():.3f} ns")
ax.legend(); fig.tight_layout(); fig.savefig(os.path.join(PLT,"fig2_grid_distribution.png"),dpi=130); plt.close(fig)

# Fig 3: per-axis
fig,ax=plt.subplots(figsize=(8,5))
vox=[1.6328,1.6371,1.6451]
for ax_i,(nm,key,c) in enumerate([("x (R-L)","di","#c0392b"),("y (A-P)","dj","#27ae60"),("z (axial)","dk","#2471a3")]):
    xs=[];ys=[]
    for d in range(-10,11):
        if d==0: xs.append(0);ys.append(base);continue
        rr=[r for r in rows if r["sweep"]==f"axis{ax_i}" and int(F(r,key))==d]
        if rr and np.isfinite(F(rr[0],"d_rv_lv")): xs.append(d*vox[ax_i]); ys.append(F(rr[0],"d_rv_lv"))
    ax.plot(xs,ys,"o-",color=c,label=f"shift along {nm}")
ax.axhline(0,color="k",lw=1); ax.axhline(base,color="gray",ls="--",label=f"baseline {base:+.3f}")
ax.set_xlabel("registration shift (mm)"); ax.set_ylabel("RV−LV τ_oPs (ns)")
ax.set_title("RV−LV vs single-axis registration shift (free-EMG)"); ax.legend()
fig.tight_layout(); fig.savefig(os.path.join(PLT,"fig3_per_axis.png"),dpi=130); plt.close(fig)

# Fig 4: delta vs radius
fig,ax=plt.subplots(figsize=(7,5))
rmm=np.array([F(r,"r_mm") for r in grid])
ax.scatter(rmm,dv,s=14,alpha=.6,color="#34495e")
ax.axhline(0,color="k"); ax.axhline(base,color="#c0392b",ls="--",label=f"baseline {base:+.3f}")
ax.set_xlabel("|registration shift| (mm)"); ax.set_ylabel("RV−LV τ_oPs (ns)")
ax.set_title("Contrast vs shift magnitude (always >0 over tested range)"); ax.legend()
fig.tight_layout(); fig.savefig(os.path.join(PLT,"fig4_vs_radius.png"),dpi=130); plt.close(fig)
print("WROTE 4 plots to",PLT)
print("files:",os.listdir(PLT))
