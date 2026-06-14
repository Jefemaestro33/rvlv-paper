#!/usr/bin/env python3
import os as _o, sys as _s
_s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), '..'))
import config as cfg
"""Clean 82Rb material characterization: isolate the dense disk source (tight box
around its centroid) from the FOV-wide scatter/randoms. Fit the clean disk core
for IRF (t0, sigma) and tau; fit the scatter separately; check tau consistency
within the core (local positional bias). Compare IRF to the human-fit assumptions.
"""
import os, numpy as np
from scipy.special import erfc
from scipy.optimize import curve_fit
LM=f"{cfg.RAW}/zenodo_12636019/Rb82_coins_Histo_Out.l"
n=os.path.getsize(LM)//32
mm=np.memmap(LM,dtype=np.float64,mode="r",shape=(n,4))
offs=np.linspace(0,n-200000,200).astype(np.int64)
samp=np.concatenate([np.array(mm[o:o+200000]) for o in offs],axis=0)
xyz=samp[:,:3]; t=samp[:,3]; phys=np.abs(t)<50
c=np.median(xyz[phys],axis=0)
print("disk centroid:",[round(x,1) for x in c])
HW=40.0  # half-width box (mm)
core=phys & (np.abs(xyz[:,0]-c[0])<HW)&(np.abs(xyz[:,1]-c[1])<HW)&(np.abs(xyz[:,2]-c[2])<HW)
scat=phys & ~core & (np.linalg.norm(xyz-c,axis=1)>150)
print(f"core (box +-{HW:.0f}mm): {core.sum():,} events ({100*core.sum()/phys.sum():.1f}% of physical)")
print(f"scatter (>150mm): {scat.sum():,} events")

BW=0.1; edges=np.arange(-15,15+BW,BW); ctr=0.5*(edges[:-1]+edges[1:])
neg=(ctr>-12)&(ctr<-6); FW=(ctr>=-3)&(ctr<=12); xw=ctr[FW]
def emg(tt,t0,s,ta):
    x=tt-t0; z=(s/ta-x/s)/np.sqrt(2); return (1/(2*ta))*np.exp((s**2)/(2*ta**2)-x/ta)*erfc(z)
def model(tt,B,A,f1,t0,s,t1,t2): return B+A*(f1*emg(tt,t0,s,t1)+(1-f1)*emg(tt,t0,s,t2))
def fit(tt,label=""):
    if tt.size<5000: return None
    h,_=np.histogram(tt,bins=edges); h=h.astype(float); y=h[FW]; w=1/np.sqrt(np.maximum(y,1))
    B0=np.median(h[neg]); A0=max((y.sum()-B0*len(y))*BW,1)
    try:
        p,co=curve_fit(model,xw,y,p0=[B0,A0,0.5,0.08,0.12,0.4,1.6],sigma=1/w,absolute_sigma=True,
                       bounds=([0,0,0.02,-1,0.05,0.1,1.0],[np.inf,np.inf,0.98,2,0.6,0.9,5]),maxfev=200000)
        chi=np.sum(((y-model(xw,*p))*w)**2)/(len(xw)-7)
        return dict(t0=p[3],t0e=np.sqrt(co[3,3]),sig=p[4],sige=np.sqrt(co[4,4]),
                    tau1=p[5],tau2=p[6],tau2e=np.sqrt(co[6,6]),Iops=1-p[2],chi=chi,n=int(tt.size),bg=p[0])
    except Exception as e: print("fail",label,e); return None

print("\n=== CLEAN DISK CORE (known 82Rb material) ===")
g=fit(t[core],"core")
print(f"  t0={g['t0']:.4f}+-{g['t0e']:.4f}  sigma(IRF)={g['sig']:.4f}+-{g['sige']:.4f}  "
      f"tau_oPs={g['tau2']:.4f}+-{g['tau2e']:.4f}  tau_short={g['tau1']:.3f}  I_oPs={g['Iops']:.3f}  chi2/ndf={g['chi']:.1f}  n={g['n']:,}")
print(f"  vs HUMAN-FIT ASSUMPTIONS: t0=0.066, sigma=0.108")
print(f"  -> sigma_material={g['sig']:.3f} (human 0.108; diff {g['sig']-0.108:+.3f})  t0_material={g['t0']:.3f} (human 0.066; diff {g['t0']-0.066:+.3f})")

s=fit(t[scat],"scatter")
if s: print(f"\n  [scatter pop. for contrast]: tau_oPs={s['tau2']:.3f}  bg-dominated; chi2/ndf={s['chi']:.1f}  n={s['n']:,}")

print("\n=== within-core tau consistency (local positional bias) ===")
xc,yc,zc=xyz[core,0],xyz[core,1],xyz[core,2]; tc=t[core]
taus=[]
for ax,nm,arr in [(0,'x',xc),(1,'y',yc),(2,'z',zc)]:
    med=np.median(arr)
    for half,sel in [("lo",arr<med),("hi",arr>=med)]:
        r=fit(tc[sel])
        if r and r['tau2e']<0.1: taus.append((f"{nm}-{half}",r['tau2'],r['tau2e']))
for nm,v,e in taus: print(f"  {nm:6s} tau={v:.3f}+-{e:.3f}")
if taus:
    vv=[v for _,v,_ in taus]
    print(f"  -> within-core tau spread = {max(vv)-min(vv):.3f} ns  (this is the local 82Rb positional/timing wobble)")
print("\nReads: sigma_material vs 0.108 = IRF validation; within-core spread = isotope-matched local instrument bias.")
