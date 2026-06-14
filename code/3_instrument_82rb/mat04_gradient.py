#!/usr/bin/env python3
import os as _o, sys as _s
_s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), '..'))
import config as cfg
"""Verify the isotope-matched instrumental positional tau-gradient in the uniform
82Rb disk material with a TIGHT core (reject edge contamination) + continuous
weighted-linear regression of tau vs position along each axis. Compare the
gradient to the human inter-chamber gradient (6.96 milli-ns/mm) to estimate how
much of the RV-LV contrast could be instrumental.
"""
import os, numpy as np
from scipy.special import erfc
from scipy.optimize import curve_fit
LM=f"{cfg.RAW}/zenodo_12636019/Rb82_coins_Histo_Out.l"
n=os.path.getsize(LM)//32
mm=np.memmap(LM,dtype=np.float64,mode="r",shape=(n,4))
offs=np.linspace(0,n-200000,300).astype(np.int64)
samp=np.concatenate([np.array(mm[o:o+200000]) for o in offs],axis=0)
xyz=samp[:,:3]; t=samp[:,3]; phys=np.abs(t)<50
c=np.median(xyz[phys],axis=0)
BW=0.1; edges=np.arange(-15,15+BW,BW); ctr=0.5*(edges[:-1]+edges[1:])
neg=(ctr>-12)&(ctr<-6); FW=(ctr>=-3)&(ctr<=12); xw=ctr[FW]
def emg(tt,t0,s,ta):
    x=tt-t0; z=(s/ta-x/s)/np.sqrt(2); return (1/(2*ta))*np.exp((s**2)/(2*ta**2)-x/ta)*erfc(z)
def model(tt,B,A,f1,t0,s,t1,t2): return B+A*(f1*emg(tt,t0,s,t1)+(1-f1)*emg(tt,t0,s,t2))
def fit(tt):
    if tt.size<5000: return None
    h,_=np.histogram(tt,bins=edges); h=h.astype(float); y=h[FW]; w=1/np.sqrt(np.maximum(y,1))
    B0=np.median(h[neg]); A0=max((y.sum()-B0*len(y))*BW,1)
    try:
        p,co=curve_fit(model,xw,y,p0=[B0,A0,0.4,0.08,0.10,0.4,1.7],sigma=1/w,absolute_sigma=True,
                       bounds=([0,0,0.02,-1,0.05,0.1,1.0],[np.inf,np.inf,0.98,2,0.6,0.9,5]),maxfev=200000)
        return p[6], float(np.sqrt(co[6,6])), p[3], p[4], int(tt.size)
    except Exception: return None

for HW in (25.0, 15.0):
    core=phys&(np.abs(xyz[:,0]-c[0])<HW)&(np.abs(xyz[:,1]-c[1])<HW)&(np.abs(xyz[:,2]-c[2])<HW)
    g=fit(t[core])
    print(f"\n=== TIGHT CORE box +-{HW:.0f}mm  ({core.sum():,} ev) ===")
    print(f"  IRF: t0={g[2]:.4f} sigma={g[3]:.4f}  tau_oPs(core)={g[0]:.4f}+-{g[1]:.4f}")
    xc=xyz[core]; tc=t[core]
    maxgrad=0; maxax=''
    for ax,nm in [(0,'x'),(1,'y'),(2,'z')]:
        a=xc[:,ax]; qs=np.quantile(a,np.linspace(0,1,6)); pos=[]; tau=[]; te=[]
        for i in range(5):
            sel=(a>=qs[i])&(a<=qs[i+1] if i==4 else a<qs[i+1])
            r=fit(tc[sel])
            if r and r[1]<0.12: pos.append(a[sel].mean()); tau.append(r[0]); te.append(r[1])
        if len(pos)>=3:
            pos=np.array(pos); tau=np.array(tau); w=1/np.array(te)**2
            # weighted linear slope
            X=np.vstack([np.ones_like(pos),pos]).T
            beta=np.linalg.solve((X*w[:,None]).T@X,(X*w[:,None]).T@tau)
            slope=beta[1]*1000  # milli-ns/mm
            print(f"  {nm}: tau over bins {[round(x,3) for x in tau]} @ pos {[round(x,0) for x in pos]} -> grad {slope:+.2f} mns/mm")
            if abs(slope)>abs(maxgrad): maxgrad=slope; maxax=nm
    print(f"  >> max |gradient| along {maxax} = {maxgrad:.2f} mns/mm")
    print(f"     human inter-chamber gradient = 6.96 mns/mm  -> instrument is {100*abs(maxgrad)/6.96:.0f}% of it")
    print(f"     over a 43mm chamber separation: instrumental dtau ~ {abs(maxgrad)*43/1000:.3f} ns (vs RV-LV +0.30 ns)")
