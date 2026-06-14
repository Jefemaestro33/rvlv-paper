#!/usr/bin/env python3
import os as _o, sys as _s
_s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), '..'))
import config as cfg
"""Is the RV-LV contrast from BLOOD or from the chamber WALL (partial volume)?
Erode each chamber mask to its geometric CORE (pure blood pool, away from walls)
and compare to the near-wall RIM. If the contrast is real blood oxygenation it
should survive/strengthen in the eroded core; if it's wall/myocardium partial
volume it should collapse in the core and live in the rim.
"""
import os, numpy as np, nibabel as nib
from scipy.special import erfc
from scipy.optimize import curve_fit
from scipy.ndimage import binary_erosion

WRK=f"{cfg.OUT}"
RV=f"{cfg.MASKS}/totalseg_chambers_11243763/heart_ventricle_right.nii.gz"
LV=f"{cfg.MASKS}/totalseg_chambers_11243763/heart_ventricle_left.nii.gz"
HIST=f"{cfg.RAW}/zenodo_11243763/20230606_histoimage3d_lm.nii.gz"
cache=np.load(os.path.join(WRK,"heart_region_events.npy"))
Nx,Ny,Nz=nib.load(HIST).shape
rvm=np.asanyarray(nib.load(RV).dataobj)>0; lvm=np.asanyarray(nib.load(LV).dataobj)>0
ri=np.round(cache[:,0]).astype(int);rj=np.round(cache[:,1]).astype(int);rk=np.round(cache[:,2]).astype(int)
ib=(ri>=0)&(ri<Nx)&(rj>=0)&(rj<Ny)&(rk>=0)&(rk<Nz); t=cache[:,3]
def evt(mask):
    s=np.zeros(cache.shape[0],bool); s[ib]=mask[ri[ib],rj[ib],rk[ib]]; return t[s]
BW=0.1; edges=np.arange(-15,15+BW,BW); ctr=0.5*(edges[:-1]+edges[1:])
neg=(ctr>-12)&(ctr<-6); FW=(ctr>=-3)&(ctr<=12); xw=ctr[FW]
def emg(tt,t0,s,ta):
    x=tt-t0; z=(s/ta-x/s)/np.sqrt(2); return (1/(2*ta))*np.exp((s**2)/(2*ta**2)-x/ta)*erfc(z)
def model(tt,B,A,f1,t0,s,t1,t2): return B+A*(f1*emg(tt,t0,s,t1)+(1-f1)*emg(tt,t0,s,t2))
def fit(tarr):
    if tarr.size<4000: return np.nan,np.nan,tarr.size
    h,_=np.histogram(tarr,bins=edges);h=h.astype(float);y=h[FW];w=1/np.sqrt(np.maximum(y,1))
    B0=np.median(h[neg]);A0=max((y.sum()-B0*len(y))*BW,1)
    try:
        p,c=curve_fit(model,xw,y,p0=[B0,A0,0.55,0.066,0.108,0.39,1.8],sigma=1/w,absolute_sigma=True,
                      bounds=([0,0,0.05,-1,0.05,0.1,1],[np.inf,np.inf,0.98,1,0.6,0.9,5]),maxfev=120000)
        return float(p[6]),float(np.sqrt(c[6,6])),int(tarr.size)
    except Exception: return np.nan,np.nan,int(tarr.size)

print(f"{'region':<18}{'RV tau':>16}{'LV tau':>16}{'RV-LV':>10}")
for label,erode in [("full mask",0),("core erode-1",1),("core erode-2",2),("core erode-3",3)]:
    rv = rvm if erode==0 else binary_erosion(rvm,iterations=erode)
    lv = lvm if erode==0 else binary_erosion(lvm,iterations=erode)
    tr=evt(rv); tl=evt(lv); ar,er,nr=fit(tr); al,el,nl=fit(tl)
    d=ar-al
    print(f"{label:<18}{ar:7.3f}+-{er:.3f}({nr//1000:>4}k){al:7.3f}+-{el:.3f}({nl//1000:>4}k){d:+8.3f}")
# rim = full minus erode-2 (near-wall shell)
for nm,full,core in [("RV",rvm,binary_erosion(rvm,iterations=2)),("LV",lvm,binary_erosion(lvm,iterations=2))]:
    rim=full&~core; tr=evt(rim); a,e,n=fit(tr)
    print(f"  {nm} RIM (near-wall shell): tau={a:.3f}+-{e:.3f}  n={n//1000}k")
print("\nRead: contrast strengthening in the eroded core -> blood-borne (oxygenation-consistent).")
print("      contrast collapsing in core / living in rim -> wall/myocardium partial volume artifact.")
