#!/usr/bin/env python3
import os as _o, sys as _s
_s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), '..'))
import config as cfg
"""Decisive confound test: is the registration-stable RV>LV contrast biological
(chamber oxygenation) or an instrumental position-dependent timing gradient?
RV is entirely deoxygenated blood -> tau_oPs should be UNIFORM within RV. If tau
slopes across RV along the RV->LV axis (toward the LV value), the 'contrast' is a
smooth spatial gradient, not chamber biology. Also test an orthogonal axis.
"""
import os, numpy as np, nibabel as nib
from scipy.special import erfc
from scipy.optimize import curve_fit

WRK=f"{cfg.OUT}"
RV=f"{cfg.MASKS}/totalseg_chambers_11243763/heart_ventricle_right.nii.gz"
LV=f"{cfg.MASKS}/totalseg_chambers_11243763/heart_ventricle_left.nii.gz"
HIST=f"{cfg.RAW}/zenodo_11243763/20230606_histoimage3d_lm.nii.gz"
cache=np.load(os.path.join(WRK,"heart_region_events.npy"))
Nx,Ny,Nz=nib.load(HIST).shape
rvm=np.asanyarray(nib.load(RV).dataobj)>0; lvm=np.asanyarray(nib.load(LV).dataobj)>0
vox_mm=np.array([1.63282919,1.63713944,1.64513969])

ri=np.round(cache[:,0]).astype(int);rj=np.round(cache[:,1]).astype(int);rk=np.round(cache[:,2]).astype(int)
ib=(ri>=0)&(ri<Nx)&(rj>=0)&(rj<Ny)&(rk>=0)&(rk<Nz)
pos=cache[:,:3]; t=cache[:,3]
def inmask(m):
    s=np.zeros(cache.shape[0],bool); s[ib]=m[ri[ib],rj[ib],rk[ib]]; return s
sel_rv=inmask(rvm); sel_lv=inmask(lvm)

# centroids (mm) and axes
rv_c=(np.argwhere(rvm).mean(0))*vox_mm; lv_c=(np.argwhere(lvm).mean(0))*vox_mm
u=(lv_c-rv_c); d_cen=np.linalg.norm(u); u/=d_cen
# orthogonal axis v (any vector perp to u)
tmp=np.array([0,0,1.0]); v=tmp-u*(tmp@u); v/=np.linalg.norm(v)
print(f"RV centroid(mm) {rv_c.round(1)}  LV centroid(mm) {lv_c.round(1)}  |RV-LV|={d_cen:.1f} mm")
print(f"RV->LV axis u={u.round(3)}   orthogonal v={v.round(3)}")

BW=0.1; edges=np.arange(-15,15+BW,BW); ctr=0.5*(edges[:-1]+edges[1:])
neg=(ctr>-12)&(ctr<-6); FW=(ctr>=-3)&(ctr<=12); xw=ctr[FW]
def emg(tt,t0,s,ta):
    x=tt-t0; z=(s/ta-x/s)/np.sqrt(2); return (1/(2*ta))*np.exp((s**2)/(2*ta**2)-x/ta)*erfc(z)
def model(tt,B,A,f1,t0,s,t1,t2): return B+A*(f1*emg(tt,t0,s,t1)+(1-f1)*emg(tt,t0,s,t2))
def fit(tarr):
    if tarr.size<4000: return np.nan,np.nan,tarr.size
    h,_=np.histogram(tarr,bins=edges); h=h.astype(float); y=h[FW]; w=1/np.sqrt(np.maximum(y,1))
    B0=np.median(h[neg]); A0=max((y.sum()-B0*len(y))*BW,1)
    try:
        p,c=curve_fit(model,xw,y,p0=[B0,A0,0.55,0.066,0.108,0.39,1.8],sigma=1/w,absolute_sigma=True,
                      bounds=([0,0,0.05,-1,0.05,0.1,1],[np.inf,np.inf,0.98,1,0.6,0.9,5]),maxfev=120000)
        return float(p[6]),float(np.sqrt(c[6,6])),int(tarr.size)
    except Exception: return np.nan,np.nan,int(tarr.size)

def proj(sel,axis): return (pos[sel]*vox_mm)@axis
def within(sel,axis,name,nbin=3):
    s=proj(sel,axis); tt=t[sel]; qs=np.quantile(s,np.linspace(0,1,nbin+1))
    print(f"\n  {name}: tau along {('u(RV->LV)' if axis is u else 'v(orth)')} in {nbin} bins:")
    out=[]
    for i in range(nbin):
        m=(s>=qs[i])&(s<=qs[i+1]) if i==nbin-1 else (s>=qs[i])&(s<qs[i+1])
        ta,te,n=fit(tt[m]); out.append((s[m].mean(),ta,te,n))
        print(f"    bin{i} <s>={s[m].mean():6.1f}mm  tau={ta:.3f}+-{te:.3f}  n={n}")
    # slope tau vs s
    ss=np.array([o[0] for o in out]); tta=np.array([o[1] for o in out])
    sl=np.polyfit(ss,tta,1)[0]
    print(f"    -> within-{name} slope = {sl*1000:.2f} milli-ns/mm  (tau change over chamber span {sl*(ss.max()-ss.min()):+.3f} ns)")
    return sl

# baseline chamber taus
trv,erv,_=fit(t[sel_rv]); tlv,elv,_=fit(t[sel_lv])
print(f"\nRV tau={trv:.3f}+-{erv:.3f}  LV tau={tlv:.3f}+-{elv:.3f}  RV-LV={trv-tlv:+.3f} ns")
print(f"Implied gradient if purely spatial: {(trv-tlv)/d_cen*1000:.2f} milli-ns/mm over {d_cen:.0f} mm")

print("\n=== WITHIN-CHAMBER tau gradients (RV is uniform deoxy blood; LV uniform oxy) ===")
sl_rv_u=within(sel_rv,u,"RV",3)
sl_lv_u=within(sel_lv,u,"LV",3)
print("\n  --- orthogonal-axis control ---")
sl_rv_v=within(sel_rv,v,"RV",3)

print("\n=== whole-heart blood pool: tau vs position along u (6 bins) ===")
sel_h=sel_rv|sel_lv
within(sel_h,u,"HEART",6)

print("\nINTERPRETATION:")
print(f"  RV->LV inter-chamber slope (the 'effect'): {(trv-tlv)/d_cen*1000:.2f} milli-ns/mm")
print(f"  within-RV slope along u: {sl_rv_u*1000:.2f} ; within-LV slope along u: {sl_lv_u*1000:.2f}")
print(f"  within-RV slope along orthogonal v: {sl_rv_v*1000:.2f}")
print("  If within-chamber u-slopes ~ inter-chamber slope -> spatial gradient (instrumental).")
print("  If within-chamber u-slopes ~ 0 and << inter-chamber -> chamber-specific (biological-consistent).")
