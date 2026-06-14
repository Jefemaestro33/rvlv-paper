#!/usr/bin/env python3
import os as _o, sys as _s
_s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), '..'))
import config as cfg
"""Capstone: does oxygenation survive AFTER controlling for position & density?
OLS tau ~ SO2 + z + HU on reliable blood compartments + partial correlations.
Plus a ladder plot. Decides whether any oxygenation signal hides under confounds.
"""
import os, csv, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
OUT=f"{cfg.OUT}"
rows=list(csv.DictReader(open(os.path.join(OUT,"oxygenation_ladder.csv"))))
def f(r,k):
    try: return float(r[k])
    except: return np.nan
# reliable blood compartments: good stats, exclude muscle & noisy fits
blood=[r for r in rows if r["cls"] in ("venous","arterial")
       and f(r,"n_core")>100000 and f(r,"tau_core_err")<0.6]
print("reliable blood compartments (n_core>100k, err<0.6):")
for r in blood: print(f"  {r['name']:14s} {r['cls']:8s} SO2={f(r,'so2'):.2f} tau={f(r,'tau_core'):.3f}+-{f(r,'tau_core_err'):.3f} z={f(r,'z_mm'):.0f} HU={f(r,'HU'):.0f} n={int(f(r,'n_core'))}")
tau=np.array([f(r,"tau_core") for r in blood])
so2=np.array([f(r,"so2") for r in blood]); z=np.array([f(r,"z_mm") for r in blood])
hu=np.array([f(r,"HU") for r in blood]); rr=np.array([f(r,"r_mm") for r in blood])
err=np.array([f(r,"tau_core_err") for r in blood])

def pcorr(a,b):
    a=a-a.mean(); b=b-b.mean(); return (a*b).sum()/np.sqrt((a*a).sum()*(b*b).sum())
def partial(x,y,covs):
    # residualize x and y on covs (+intercept), correlate residuals
    C=np.column_stack([np.ones_like(x)]+covs)
    bx=np.linalg.lstsq(C,x,rcond=None)[0]; by=np.linalg.lstsq(C,y,rcond=None)[0]
    return pcorr(x-C@bx, y-C@by)

print(f"\nn={len(blood)} compartments")
print(f"raw corr(tau, SO2)              = {pcorr(tau,so2):+.3f}")
print(f"partial corr(tau,SO2 | z)       = {partial(so2,tau,[z]):+.3f}")
print(f"partial corr(tau,SO2 | z,HU)    = {partial(so2,tau,[z,hu]):+.3f}")
print(f"partial corr(tau,SO2 | z,HU,r)  = {partial(so2,tau,[z,hu,rr]):+.3f}")
print(f"(reference) corr(tau, z)        = {pcorr(tau,z):+.3f}")
print(f"(reference) corr(tau, HU)       = {pcorr(tau,hu):+.3f}")

# weighted OLS tau ~ SO2 + z + HU
X=np.column_stack([np.ones_like(tau),(so2-so2.mean())/so2.std(),(z-z.mean())/z.std(),(hu-hu.mean())/hu.std()])
W=1/err**2
beta=np.linalg.solve((X*W[:,None]).T@X,(X*W[:,None]).T@(tau*W)/1) if False else np.linalg.lstsq(X*np.sqrt(W)[:,None],tau*np.sqrt(W),rcond=None)[0]
resid=tau-X@beta; dof=max(len(tau)-X.shape[1],1)
s2=np.sum(W*resid**2)/dof
cov=s2*np.linalg.inv((X*W[:,None]).T@X)
se=np.sqrt(np.diag(cov))
print("\nweighted OLS tau ~ 1 + SO2_z + z_z + HU_z (standardized predictors, ns per 1sd):")
for nm,b,s in zip(["intercept","SO2","axial_z","HU"],beta,se):
    print(f"  {nm:10s} {b:+.3f} +- {s:.3f}   {'<-- oxygenation' if nm=='SO2' else ''}")
print(f"  => SO2 coefficient is {abs(beta[1]/se[1]):.2f} sigma "
      f"({'NOT significant' if abs(beta[1]/se[1])<2 else 'significant'}); "
      f"sign {'NEGATIVE(oxy-consistent)' if beta[1]<0 else 'POSITIVE(oxy-INCONSISTENT)'}")

# plot ladder
allb=[r for r in rows if r["cls"] in ("venous","arterial") and f(r,"n_core")>50000 and f(r,"tau_core_err")<0.6]
allb=sorted(allb,key=lambda r:-f(r,"tau_core"))
fig,ax=plt.subplots(1,2,figsize=(13,5))
cols={"venous":"#c0392b","arterial":"#2471a3"}
for i,r in enumerate(allb):
    ax[0].errorbar(f(r,"tau_core"),i,xerr=f(r,"tau_core_err"),fmt="o",color=cols[r["cls"]],ms=7)
    lab=r["name"]+(" *" if r["name"] in("pulm_ARTERY","pulm_VEIN") else "")
    ax[0].text(f(r,"tau_core")+f(r,"tau_core_err")+0.03,i,lab,va="center",fontsize=9,
               fontweight="bold" if "*" in lab else "normal")
ax[0].set_yticks([]); ax[0].set_xlabel("tau_oPs core (ns)")
ax[0].set_title("Oxygenation ladder (red=deoxy/venous, blue=oxy/arterial)\n* = pulmonary crossover controls")
ax[0].invert_yaxis()
import matplotlib.patches as mp
ax[0].legend(handles=[mp.Patch(color="#c0392b",label="deoxygenated"),mp.Patch(color="#2471a3",label="oxygenated")],loc="lower right")
# tau vs SO2 (jittered) and vs z
ax[1].scatter(so2,tau,c=[cols[r["cls"]] for r in blood],s=60)
for r,x,y in zip(blood,so2,tau): ax[1].annotate(r["name"],(x,y),fontsize=7,xytext=(3,3),textcoords="offset points")
ax[1].set_xlabel("nominal SO2 (oxygenation)"); ax[1].set_ylabel("tau_oPs core (ns)")
ax[1].set_title(f"tau vs oxygenation: r={pcorr(tau,so2):+.2f} (flat = NOT oxygenation)")
fig.tight_layout(); fig.savefig(os.path.join(OUT,"plots","fig5_oxygenation_ladder.png"),dpi=130)
print("\nWROTE plots/fig5_oxygenation_ladder.png")
