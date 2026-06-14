#!/usr/bin/env python3
import os as _o, sys as _s
_s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), '..'))
import config as cfg
"""Formal orientation determination (Codex-approved, NO tau).
STEP A (decisive, no CT masks): which listmode->histo convention reproduces the
PROVIDER's histoimage? Bin events under all 48 perm/sign orientations into the
histo grid and correlate (coarse 16-vox blocks, flip-sensitive) with the provided
histoimage. The histoimage is the authors' own reconstruction => whatever
reproduces it IS the correct convention. Since CT/histo/masks are all LPS-consistent
(no flip among them), this also fixes chamber identity.
STEP B (cross-check): for the top orientations, register to CT by a JOINT organ
panel (heart+liver+spleen+kidneys) and report per-organ enrichment + joint score.
"""
import os, itertools, numpy as np, nibabel as nib
B=cfg.RAW
LM=f"{cfg.RAW}/zenodo_11243763/20230606_positronium_patient_evaluated_data_Histo_Out.l"
HIST=f"{cfg.RAW}/zenodo_11243763/20230606_histoimage3d_lm.nii.gz"
CH=f"{cfg.MASKS}/totalseg_chambers_11243763"; TS=f"{cfg.MASKS}/totalseg_11243763"
H=np.asanyarray(nib.load(HIST).dataobj).astype(np.float32)
hox,hoy,hoz=418.00427246,419.10751343,-529.73498535
hsx,hsy,hsz=-1.63282919,-1.63713944,1.64513969
n=os.path.getsize(LM)//32
mm=np.memmap(LM,dtype=np.float64,mode="r",shape=(n,4))
offs=np.linspace(0,n-200000,60).astype(np.int64)
s=np.concatenate([np.array(mm[o:o+200000]) for o in offs],axis=0)
cols=[s[:,0],s[:,1],s[:,2]]
BF=16; NBX,NBY,NBZ=512//BF,512//BF,640//BF
def coarse_img(v):
    return v[:512,:512,:640].reshape(NBX,BF,NBY,BF,NBZ,BF).sum(axis=(1,3,5)).reshape(-1).astype(np.float64)
Hc=coarse_img(H)
def corr(a,b):
    a=a-a.mean(); b=b-b.mean(); return float((a*b).sum()/np.sqrt((a*a).sum()*(b*b).sum()+1e-9))
print("=== STEP A: histoimage reproduction (coarse 16-vox blocks, no CT masks) ===")
res=[]
for perm in itertools.permutations(range(3)):
  for sgn in itertools.product([1,-1],repeat=3):
    wx=sgn[0]*cols[perm[0]]; wy=sgn[1]*cols[perm[1]]; wz=sgn[2]*cols[perm[2]]
    i=np.round((wx-hox)/hsx); j=np.round((wy-hoy)/hsy); k=np.round((wz-hoz)/hsz)
    ib=(i>=0)&(i<512)&(j>=0)&(j<512)&(k>=0)&(k<640)
    cb=((i[ib]//BF)*NBY*NBZ+(j[ib]//BF)*NBZ+(k[ib]//BF)).astype(np.int64)
    cnt=np.bincount(cb,minlength=NBX*NBY*NBZ).astype(np.float64)
    res.append((perm,sgn,float(ib.mean()),corr(cnt,Hc)))
res.sort(key=lambda r:-r[3])
print(f"{'perm':12s}{'signs':14s}{'in_bounds':>10s}{'histo_corr':>12s}")
for perm,sgn,ibm,cc in res[:8]:
    tag=" <-- MINE" if (perm==(0,2,1) and sgn==(-1,1,-1)) else ("  (OFFICIAL)" if (perm==(0,2,1) and sgn==(1,-1,-1)) else "")
    print(f"  {str(perm):10s}{str(sgn):14s}{ibm:10.3f}{cc:12.3f}{tag}")
best=res[0]
print(f"\nBEST (reproduces histoimage) = perm {best[0]} signs {best[1]}  corr {best[3]:.3f}")

# ---- STEP B: joint organ panel via CT registration for top orientations ----
def m(p):
    try: return np.asanyarray(nib.load(p).dataobj)>0
    except Exception: return None
organs={"heart":m(f"{TS}/heart.nii.gz"),"liver":m(f"{TS}/liver.nii.gz"),"spleen":m(f"{TS}/spleen.nii.gz"),
        "kidney_L":m(f"{TS}/kidney_left.nii.gz"),"kidney_R":m(f"{TS}/kidney_right.nii.gz")}
organs={k:v for k,v in organs.items() if v is not None}
vf={k:v.mean() for k,v in organs.items()}
rv=m(f"{CH}/heart_ventricle_right.nii.gz"); lv=m(f"{CH}/heart_ventricle_left.nii.gz")
CSX,CSY,CSZ=-1.52343750,-1.52343750,1.64999998; COX,COY,COZ=389.238281,534.738281,-1200.94995
ss=s[:6_000_000]; c0,c1,c2=ss[:,0],ss[:,1],ss[:,2]
def enrich(wx,wy,wz,T,mask,mvf):
    i=np.round((wx+T[0]-COX)/CSX).astype(np.int64); j=np.round((wy+T[1]-COY)/CSY).astype(np.int64); k=np.round((wz+T[2]-COZ)/CSZ).astype(np.int64)
    ib=(i>=0)&(i<512)&(j>=0)&(j<512)&(k>=0)&(k<644)
    o=np.zeros(i.shape,bool); o[ib]=mask[i[ib],j[ib],k[ib]]; return o.mean()/mvf
def cand(perm,sgn):
    wx=sgn[0]*ss[:,perm[0]]; wy=sgn[1]*ss[:,perm[1]]; wz=sgn[2]*ss[:,perm[2]]
    # search T maximizing JOINT organ enrichment (sum of log-enrich over all organs)
    best=(-1e9,None)
    for Tx in range(-70,71,20):
      for Ty in range(80,211,20):
        for Tz in range(-740,-599,20):
            js=sum(np.log(max(enrich(wx,wy,wz,(Tx,Ty,Tz),mask,vf[nm]),1e-3)) for nm,mask in organs.items())
            if js>best[0]: best=(js,(Tx,Ty,Tz))
    T=best[1]
    ens={nm:enrich(wx,wy,wz,T,mask,vf[nm]) for nm,mask in organs.items()}
    return T,best[0],ens
print("\n=== STEP B: joint organ panel (optimize T over ALL organs jointly, no tau) ===")
cands=[(best[0],best[1])]
# add the transaxial flip of the best (flip x,y signs) for contrast
fp,fs=best[0],list(best[1]);
# flip the two non-axial signs (axial axis = the one mapping to k/z); approximate by flipping all-but-largest-range
cands.append((best[0],tuple(-x if a!=2 else x for a,x in enumerate(best[1]))))
# also explicitly mine and official
for o in [((0,2,1),(-1,1,-1)),((0,2,1),(1,-1,-1))]:
    if o not in cands: cands.append(o)
print(f"{'perm/signs':26s}{'jointJ':>8s}  organ enrichments")
for perm,sgn in cands:
    T,J,ens=cand(perm,sgn)
    es="  ".join(f"{k}:{v:.1f}x" for k,v in ens.items())
    tag=" <-MINE" if (perm,sgn)==((0,2,1),(-1,1,-1)) else (" <-OFFICIAL" if (perm,sgn)==((0,2,1),(1,-1,-1)) else "")
    print(f"  {str(perm)}{str(sgn):14s}{J:8.2f}  T={T}  {es}{tag}")
print("\nDecision: orientation with highest histo_corr (Step A) AND highest jointJ + all-organ enrichment (Step B) is correct.")
print("Chamber identity then follows from the (correctly-oriented) events + the labeled RV/LV masks.")
