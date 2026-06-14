#!/usr/bin/env python3
import os as _o, sys as _s
_s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), '..'))
import config as cfg
"""Inspect the Bern 82Rb reference-material listmode (Zenodo 12636019,
Rb82_coins_Histo_Out.l). Same isotope as the human (11243763) -> use it to
characterize the instrument IRF / prompt-offset and the positional dependence of
tau in a KNOWN material. This script: verify format, column stats, t-spectrum,
and locate the disk source(s) via spatial occupancy.
"""
import os, numpy as np
LM = f"{cfg.RAW}/zenodo_12636019/Rb82_coins_Histo_Out.l"
sz = os.path.getsize(LM); n = sz // 32
print(f"size_bytes {sz}  size_mod_32 {sz%32}  event_count {n:,}")
mm = np.memmap(LM, dtype=np.float64, mode="r", shape=(n, 4))
print("first 5 rows (x,y,z,t):"); print(np.array(mm[:5]))

# sample ~6M events evenly
nblk, blk = 30, 200_000
offs = np.linspace(0, max(n-blk,0), nblk).astype(np.int64)
samp = np.concatenate([np.array(mm[o:o+blk]) for o in offs], axis=0)
for i,c in enumerate(["x","y","z","t"]):
    v = samp[:,i]
    print(f"  {c}: min {v.min():.2f} max {v.max():.2f} mean {v.mean():.2f} "
          f"p1 {np.percentile(v,1):.2f} p50 {np.percentile(v,50):.2f} p99 {np.percentile(v,99):.2f}")

# t-spectrum (coarse) over physical window
t = samp[:,3]; tl,th = np.percentile(t,[0.01,99.99])
cb = np.linspace(max(tl,-15), min(th,15), 41); ch,_ = np.histogram(t,bins=cb); mx=ch.max()
print(f"t-spectrum [{cb[0]:.2f},{cb[-1]:.2f}] (prompt peak locates t0):")
for i in range(40):
    print(f"  [{cb[i]:6.2f},{cb[i+1]:6.2f}) {ch[i]:8d} {'#'*int(50*ch[i]/mx)}")

# spatial occupancy to find disks: 2D histogram in x-z (axial) and x-y
def occ(a,b,na,nb,nm):
    H,xe,ye = np.histogram2d(a,b,bins=[na,nb])
    print(f"\n{nm} occupancy ({na}x{nb}), '#'=dense:")
    mxh=H.max()
    for r in range(na):
        row="".join("#" if H[r,c]>0.15*mxh else ("." if H[r,c]>0.02*mxh else " ") for c in range(nb))
        print(f"  {xe[r]:7.0f} |{row}|")
occ(samp[:,0],samp[:,2],24,24,"x(rows) vs z(cols)")
occ(samp[:,0],samp[:,1],24,24,"x(rows) vs y(cols)")
print("\nz spans the long body axis; clusters = disk source positions.")
