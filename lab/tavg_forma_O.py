import sys,os; sys.path.insert(0,"/home/claude/particle-sim/lab")
import numpy as np, tavg_radial as tr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
a=tr.Atom(8,tr.config_from_string("1s2 2s2 2p4")); a.scf(); r=a.g.r
print("E",a.E)
rng=np.random.default_rng(1)
def rad(i,n):
    w=a.P[i]**2*r; c=np.cumsum(w); c/=c[-1]   # P=rR, measure P^2 dr, dr=r dx
    return np.interp(rng.random(n),c,r)
pts=[];col=[]
N=12000
cfg=a.config
for i,(n,l,q) in enumerate(cfg):
    if l==0:
        m=N*q; rr=rad(i,m); v=rng.normal(size=(m,3)); v/=np.linalg.norm(v,axis=1)[:,None]
        pts.append(rr[:,None]*v); col+= [0 if n==1 else 1]*m
    else:  # 2p4 of the 3P term: pz doubly, px and py singly (real orbitals)
        for ax,occ in ((2,2),(0,1),(1,1)):
            m=N*occ; rr=rad(i,m*4); v=rng.normal(size=(m*4,3)); v/=np.linalg.norm(v,axis=1)[:,None]
            keep=rng.random(m*4)<v[:,ax]**2; p=(rr[:,None]*v)[keep][:m]
            pts.append(p); col+=[2]*len(p)
P=np.vstack(pts); col=np.array(col)
# thin slab view from slightly above
th=0.0; R=np.array([[1,0,0],[0,np.cos(th),-np.sin(th)],[0,np.sin(th),np.cos(th)]]); Q=P@R.T
sl=np.abs(Q[:,1])<0.45; Q=Q[sl]; col=col[sl]
fig=plt.figure(figsize=(9,9),facecolor="black"); ax=fig.add_axes([0,0,1,1]); ax.set_facecolor("black")
cols={0:(1.0,0.85,0.5),1:(0.45,0.75,1.0),2:(0.75,0.5,1.0)}
for k,s,al in ((2,1,0.10),(1,1,0.07),(0,1,0.02)):
    m=col==k; c=cols[k]
    for size,alpha in ((60,al*0.15),(14,al*0.5),(2,al*2)):
        ax.scatter(Q[m,0],Q[m,2],s=size,color=c,alpha=min(alpha,1),linewidths=0)
ax.scatter([0],[0],s=30,color="white")
ax.set_xlim(-2.6,2.6); ax.set_ylim(-2.6,2.6); ax.set_aspect("equal"); ax.axis("off")
ax.text(0.03,0.97,"Ossigeno (³P) dalla regola: Hartree–Fock del livello fine",color="white",fontsize=14,transform=ax.transAxes,va="top")
ax.text(0.03,0.93,"oro: 1s²   azzurro: 2s²   viola: 2p⁴ (pz doppio, px e py singoli)",color="#bbbbbb",fontsize=11,transform=ax.transAxes,va="top")
ax.text(0.03,0.04,"APPROSSIMATO, NON ANCORA ACCESO (controllo 1s aperto) · E_c CCSD(T) −256,2 mEh (−0,65 %) · IP 13,601 eV",color="#ffb347",fontsize=11,transform=ax.transAxes)
ax.text(0.97,0.08,"fetta sottile, vista di lato · 1 a₀ ≈ 0,53 Å",color="#888888",fontsize=10,transform=ax.transAxes,ha="right")
fig.savefig(sys.argv[1],dpi=110,facecolor="black")
