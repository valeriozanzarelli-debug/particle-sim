"""Diagnosis: guide the photon ball with the EXACT field. If this gives ≈0.99
the guidance law is right and only the neighbour estimate needs work."""
import numpy as np, time
from rules import K, LAM, fresnel_double_slit
W0,D,A,L,BINS,HALF = 7.0,8.0,2.0,160.0,240,60.0
m=600
s=np.concatenate([np.linspace(c-A/2,c+A/2,m) for c in (-D/2,D/2)]); ds=A/m
inc=np.exp(-s**2/W0**2)
def field(x,z):
    r=np.sqrt((x[:,None]-s[None,:])**2+z**2)
    e=inc[None,:]*np.exp(1j*K*r)/np.sqrt(r)*(0.5*(1+z/r))
    psi=e.sum(1); dpsi=(e*(1j*K*(x[:,None]-s[None,:])/r)).sum(1)
    return psi,dpsi
def vel(x,z):
    p,dp=field(x,z); return np.imag(np.conj(p)*dp)/(K*np.abs(p)**2)
rng=np.random.default_rng(1)
# photon ball starts at a ball position inside the slits, weighted by |U_inc|^2
xs0=rng.normal(0,W0/2,200000); xs0=xs0[(np.abs(xs0-D/2)<=A/2)|(np.abs(xs0+D/2)<=A/2)][:3000]
edges=np.linspace(-HALF,HALF,BINS+1); ref=fresnel_double_slit(0.5*(edges[1:]+edges[:-1]),W0,D,A,L)
for dz in (0.5,0.1):
    t=time.time(); x=xs0.copy(); z=0.05
    while z<L:                      # RK2 midpoint
        k1=vel(x,z); k2=vel(x+0.5*dz*k1,z+0.5*dz); x=x+dz*k2; z+=dz
    h,_=np.histogram(x,edges); c=np.corrcoef(h.reshape(-1,8).sum(1),ref.reshape(-1,8).sum(1))[0,1]
    print(f'exact guidance, RK2 dz={dz}: {x.size} clicks, corr {c:.3f} [{time.time()-t:.0f}s]',flush=True)
