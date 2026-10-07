"""
"""
import numpy as np
from scipy.optimize import brentq

def _basis(x,a,b,L):
    # w-basis: exp(-a x), exp(-a(L-x)), cos(b x), sin(b x) and derivatives 0..3
    e1=np.exp(-a*x); e2=np.exp(-a*(L-x)); c=np.cos(b*x); s=np.sin(b*x)
    W=np.array([[e1,e2,c,s],
                [-a*e1,a*e2,-b*s,b*c],
                [a*a*e1,a*a*e2,-b*b*c,-b*b*s],
                [-a**3*e1,a**3*e2,b**3*s,-b**3*c]])
    return W  # W[d,j] = d-th derivative of basis j

def _ab(p,c):
    s=np.sqrt(1+4*c*c*p)
    a=np.sqrt((1+s)/(2*c*c)); b=np.sqrt((s-1)/(2*c*c))
    return a,b

def det_sd(p,L,c,bc):
    a,b=_ab(p,c)
    W0=_basis(0.0,a,b,L); WL=_basis(L,a,b,L)
    if bc=='CF':
        # unknowns C1..C4, delta ; w = delta + sum Cj phi_j
        M=np.zeros((5,5))
        M[0,:4]=W0[0]; M[0,4]=1          # w(0)=0
        M[1,:4]=W0[1]                     # w'(0)=0
        M[2,:4]=WL[0]                     # w(L)=delta
        M[3,:4]=W0[3]-W0[2]/c             # kappa'(0)-kappa(0)/c=0
        M[4,:4]=WL[3]+WL[2]/c             # kappa'(L)+kappa(L)/c=0
    elif bc=='CC':
        # w = (M0 + V x)/P + sum Cj phi_j ; unknowns C1..C4, m0=M0/P, v=V/P
        M=np.zeros((6,6))
        M[0,:4]=W0[0]; M[0,4]=1
        M[1,:4]=W0[1]; M[1,5]=1
        M[2,:4]=WL[0]; M[2,4]=1; M[2,5]=L
        M[3,:4]=WL[1]; M[3,5]=1
        M[4,:4]=W0[3]-W0[2]/c
        M[5,:4]=WL[3]+WL[2]/c
    elif bc=='SS':
        # pinned-pinned: w(0)=w(L)=0, M(0)=M(L)=0 -> w = sum Cj phi_j (M=-Pw)
        M=np.zeros((4,4))
        M[0]=W0[0]; M[1]=WL[0]
        M[2]=W0[3]-W0[2]/c; M[3]=WL[3]+WL[2]/c
    return np.linalg.det(M)

def pcr_sd(L,EI,c,bc='CF'):
    if c<=1e-9*L: return pcr_local(L,EI,bc)
    P0=pcr_local(L,EI,bc)
    ps=np.linspace(0.5,40,4000)*P0/EI
    d=np.array([det_sd(p,L,c,bc) for p in ps])
    for i in range(len(ps)-1):
        if np.sign(d[i])!=np.sign(d[i+1]):
            return brentq(det_sd,ps[i],ps[i+1],args=(L,c,bc),xtol=1e-14)*EI
    return np.nan

def pcr_local(L,EI,bc='CF'):
    K={'CF':2.0,'CC':0.5,'SS':1.0}[bc]
    return np.pi**2*EI/(K*L)**2

def pcr_eringen(L,EI,c,bc='CF'):
    P=pcr_local(L,EI,bc); return P/(1+c*c*P/EI)

if __name__=="__main__":
    lam=[0.0001,0.001,0.01,0.1,0.25,0.5]
    paper=[0.250025,0.250251,0.252580,0.282722,0.353116,0.491844]
    print("CF check vs Barretta et al. table (alpha = P L^2/(pi^2 EI))")
    for l_,pp in zip(lam,paper):
        a=pcr_sd(1.0,1.0,l_,'CF')/np.pi**2
        print(f"  lambda={l_:7.4f}  alpha={a:.6f}  paper={pp:.6f}")
    for bc in ('CC','SS'):
        print(bc,[round(pcr_sd(1.0,1.0,l_,bc)/np.pi**2,4) for l_ in (0.01,0.1,0.25)])

# ======================================================================= vibration
def _cubic_roots_real(c, b4):
    s = np.roots([c*c, -1.0, 0.0, b4])
    if np.max(np.abs(s.imag)) > 1e-9*np.max(np.abs(s)): return None
    return np.sort(s.real)

def _basis6(x, roots, L):
    """columns: for s>0: exp(-r x), exp(-r(L-x)); for s<0: cos(q x), sin(q x); rows: derivative 0..5"""
    cols = []
    for s in roots:
        if s > 0:
            r = np.sqrt(s); e1 = np.exp(-r*x); e2 = np.exp(-r*(L-x))
            cols.append([e1*(-r)**d for d in range(6)]); cols.append([e2*r**d for d in range(6)])
        else:
            q = np.sqrt(-s); c_, s_ = np.cos(q*x), np.sin(q*x)
            cc = [c_, -q*s_, -q*q*c_, q**3*s_, q**4*c_, -q**5*s_]
            ss = [s_, q*c_, -q*q*s_, -q**3*c_, q**4*s_, q**5*c_]
            cols.append(cc); cols.append(ss)
    return np.array(cols).T            # (6 derivatives, 6 functions)

def det_sd_vib(b4, L, c, bc):
    roots = _cubic_roots_real(c, b4)
    if roots is None: return np.nan
    W0 = _basis6(0.0, roots, L); WL = _basis6(L, roots, L)
    M = lambda W: W[2]-c*c*W[4]       # moment / EI  (kappa - c^2 kappa'')
    V = lambda W: W[3]-c*c*W[5]       # shear  / EI
    cb0 = W0[3]-W0[2]/c; cbL = WL[3]+WL[2]/c     # constitutive BCs on kappa = w''
    if bc == 'CF':   rows = [W0[0], W0[1], cb0, M(WL), V(WL), cbL]
    elif bc == 'CC': rows = [W0[0], W0[1], WL[0], WL[1], cb0, cbL]
    elif bc == 'SS': rows = [W0[0], M(W0), WL[0], M(WL), cb0, cbL]
    return np.linalg.det(np.array(rows))

def beta_local(bc, mode=1):
    t = {'CF': [1.875104069, 4.694091133, 7.854757438], 'CC': [4.730040745, 7.853204624, 10.99560784],
         'SS': [np.pi, 2*np.pi, 3*np.pi]}
    return t[bc][mode-1]

def omega_sd_beam(L, EI, rhoA, c, bc='CF', mode=1):
    """natural circular frequency of a stress-driven Euler-Bernoulli beam (units consistent with EI, rhoA)"""
    bl = beta_local(bc, mode)/L
    if c <= 1e-9*L: return bl**2*np.sqrt(EI/rhoA)
    b1 = beta_local(bc, 1)/L
    bs = np.linspace(0.6*b1, 6.0*bl, 3000*mode)
    d = np.array([det_sd_vib(b**4, L, c, bc) for b in bs])
    found = 0
    for i in range(len(bs)-1):
        if np.isfinite(d[i]) and np.isfinite(d[i+1]) and np.sign(d[i]) != np.sign(d[i+1]):
            found += 1
            if found == mode:
                b = brentq(lambda bb: det_sd_vib(bb**4, L, c, bc), bs[i], bs[i+1], xtol=1e-14)
                return b**2*np.sqrt(EI/rhoA)
    return np.nan

def omega_sd_rod(L, EA, rhoA, c, bc='CF', mode=1):
    """axial natural frequency of a stress-driven rod; CF = fixed-free, CC = fixed-fixed"""
    gl = ((2*mode-1)*np.pi/2 if bc == 'CF' else mode*np.pi)/L
    if c <= 1e-9*L: return gl*np.sqrt(EA/rhoA)
    def det(g):
        s1 = (1+np.sqrt(1+4*c*c*g*g))/(2*c*c); s2 = (1-np.sqrt(1+4*c*c*g*g))/(2*c*c)
        r = np.sqrt(s1); q = np.sqrt(-s2)
        def B(x):
            e1, e2 = np.exp(-r*x), np.exp(-r*(L-x)); cc, ss = np.cos(q*x), np.sin(q*x)
            return np.array([[e1, e2, cc, ss], [-r*e1, r*e2, -q*ss, q*cc], [r*r*e1, r*r*e2, -q*q*cc, -q*q*ss],
                             [-r**3*e1, r**3*e2, q**3*ss, -q**3*cc]])
        B0, BL = B(0.0), B(L)
        cb0 = B0[2]-B0[1]/c; cbL = BL[2]+BL[1]/c          # CBCs on strain u'
        rows = [B0[0], cb0, cbL, (BL[1]-c*c*BL[3]) if bc == 'CF' else BL[0]]
        return np.linalg.det(np.array(rows))
    gs = np.linspace(0.6, 6.0, 3000)*gl; d = [det(g) for g in gs]
    for i in range(len(gs)-1):
        if np.sign(d[i]) != np.sign(d[i+1]):
            return brentq(det, gs[i], gs[i+1], xtol=1e-14)*np.sqrt(EA/rhoA)
    return np.nan

def omega_disp(k, EI0, rhoA, c, model='SD', kGA=None, rhoI=None):
    """flexural dispersion: Euler-Bernoulli, or Timoshenko (kGA, rotary inertia rhoI) with
       wavenumber-dependent bending stiffness EI(k) = EI0(1+c^2k^2) [SD] or EI0/(1+c^2k^2) [Eringen]"""
    k = np.asarray(k, float)
    EIk = EI0*(1+c*c*k*k) if model == 'SD' else (EI0/(1+c*c*k*k) if model == 'ER' else EI0+0*k)
    if kGA is None: return np.sqrt(EIk*k**4/rhoA)
    out = []
    for kk, ei in zip(np.atleast_1d(k), np.atleast_1d(EIk)):
        # det [[kGA k^2 - rhoA w^2, i kGA k],[-i kGA k, ei k^2 + kGA - rhoI w^2]] = 0 -> quadratic in w^2
        a = rhoA*rhoI; b = -(rhoA*(ei*kk*kk+kGA)+rhoI*kGA*kk*kk); cc = kGA*kk*kk*(ei*kk*kk+kGA)-(kGA*kk)**2
        w2 = (-b-np.sqrt(b*b-4*a*cc))/(2*a); out.append(np.sqrt(w2))
    return np.array(out)

def omega_disp_rod(k, EA0, rhoA, c, model='SD', lateral_inertia=0.0):
    """longitudinal dispersion; lateral_inertia = nu_x^2 rhoI_x + nu_y^2 rhoI_y adds the Rayleigh-Love
       (Poisson lateral inertia) correction  rhoA -> rhoA + lateral_inertia*k^2"""
    k = np.asarray(k, float)
    EAk = EA0*(1+c*c*k*k) if model == 'SD' else (EA0/(1+c*c*k*k) if model == 'ER' else EA0+0*k)
    return np.sqrt(EAk*k*k/(rhoA+lateral_inertia*k*k))

def omega_sd_beam_integral(L, EI, rhoA, c, n=600):
    """independent check (CF only): integral SD kernel, collocation; w = (rhoA w^2/EI) T Phi S w"""
    x = np.linspace(0, L, n+1); h = L/n; wq = np.full(n+1, h); wq[[0, -1]] = h/2
    Phi = np.exp(-np.abs(x[:, None]-x[None, :])/c)/(2*c)*wq[None, :]
    S = np.array([[(xi-xx)*wq[j] if xi >= xx else 0.0 for j, xi in enumerate(x)] for xx in x])   # M(x)=int_x^L (xi-x) q
    T = np.array([[(xx-xi)*wq[j] if xi <= xx else 0.0 for j, xi in enumerate(x)] for xx in x])   # w(x)=int_0^x (x-xi) kappa
    G = T@Phi@S
    ev = np.linalg.eigvals(G); ev = ev[np.abs(ev.imag) < 1e-9*np.abs(ev).max()].real
    b4 = 1/ev[ev > 0].max()
    return np.sqrt(b4)*np.sqrt(EI/rhoA)

def omega_sd_rod_integral(L, EA, rhoA, c, n=600):
    """independent check (fixed-free rod): integral SD kernel, collocation"""
    x = np.linspace(0, L, n+1); h = L/n; wq = np.full(n+1, h); wq[[0, -1]] = h/2
    Phi = np.exp(-np.abs(x[:, None]-x[None, :])/c)/(2*c)*wq[None, :]
    S = np.array([[wq[j] if xi >= xx else 0.0 for j, xi in enumerate(x)] for xx in x])    # N(x) = int_x^L q
    T = np.array([[wq[j] if xi <= xx else 0.0 for j, xi in enumerate(x)] for xx in x])    # u(x) = int_0^x eps
    ev = np.linalg.eigvals(T@Phi@S); ev = ev[np.abs(ev.imag) < 1e-9*np.abs(ev).max()].real
    return np.sqrt(1/ev[ev > 0].max())*np.sqrt(EA/rhoA)
