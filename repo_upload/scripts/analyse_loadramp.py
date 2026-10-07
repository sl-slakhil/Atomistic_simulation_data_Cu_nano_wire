"""Analyse load-ramp buckling runs (s2_loadramp.py JSON files) -> CSVs for Paper A.
P_cr is found from the stable steps by fitting the tip deflection with the amplification law
    d(P) = a/(1 - P/P_cr) + b      (same form as the flexibility method)
The first step that jumps to a large (post-buckled) deflection, if any, brackets P_cr from above.
usage: python3 analyse_loadramp.py <folder with ramp_*.json> <output data folder> [csv folder with params_*.csv] [orientations, default 100]"""
import sys, os, glob, json, numpy as np
from scipy.optimize import least_squares
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'src'))
import nlbeam as nb

# frozen-wave parameters per wire, read from params_<o>_<width>.csv (run_all.py --analyze);
# fall back to the 2.2 nm values below if that file is missing.
PAR = {('100', 2.2, 'x'): (4547.18, 5.98496, 311.06, 23.8, 134.24)}
def par(o, w, dr, csvdir):
    import csv
    fn = os.path.join(csvdir, f'params_{o}_{w}.csv')
    if os.path.exists(fn):
        p = next(csv.DictReader(open(fn))); s = dr
        return float(p['EI0_'+s]), float(p['c_'+s]), float(p['kGA_'+s]), float(p['d']), float(p['EA0'])
    return PAR[(o, w, dr)]

def fit_pcr(P, d):
    f = lambda q: (q[0]/(1-P/q[1]) + q[2] - d)/d
    best = None
    for Pc0 in (1.05, 1.2, 1.5, 2.0):
        r = least_squares(f, [d[0]*(1-P[0]/(Pc0*P.max())), Pc0*P.max(), 0.0],
                          bounds=([0, P.max()*1.0001, -np.inf], [np.inf, 10*P.max(), np.inf]))
        if best is None or r.cost < best.cost: best = r
    return best.x

src, out = sys.argv[1], sys.argv[2]; os.makedirs(out, exist_ok=True)
csvdir = sys.argv[3] if len(sys.argv) > 3 else os.path.join(os.path.dirname(os.path.abspath(src)), 'csv')
only = sys.argv[4].split(',') if len(sys.argv) > 4 else ['100']   # orientations to analyse
rows, raw = [], []
for fn in sorted(glob.glob(os.path.join(src, 'ramp_*.json'))):
    r = json.load(open(fn)); o = r['o']; dr = r.get('direction', 'x'); w = float(r['width'])
    if o not in only: continue
    EI0, c, kGA, dw, EA = par(o, w, dr, csvdir)
    P, d, L = np.array(r['P']), np.array(r['d']), r['L']
    jump = np.where(d > 0.06*L)[0]
    stable = np.arange(len(P)) if len(jump) == 0 else np.arange(jump[0])
    if len(stable) > 3 and d[stable[0]] < 0.3*d[stable[1]]: stable = stable[1:]   # first step not converged (outlier)
    a, Pc, b = fit_pcr(P[stable], d[stable])
    Pup = P[jump[0]] if len(jump) else np.nan
    Ploc = 1/(1/nb.pcr_local(L, EI0, 'CF') + 1/kGA); Psd = 1/(1/nb.pcr_sd(L, EI0, c, 'CF') + 1/kGA)
    rows.append([o, w, dr, L, L/dw, c, c/L, Pc, P[stable].max(), Pup, Ploc, Psd, Pc/Ploc, Psd/Ploc, 100*Pc/EA, a, b])
    for Pi, di in zip(P, d): raw.append([o, w, dr, round(L, 2), Pi, di, Pi/Ploc, int(di > 0.06*L)])
    print(f"{o} {w} {dr} L={L:.1f} L/d={L/dw:.2f}  Pcr_MD={Pc:.4f} (stable to {P[stable].max():.4f}, jump at {Pup:.4f})"
          f"  local={Ploc:.4f} SD={Psd:.4f}  MD/loc={Pc/Ploc:.4f} SD/loc={Psd/Ploc:.4f}")
rows.sort(key=lambda q: (q[0], q[1], q[3]))
hdr = 'axis,width_nm,direction,L_A,L_over_d,c_A,lambda,Pcr_MD,P_last_stable,P_first_buckled,Pcr_local_Timoshenko,Pcr_SD,ratio_MD_local,ratio_SD_local,eps_cr_percent,fit_a,fit_b'
np.savetxt(os.path.join(out, 'S2_loadramp_buckling.csv'), np.array(rows, dtype=object), fmt='%s', delimiter=',', header=hdr, comments='')
np.savetxt(os.path.join(out, 'S2_loadramp_raw.csv'), np.array(raw, dtype=object), fmt='%s', delimiter=',',
           header='axis,width_nm,direction,L_A,P_eV_per_A,tip_deflection_A,P_over_Pcr_local,buckled', comments='')
