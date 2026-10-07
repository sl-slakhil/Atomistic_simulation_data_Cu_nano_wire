"""
"""
import sys, os, csv, numpy as np
from scipy.optimize import brentq
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(here, 'src')); import nlbeam as nb
D, C = sys.argv[1], sys.argv[2]
R = list(csv.DictReader(open(os.path.join(D, 'S2_loadramp_buckling.csv'))))
W = ['1.8', '2.2', '2.9']; tag = {w: w.replace('.', 'p') for w in W}
P = {w: next(csv.DictReader(open(os.path.join(C, f'params_100_{w}.csv')))) for w in W}
slope = {}
for w in W:
    e = np.array([float(r['eps_cr_percent']) for r in R if r['width_nm'] == w])
    y = np.array([float(r['Pcr_MD'])/float(r['Pcr_SD']) for r in R if r['width_nm'] == w])
    p = np.polyfit(e, y, 1); slope[w] = (p[0]/p[1], p[1])
psd = lambda L, EI, k, c: 1/(1/nb.pcr_sd(L, EI, c) + 1/k)
ploc = lambda L, EI, k: 1/(1/nb.pcr_local(L, EI) + 1/k)
def invert(L, EI, k, Pt):
    if ploc(L, EI, k) >= Pt: return np.nan          # MD below the local load: no positive c
    return brentq(lambda c: psd(L, EI, k, c) - Pt, 0.02, 40)
rows = []
for r in R:
    w = r['width_nm']; p = P[w]; EI, k, d, cfw = float(p['EI0_x']), float(p['kGA_x']), float(p['d']), float(p['c_x'])
    L, Pm, e = float(r['L_A']), float(r['Pcr_MD']), float(r['eps_cr_percent'])
    Pc = Pm/(1 + slope[w][0]*e)
    rows.append(dict(width_nm=float(w), L_A=L, L_over_d=L/d, eps_cr_percent=e, Pcr_MD=Pm, Pcr_MD_zero_strain=Pc,
                     Pcr_local=ploc(L, EI, k), Pcr_SD=psd(L, EI, k, cfw), c_frozenwave_A=cfw,
                     c_from_raw_MD_A=invert(L, EI, k, Pm), c_from_zero_strain_MD_A=invert(L, EI, k, Pc),
                     ratio_MD_local=Pm/ploc(L, EI, k), ratio_MDzero_local=Pc/ploc(L, EI, k), ratio_SD_local=psd(L, EI, k, cfw)/ploc(L, EI, k),
                     zero_strain_intercept=slope[w][1]))
fmt = lambda v: '' if (isinstance(v, float) and np.isnan(v)) else (f'{v:.6g}' if isinstance(v, float) else str(v))
keys = list(rows[0])
with open(os.path.join(D, 'S2_c_from_buckling.csv'), 'w', newline='') as f:
    wr = csv.writer(f); wr.writerow(keys); [wr.writerow([fmt(r[k]) for k in keys]) for r in rows]
# wide Veusz file: L_over_d_<w>, c_raw_<w>, c_zero_<w>, ratio_MD_local_<w>, ratio_SD_local_<w>, plus frozen-wave c lines
cols = {}
for w in W:
    rr = [r for r in rows if r['width_nm'] == float(w)]; t = tag[w]
    cols[f'L_over_d_{t}'] = [r['L_over_d'] for r in rr]
    cols[f'c_raw_{t}'] = [r['c_from_raw_MD_A'] for r in rr]
    cols[f'c_zero_{t}'] = [r['c_from_zero_strain_MD_A'] for r in rr]
    cols[f'ratio_MD_local_{t}'] = [r['ratio_MD_local'] for r in rr]
    cols[f'ratio_MDzero_local_{t}'] = [r['ratio_MDzero_local'] for r in rr]
    cols[f'ratio_SD_local_{t}'] = [r['ratio_SD_local'] for r in rr]
    cols[f'cFW_line_x_{t}'] = [3.5, 8.5]; cols[f'cFW_line_y_{t}'] = [float(P[w]['c_x'])]*2
n = max(len(v) for v in cols.values())
with open(os.path.join(D, 'S2_c_from_buckling_veusz.csv'), 'w', newline='') as f:
    wr = csv.writer(f); wr.writerow(list(cols))
    for i in range(n): wr.writerow([fmt(v[i]) if i < len(v) else '' for v in cols.values()])
# inversion curves P_SD(c)/P_local vs c for every column
cg = np.linspace(0, 14, 141); sc = {'c_A': list(cg)}
for r in rows:
    w = str(r['width_nm']); p = P[w]; EI, k = float(p['EI0_x']), float(p['kGA_x']); L = r['L_A']; t = f"{tag[w]}_Ld{r['L_over_d']:.1f}".replace('.', 'p')
    sc[f'SDratio_{t}'] = [1.0] + [psd(L, EI, k, c)/ploc(L, EI, k) for c in cg[1:]]
    sc[f'MDraw_{t}'] = [r['ratio_MD_local']]*len(cg); sc[f'MDzero_{t}'] = [r['ratio_MDzero_local']]*len(cg)
# vertical lines at the frozen-wave c (two points each, for Veusz xy plots)
sc['vline_y'] = [0.99, 1.10] + [np.nan]*(len(cg)-2)
for w in W: sc[f'vline_cFW_{tag[w]}'] = [float(P[w]['c_x'])]*2 + [np.nan]*(len(cg)-2)
with open(os.path.join(D, 'S2_c_scan_veusz.csv'), 'w', newline='') as f:
    wr = csv.writer(f); wr.writerow(list(sc))
    for i in range(len(cg)): wr.writerow([fmt(float(v[i])) for v in sc.values()])  # blank = missing
for r in rows: print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items() if k in ('width_nm', 'L_over_d', 'c_frozenwave_A', 'c_from_raw_MD_A', 'c_from_zero_strain_MD_A')})
