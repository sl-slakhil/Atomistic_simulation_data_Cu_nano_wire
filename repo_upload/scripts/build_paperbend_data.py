"""Build the Paper A data files (frozen-wave test and load-ramp buckling) from the full LAMMPS run.
usage: python3 build_paperA_data.py <results_full folder> <output data folder>
   needs analyse_loadramp.py in the same folder.
Writes
  S3_frozenwave_<o>_<w>.csv : k, wavelength, MD stiffness, local-3D ratio, MD/X0, SD ratio 1+sign(c)c^2k^2, SD stiffness
  S3_frozenwave_fits.csv    : X0, c (SD), rms, Eringen c^2, excess c for every wire and stiffness
  S3_frozenwave_fitcurves_SD.csv : smooth SD curves (k = 0 ... 0.25 1/A) for every wire and stiffness
  S2_loadramp_buckling.csv, S2_loadramp_raw.csv : [100] load-ramp columns (via analyse_loadramp.py)
"""
import sys, os, csv, glob, math, subprocess
src, out = sys.argv[1], sys.argv[2]; os.makedirs(out, exist_ok=True)
C = os.path.join(src, 'csv')
fits = list(csv.DictReader(open(os.path.join(C, 'frozenwave_fits.csv'))))
F = {(r['orientation'], r['width_nm'], r['quantity']): r for r in fits}
sd = lambda c, k: 1 + math.copysign(1, c)*c*c*k*k
with open(os.path.join(out, 'S3_frozenwave_fits.csv'), 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['axis', 'width_nm', 'quantity', 'X0', 'c_SD_A', 'rms_SD', 'c2_Eringen_A2', 'rms_Eringen', 'c_excess_A', 'rms_excess', 'kGA_longwave'])
    for r in fits: w.writerow([r['orientation'], r['width_nm'], r['quantity'], r['X0'], r['c_SD_A'], r['rms_SD'], r['c2_Eringen_A2'], r['rms_Eringen'], r['c_excess_A'], r['rms_excess'], r['kGA_longwave']])
curves = {}
for fn in sorted(glob.glob(os.path.join(C, 'md_frozenwave_*.csv'))):
    o, wd = os.path.basename(fn)[:-4].split('_')[2:4]
    rows = list(csv.DictReader(open(fn))); qs = [q for q in ('EA', 'EIx', 'EIy') if q in rows[0]]
    hdr = ['k_invA', 'wavelength_A']
    for q in qs: hdr += [q, q+'_Rlocal3D', q+'_over_X0_MD', q+'_over_X0_SD', q+'_SD']
    with open(os.path.join(out, f'S3_frozenwave_{o}_{wd}.csv'), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(hdr)
        for r in rows:
            k = float(r['k_invA']); line = [k, 2*math.pi/k]
            for q in qs:
                X0, c = float(F[(o, wd, q)]['X0']), float(F[(o, wd, q)]['c_SD_A'])
                line += [float(r[q]), float(r[q+'_Rlocal3D']), float(r[q])/X0, sd(c, k), X0*sd(c, k)]
            w.writerow(line)
    for q in qs: curves[f'{o}_{wd}_{q}'] = float(F[(o, wd, q)]['c_SD_A'])
with open(os.path.join(out, 'S3_frozenwave_fitcurves_SD.csv'), 'w', newline='') as f:
    w = csv.writer(f); keys = list(curves); w.writerow(['k_invA'] + [k_+'_over_X0_SD' for k_ in keys])
    for i in range(101):
        k = 0.25*i/100; w.writerow([round(k, 5)] + [sd(curves[k_], k) for k_ in keys])
here = os.path.dirname(os.path.abspath(__file__))
subprocess.run([sys.executable, os.path.join(here, 'analyse_loadramp.py'), os.path.join(src, 's2'), out, C, '100'], check=True)
print('Paper A data written to', out)
