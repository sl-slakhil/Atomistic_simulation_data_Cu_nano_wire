# Nonlocal length of Cu nanowires from atomistic simulations — figure data

This folder holds the data behind two figures:

- **Fig. 1:** the frozen-wave (sinusoidal-load) test.
- **Fig. 2:** the nonlocal length recovered from buckling.

All files are plain CSV: comma separated, with one header row. A blank cell means "no value".

## What was done (short)

1. **Atomistic model**
   - Molecular statics at 0 K with LAMMPS (Thompson et al., *Comput. Phys. Commun.* 271 (2022) 108171).
   - Embedded-atom potential `Cu_u3.eam` (Foiles, Baskes, Daw, *Phys. Rev. B* 33 (1986) 7983).
   - Square [100] wires of 1.8, 2.2 and 2.9 nm width, and [110] wires of 2.2 and 2.9 nm.
   - Section properties come from the atoms: A = NΩ/ℓ and I = (Ω/ℓ)Σx², with Ω = a0³/4 and a0 = 3.615 Å.
   - d = √A.

2. **Frozen-wave test (Fig. 1)**
   - A periodic wire, exactly one wavelength long (λ = 2π/k), is loaded by a sinusoidal force on every atom.
   - The exact linear response u = K⁻¹f is computed matrix-free, by finite-difference Hessian–vector products and conjugate gradients.
   - The stiffness at that wavelength is load divided by deformation:
     - EI(k) = q0/(k³ψ̂), from the cross-section rotation ψ̂
     - EA(k) = q0/(k²û), from the axial displacement û
   - X(k)/X0 is this stiffness divided by its long-wave value.
   - The stress-driven (SD, Romano–Barretta) beam predicts X/X0 = 1 + c²k². Fitting this to the MD points gives the nonlocal length c.
   - Only points with λ ≥ 2d are used in the fit.
   - The strain-driven (Eringen) form 1/(1 + c²k²) can only fit with c² < 0, so it predicts the wrong trend.
   - A local 3D anisotropic continuum model (semi-analytical FE) shows how much of the stiffening is ordinary 3D elasticity.

3. **Buckling (Fig. 2)**
   - Clamped–free [100] columns: 7 columns, L/d = 4.1–7.8.
   - Load ramp: an axial force on the top plane plus a small lateral force. P_cr is found by fitting the tip deflection with δ = a/(1 − P/P_cr) + b.
   - The MD P_cr is compared with:
     - the local Timoshenko–Engesser load
     - the SD load: Euler–Bernoulli SD with constitutive boundary conditions, plus the same shear correction
   - Inverse problem: for each column, find the c that makes the SD load equal to
     - (i) the raw MD load
     - (ii) the MD load corrected to zero axial strain. The correction is P_MD/(1 + s·ε_cr). Here s = slope/intercept of a linear fit of P_MD/P_SD against ε_cr, done for each wire.
   - Why the correction is needed: the compressed core of a [100] wire softens, and this hides the nonlocal stiffening.

**Main result**
- Frozen wave: c ≈ 6.0, 6.0 and 6.4 Å for the 1.8, 2.2 and 2.9 nm [100] wires.
- Buckling, after the zero-strain correction (1.8 and 2.2 nm wires): c = 6.2–6.5 Å.
- 2.9 nm wire, raw load at L/d = 6.6: c = 6.6 Å.
- Raw buckling loads alone give no c or a far too small one.

**Units:** energy in eV, length in Å.
- Force: eV/Å (1 eV/Å = 1.602 nN)
- EA and κGA: eV/Å
- EI: eV·Å
- k: 1/Å

---

## fig1_frozen_wave/

### `S3_frozenwave_<axis>_<width>.csv`

One file per wire:

| File | Wire | d (Å) |
|---|---|---|
| `100_1.8` | [100], 1.8 nm | 20.2 |
| `100_2.2` | [100], 2.2 nm | 23.8 |
| `100_2.9` | [100], 2.9 nm | 31.0 |
| `110_2.2` | [110], 2.2 nm | 22.7 |
| `110_2.9` | [110], 2.9 nm | 31.5 |

One row per wavelength. Points with k > π/d (λ < 2d) are shown in the figure but not used in the fit.

| Header | Meaning |
|---|---|
| `k_invA` | wavenumber k = 2π/λ (1/Å) — x-axis of Fig. 1 |
| `wavelength_A` | wavelength λ (Å) |
| `EA` | MD axial stiffness EA(k) (eV/Å) |
| `EA_Rlocal3D` | local 3D continuum stiffness ratio EA(k)/EA0 for the same cross-section (no nonlocality) |
| `EA_over_X0_MD` | MD ratio EA(k)/EA0 — plotted points |
| `EA_over_X0_SD` | SD fit 1 + sign(c)·c²k² at the same k — fitted line |
| `EA_SD` | SD stiffness EA0·(1 + sign(c)·c²k²) (eV/Å) |
| `EIx`, `EIx_Rlocal3D`, `EIx_over_X0_MD`, `EIx_over_X0_SD`, `EIx_SD` | same as above for bending about the x axis (EI in eV·Å) |
| `EIy`, … | same for bending about the y axis ([110] files only; [100] is symmetric) |

### `S3_frozenwave_fits.csv`

One row per wire and stiffness.

| Header | Meaning |
|---|---|
| `axis` | wire axis (100 or 110) |
| `width_nm` | nominal width (nm) |
| `quantity` | `EA`, `EIx` or `EIy` |
| `X0` | long-wave stiffness EA0 (eV/Å) or EI0 (eV·Å) |
| `c_SD_A` | stress-driven nonlocal length c (Å). Negative means the fit needs softening (1 − c²k²). |
| `rms_SD` | rms relative error of the SD fit (fraction, 0.006 = 0.6%) |
| `c2_Eringen_A2` | c² of the strain-driven (Eringen) fit (Å²). Negative values are unphysical. |
| `rms_Eringen` | rms relative error of the Eringen fit (fraction) |
| `c_excess_A` | SD length fitted to MD / local-3D, i.e. after removing ordinary 3D elasticity (Å) |
| `rms_excess` | rms relative error of that fit (fraction) |
| `kGA_longwave` | shear stiffness κGA at the longest wavelength (eV/Å), EI rows only |

### `S3_frozenwave_fitcurves_SD.csv`

Smooth SD curves for plotting the lines.

| Header | Meaning |
|---|---|
| `k_invA` | k from 0 to 0.25 1/Å (101 points) |
| `<axis>_<width>_<quantity>_over_X0_SD` | 1 + sign(c)·c²k² for that wire and stiffness, e.g. `100_2.2_EIx_over_X0_SD` |

---

## fig2_c_from_buckling/

### `S2_loadramp_buckling.csv` (source buckling results)

One row per MD column.

| Header | Meaning |
|---|---|
| `axis`, `width_nm`, `direction` | wire axis, width (nm), bending direction |
| `L_A` | free length L from the last clamped plane to the top plane (Å) |
| `L_over_d` | slenderness L/d |
| `c_A` | frozen-wave c of that wire (Å) |
| `lambda` | nonlocal ratio λ = c/L |
| `Pcr_MD` | MD critical load from the amplification fit (eV/Å) |
| `P_last_stable` | largest load step that stayed stable (eV/Å) |
| `P_first_buckled` | first load step that jumped to a large deflection (eV/Å). Blank means it did not buckle within the ramp. |
| `Pcr_local_Timoshenko` | local critical load P_E/(1 + P_E/κGA), with P_E = π²EI0/(4L²) (eV/Å) |
| `Pcr_SD` | stress-driven critical load with the same shear correction (eV/Å) |
| `ratio_MD_local`, `ratio_SD_local` | P_MD/P_local and P_SD/P_local |
| `eps_cr_percent` | axial strain at buckling, 100·P_MD/EA0 (%) |
| `fit_a`, `fit_b` | parameters a and b of δ = a/(1 − P/P_cr) + b (Å) |

### `S2_c_from_buckling.csv`

c recovered from each column (all data, one row per column).

| Header | Meaning |
|---|---|
| `width_nm` | wire width (nm) |
| `L_A`, `L_over_d` | column length (Å) and L/d |
| `eps_cr_percent` | axial strain at buckling (%) |
| `Pcr_MD` | raw MD critical load (eV/Å) |
| `Pcr_MD_zero_strain` | MD load corrected to zero axial strain (eV/Å) |
| `Pcr_local` | local Timoshenko load (eV/Å) |
| `Pcr_SD` | SD load with the frozen-wave c (eV/Å) |
| `c_frozenwave_A` | c from the frozen-wave test (Å) |
| `c_from_raw_MD_A` | c that makes P_SD = raw P_MD (Å). Blank where MD is below the local load (no positive c exists). |
| `c_from_zero_strain_MD_A` | c that makes P_SD = zero-strain P_MD (Å) |
| `ratio_MD_local`, `ratio_MDzero_local`, `ratio_SD_local` | raw MD, zero-strain MD and SD loads divided by the local load |
| `zero_strain_intercept` | intercept at ε_cr = 0 of the linear fit of P_MD/P_SD against ε_cr for that wire |

### `S2_c_from_buckling_veusz.csv`

Same data in wide format, one column per wire (`<w>` = `1p8`, `2p2`, `2p9`), for Fig. 2a.

| Header | Meaning |
|---|---|
| `L_over_d_<w>` | x-axis: L/d of each column |
| `c_raw_<w>` | c from the raw MD load (Å). Blank = no solution. |
| `c_zero_<w>` | c from the zero-strain MD load (Å) |
| `ratio_MD_local_<w>`, `ratio_MDzero_local_<w>`, `ratio_SD_local_<w>` | load ratios as above |
| `cFW_line_x_<w>`, `cFW_line_y_<w>` | two points (L/d = 3.5 and 8.5) for the horizontal dashed line at the frozen-wave c |

### `S2_c_scan_veusz.csv`

Inversion curves for Fig. 2b. Tags `Ld4p8` and similar mean L/d = 4.8.

| Header | Meaning |
|---|---|
| `c_A` | x-axis: trial nonlocal length c from 0 to 14 Å |
| `SDratio_<w>_Ld<L/d>` | P_SD(c)/P_local for that column (solid curve) |
| `MDraw_<w>_Ld<L/d>` | raw MD level P_MD/P_local, constant (dashed line) |
| `MDzero_<w>_Ld<L/d>` | zero-strain MD level, constant (dotted line) |
| `vline_y` | y range (0.99, 1.10) for the vertical lines |
| `vline_cFW_<w>` | x position of the vertical line at the frozen-wave c (two equal values) |

---

## scripts/

| File | Purpose |
|---|---|
| `nlbeam.py` | SD, local and Eringen beam solvers: buckling (`pcr_sd`, `pcr_local`), frequencies, dispersion |
| `analyse_loadramp.py` | MD load-ramp output (JSON) → `S2_loadramp_buckling.csv` |
| `build_paperA_data.py` | full LAMMPS results → all `S3_*` and `S2_loadramp_*` CSVs |
| `c_from_buckling.py` | `S2_loadramp_buckling.csv` → `S2_c_from_buckling*.csv` and `S2_c_scan_veusz.csv` |

`fig1_frozen_wave.png` and `fig2_c_from_buckling.png` are the figures made from these files with Veusz.
