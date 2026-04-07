# Kinetic Sunyaev-Zel'dovich (kSZ) Analysis — ACT DR6-like Synthetic Data

**Scientist:** denario-3 (Denario AI Research Scientist)
**Date:** 2026-04-07
**Status:** Project initialized — ready for idea generation

## Dataset

Synthetic kSZ temperature map at ACT DR6-like sensitivity:
- 10×10 degree flat-sky patch, 0.5 arcmin pixels (1200×1200)
- Beam FWHM: 1.4 arcmin | Noise: 10 μK-arcmin
- 5000 halos with masses 10^13–10^15 M_sun, Gaussian v_los (σ=300 km/s)
- kSZ signal: ΔT = -τ v_r/c T_CMB, painted + beam-convolved

## Files

- `observed_map.npy` — total observed map (CMB + kSZ + noise), μK
- `ksz_map_truth.npy` — ground-truth kSZ signal
- `cmb_map_truth.npy` — ground-truth CMB realization
- `noise_map.npy` — noise realization
- `halo_catalog.npz` — positions, masses, velocities, optical depths

---


# Synthetic Kinetic Sunyaev-Zel'dovich (kSZ) Dataset — ACT DR6-like

## Overview

This dataset simulates a kinetic Sunyaev-Zel'dovich (kSZ) temperature map at ACT DR6-like
sensitivity, combined with a galaxy/cluster halo catalog. The kSZ effect arises from inverse
Compton scattering of CMB photons off free electrons in moving galaxy clusters; it imprints
a temperature anisotropy ΔT/T_CMB = -τ v_r/c, where τ is the Thomson optical depth and v_r
is the line-of-sight peculiar velocity.

## Map parameters

- Sky patch: 10 × 10 degrees (flat-sky approximation)
- Pixel size: 0.5 arcmin → 1200 × 1200 pixels
- Beam FWHM: 1.4 arcmin (Gaussian, representative of ACT DR6 at 150 GHz)
- White noise level: 10 μK-arcmin → 20 μK RMS per 0.5-arcmin pixel
- Temperature unit: μK (micro-Kelvin)

## File inventory (all absolute paths)

### Temperature maps — shape (1200, 1200), dtype float32, units μK

- `/home/node/work/projects/ksz_v1/observed_map.npy`
  The total observed temperature map: CMB + beam-convolved kSZ signal + white noise.
  This is the only quantity available to an observer. RMS ≈ 60 μK.

- `/home/node/work/projects/ksz_v1/cmb_map_truth.npy`
  The ground-truth beam-convolved primary CMB realization (Gaussian random field,
  approximate ΛCDM TT power spectrum, Silk-damped above ℓ ~ 4500). RMS ≈ 57 μK.

- `/home/node/work/projects/ksz_v1/ksz_map_truth.npy`
  The ground-truth beam-convolved kSZ temperature map (sum of point-source contributions
  from all halos, then beam-convolved). RMS ≈ 0.14 μK. The kSZ signal is deeply buried
  under CMB and noise on a per-pixel basis.

- `/home/node/work/projects/ksz_v1/noise_map.npy`
  The white noise realization. RMS = 20 μK exactly (10 μK-arcmin / 0.5 arcmin pixel).

### Halo catalog — `/home/node/work/projects/ksz_v1/halo_catalog.npz`

Load with: `catalog = np.load('/home/node/work/projects/ksz_v1/halo_catalog.npz')`

Arrays (all length 5000):
- `ra` (float64): right ascension in degrees, uniform ∈ [0, 10]
- `dec` (float64): declination in degrees, uniform ∈ [0, 10]
- `mass_Msun` (float64): halo mass M_500 in solar masses, range 10^13 – 10^15 M_sun, log-uniform
- `log_mass` (float64): log10(M_500 / M_sun), uniform ∈ [13, 15]
- `v_los_km_s` (float64): line-of-sight peculiar velocity in km/s, Gaussian with σ = 300 km/s
- `tau` (float64): Thomson optical depth, τ = 2×10^-3 × (M/10^14)^{2/3}, range [4.3×10^-4, 9.3×10^-3]
- `delta_T_ksz_uK` (float64): kSZ temperature decrement ΔT = -τ v_r/c × T_CMB in μK, range [-71, +76] μK per halo
- `ix_pixel` (int): map x-pixel index (column) of halo position
- `iy_pixel` (int): map y-pixel index (row) of halo position

## Data generating process

1. **CMB map**: Gaussian random field in Fourier space with approximate ΛCDM TT power spectrum
   D_ℓ = 6000 × (ℓ/200)² / (1 + (ℓ/200)^1.2)^2.5 × exp(−(ℓ/4500)²) μK²,
   then beam-convolved (Gaussian, FWHM 1.4 arcmin).

2. **Halo catalog**: 5000 halos with positions uniform on the 10×10 deg patch, masses
   log-uniform in [10^13, 10^15] M_sun, v_los ~ N(0, 300² km²/s²).

3. **kSZ map**: each halo is painted as a delta function at its pixel with amplitude
   ΔT_kSZ, then the full map is Gaussian beam-convolved.

4. **Noise**: i.i.d. Gaussian white noise, σ = 20 μK per pixel.

5. **Observed map** = CMB + kSZ + noise.

Random seed: 42.

## Physical context and key scales

- The kSZ signal is sub-dominant: kSZ map RMS (0.14 μK) ≪ noise (20 μK) ≪ CMB (57 μK).
- Individual kSZ detections per halo are not possible (per-halo SNR ≪ 1).
- Statistical detection requires cross-correlating the map with the halo catalog, exploiting
  the known positions and (in real surveys) photometric/spectroscopic redshifts.
- The canonical statistic is the *pairwise kSZ estimator*: Σ_{i≠j} (T_i − T_j) / v̂_{ij,r}
  weighted by separation, which measures ⟨τ⟩ × pairwise velocity correlation.
- The optical depth × velocity product ⟨τ v_r⟩ encodes both the baryon content of halos
  and the large-scale velocity field (cosmological information).

## Suggested analyses

1. **Pairwise kSZ estimator**: measure ΔT̄(r) = Σ w_ij (T_i − T_j) as a function of
   comoving separation r_ij; fit amplitude to extract ⟨τ_e⟩.
2. **Aperture photometry / filtered maps**: apply a CMB-cleaning filter (e.g. matched
   filter or simple high-pass) at halo positions; stack to measure mean kSZ signal as
   a function of halo mass.
3. **kSZ power spectrum**: estimate the kSZ contribution to the CMB temperature power
   spectrum by cross-correlating the observed map with the halo density field.
4. **Mass dependence of τ**: bin halos by mass and measure the τ–mass scaling relation
   from stacked kSZ signal.
5. **Velocity reconstruction**: compare the kSZ-inferred velocities against the known
   true v_los to quantify reconstruction fidelity.
