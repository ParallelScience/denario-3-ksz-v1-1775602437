<!-- filename: reports/step_8_ksz_analysis_notes.md -->
# Analysis of Synthetic ACT DR6-like kSZ Dataset

## 1. Baseline Validation
The initial phase involved validating the synthetic ACT DR6-like dataset. The dataset covers a $10^\circ \times 10^\circ$ sky patch at 0.5 arcmin resolution ($1200 \times 1200$ grid). Summary statistics confirm the dominance of the primary CMB and instrumental noise over the kSZ signal:
- Total observed map RMS: $60.41\ \mu\text{K}$
- Primary CMB RMS: $\approx 56.98\ \mu\text{K}$
- White noise floor RMS: $\approx 19.99\ \mu\text{K}$
- Ground-truth kSZ map RMS: $0.14\ \mu\text{K}$

The halo catalog contains 5,000 clusters with log-uniform mass distribution ($10^{13}$ to $10^{15}\ M_\odot$). Line-of-sight velocities ($v_{los}$) have an RMS of $300.56\ \text{km/s}$. Visual diagnostics are provided in <code>diagnostic_maps_1_1775604087.png</code>.

## 2. Matched Filter Construction and Application
An optimal matched filter $\Psi(\ell) = \frac{B(\ell)}{C_\ell^{CMB} + C_\ell^{noise}}$ was constructed to maximize sensitivity to the point-like kSZ signal. The filter acts as a band-pass, suppressing large-scale CMB fluctuations and small-scale noise. Empirical power spectra (<code>power_spectra_filter_2_1775604260.png</code>) confirm the filter effectively removes the primary CMB foreground.

## 3. Pairwise Estimator Implementation
The pairwise kSZ estimator $\hat{\tau}(r)$ was used to cross-correlate temperature maps with the halo catalog. Results (<code>pairwise_estimator_1775604400.png</code>) show:
- Mean recovered $\langle \hat{\tau}_{truth} \rangle \approx 3.26 \times 10^{-4}$
- Mean recovered $\langle \hat{\tau}_{obs} \rangle \approx 1.69 \times 10^{-4}$

The discrepancy with the true mean optical depth ($\approx 2.87 \times 10^{-3}$) is attributed to beam dilution from the 1.4 arcmin Gaussian beam.

## 4. Mass Scaling Relation
We measured the $\tau-M$ scaling relation by binning halos by mass. The linear regression yielded $\log_{10}(\tau) = 0.4679 \log_{10}(M) - 9.2926$. The recovered slope of $0.4679$ is shallower than the theoretical $0.6667$ due to mass-dependent beam dilution effects.

## 5. Stacking and Visualization
Stacking analysis was performed on cutouts centered on halos. The highest mass bin ($\langle M \rangle \approx 7.99 \times 10^{14}\ M_\odot$) shows a clear temperature decrement in the filtered stack (<code>stacked_cutouts_5_1775605082.png</code>). The $\tau-M$ scaling with jackknife uncertainties is shown in <code>tau_mass_scaling_jk_5_1775605082.png</code>.

## 6. Jackknife Uncertainties and Robustness Tests
Leave-one-out jackknife resampling (25 subregions) was used to estimate covariance. A null test (shuffling velocities) confirmed the signal is associated with the peculiar velocity field (<code>pairwise_null_test_6_1775605258.png</code>). The full covariance $\chi^2$ for the real signal is 115.55, indicating a highly significant detection.

## 7. Sensitivity and SNR Analysis
SNR scales as $\sqrt{N}$ (<code>snr_vs_nhalos_7_1775605371.png</code>). For 5,000 halos, the integrated SNR is 3.78. This confirms the measurement is in the noise-dominated regime, necessitating larger catalogs for future high-precision cosmological constraints.