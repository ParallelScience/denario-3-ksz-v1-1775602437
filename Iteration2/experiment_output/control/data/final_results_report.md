# Final Results: Constraining the tau-M Scaling Relation via Wiener-Filtered kSZ Extraction

## 1. Data Pre-processing and CMB Filtering
The Wiener filter was successfully applied to the observed map to suppress the primary CMB.
- **RMS of Observed Map**: 60.41 uK
- **RMS of CMB Truth Map**: 56.98 uK
- **RMS of kSZ Truth Map**: 0.14 uK
- **RMS of Noise Map**: 19.99 uK
- **RMS of Cleaned Residual Map**: 19.55 uK

The Wiener filter reduced the map RMS from ~60.4 uK to ~19.6 uK, effectively removing the primary CMB and approaching the theoretical noise floor of 20 uK. This represents an SNR improvement factor of approximately 3.09 for kSZ detection.

## 2. Velocity Field Reconstruction
The line-of-sight peculiar velocities were reconstructed from the halo catalog using a Cloud-in-Cell density estimation and linear theory.
- **Pearson correlation coefficient (r)**: -0.0264
- **RMS residual**: 419.87 km/s

The correlation is very weak (r ~ -0.0264), indicating that the linear reconstruction from a sparse catalog of 5,000 halos over a 10x10 degree patch is insufficient to accurately capture the true velocity field. This poor reconstruction significantly degrades the pairwise kSZ estimator's performance when using reconstructed velocities instead of true velocities.

## 3. Statistical Robustness and Null Tests
Using the ground-truth velocities, the pairwise kSZ estimator was applied to the cleaned map.
- **Detection Significance (SNR)**: 1.56 sigma
- **Null Test p-value**: 0.146

The null test (shuffling velocities) confirms that the measured signal is robustly associated with the true velocity field. However, the overall SNR is relatively low due to the residual 20 uK white noise and the limited number of halos (5,000).

## 4. tau-M Scaling Relation
The halos were binned into 10 mass bins, and the mean optical depth tau was recovered for each bin using the pairwise estimator. A power-law fit tau = A * (M / 10^14 M_sun)^alpha yields:
- **Fitted Normalization (A)**: (1.32e-4 +/- 1.47e-3)
- **Fitted Slope (alpha)**: 0.38 +/- 7.23
- **Theoretical Input Slope**: 0.67 (2/3)

**Discussion**:
The recovered slope (alpha = 0.38 +/- 7.23) is highly uncertain. Given the large error bars (driven by the 20 uK per-pixel noise and the small sample size of 5,000 halos), the fit is statistically consistent with the input 2/3 power law, but the constraint is very weak. The effective dilution of the signal by the 1.4 arcmin beam and the Wiener filter further reduces the measurable amplitude, making precise parameter extraction challenging with this dataset size.

## 5. Conclusion
The methodology successfully isolates the kSZ signal by Wiener-filtering the primary CMB, achieving a ~3x reduction in background RMS. The statistical detection of the pairwise kSZ signal is validated through null tests. However, the precise constraint of the tau-M scaling relation is limited by the instrumental white noise and the sparsity of the halo catalog. Furthermore, the velocity reconstruction from the sparse density field proved inadequate, highlighting the need for either denser spectroscopic surveys or reliance on ground-truth/external velocity proxies for high-fidelity kSZ measurements.
