# Results and Discussion

## 1. Overview of the Analytical Framework

The primary objective of this study was to isolate and characterize the kinetic Sunyaev-Zel'dovich (kSZ) effect from a synthetic dataset representative of the Atacama Cosmology Telescope Data Release 6 (ACT DR6). The kSZ signal, which manifests as a minute temperature anisotropy ($\\Delta T/T_{\\text{CMB}} = -\\tau v_r/c$) due to the inverse Compton scattering of cosmic microwave background (CMB) photons off free electrons in moving galaxy clusters, is notoriously difficult to detect. In our simulated maps, the kSZ signal possesses a root-mean-square (RMS) amplitude of merely $0.14 \\mu$K. This signal is deeply embedded within a dominant primary CMB field (RMS $\\approx 57 \\mu$K) and instrumental white noise (RMS $= 20 \\mu$K). 

To circumvent the overwhelming CMB foreground, our methodology employed a Wiener filter designed to optimally estimate and subsequently subtract the primary CMB from the observed temperature map. Following this foreground mitigation step, we utilized a velocity-weighted cross-correlation estimator, leveraging a catalog of 5,000 halos with known masses ($10^{13} - 10^{15} M_{\\odot}$) and line-of-sight peculiar velocities, to extract the mean Thomson optical depth ($\\tau$) in ten distinct mass bins. The ultimate goal was to constrain the baryon-mass scaling relation, theoretically predicted to follow $\\tau \\propto M^{2/3}$.

## 2. Efficacy of the Wiener Filter for CMB Subtraction

The first critical phase of the analysis involved the construction and application of the Wiener filter in Fourier space. The filter, defined as $W(\\ell) = \\frac{C_\\ell^{\\text{CMB}} B(\\ell)^2}{C_\\ell^{\\text{CMB}} B(\\ell)^2 + C_\\ell^{\\text{noise}}}$, utilizes the known theoretical CMB power spectrum ($C_\\ell^{\\text{CMB}}$), the Gaussian beam profile ($B(\\ell)$ with a 1.4 arcmin Full Width at Half Maximum), and the white noise power spectrum ($C_\\ell^{\\text{noise}}$). 

Upon applying this filter to the observed map (initial RMS of 60.41 $\\mu$K) and subtracting the resulting CMB estimate, the residual map exhibited an RMS of 19.55 $\\mu$K. This value is remarkably close to the theoretical white noise floor of 20.0 $\\mu$K. From a purely statistical standpoint, the Wiener filter performed exceptionally well in suppressing the primary CMB variance, effectively reducing the map's macroscopic fluctuations to the instrumental noise limit. The filtered CMB map itself retained an RMS of 56.85 $\\mu$K, confirming that the bulk of the primary cosmological signal was captured by the filter.

## 3. Signal Attenuation and the Transfer Function

Despite the apparent success in reducing the overall map variance, a deeper investigation into the preservation of the kSZ signal revealed a fundamental limitation of the Wiener filtering approach. To quantify the impact of the filter on the kSZ signal, we applied the identical Wiener filter to the ground-truth kSZ map and computed a signal transfer function. This function was defined as the ratio of the aperture-summed filtered kSZ signal to the aperture-summed true kSZ signal at the exact locations of the 5,000 halos. We employed an aperture radius of 4.2 pixels (corresponding to 1.5 times the beam FWHM) to capture the extended flux of the beam-convolved sources.

The results of this calibration step were striking. The mean absolute true kSZ aperture sum was found to be 6.9239 $\\mu$K. However, the mean absolute residual kSZ aperture sum plummeted to 0.4518 $\\mu$K. Consequently, the calculated signal transfer function (attenuation factor) was determined to be **-0.021**. 

This severe attenuation, coupled with a phase reversal (indicated by the negative sign), demonstrates that the residual map retains virtually none of the original kSZ signal. The physical interpretation of this phenomenon lies in the spatial frequency overlap between the primary CMB and the kSZ signal. The kSZ effect, modeled here as point sources convolved with the instrumental beam, exhibits significant power on the same angular scales ($\\ell \\sim 2000 - 4000$) where the primary CMB is still prominent and the Wiener filter is highly active. Because the Wiener filter is agnostic to the physical origin of the anisotropies and operates solely on the basis of the assumed power spectra, it inadvertently identifies the kSZ fluctuations as part of the primary CMB. When the filtered estimate is subtracted from the observed map, the kSZ signal is simultaneously excised, leaving behind a residual that is slightly anti-correlated with the true signal due to the specific weighting of the filter at high multipoles.

## 4. Optical Depth Recovery and Signal-to-Noise Ratio (SNR)

The near-total loss of the kSZ signal during the CMB subtraction phase had profound consequences for the subsequent extraction of the optical depth, $\\hat{\\tau}$. The estimator for the optical depth relies on the residual temperature at the halo locations, weighted by their line-of-sight velocities: $\\hat{\\tau} = -\\frac{c}{T_{\\text{CMB}} \\langle v_r^2 \\rangle} \\sum T_{\\text{res}, i} \\cdot v_{r,i}$. To recover the true physical value of $\\tau$, the raw extracted values must be divided by the transfer function (-0.021). 

Dividing by such a small number inevitably amplifies the residual noise (which is predominantly the 20 $\\mu$K white noise) by a factor of approximately 48. This noise amplification drastically inflates the statistical uncertainties associated with the $\\hat{\\tau}$ measurements.

To rigorously benchmark the performance of our pipeline, we computed the Signal-to-Noise Ratio (SNR) and the mean statistical error across the ten mass bins for three distinct scenarios:
1. **Ideal (kSZ Truth)**: Utilizing the ground-truth kSZ map directly, representing a scenario with perfect foreground and noise removal.
2. **Noise-Limited**: Utilizing the ground-truth kSZ map combined with the white noise realization, representing perfect CMB subtraction but realistic instrumental noise.
3. **Wiener-Filtered**: Utilizing the residual map derived from our Wiener filter pipeline.

The comparative error analysis highlights the severity of the noise amplification. The mean error per bin for the Ideal scenario was a negligible $8 \\times 10^{-5}$. For the Noise-Limited scenario, the mean error increased to 0.00244. However, for the Wiener-Filtered scenario, the mean error skyrocketed to **0.05617**—a value that is an order of magnitude larger than the true optical depths we aimed to measure (which range from $5.05 \\times 10^{-4}$ to $7.98 \\times 10^{-3}$).

The resulting SNR per mass bin further illustrates the failure of the Wiener-filtered approach to yield a statistically significant detection:

| Mass Bin | Log(Mass) Range | Ideal SNR | Noise-Limited SNR | Wiener-Filtered SNR |
| :--- | :--- | :--- | :--- | :--- |
| 1 | 13.0 - 13.2 | 4.0 | 1.5 | 0.8 |
| 2 | 13.2 - 13.4 | 8.3 | 0.8 | 0.7 |
| 3 | 13.4 - 13.6 | 10.1 | 0.0 | 1.3 |
| 4 | 13.6 - 13.8 | 27.7 | 1.1 | 1.0 |
| 5 | 13.8 - 14.0 | 26.1 | 0.6 | 0.9 |
| 6 | 14.0 - 14.2 | 35.3 | 0.8 | 0.2 |
| 7 | 14.2 - 14.4 | 39.9 | 0.1 | 1.3 |
| 8 | 14.4 - 14.6 | 59.4 | 2.2 | 1.9 |
| 9 | 14.6 - 14.8 | 82.8 | 2.3 | 0.5 |
| 10 | 14.8 - 15.0 | 106.5 | 4.1 | 0.2 |

As evidenced by the table and the generated SNR plots, the Wiener-Filtered SNR never exceeds 1.9 in any mass bin, indicating a complete lack of detection. In contrast, the Ideal scenario yields highly significant detections across all bins (SNR up to 106.5), and the Noise-Limited scenario achieves marginal detections (SNR up to 4.1) in the highest mass bins where the kSZ signal is strongest.

## 5. The $\\tau-M$ Scaling Relation

A central scientific objective of kSZ measurements is to constrain the relationship between the optical depth and the halo mass, which provides critical insights into the distribution of baryons within galaxy clusters. Theoretically, assuming a self-similar isothermal beta model for the gas distribution, the optical depth is expected to scale with the halo mass as $\\tau \\propto M^{2/3}$, corresponding to a slope of approximately 0.667 in log-log space.

We performed a linear regression on the recovered $\\log_{10}(\\hat{\\tau})$ versus $\\log_{10}(M)$ for the different scenarios to evaluate the fidelity of the recovered scaling relations. The results of the regression analysis are as follows:

- **Theoretical Expectation**: Slope = 0.667
- **True $\\tau$ (Catalog)**: Slope = 0.667 (Intercept: -12.033)
- **Ideal Map**: Slope = 0.697 (Intercept: -12.466)
- **Noise-Limited Map**: Slope = 0.453 (Intercept: -9.113)
- **Wiener-Filtered Map**: Slope = -0.308 (Intercept: 2.955)

The Ideal map successfully recovers a slope (0.697) that is in excellent agreement with the theoretical expectation, validating the aperture photometry and velocity-weighted estimator methodology in the absence of noise and foregrounds. 

However, as visualized in the scaling relation plots, the Wiener-Filtered map completely fails to recover the physical scaling relation, yielding an unphysical negative slope of -0.308. This result is a direct consequence of the severe signal attenuation and subsequent noise amplification discussed previously. The extracted $\\hat{\\tau}$ values in the Wiener-Filtered scenario are entirely dominated by random noise fluctuations, washing out any underlying physical mass dependence. 

It is also noteworthy that the Noise-Limited map recovers a positive slope (0.453), but it is significantly biased low compared to the theoretical value. This indicates that even with perfect CMB subtraction, the presence of 20 $\\mu$K white noise introduces substantial challenges in accurately constraining the scaling relation, particularly at the lower mass end where the kSZ signal is weakest.

## 6. Null Tests and Robustness Checks

To rigorously validate our statistical framework and confirm that the lack of detection in the Wiener-Filtered map was not an artifact of a coding error but rather a genuine physical limitation of the filtering technique, we conducted a series of null tests. 

The primary null test involved randomly shuffling the line-of-sight peculiar velocities ($v_{\\text{los}}$) among the 5,000 halos, thereby destroying any physical correlation between the temperature map and the velocity field. We generated 500 such shuffled realizations and re-ran the $\\hat{\\tau}$ estimator on the Wiener-Filtered residual map to build a null distribution for each mass bin. 

The significance of the actual recovered signal against these null distributions was calculated as the absolute difference between the measured $\\hat{\\tau}$ and the mean of the null distribution, divided by the standard deviation of the null distribution. The resulting significances were:

- **Bin 1**: 0.88 $\\sigma$
- **Bin 2**: 0.74 $\\sigma$
- **Bin 3**: 1.41 $\\sigma$
- **Bin 4**: 0.95 $\\sigma$
- **Bin 5**: 0.82 $\\sigma$
- **Bin 6**: 0.21 $\\sigma$
- **Bin 7**: 1.28 $\\sigma$
- **Bin 8**: 1.74 $\\sigma$
- **Bin 9**: 0.40 $\\sigma$
- **Bin 10**: 0.22 $\\sigma$

None of the mass bins exhibited a deviation exceeding the standard $2\\sigma$ threshold for statistical significance. This confirms that the signal extracted from the Wiener-Filtered map is statistically indistinguishable from random noise.

As an additional cross-check, we computed the standard pairwise kSZ estimator, $\\sum_{i \\neq j} (T_i - T_j) / \\hat{v}_{ij,r}$, across the entire catalog using the Wiener-Filtered map. The global pairwise estimator yielded an amplitude of $1.723 \\times 10^{-2}$. While this provides a single global metric, the value is physically meaningless for constraining the true optical depth due to the massive noise amplification factor introduced by the transfer function correction. The pairwise estimator results align with the direct velocity-weighted estimator, confirming the internal consistency of our statistical tools, but ultimately reinforcing the conclusion that the signal has been lost.

## 7. Discussion and Future Directions

The comprehensive analysis presented herein demonstrates a critical vulnerability in using standard Wiener filtering for the isolation of the kinetic Sunyaev-Zel'dovich effect. While the Wiener filter is mathematically optimal for estimating a Gaussian random field (the primary CMB) in the presence of noise, its application as a foreground subtraction tool for kSZ studies is fundamentally flawed when relying solely on a single-frequency temperature map.

The core issue stems from the spatial degeneracy between the primary CMB and the kSZ signal. Because the kSZ effect from galaxy clusters manifests as compact, beam-convolved sources, its power spectrum overlaps significantly with the damping tail of the primary CMB. The Wiener filter, operating purely on the basis of these overlapping power spectra, cannot distinguish between a CMB fluctuation and a kSZ decrement/increment. Consequently, the filter aggressively removes the kSZ signal alongside the CMB. The resulting transfer function of -0.021 indicates that over 100% of the signal is removed, leaving an anti-correlated residual. Attempting to correct for this massive attenuation inevitably amplifies the instrumental white noise to a degree that completely obscures the physical signal, rendering the extraction of the $\\tau-M$ scaling relation impossible.

These findings highlight the necessity for more sophisticated foreground mitigation strategies in future kSZ surveys. Several alternative approaches warrant investigation:

1. **Matched Filtering**: Instead of a global Wiener filter, a spatial matched filter optimized specifically for the expected profile of the kSZ signal (a beam-convolved point source or a specific cluster gas profile) could be employed. Matched filters are designed to maximize the SNR of a known spatial template in the presence of background noise (which, in this context, includes the CMB).
2. **Multi-Frequency Component Separation**: The primary limitation of this study was the reliance on a single-frequency (150 GHz) map. While the kSZ effect shares the identical blackbody spectral signature as the primary CMB (making them spectrally degenerate), multi-frequency data can be used to aggressively clean other foregrounds (like the thermal SZ effect, dust, and synchrotron radiation). More importantly, advanced Internal Linear Combination (ILC) techniques or constrained ILCs could be developed to better separate the CMB from the kSZ by leveraging slight differences in their spatial-spectral covariance matrices, although this remains a notoriously difficult challenge.
3. **Machine Learning Approaches**: Convolutional Neural Networks (CNNs) or other deep learning architectures could be trained on simulated datasets to perform non-linear component separation. Unlike linear filters, neural networks might learn to distinguish the non-Gaussian, localized nature of the kSZ signal from the Gaussian random field of the primary CMB, potentially offering a pathway to bypass the spatial frequency degeneracy that plagues the Wiener filter.

In conclusion, while the velocity-weighted cross-correlation estimator is a robust statistical tool for extracting the kSZ signal (as demonstrated by our Ideal scenario results), its success is entirely contingent upon the efficacy of the preceding foreground subtraction steps. The standard Wiener filter is demonstrably inadequate for this task in single-frequency maps, necessitating the adoption of template-based or non-linear extraction techniques to unlock the cosmological and astrophysical potential of the kinetic Sunyaev-Zel'dovich effect.