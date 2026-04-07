# Results Interpretation and Reporting

## 1. Effectiveness of Wiener-Filtered CMB Subtraction

The primary goal of this study was to isolate the kinetic Sunyaev-Zel'dovich (kSZ) signal by subtracting a Wiener-filtered estimate of the primary CMB from the observed map. The Wiener filter was constructed to optimally estimate the CMB given the known CMB power spectrum, beam, and white noise levels. The residual map was then used to extract the kSZ signal via aperture photometry at the known halo locations.

However, the analysis reveals that the Wiener filter approach, as implemented, severely attenuates the kSZ signal. The signal transfer function (attenuation factor) was calculated to be **-0.021**. This indicates that the residual map retains almost none of the original kSZ signal (and is slightly anti-correlated). Because the kSZ signal is a point source convolved with the beam, its spatial frequencies overlap significantly with the primary CMB. The Wiener filter, designed to extract the CMB, inadvertently captures the kSZ signal as well, leading to its removal when the filtered CMB is subtracted from the observed map.

Consequently, correcting the extracted signal by dividing by this small transfer function amplifies the residual noise by a factor of ~50. This results in very large statistical uncertainties in the recovered optical depths (\hat{\tau}).

## 2. Signal-to-Noise Ratio (SNR) and Error Analysis

We compared the SNR across ten mass bins for three scenarios:
1. **Ideal (kSZ Truth)**: Using the ground-truth kSZ map.
2. **Noise-Limited**: Using the kSZ map + white noise (perfect CMB subtraction).
3. **Wiener-Filtered**: Using the residual map from the Wiener filter pipeline.

The mean error per bin for the Wiener-filtered map is **0.056**, which is orders of magnitude larger than the true \tau values (\sim 10^{-4} - 10^{-3}). In contrast, the noise-limited scenario has a mean error of **0.0024**, and the ideal scenario has **0.00008**.

The SNR per bin for the Wiener-filtered map ranges from **0.2 to 1.9**, meaning no statistically significant detection of the kSZ signal was achieved in any individual mass bin. The noise-limited scenario achieved marginal detections (SNR up to 4.1 in the highest mass bin), while the ideal scenario achieved highly significant detections (SNR 4.0 to 106.5).

**SNR per Mass Bin:**
- Bin 1: Ideal 4.0, Noise-Limited 1.5, Wiener-Filtered 0.8
- Bin 2: Ideal 8.3, Noise-Limited 0.8, Wiener-Filtered 0.7
- Bin 3: Ideal 10.1, Noise-Limited 0.0, Wiener-Filtered 1.3
- Bin 4: Ideal 27.7, Noise-Limited 1.1, Wiener-Filtered 1.0
- Bin 5: Ideal 26.1, Noise-Limited 0.6, Wiener-Filtered 0.9
- Bin 6: Ideal 35.3, Noise-Limited 0.8, Wiener-Filtered 0.2
- Bin 7: Ideal 39.9, Noise-Limited 0.1, Wiener-Filtered 1.3
- Bin 8: Ideal 59.4, Noise-Limited 2.2, Wiener-Filtered 1.9
- Bin 9: Ideal 82.8, Noise-Limited 2.3, Wiener-Filtered 0.5
- Bin 10: Ideal 106.5, Noise-Limited 4.1, Wiener-Filtered 0.2

## 3. Scaling Relation (\tau - M)

We performed a log-log linear regression to determine the scaling relation between the recovered optical depth and halo mass. The theoretical expectation is \tau \propto M^{2/3} (slope \approx 0.667).

- **Theoretical Slope**: 0.667
- **True \tau Slope**: 0.667 (Intercept: -12.033)
- **Ideal Map Slope**: 0.697 (Intercept: -12.466)
- **Noise-Limited Map Slope**: 0.453 (Intercept: -9.113)
- **Wiener-Filtered Map Slope**: -0.308 (Intercept: 2.955)

The Wiener-filtered map fails to recover the theoretical scaling relation, yielding a negative slope (-0.308). This is a direct consequence of the severe signal attenuation and noise amplification, which completely washes out the physical mass dependence. The noise-limited map recovers a positive slope (0.453), though biased low due to noise. The ideal map successfully recovers a slope (0.697) very close to the theoretical value.

## 4. Null Tests and Pairwise kSZ Estimator

To validate the pipeline, we performed a null test by shuffling the line-of-sight peculiar velocities 500 times and re-running the estimator on the Wiener-filtered map. The significance of the recovered signal against the null distributions is as follows:
- Bin 1: 0.88 \sigma
- Bin 2: 0.74 \sigma
- Bin 3: 1.41 \sigma
- Bin 4: 0.95 \sigma
- Bin 5: 0.82 \sigma
- Bin 6: 0.21 \sigma
- Bin 7: 1.28 \sigma
- Bin 8: 1.74 \sigma
- Bin 9: 0.40 \sigma
- Bin 10: 0.22 \sigma

This confirms that the recovered signal in the Wiener-filtered map is consistent with noise, as no bin exceeds a 2\sigma detection threshold.

As a cross-check, we computed the global pairwise kSZ estimator on the Wiener-filtered map, yielding an amplitude of **1.723 \times 10^{-2}**. Given the large noise amplification, this value is not physically meaningful for constraining the true optical depth.

## 5. Conclusion

The Wiener-filtered CMB subtraction method, while effective at removing the primary CMB variance, simultaneously removes the vast majority of the kSZ signal due to the overlap in their spatial power spectra. The resulting signal attenuation (-0.021) makes it impossible to recover the \tau-M scaling relation in the presence of 20 \mu K white noise. Future work should explore alternative foreground cleaning techniques, such as matched filtering optimized specifically for the kSZ spatial profile, or utilizing multi-frequency data to separate the CMB and kSZ signals based on their spectral signatures.
