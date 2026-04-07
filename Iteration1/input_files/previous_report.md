

Iteration 0:
### Project Summary: kSZ Extraction via Wiener Filtering

**Objective:** Constrain the $\tau-M$ scaling relation ($\tau \propto M^{2/3}$) using a velocity-weighted cross-correlation estimator on ACT DR6-like synthetic maps, employing Wiener filtering for primary CMB foreground mitigation.

**Methodology:**
1. **Foreground Mitigation:** Applied a Wiener filter $W(\ell)$ to the observed map (60.4 $\mu$K RMS) to estimate and subtract the primary CMB.
2. **Signal Extraction:** Used a velocity-weighted estimator $\hat{\tau} = -\frac{c}{T_{CMB} \langle v_r^2 \rangle} \sum (T_{res, i}) \cdot v_{r,i}$ at halo positions.
3. **Calibration:** Calculated a signal transfer function by applying the Wiener filter to the ground-truth kSZ map to quantify signal attenuation.
4. **Validation:** Performed mass binning (10 bins), bootstrapping for uncertainty, and null tests (velocity shuffling).

**Key Findings:**
* **Filter Failure:** The Wiener filter is spatially degenerate with the kSZ signal (both occupy $\ell \sim 2000-4000$). The filter effectively excised the kSZ signal, resulting in a transfer function of -0.021 (near-total signal loss and phase reversal).
* **Noise Amplification:** Correcting for this attenuation amplified the 20 $\mu$K white noise by a factor of ~48, rendering the residual map dominated by noise.
* **Detection Status:** No statistically significant detection of the kSZ signal was achieved (SNR < 1.9 in all bins). The recovered $\tau-M$ slope was -0.308, contradicting the theoretical 0.667.
* **Benchmark Success:** The "Ideal" scenario (using kSZ truth) confirmed the estimator's validity, yielding a slope of 0.697 and high SNR.

**Constraints & Future Directions:**
* **Constraint:** Single-frequency Wiener filtering is fundamentally inadequate for kSZ extraction due to the spatial frequency overlap with the primary CMB.
* **Decision:** Abandon global Wiener filtering for this task.
* **Future Work:** 
    * Implement **Matched Filtering** optimized for the beam-convolved cluster profile.
    * Explore **non-linear component separation** (e.g., CNNs) to exploit the non-Gaussian, localized nature of kSZ signals.
    * Investigate multi-frequency ILC methods if additional spectral data becomes available.
        