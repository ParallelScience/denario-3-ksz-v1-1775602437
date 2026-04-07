# filename: codebase/step_5.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import datetime
import warnings
import time

warnings.filterwarnings('ignore', category=UserWarning)
plt.rcParams['text.usetex'] = False

def aperture_sum(map_data, ix_arr, iy_arr, r_pix):
    ny, nx = map_data.shape
    r_int = int(np.ceil(r_pix))
    y, x = np.ogrid[-r_int:r_int+1, -r_int:r_int+1]
    mask_circle = x**2 + y**2 <= r_pix**2
    sums = np.zeros(len(ix_arr))
    for i in range(len(ix_arr)):
        cx, cy = ix_arr[i], iy_arr[i]
        y_indices = cy + np.arange(-r_int, r_int+1)
        x_indices = cx + np.arange(-r_int, r_int+1)
        valid_y = (y_indices >= 0) & (y_indices < ny)
        valid_x = (x_indices >= 0) & (x_indices < nx)
        if not (valid_y.any() and valid_x.any()):
            continue
        y_min_sub = np.argmax(valid_y)
        y_max_sub = len(valid_y) - np.argmax(valid_y[::-1])
        x_min_sub = np.argmax(valid_x)
        x_max_sub = len(valid_x) - np.argmax(valid_x[::-1])
        map_y_min = y_indices[y_min_sub]
        map_y_max = y_indices[y_max_sub-1] + 1
        map_x_min = x_indices[x_min_sub]
        map_x_max = x_indices[x_max_sub-1] + 1
        sub_map = np.zeros((2*r_int+1, 2*r_int+1))
        sub_map[y_min_sub:y_max_sub, x_min_sub:x_max_sub] = map_data[map_y_min:map_y_max, map_x_min:map_x_max]
        sums[i] = np.sum(sub_map[mask_circle])
    return sums

if __name__ == '__main__':
    catalog = np.load('/home/node/work/projects/ksz_v1/halo_catalog.npz')
    ksz_map = np.load('/home/node/work/projects/ksz_v1/ksz_map_truth.npy')
    noise_map = np.load('/home/node/work/projects/ksz_v1/noise_map.npy')
    noise_limited_map = ksz_map + noise_map
    ix = catalog['ix_pixel']
    iy = catalog['iy_pixel']
    v_los = catalog['v_los_km_s']
    log_mass = catalog['log_mass']
    mass = catalog['mass_Msun']
    true_tau = catalog['tau']
    delta_T = catalog['delta_T_ksz_uK']
    radius_pix = 1.5 * 1.4 / 0.5
    T_ideal = aperture_sum(ksz_map, ix, iy, radius_pix)
    T_nl = aperture_sum(noise_limited_map, ix, iy, radius_pix)
    calib_ideal = np.sum(T_ideal * delta_T) / np.sum(delta_T**2)
    print("--- Calibration ---")
    print("Ideal map calibration factor (aperture sum / true delta T): " + str(np.round(calib_ideal, 4)))
    n_bins = 10
    bins = np.linspace(13, 15, n_bins + 1)
    bin_indices = np.digitize(log_mass, bins) - 1
    bin_indices[bin_indices == n_bins] = n_bins - 1
    c_km_s = 299792.458
    T_CMB_uK = 2.7255e6
    def get_tau_stats(T_ext, calib, n_boot=500):
        tau_b = np.zeros(n_bins)
        tau_err_b = np.zeros(n_bins)
        for b in range(n_bins):
            mask = bin_indices == b
            n_halos = np.sum(mask)
            if n_halos == 0:
                continue
            T_b = T_ext[mask]
            v_b = v_los[mask]
            num = np.sum(T_b * v_b)
            den = np.sum(v_b**2)
            tau_raw = - (c_km_s / T_CMB_uK) * (num / den) if den != 0 else 0.0
            tau_b[b] = tau_raw / calib
            boot_tau = np.zeros(n_boot)
            for i in range(n_boot):
                idx = np.random.randint(0, n_halos, n_halos)
                num_boot = np.sum(T_b[idx] * v_b[idx])
                den_boot = np.sum(v_b[idx]**2)
                boot_tau[i] = - (c_km_s / T_CMB_uK) * (num_boot / den_boot) / calib if den_boot != 0 else 0.0
            tau_err_b[b] = np.std(boot_tau)
        return tau_b, tau_err_b
    tau_ideal, err_ideal = get_tau_stats(T_ideal, calib_ideal)
    tau_nl, err_nl = get_tau_stats(T_nl, calib_ideal)
    res_step4 = np.load('data/optical_depth_results.npz')
    mass_mean_bins = res_step4['mass_mean_bins']
    tau_true_mean_bins = res_step4['tau_true_mean_bins']
    tau_hat_bins = res_step4['tau_hat_bins']
    tau_hat_err_bins = res_step4['tau_hat_err_bins']
    null_tau_distributions = res_step4['null_tau_distributions']
    tau_corrected_pw_global = res_step4['tau_corrected_pw_global'][0]
    transfer_function = np.load('data/transfer_function.npy')[0]
    snr_wf = np.abs(tau_hat_bins) / (tau_hat_err_bins + 1e-10)
    snr_ideal = np.abs(tau_ideal) / (err_ideal + 1e-10)
    snr_nl = np.abs(tau_nl) / (err_nl + 1e-10)
    print("\n--- Signal-to-Noise Ratio (SNR) per Mass Bin ---")
    print("Bin | Ideal SNR | Noise-Limited SNR | Wiener-Filtered SNR")
    for b in range(n_bins):
        print(" " + str(b+1).ljust(2) + " | " + str(np.round(snr_ideal[b], 1)).rjust(9) + " | " + str(np.round(snr_nl[b], 1)).rjust(17) + " | " + str(np.round(snr_wf[b], 1)).rjust(19))
    print("\n--- Error Analysis ---")
    print("Mean error per bin (Ideal): " + str(np.round(np.mean(err_ideal), 5)))
    print("Mean error per bin (Noise-Limited): " + str(np.round(np.mean(err_nl), 5)))
    print("Mean error per bin (Wiener-Filtered): " + str(np.round(np.mean(tau_hat_err_bins), 5)))
    def do_regression(x, y, name):
        valid = y > 0
        if np.sum(valid) >= 2:
            res = linregress(np.log10(x[valid]), np.log10(y[valid]))
            print(name + " slope: " + str(np.round(res.slope, 3)) + " (Intercept: " + str(np.round(res.intercept, 3)) + ")")
            return res
        else:
            print(name + " slope: NaN (not enough positive bins)")
            return None
    print("\n--- Scaling Relation (log-log linear regression) ---")
    print("Theoretical slope: 0.667")
    res_true = do_regression(mass_mean_bins, tau_true_mean_bins, "True tau")
    res_ideal = do_regression(mass_mean_bins, tau_ideal, "Ideal map")
    res_nl = do_regression(mass_mean_bins, tau_nl, "Noise-limited map")
    res_wf = do_regression(mass_mean_bins, tau_hat_bins, "Wiener-filtered map")
    print("\n--- Global Pairwise kSZ Estimator ---")
    print("Amplitude (tau): " + str(np.format_float_scientific(tau_corrected_pw_global, precision=3)))
    print("\n--- Transfer Function ---")
    print("Correction factor (Attenuation): " + str(np.round(transfer_function, 4)))
    print("\n--- Null Test Significance (Wiener-Filtered) ---")
    null_means = np.mean(null_tau_distributions, axis=1)
    null_stds = np.std(null_tau_distributions, axis=1)
    for b in range(n_bins):
        sig = np.abs(tau_hat_bins[b] - null_means[b]) / null_stds[b] if null_stds[b] > 0 else 0.0
        print("Bin " + str(b+1) + ": " + str(np.round(sig, 2)) + " sigma")
    print("\nNote: Negative tau values in the Wiener-filtered map are omitted from the log-log plot.")
    fig, axs = plt.subplots(1, 3, figsize=(18, 5))
    axs[0].errorbar(mass_mean_bins, tau_ideal, yerr=err_ideal, fmt='o-', label='Ideal (kSZ Truth)')
    axs[0].errorbar(mass_mean_bins, tau_nl, yerr=err_nl, fmt='s-', label='Noise-Limited')
    axs[0].errorbar(mass_mean_bins, tau_hat_bins, yerr=tau_hat_err_bins, fmt='^-', label='Wiener-Filtered')
    axs[0].plot(mass_mean_bins, tau_true_mean_bins, 'k--', label='True tau proportional to M^(2/3)')
    axs[0].set_xscale('log')
    axs[0].set_yscale('log')
    axs[0].set_xlabel('Halo Mass M_500 [M_sun]')
    axs[0].set_ylabel('Optical Depth tau')
    axs[0].set_title('tau - M Scaling Relation')
    axs[0].legend()
    axs[0].grid(True, which='both', ls='--', alpha=0.5)
    axs[0].set_ylim(1e-4, 1e-1)
    axs[1].plot(mass_mean_bins, snr_ideal, 'o-', label='Ideal')
    axs[1].plot(mass_mean_bins, snr_nl, 's-', label='Noise-Limited')
    axs[1].plot(mass_mean_bins, snr_wf, '^-', label='Wiener-Filtered')
    axs[1].set_xscale('log')
    axs[1].set_yscale('log')
    axs[1].set_xlabel('Halo Mass M_500 [M_sun]')
    axs[1].set_ylabel('Signal-to-Noise Ratio (SNR)')
    axs[1].set_title('SNR per Mass Bin')
    axs[1].legend()
    axs[1].grid(True, which='both', ls='--', alpha=0.5)
    axs[2].errorbar(mass_mean_bins, null_means, yerr=null_stds, fmt='o-', color='gray', label='Null Test (Shuffled v_los)')
    axs[2].errorbar(mass_mean_bins, tau_hat_bins, yerr=tau_hat_err_bins, fmt='^-', color='red', label='Wiener-Filtered Signal')
    axs[2].set_xscale('log')
    axs[2].set_xlabel('Halo Mass M_500 [M_sun]')
    axs[2].set_ylabel('Estimated tau')
    axs[2].set_title('Null Test vs Recovered Signal')
    axs[2].legend()
    axs[2].grid(True, which='both', ls='--', alpha=0.5)
    plt.tight_layout()
    timestamp = int(time.time())
    plot_filename = "data/scaling_relation_analysis_" + str(timestamp) + ".png"
    plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
    print("\nPlot saved to " + plot_filename)