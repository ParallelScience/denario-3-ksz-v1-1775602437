# filename: codebase/step_7.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

def sensitivity_and_snr(data_dir, filtered_map_path, ksz_truth_path, catalog_path):
    filt_map = np.load(filtered_map_path)
    ksz_truth = np.load(ksz_truth_path)
    catalog = pd.read_csv(catalog_path)
    iy = catalog['iy_pixel'].values
    ix = catalog['ix_pixel'].values
    v_los = catalog['v_los_km_s'].values
    T_filt = filt_map[iy, ix]
    T_truth = ksz_truth[iy, ix]
    x_deg = catalog['ra'].values
    y_deg = catalog['dec'].values
    x_idx = np.floor(x_deg / 2.0).astype(int)
    y_idx = np.floor(y_deg / 2.0).astype(int)
    x_idx[x_idx == 5] = 4
    y_idx[y_idx == 5] = 4
    region_id = x_idx + 5 * y_idx
    N_halos_list = np.array([500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000])
    n_iterations = 20
    snr_obs_mean = np.zeros(len(N_halos_list))
    snr_obs_std = np.zeros(len(N_halos_list))
    snr_truth_mean = np.zeros(len(N_halos_list))
    snr_truth_std = np.zeros(len(N_halos_list))
    c = 299792.458
    T_cmb = 2.7255e6
    def compute_tau_and_err(T, v, regions):
        num = np.sum(T * v)
        den = np.sum(v**2)
        tau = - (c / T_cmb) * (num / den) if den > 0 else np.nan
        unique_regions = np.unique(regions)
        N_jk = len(unique_regions)
        tau_jk = np.zeros(N_jk)
        for k, reg in enumerate(unique_regions):
            mask = (regions != reg)
            num_jk = np.sum(T[mask] * v[mask])
            den_jk = np.sum(v[mask]**2)
            tau_jk[k] = - (c / T_cmb) * (num_jk / den_jk) if den_jk > 0 else np.nan
        valid = ~np.isnan(tau_jk)
        N_valid = np.sum(valid)
        if N_valid > 1:
            mean_tau_jk = np.mean(tau_jk[valid])
            var_jk = (N_valid - 1) / N_valid * np.sum((tau_jk[valid] - mean_tau_jk)**2)
            err = np.sqrt(var_jk)
        else:
            err = np.nan
        return tau, err
    np.random.seed(42)
    for i, N in enumerate(N_halos_list):
        snr_obs_iter = []
        snr_truth_iter = []
        iters = 1 if N == len(catalog) else n_iterations
        for it in range(iters):
            idx = np.random.choice(len(catalog), N, replace=False)
            T_filt_sub = T_filt[idx]
            T_truth_sub = T_truth[idx]
            v_sub = v_los[idx]
            reg_sub = region_id[idx]
            T_filt_sub = T_filt_sub - np.mean(T_filt_sub)
            T_truth_sub = T_truth_sub - np.mean(T_truth_sub)
            v_sub = v_sub - np.mean(v_sub)
            tau_obs, err_obs = compute_tau_and_err(T_filt_sub, v_sub, reg_sub)
            tau_truth, err_truth = compute_tau_and_err(T_truth_sub, v_sub, reg_sub)
            if err_obs > 0:
                snr_obs_iter.append(tau_obs / err_obs)
            if err_truth > 0:
                snr_truth_iter.append(tau_truth / err_truth)
        snr_obs_mean[i] = np.mean(snr_obs_iter)
        snr_obs_std[i] = np.std(snr_obs_iter) if iters > 1 else 0.0
        snr_truth_mean[i] = np.mean(snr_truth_iter)
        snr_truth_std[i] = np.std(snr_truth_iter) if iters > 1 else 0.0
    A_obs = snr_obs_mean[-1] / np.sqrt(N_halos_list[-1])
    snr_obs_theory = A_obs * np.sqrt(N_halos_list)
    A_truth = snr_truth_mean[-1] / np.sqrt(N_halos_list[-1])
    snr_truth_theory = A_truth * np.sqrt(N_halos_list)
    print('--- Sensitivity and SNR Analysis ---')
    print('Formulas used: tau_hat = - (c / T_cmb) * sum(dT_i * dv_i) / sum(dv_i^2)')
    print('SNR = tau_hat / sigma_tau')
    print('\nFinal Integrated SNR (N=5000):')
    print('  Observed Map (Filtered): ' + str(round(snr_obs_mean[-1], 2)))
    print('  Truth kSZ Map:           ' + str(round(snr_truth_mean[-1], 2)))
    mpl.rcParams['text.usetex'] = False
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    axes[0].errorbar(N_halos_list, snr_obs_mean, yerr=snr_obs_std, fmt='ro', label='Observed Map (Filtered)', capsize=4)
    axes[0].plot(N_halos_list, snr_obs_theory, 'r--', label='Theory sqrt(N)')
    axes[0].set_xlabel('Number of Halos (N)')
    axes[0].set_ylabel('Signal-to-Noise Ratio (SNR)')
    axes[0].set_title('Observed Map (Filtered) SNR vs N')
    axes[0].legend()
    axes[0].grid(True, linestyle='--', alpha=0.7)
    axes[1].errorbar(N_halos_list, snr_truth_mean, yerr=snr_truth_std, fmt='bs', label='Truth kSZ Map', capsize=4)
    axes[1].plot(N_halos_list, snr_truth_theory, 'b--', label='Theory sqrt(N)')
    axes[1].set_xlabel('Number of Halos (N)')
    axes[1].set_ylabel('Signal-to-Noise Ratio (SNR)')
    axes[1].set_title('Truth kSZ Map SNR vs N')
    axes[1].legend()
    axes[1].grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    timestamp = str(int(time.time()))
    plot_filename = 'snr_vs_nhalos_7_' + timestamp + '.png'
    plot_filepath = os.path.join(data_dir, plot_filename)
    plt.savefig(plot_filepath, dpi=300)
    plt.close()
    print('\nSNR plot saved to ' + plot_filepath)
    results_df = pd.DataFrame({'N_halos': N_halos_list, 'snr_obs_mean': snr_obs_mean, 'snr_obs_std': snr_obs_std, 'snr_truth_mean': snr_truth_mean, 'snr_truth_std': snr_truth_std})
    csv_path = os.path.join(data_dir, 'snr_analysis_results.csv')
    results_df.to_csv(csv_path, index=False)
    print('SNR results saved to ' + csv_path)

if __name__ == '__main__':
    data_directory = 'data/'
    filtered_map_path = os.path.join(data_directory, 'filtered_observed_map.npy')
    ksz_truth_path = '/home/node/work/projects/ksz_v1/ksz_map_truth.npy'
    catalog_path = os.path.join(data_directory, 'processed_halo_catalog.csv')
    sensitivity_and_snr(data_directory, filtered_map_path, ksz_truth_path, catalog_path)