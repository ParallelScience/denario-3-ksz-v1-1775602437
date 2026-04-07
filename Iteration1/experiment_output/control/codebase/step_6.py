# filename: codebase/step_6.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

def statistical_uncertainty_and_robustness(data_dir, obs_map_path, catalog_path):
    obs_map = np.load(obs_map_path)
    catalog = pd.read_csv(catalog_path)
    iy = catalog['iy_pixel'].values
    ix = catalog['ix_pixel'].values
    T_obs = obs_map[iy, ix]
    v_los = catalog['v_los_km_s'].values
    x_deg = catalog['ra'].values
    y_deg = catalog['dec'].values
    N = len(catalog)
    x_idx = np.floor(x_deg / 2.0).astype(int)
    y_idx = np.floor(y_deg / 2.0).astype(int)
    x_idx[x_idx == 5] = 4
    y_idx[y_idx == 5] = 4
    region_id = x_idx + 5 * y_idx
    N_jk = 25
    i, j = np.triu_indices(N, k=1)
    dist_deg = np.sqrt((x_deg[i] - x_deg[j])**2 + (y_deg[i] - y_deg[j])**2)
    dist_arcmin = dist_deg * 60.0
    dv = v_los[i] - v_los[j]
    dT_obs = T_obs[i] - T_obs[j]
    reg_i = region_id[i]
    reg_j = region_id[j]
    bins = np.arange(0, 601, 30)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    num_bins = len(bins) - 1
    c = 299792.458
    T_cmb = 2.7255e6
    bin_masks = []
    for k in range(num_bins):
        bin_masks.append((dist_arcmin >= bins[k]) & (dist_arcmin < bins[k+1]))
    tau_obs = np.zeros(num_bins)
    for k in range(num_bins):
        mask = bin_masks[k]
        if np.sum(mask) > 0:
            sum_dv2 = np.sum(dv[mask]**2)
            if sum_dv2 > 0:
                tau_obs[k] = - (c / T_cmb) * np.sum(dT_obs[mask] * dv[mask]) / sum_dv2
            else:
                tau_obs[k] = np.nan
        else:
            tau_obs[k] = np.nan
    tau_jk = np.zeros((num_bins, N_jk))
    for jk in range(N_jk):
        valid_pairs = (reg_i != jk) & (reg_j != jk)
        for k in range(num_bins):
            mask = bin_masks[k] & valid_pairs
            if np.sum(mask) > 0:
                sum_dv2 = np.sum(dv[mask]**2)
                if sum_dv2 > 0:
                    tau_jk[k, jk] = - (c / T_cmb) * np.sum(dT_obs[mask] * dv[mask]) / sum_dv2
                else:
                    tau_jk[k, jk] = np.nan
            else:
                tau_jk[k, jk] = np.nan
    cov_matrix = np.zeros((num_bins, num_bins))
    mean_tau_jk = np.nanmean(tau_jk, axis=1)
    for k1 in range(num_bins):
        for k2 in range(num_bins):
            valid_jk = ~np.isnan(tau_jk[k1]) & ~np.isnan(tau_jk[k2])
            N_valid = np.sum(valid_jk)
            if N_valid > 1:
                diff1 = tau_jk[k1, valid_jk] - mean_tau_jk[k1]
                diff2 = tau_jk[k2, valid_jk] - mean_tau_jk[k2]
                cov_matrix[k1, k2] = (N_valid - 1) / N_valid * np.sum(diff1 * diff2)
            else:
                cov_matrix[k1, k2] = np.nan
    tau_err = np.sqrt(np.maximum(np.diag(cov_matrix), 0))
    np.random.seed(42)
    v_los_shuffled = np.random.permutation(v_los)
    dv_shuffled = v_los_shuffled[i] - v_los_shuffled[j]
    tau_null = np.zeros(num_bins)
    for k in range(num_bins):
        mask = bin_masks[k]
        if np.sum(mask) > 0:
            sum_dv2 = np.sum(dv_shuffled[mask]**2)
            if sum_dv2 > 0:
                tau_null[k] = - (c / T_cmb) * np.sum(dT_obs[mask] * dv_shuffled[mask]) / sum_dv2
            else:
                tau_null[k] = np.nan
        else:
            tau_null[k] = np.nan
    results_df = pd.DataFrame({'r_arcmin': bin_centers, 'tau_obs': tau_obs, 'tau_err': tau_err, 'tau_null': tau_null})
    results_path = os.path.join(data_dir, 'pairwise_jackknife_null_results.csv')
    results_df.to_csv(results_path, index=False)
    cov_path = os.path.join(data_dir, 'pairwise_covariance_matrix.npy')
    np.save(cov_path, cov_matrix)
    print("--- Statistical Uncertainty and Robustness Tests ---")
    print("Jackknife covariance matrix saved to " + cov_path)
    print("Results saved to " + results_path)
    print("\nPairwise Estimator Results (Real vs Null):")
    for k in range(num_bins):
        sig_real = tau_obs[k] / tau_err[k] if tau_err[k] > 0 else np.nan
        sig_null = tau_null[k] / tau_err[k] if tau_err[k] > 0 else np.nan
        print("Bin " + str(k+1) + " (r=" + str(round(bin_centers[k], 1)) + " arcmin):")
        print("  Real tau: " + str(tau_obs[k]) + " +/- " + str(tau_err[k]) + " (SNR: " + str(round(sig_real, 2)) + ")")
        print("  Null tau: " + str(tau_null[k]) + " (SNR: " + str(round(sig_null, 2)) + ")")
        print("")
    valid_bins = ~np.isnan(tau_obs) & ~np.isnan(tau_err) & (tau_err > 0)
    if np.sum(valid_bins) > 0:
        chi2_diag_real = np.sum((tau_obs[valid_bins] / tau_err[valid_bins])**2)
        chi2_diag_null = np.sum((tau_null[valid_bins] / tau_err[valid_bins])**2)
        dof = np.sum(valid_bins)
        print("Diagonal Chi^2 (Real Signal): " + str(round(chi2_diag_real, 2)) + " for " + str(dof) + " degrees of freedom")
        print("Diagonal Chi^2 (Null Test): " + str(round(chi2_diag_null, 2)) + " for " + str(dof) + " degrees of freedom")
        cov_valid = cov_matrix[np.ix_(valid_bins, valid_bins)]
        try:
            inv_cov = np.linalg.pinv(cov_valid, rcond=1e-15)
            chi2_real = np.dot(tau_obs[valid_bins], np.dot(inv_cov, tau_obs[valid_bins]))
            chi2_null = np.dot(tau_null[valid_bins], np.dot(inv_cov, tau_null[valid_bins]))
            print("Full Covariance Chi^2 (Real Signal): " + str(round(chi2_real, 2)) + " for " + str(dof) + " degrees of freedom")
            print("Full Covariance Chi^2 (Null Test): " + str(round(chi2_null, 2)) + " for " + str(dof) + " degrees of freedom")
        except Exception as e:
            print("Could not compute full covariance Chi^2: " + str(e))
    mpl.rcParams['text.usetex'] = False
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.errorbar(bin_centers, tau_obs, yerr=tau_err, fmt='ro-', capsize=4, label='Real Signal (Observed Map)')
    ax.plot(bin_centers, tau_null, 'bs--', label='Null Test (Shuffled Velocities)')
    ax.axhline(0, color='black', linestyle=':', alpha=0.7)
    ax.set_xlabel('Angular Separation (arcmin)')
    ax.set_ylabel('Recovered Pairwise Estimator tau_hat(r)')
    ax.set_title('Pairwise kSZ Estimator: Real Signal vs Null Test')
    ax.legend()
    ax.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    timestamp = str(int(time.time()))
    plot_filename = 'pairwise_null_test_6_' + timestamp + '.png'
    plot_filepath = os.path.join(data_dir, plot_filename)
    plt.savefig(plot_filepath, dpi=300)
    plt.close()
    print("\nNull test plot saved to " + plot_filepath)

if __name__ == '__main__':
    data_directory = 'data/'
    obs_map_path = '/home/node/work/projects/ksz_v1/observed_map.npy'
    catalog_path = os.path.join(data_directory, 'processed_halo_catalog.csv')
    statistical_uncertainty_and_robustness(data_directory, obs_map_path, catalog_path)