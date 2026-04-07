# filename: codebase/step_3.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

def compute_pairwise_estimator(data_dir, obs_map_path, ksz_truth_path, catalog_path):
    obs_map = np.load(obs_map_path)
    ksz_truth = np.load(ksz_truth_path)
    catalog = pd.read_csv(catalog_path)
    iy = catalog['iy_pixel'].values
    ix = catalog['ix_pixel'].values
    T_obs = obs_map[iy, ix]
    T_truth = ksz_truth[iy, ix]
    v_los = catalog['v_los_km_s'].values
    x_deg = catalog['ra'].values
    y_deg = catalog['dec'].values
    N = len(catalog)
    i, j = np.triu_indices(N, k=1)
    dist_deg = np.sqrt((x_deg[i] - x_deg[j])**2 + (y_deg[i] - y_deg[j])**2)
    dist_arcmin = dist_deg * 60.0
    dv = v_los[i] - v_los[j]
    dT_obs = T_obs[i] - T_obs[j]
    dT_truth = T_truth[i] - T_truth[j]
    bins = np.arange(0, 601, 30)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    tau_obs = np.zeros(len(bins)-1)
    tau_truth = np.zeros(len(bins)-1)
    c = 299792.458
    T_cmb = 2.7255e6
    for k in range(len(bins)-1):
        mask = (dist_arcmin >= bins[k]) & (dist_arcmin < bins[k+1])
        if np.sum(mask) > 0:
            sum_dv2 = np.sum(dv[mask]**2)
            if sum_dv2 > 0:
                tau_obs[k] = - (c / T_cmb) * np.sum(dT_obs[mask] * dv[mask]) / sum_dv2
                tau_truth[k] = - (c / T_cmb) * np.sum(dT_truth[mask] * dv[mask]) / sum_dv2
            else:
                tau_obs[k] = np.nan
                tau_truth[k] = np.nan
        else:
            tau_obs[k] = np.nan
            tau_truth[k] = np.nan
    results_df = pd.DataFrame({'r_arcmin': bin_centers, 'tau_obs': tau_obs, 'tau_truth': tau_truth})
    results_path = os.path.join(data_dir, 'pairwise_estimator_results.csv')
    results_df.to_csv(results_path, index=False)
    print('--- Pairwise Estimator Summary ---')
    mean_tau_truth = np.nanmean(tau_truth)
    mean_tau_obs = np.nanmean(tau_obs)
    true_mean_tau = catalog['tau'].mean()
    print('Mean recovered tau (Truth Map, all bins): ' + str(mean_tau_truth))
    print('Mean recovered tau (Observed Map, all bins): ' + str(mean_tau_obs))
    print('True mean tau from catalog: ' + str(true_mean_tau))
    mpl.rcParams['text.usetex'] = False
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(bin_centers, tau_truth, marker='o', linestyle='-', color='blue', label='kSZ Truth Map')
    ax.plot(bin_centers, tau_obs, marker='s', linestyle='--', color='red', label='Observed Map (Raw)')
    ax.axhline(true_mean_tau, color='black', linestyle=':', label='True Mean tau (Catalog)')
    dilution_factor = 0.25 / (2.0 * np.pi * (1.4 / np.sqrt(8.0 * np.log(2.0)))**2)
    expected_tau = true_mean_tau * dilution_factor
    ax.axhline(expected_tau, color='green', linestyle='-.', label='Expected Beam-Diluted tau')
    ax.set_xlabel('Angular Separation (arcmin)')
    ax.set_ylabel('Recovered Pairwise Estimator tau_hat(r)')
    ax.set_title('Pairwise kSZ Estimator: Observed vs Truth Map')
    ax.legend()
    ax.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plot_filename = 'pairwise_estimator_' + str(int(time.time())) + '.png'
    plot_filepath = os.path.join(data_dir, plot_filename)
    plt.savefig(plot_filepath, dpi=300)
    plt.close()
    print('Pairwise estimator plot saved to ' + plot_filepath)

if __name__ == '__main__':
    data_directory = 'data/'
    obs_map_path = '/home/node/work/projects/ksz_v1/observed_map.npy'
    ksz_truth_path = '/home/node/work/projects/ksz_v1/ksz_map_truth.npy'
    catalog_path = os.path.join(data_directory, 'halo_catalog.csv')
    compute_pairwise_estimator(data_directory, obs_map_path, ksz_truth_path, catalog_path)