# filename: codebase/step_4.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
import numpy as np
import pandas as pd
from scipy.stats import linregress

def mass_binning_and_scaling(data_dir, filtered_map_path, catalog_path):
    filtered_map = np.load(filtered_map_path)
    catalog = pd.read_csv(catalog_path)
    c = 299792.458
    T_cmb = 2.7255e6
    log_mass = catalog['log_mass'].values
    mass = catalog['mass_Msun'].values
    v_los = catalog['v_los_km_s'].values
    tau_true = catalog['tau'].values
    iy = catalog['iy_pixel'].values
    ix = catalog['ix_pixel'].values
    T_i = filtered_map[iy, ix]
    num_bins = 10
    bins = np.linspace(log_mass.min(), log_mass.max(), num_bins + 1)
    mean_log_mass = np.zeros(num_bins)
    mean_mass = np.zeros(num_bins)
    tau_recovered = np.zeros(num_bins)
    tau_true_mean = np.zeros(num_bins)
    for i in range(num_bins):
        if i == num_bins - 1:
            mask = (log_mass >= bins[i]) & (log_mass <= bins[i+1])
        else:
            mask = (log_mass >= bins[i]) & (log_mass < bins[i+1])
        mean_log_mass[i] = np.mean(log_mass[mask])
        mean_mass[i] = np.mean(mass[mask])
        tau_true_mean[i] = np.mean(tau_true[mask])
        T_i_bin = T_i[mask]
        v_los_bin = v_los[mask]
        numerator = np.sum(T_i_bin * v_los_bin)
        denominator = np.sum(v_los_bin**2)
        if denominator > 0:
            tau_recovered[i] = - (c / T_cmb) * (numerator / denominator)
        else:
            tau_recovered[i] = np.nan
    valid = ~np.isnan(tau_recovered) & (tau_recovered > 0)
    log_M_valid = np.log10(mean_mass[valid])
    log_tau_valid = np.log10(tau_recovered[valid])
    slope, intercept, r_value, p_value, std_err = linregress(log_M_valid, log_tau_valid)
    log_tau_true = np.log10(tau_true_mean)
    slope_true, intercept_true, _, _, _ = linregress(np.log10(mean_mass), log_tau_true)
    results_df = pd.DataFrame({'log_mass_bin_center': mean_log_mass, 'mean_mass': mean_mass, 'tau_recovered': tau_recovered, 'tau_true_mean': tau_true_mean})
    results_path = os.path.join(data_dir, 'tau_mass_scaling_results.csv')
    results_df.to_csv(results_path, index=False)
    print('--- Mass Binning and Scaling Relation Results ---')
    print('Recovered tau-M relation: log10(tau) = ' + str(round(slope, 4)) + ' * log10(M) + ' + str(round(intercept, 4)))
    print('True tau-M relation:      log10(tau) = ' + str(round(slope_true, 4)) + ' * log10(M) + ' + str(round(intercept_true, 4)))
    print('Theoretical slope:        0.6667')
    print('\nBin-by-bin results:')
    for i in range(num_bins):
        print('Bin ' + str(i+1) + ': logM = ' + str(round(mean_log_mass[i], 2)) + ', True tau = ' + str(tau_true_mean[i]) + ', Recovered tau = ' + str(tau_recovered[i]))
    params_df = pd.DataFrame({'parameter': ['slope_recovered', 'intercept_recovered', 'slope_true', 'intercept_true'], 'value': [slope, intercept, slope_true, intercept_true]})
    params_path = os.path.join(data_dir, 'tau_mass_scaling_parameters.csv')
    params_df.to_csv(params_path, index=False)
    print('\nResults saved to ' + results_path)
    print('Parameters saved to ' + params_path)

if __name__ == '__main__':
    data_directory = 'data/'
    filtered_map_path = os.path.join(data_directory, 'filtered_observed_map.npy')
    catalog_path = os.path.join(data_directory, 'processed_halo_catalog.csv')
    mass_binning_and_scaling(data_directory, filtered_map_path, catalog_path)