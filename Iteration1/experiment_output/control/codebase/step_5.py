# filename: codebase/step_5.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

def remove_2d_plane(cutout):
    ny, nx = cutout.shape
    y, x = np.mgrid[:ny, :nx]
    A = np.c_[x.ravel(), y.ravel(), np.ones(nx*ny)]
    C, _, _, _ = np.linalg.lstsq(A, cutout.ravel(), rcond=None)
    plane = C[0]*x + C[1]*y + C[2]
    return cutout - plane

def stacking_and_visualization(data_dir, raw_map_path, filtered_map_path, catalog_path):
    raw_map = np.load(raw_map_path)
    filt_map = np.load(filtered_map_path)
    catalog = pd.read_csv(catalog_path)
    c = 299792.458
    T_cmb = 2.7255e6
    log_mass = catalog['log_mass'].values
    mass = catalog['mass_Msun'].values
    v_los = catalog['v_los_km_s'].values
    tau_true = catalog['tau'].values
    iy = catalog['iy_pixel'].values
    ix = catalog['ix_pixel'].values
    num_bins = 10
    bins = np.linspace(log_mass.min(), log_mass.max(), num_bins + 1)
    hw = 15
    grid_size = 2 * hw + 1
    stacks_raw = np.zeros((num_bins, grid_size, grid_size))
    stacks_filt = np.zeros((num_bins, grid_size, grid_size))
    mean_mass = np.zeros(num_bins)
    tau_true_mean = np.zeros(num_bins)
    tau_recovered = np.zeros(num_bins)
    x_deg = catalog['ra'].values
    y_deg = catalog['dec'].values
    x_idx = np.floor(x_deg / 2.0).astype(int)
    y_idx = np.floor(y_deg / 2.0).astype(int)
    x_idx[x_idx == 5] = 4
    y_idx[y_idx == 5] = 4
    region_id = x_idx + 5 * y_idx
    N_jk = 25
    tau_jk = np.zeros((num_bins, N_jk))
    T_i_filt = filt_map[iy, ix]
    for b in range(num_bins):
        if b == num_bins - 1:
            mask_b = (log_mass >= bins[b]) & (log_mass <= bins[b+1])
        else:
            mask_b = (log_mass >= bins[b]) & (log_mass < bins[b+1])
        mean_mass[b] = np.mean(mass[mask_b])
        tau_true_mean[b] = np.mean(tau_true[mask_b])
        v_los_b = v_los[mask_b]
        T_i_b = T_i_filt[mask_b]
        dv = v_los_b - np.mean(v_los_b)
        dT = T_i_b - np.mean(T_i_b)
        num = np.sum(dT * dv)
        den = np.sum(dv**2)
        if den > 0:
            tau_recovered[b] = - (c / T_cmb) * (num / den)
        else:
            tau_recovered[b] = np.nan
        for k in range(N_jk):
            mask_jk = mask_b & (region_id != k)
            v_los_jk = v_los[mask_jk]
            T_i_jk = T_i_filt[mask_jk]
            dv_jk = v_los_jk - np.mean(v_los_jk)
            dT_jk = T_i_jk - np.mean(T_i_jk)
            num_jk = np.sum(dT_jk * dv_jk)
            den_jk = np.sum(dv_jk**2)
            if den_jk > 0:
                tau_jk[b, k] = - (c / T_cmb) * (num_jk / den_jk)
            else:
                tau_jk[b, k] = np.nan
        idx_b = np.where(mask_b)[0]
        valid_idx = [i for i in idx_b if hw <= ix[i] < 1200 - hw and hw <= iy[i] < 1200 - hw]
        if len(valid_idx) > 0:
            v_valid = v_los[valid_idx]
            v_valid_mean = np.mean(v_valid)
            dv_valid = v_valid - v_valid_mean
            sigma_v = np.std(v_valid)
            sum_dv2 = np.sum(dv_valid**2)
            if sum_dv2 > 0:
                for i, idx_i in enumerate(valid_idx):
                    w_i = - dv_valid[i] / sum_dv2 * sigma_v
                    raw_cutout = raw_map[iy[idx_i]-hw:iy[idx_i]+hw+1, ix[idx_i]-hw:ix[idx_i]+hw+1]
                    raw_cutout_clean = remove_2d_plane(raw_cutout)
                    stacks_raw[b] += w_i * raw_cutout_clean
                    filt_cutout = filt_map[iy[idx_i]-hw:iy[idx_i]+hw+1, ix[idx_i]-hw:ix[idx_i]+hw+1]
                    filt_cutout_clean = filt_cutout - np.mean(filt_cutout)
                    stacks_filt[b] += w_i * filt_cutout_clean
    tau_err = np.zeros(num_bins)
    for b in range(num_bins):
        valid_jk = ~np.isnan(tau_jk[b])
        N_valid = np.sum(valid_jk)
        if N_valid > 1:
            mean_tau_jk = np.mean(tau_jk[b][valid_jk])
            var_jk = (N_valid - 1) / N_valid * np.sum((tau_jk[b][valid_jk] - mean_tau_jk)**2)
            tau_err[b] = np.sqrt(var_jk)
        else:
            tau_err[b] = np.nan
    print('--- Stacking and Jackknife Results ---')
    print('Note: Stacking weights are designed to remove monopole leakage and yield maps in uK.')
    print('Note: Raw map cutouts have a 2D plane removed to suppress large-scale CMB gradients.')
    print('Note: Filtered map stacks may show a \'Mexican hat\' profile due to the matched filter suppressing low-ell modes.')
    for b in range(num_bins):
        print('Bin ' + str(b+1) + ' (Mean Mass: ' + str(round(mean_mass[b], 2)) + ' M_sun):')
        print('  Recovered tau: ' + str(tau_recovered[b]))
        print('  Jackknife error: ' + str(tau_err[b]))
        print('  True tau: ' + str(tau_true_mean[b]))
        print('  SNR (Recovered / Error): ' + str(tau_recovered[b] / tau_err[b] if tau_err[b] > 0 else np.nan))
        print('')
    mpl.rcParams['text.usetex'] = False
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    extent = [-(hw + 0.5) * 0.5, (hw + 0.5) * 0.5, -(hw + 0.5) * 0.5, (hw + 0.5) * 0.5]
    im00 = axes[0, 0].imshow(stacks_raw[0], origin='lower', extent=extent, cmap='viridis')
    axes[0, 0].set_title('Lowest Mass Bin - Raw Map (Plane Removed)\nMean Mass: ' + str(round(mean_mass[0], 2)) + ' M_sun')
    fig.colorbar(im00, ax=axes[0, 0], label='Temperature (uK)')
    im01 = axes[0, 1].imshow(stacks_filt[0], origin='lower', extent=extent, cmap='viridis')
    axes[0, 1].set_title('Lowest Mass Bin - Filtered Map\nMean Mass: ' + str(round(mean_mass[0], 2)) + ' M_sun')
    fig.colorbar(im01, ax=axes[0, 1], label='Temperature (uK)')
    im10 = axes[1, 0].imshow(stacks_raw[-1], origin='lower', extent=extent, cmap='viridis')
    axes[1, 0].set_title('Highest Mass Bin - Raw Map (Plane Removed)\nMean Mass: ' + str(round(mean_mass[-1], 2)) + ' M_sun')
    fig.colorbar(im10, ax=axes[1, 0], label='Temperature (uK)')
    im11 = axes[1, 1].imshow(stacks_filt[-1], origin='lower', extent=extent, cmap='viridis')
    axes[1, 1].set_title('Highest Mass Bin - Filtered Map\nMean Mass: ' + str(round(mean_mass[-1], 2)) + ' M_sun')
    fig.colorbar(im11, ax=axes[1, 1], label='Temperature (uK)')
    for ax in axes.flat:
        ax.set_xlabel('dx (arcmin)')
        ax.set_ylabel('dy (arcmin)')
    plt.tight_layout()
    timestamp = str(int(time.time()))
    plot_filename_1 = 'stacked_cutouts_5_' + timestamp + '.png'
    plot_filepath_1 = os.path.join(data_dir, plot_filename_1)
    plt.savefig(plot_filepath_1, dpi=300)
    plt.close()
    print('Stacked cutouts plot saved to ' + plot_filepath_1)
    fig, ax = plt.subplots(figsize=(9, 7))
    M_range = np.logspace(13, 15, 100)
    tau_theory = 2e-3 * (M_range / 1e14)**(2/3)
    ax.plot(M_range, tau_theory, 'k-', linewidth=2, label='Theoretical M^(2/3) scaling')
    ax.plot(mean_mass, tau_true_mean, 'bs-', markersize=6, label='True Mean tau (Catalog)')
    ax.errorbar(mean_mass, tau_recovered, yerr=tau_err, fmt='ro', capsize=4, markersize=6, label='Recovered tau')
    ax.set_xscale('log')
    ax.set_yscale('symlog', linthresh=1e-4, linscale=1)
    yticks = [-1e-2, -1e-3, -1e-4, 0, 1e-4, 1e-3, 1e-2]
    ax.set_yticks(yticks)
    ax.set_yticklabels(['-0.01', '-0.001', '-0.0001', '0', '0.0001', '0.001', '0.01'])
    ax.set_xlabel('Halo Mass M_500 (M_sun)')
    ax.set_ylabel('Thomson Optical Depth tau')
    ax.set_title('tau-M Scaling Relation with Jackknife Uncertainties')
    ax.legend()
    ax.grid(True, which='both', ls='--', alpha=0.5)
    plt.tight_layout()
    plot_filename_2 = 'tau_mass_scaling_jk_5_' + timestamp + '.png'
    plot_filepath_2 = os.path.join(data_dir, plot_filename_2)
    plt.savefig(plot_filepath_2, dpi=300)
    plt.close()
    print('tau-M scaling plot saved to ' + plot_filepath_2)

if __name__ == '__main__':
    data_directory = 'data/'
    raw_map_path = '/home/node/work/projects/ksz_v1/observed_map.npy'
    filtered_map_path = os.path.join(data_directory, 'filtered_observed_map.npy')
    catalog_path = os.path.join(data_directory, 'processed_halo_catalog.csv')
    stacking_and_visualization(data_directory, raw_map_path, filtered_map_path, catalog_path)