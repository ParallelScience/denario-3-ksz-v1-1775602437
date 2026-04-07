# filename: codebase/step_3.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

if __name__ == '__main__':
    plt.rcParams['text.usetex'] = False
    catalog_path = '/home/node/work/projects/ksz_v1/halo_catalog.npz'
    residual_map_path = 'data/residual_map.npy'
    catalog = np.load(catalog_path)
    residual_map = np.load(residual_map_path)
    log_mass = catalog['log_mass']
    delta_T_ksz = catalog['delta_T_ksz_uK']
    ix = catalog['ix_pixel']
    iy = catalog['iy_pixel']
    num_bins = 10
    bins = np.linspace(13, 15, num_bins + 1)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    mean_delta_T = np.zeros(num_bins)
    se_delta_T = np.zeros(num_bins)
    mean_mass = np.zeros(num_bins)
    halo_counts = np.zeros(num_bins, dtype=int)
    r_pix = 10
    stacked_stamps = np.zeros((num_bins, 2*r_pix+1, 2*r_pix+1))
    print('--- Mass Binning and Stacking ---')
    print('Number of bins: ' + str(num_bins))
    print('Mass range: 10^' + str(bins[0]) + ' to 10^' + str(bins[-1]) + ' M_sun')
    print('Stamp size: ' + str(2*r_pix+1) + 'x' + str(2*r_pix+1) + ' pixels')
    for i in range(num_bins):
        if i == num_bins - 1:
            mask = (log_mass >= bins[i]) & (log_mass <= bins[i+1])
        else:
            mask = (log_mass >= bins[i]) & (log_mass < bins[i+1])
        idx = np.where(mask)[0]
        halo_counts[i] = len(idx)
        if len(idx) > 0:
            mean_mass[i] = np.mean(log_mass[idx])
            mean_delta_T[i] = np.mean(delta_T_ksz[idx])
            se_delta_T[i] = np.std(delta_T_ksz[idx]) / np.sqrt(len(idx))
            stamps = []
            for j in idx:
                cx, cy = ix[j], iy[j]
                if cx - r_pix >= 0 and cx + r_pix < residual_map.shape[1] and cy - r_pix >= 0 and cy + r_pix < residual_map.shape[0]:
                    stamp = residual_map[cy - r_pix : cy + r_pix + 1, cx - r_pix : cx + r_pix + 1]
                    stamps.append(stamp)
            if len(stamps) > 0:
                stamps = np.array(stamps)
                stacked_stamps[i] = np.mean(stamps, axis=0)
        print('Bin ' + str(i+1) + ' (' + str(np.round(bins[i], 1)) + '-' + str(np.round(bins[i+1], 1)) + '): ' + str(halo_counts[i]) + ' halos, ' + 'True dT = ' + str(np.round(mean_delta_T[i], 4)) + ' +/- ' + str(np.round(se_delta_T[i], 4)) + ' uK')
    output_npz = 'data/bin_statistics.npz'
    np.savez(output_npz, bins=bins, bin_centers=bin_centers, mean_mass=mean_mass, halo_counts=halo_counts, mean_delta_T=mean_delta_T, se_delta_T=se_delta_T, stacked_stamps=stacked_stamps)
    print('\nSaved bin statistics and stacked stamps to ' + output_npz)
    fig, axes = plt.subplots(2, 5, figsize=(18, 7))
    axes = axes.flatten()
    vmax = np.max(np.abs(stacked_stamps))
    if vmax == 0: vmax = 1
    for i in range(num_bins):
        ax = axes[i]
        im = ax.imshow(stacked_stamps[i], cmap='RdBu_r', vmin=-vmax, vmax=vmax, origin='lower', extent=[-r_pix*0.5, r_pix*0.5, -r_pix*0.5, r_pix*0.5])
        title_str = 'logM: ' + str(np.round(bins[i], 1)) + '-' + str(np.round(bins[i+1], 1)) + '\n'
        title_str += 'True dT: ' + str(np.round(mean_delta_T[i], 3)) + ' ± ' + str(np.round(se_delta_T[i], 3)) + ' μK'
        ax.set_title(title_str, fontsize=10)
        if i >= 5:
            ax.set_xlabel('Delta RA (arcmin)')
        if i % 5 == 0:
            ax.set_ylabel('Delta Dec (arcmin)')
    fig.tight_layout(rect=[0, 0, 0.92, 0.95])
    cbar_ax = fig.add_axes([0.93, 0.15, 0.015, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label('Mean Residual (μK)')
    fig.suptitle('Stacked Residual Stamps by Mass Bin', fontsize=16)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    plot_filename = 'data/stacked_stamps_1_' + timestamp + '.png'
    plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
    print('Stacked stamps plot saved to ' + plot_filename)