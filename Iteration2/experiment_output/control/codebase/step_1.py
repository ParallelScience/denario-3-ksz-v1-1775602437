# filename: codebase/step_1.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
import numpy as np
import matplotlib.pyplot as plt
import os
import scipy.stats
import time

plt.rcParams['text.usetex'] = False

def get_1d_ps(P2d, ell2d, pixel_size_rad, bins=100):
    ell_max = np.pi / pixel_size_rad
    ell_bins = np.linspace(100, ell_max, bins)
    ell_1d = 0.5 * (ell_bins[1:] + ell_bins[:-1])
    P_1d, _, _ = scipy.stats.binned_statistic(ell2d.ravel(), P2d.ravel(), statistic='mean', bins=ell_bins)
    return ell_1d, P_1d

if __name__ == '__main__':
    obs_map_path = '/home/node/work/projects/ksz_v1/observed_map.npy'
    cmb_map_path = '/home/node/work/projects/ksz_v1/cmb_map_truth.npy'
    ksz_map_path = '/home/node/work/projects/ksz_v1/ksz_map_truth.npy'
    noise_map_path = '/home/node/work/projects/ksz_v1/noise_map.npy'
    halo_catalog_path = '/home/node/work/projects/ksz_v1/halo_catalog.npz'
    obs_map = np.load(obs_map_path)
    cmb_map = np.load(cmb_map_path)
    ksz_map = np.load(ksz_map_path)
    noise_map = np.load(noise_map_path)
    catalog = np.load(halo_catalog_path)
    rms_obs = np.std(obs_map)
    rms_cmb = np.std(cmb_map)
    rms_ksz = np.std(ksz_map)
    rms_noise = np.std(noise_map)
    print('--- Map RMS Values ---')
    print('RMS of Observed Map: ' + str(np.round(rms_obs, 4)) + ' uK')
    print('RMS of CMB Truth Map: ' + str(np.round(rms_cmb, 4)) + ' uK')
    print('RMS of kSZ Truth Map: ' + str(np.round(rms_ksz, 4)) + ' uK')
    print('RMS of Noise Map: ' + str(np.round(rms_noise, 4)) + ' uK')
    N = 1200
    pixel_size_arcmin = 0.5
    pixel_size_rad = pixel_size_arcmin * np.pi / (180 * 60)
    omega_pix = pixel_size_rad**2
    freqs = np.fft.fftfreq(N, d=pixel_size_rad) * 2 * np.pi
    lx, ly = np.meshgrid(freqs, freqs, indexing='ij')
    ell = np.sqrt(lx**2 + ly**2)
    ell[0, 0] = 1e-5
    Dl = 6000 * (ell/200)**2 / (1 + (ell/200)**1.2)**2.5 * np.exp(-(ell/4500)**2)
    Cl_cmb_theory = 2 * np.pi * Dl / (ell * (ell + 1))
    Cl_cmb_theory[0, 0] = 0
    fwhm_rad = 1.4 * np.pi / (180 * 60)
    sigma_b = fwhm_rad / np.sqrt(8 * np.log(2))
    Bl = np.exp(-0.5 * ell**2 * sigma_b**2)
    Cl_cmb_map = Cl_cmb_theory * Bl**2
    P_noise = np.abs(np.fft.fft2(noise_map))**2 * omega_pix / (N * N)
    P_ksz = np.abs(np.fft.fft2(ksz_map))**2 * omega_pix / (N * N)
    P_cmb_empirical = np.abs(np.fft.fft2(cmb_map))**2 * omega_pix / (N * N)
    Cl_noise = 20**2 * omega_pix
    W_ell = Cl_cmb_map / (Cl_cmb_map + Cl_noise + P_ksz)
    W_ell[0, 0] = 0
    obs_fft = np.fft.fft2(obs_map)
    cmb_wiener_fft = W_ell * obs_fft
    cmb_wiener = np.real(np.fft.ifft2(cmb_wiener_fft))
    cleaned_map = obs_map - cmb_wiener
    rms_cleaned = np.std(cleaned_map)
    print('RMS of Cleaned Residual Map: ' + str(np.round(rms_cleaned, 4)) + ' uK')
    perfect_cleaned_map = obs_map - cmb_map
    rms_perfect_cleaned = np.std(perfect_cleaned_map)
    print('RMS of Perfect Cleaned Map (Obs - True CMB): ' + str(np.round(rms_perfect_cleaned, 4)) + ' uK')
    data_dir = 'data/'
    cleaned_map_path = os.path.join(data_dir, 'wiener_cleaned_map.npy')
    np.save(cleaned_map_path, cleaned_map)
    print('Wiener cleaned map saved to ' + cleaned_map_path)
    perfect_cleaned_map_path = os.path.join(data_dir, 'perfect_cleaned_map.npy')
    np.save(perfect_cleaned_map_path, perfect_cleaned_map)
    print('Perfect cleaned map saved to ' + perfect_cleaned_map_path)
    filter_path = os.path.join(data_dir, 'wiener_filter.npy')
    np.save(filter_path, W_ell)
    print('Wiener filter saved to ' + filter_path)
    beam_path = os.path.join(data_dir, 'beam_transfer_function.npy')
    np.save(beam_path, Bl)
    print('Beam transfer function saved to ' + beam_path)
    np.save(os.path.join(data_dir, 'P_obs_2d.npy'), np.abs(obs_fft)**2 * omega_pix / (N * N))
    np.save(os.path.join(data_dir, 'P_cleaned_2d.npy'), np.abs(np.fft.fft2(cleaned_map))**2 * omega_pix / (N * N))
    ell_1d, Cl_cmb_1d = get_1d_ps(Cl_cmb_map, ell, pixel_size_rad)
    _, P_cmb_emp_1d = get_1d_ps(P_cmb_empirical, ell, pixel_size_rad)
    _, Cl_noise_1d = get_1d_ps(np.full_like(ell, Cl_noise), ell, pixel_size_rad)
    _, P_ksz_1d = get_1d_ps(P_ksz, ell, pixel_size_rad)
    _, W_1d = get_1d_ps(W_ell, ell, pixel_size_rad)
    timestamp = int(time.time())
    plot_filename = os.path.join(data_dir, 'wiener_filter_diagnostic_1_' + str(timestamp) + '.png')
    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax1.plot(ell_1d, Cl_cmb_1d, label='CMB (Analyt.+Beam)', color='blue', linewidth=2)
    ax1.plot(ell_1d, P_cmb_emp_1d, label='CMB (Empirical)', color='cyan', linestyle='--', alpha=0.7)
    ax1.plot(ell_1d, Cl_noise_1d, label='Noise (Analyt.)', color='red', linestyle='--', linewidth=2)
    ax1.plot(ell_1d, P_ksz_1d, label='kSZ (Empirical)', color='green', linewidth=2)
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.set_xlabel('Multipole ell')
    ax1.set_ylabel('Power Spectrum C_ell [uK^2 sr]')
    ax1.set_title('Power Spectra and Wiener Filter Transfer Function')
    ax1.grid(True, which='both', ls='--', alpha=0.5)
    ax1.legend(loc='lower left')
    ax2 = ax1.twinx()
    ax2.plot(ell_1d, W_1d, label='Wiener Filter W(ell)', color='purple', linestyle='-.', linewidth=2)
    ax2.set_yscale('log')
    ax2.set_ylabel('Filter Transfer Function W(ell)')
    ax2.set_ylim(1e-6, 2)
    ax2.legend(loc='upper right')
    fig.tight_layout()
    fig.savefig(plot_filename, dpi=300)
    print('Diagnostic plot saved to ' + plot_filename)