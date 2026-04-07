# filename: codebase/step_2.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

def construct_and_apply_matched_filter(data_dir, obs_map_path, cmb_truth_path, noise_truth_path):
    obs_map = np.load(obs_map_path)
    cmb_map = np.load(cmb_truth_path)
    noise_map = np.load(noise_truth_path)
    N = obs_map.shape[0]
    pixel_size_arcmin = 0.5
    pixel_size_rad = pixel_size_arcmin * np.pi / (180.0 * 60.0)
    omega_pix = pixel_size_rad**2
    freq_x = np.fft.fftfreq(N, d=pixel_size_rad)
    freq_y = np.fft.fftfreq(N, d=pixel_size_rad)
    fx, fy = np.meshgrid(freq_x, freq_y)
    ell_2d = 2 * np.pi * np.sqrt(fx**2 + fy**2)
    ell_2d[0, 0] = 1e-5
    D_ell_2d = 6000.0 * (ell_2d/200.0)**2 / (1.0 + (ell_2d/200.0)**1.2)**2.5 * np.exp(-(ell_2d/4500.0)**2)
    C_ell_CMB_unconv_2d = D_ell_2d * 2.0 * np.pi / (ell_2d * (ell_2d + 1.0))
    C_ell_CMB_unconv_2d[0, 0] = 0.0
    fwhm_arcmin = 1.4
    fwhm_rad = fwhm_arcmin * np.pi / (180.0 * 60.0)
    sigma_beam_rad = fwhm_rad / np.sqrt(8.0 * np.log(2.0))
    B_ell_2d = np.exp(-0.5 * ell_2d**2 * sigma_beam_rad**2)
    C_ell_CMB_conv_2d = C_ell_CMB_unconv_2d * B_ell_2d**2
    rms_noise_uK = 20.0
    C_ell_noise_2d = np.full_like(ell_2d, rms_noise_uK**2 * omega_pix)
    P_total_2d = C_ell_CMB_conv_2d + C_ell_noise_2d
    Psi_ell_2d = B_ell_2d / P_total_2d
    norm_factor = np.mean(Psi_ell_2d * B_ell_2d)
    Psi_ell_2d_norm = Psi_ell_2d / norm_factor
    obs_fft = np.fft.fft2(obs_map)
    filtered_fft = obs_fft * Psi_ell_2d_norm
    filtered_map = np.fft.ifft2(filtered_fft).real
    filtered_map_path = os.path.join(data_dir, 'filtered_observed_map.npy')
    np.save(filtered_map_path, filtered_map)
    def get_empirical_cell(map_2d):
        map_fft = np.fft.fft2(map_2d)
        P2D = np.abs(map_fft)**2 * omega_pix / (N**2)
        ell_bins = np.arange(100, 10000, 100)
        ell_centers = 0.5 * (ell_bins[1:] + ell_bins[:-1])
        P1D = np.zeros_like(ell_centers)
        for i in range(len(ell_centers)):
            mask = (ell_2d >= ell_bins[i]) & (ell_2d < ell_bins[i+1])
            if np.any(mask):
                P1D[i] = np.mean(P2D[mask])
            else:
                P1D[i] = np.nan
        return ell_centers, P1D
    ell_emp, Cell_obs = get_empirical_cell(obs_map)
    _, Cell_cmb = get_empirical_cell(cmb_map)
    _, Cell_noise = get_empirical_cell(noise_map)
    ell_1d = ell_emp
    D_ell_1d = 6000.0 * (ell_1d/200.0)**2 / (1.0 + (ell_1d/200.0)**1.2)**2.5 * np.exp(-(ell_1d/4500.0)**2)
    C_ell_CMB_unconv_1d = D_ell_1d * 2.0 * np.pi / (ell_1d * (ell_1d + 1.0))
    B_ell_1d = np.exp(-0.5 * ell_1d**2 * sigma_beam_rad**2)
    C_ell_CMB_conv_1d = C_ell_CMB_unconv_1d * B_ell_1d**2
    C_ell_noise_1d = np.full_like(ell_1d, rms_noise_uK**2 * omega_pix)
    P_total_1d = C_ell_CMB_conv_1d + C_ell_noise_1d
    Psi_ell_1d = B_ell_1d / P_total_1d
    Psi_ell_1d_norm = Psi_ell_1d / norm_factor
    mpl.rcParams['text.usetex'] = False
    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax1.plot(ell_emp, Cell_obs, label='Observed Map (Empirical)', color='black', alpha=0.7)
    ax1.plot(ell_emp, Cell_cmb, label='CMB Truth (Empirical)', color='blue', alpha=0.7)
    ax1.plot(ell_emp, Cell_noise, label='Noise (Empirical)', color='red', alpha=0.7)
    ax1.plot(ell_1d, C_ell_CMB_conv_1d, label='CMB Convolved (Theory)', color='cyan', linestyle='--')
    ax1.plot(ell_1d, C_ell_noise_1d, label='Noise (Theory)', color='orange', linestyle='--')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel('Multipole ell')
    ax1.set_ylabel('Power Spectrum C_ell (uK^2 sr)')
    ax1.set_title('Power Spectra and Matched Filter Transfer Function')
    ax1.grid(True, which='both', ls='--', alpha=0.5)
    ax2 = ax1.twinx()
    ax2.plot(ell_1d, Psi_ell_1d_norm, label='Matched Filter Psi(ell)', color='green', linewidth=2)
    ax2.set_ylabel('Filter Transfer Function (dimensionless)')
    ax2.set_yscale('log')
    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper right')
    plt.tight_layout()
    timestamp = str(int(time.time()))
    plot_filename = 'power_spectra_filter_2_' + timestamp + '.png'
    plot_filepath = os.path.join(data_dir, plot_filename)
    plt.savefig(plot_filepath, dpi=300)
    plt.close()

if __name__ == '__main__':
    data_directory = 'data/'
    obs_map_path = '/home/node/work/projects/ksz_v1/observed_map.npy'
    cmb_truth_path = '/home/node/work/projects/ksz_v1/cmb_map_truth.npy'
    noise_truth_path = '/home/node/work/projects/ksz_v1/noise_map.npy'
    construct_and_apply_matched_filter(data_directory, obs_map_path, cmb_truth_path, noise_truth_path)