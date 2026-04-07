# filename: codebase/step_1.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
import numpy as np
import os

if __name__ == '__main__':
    observed_map_path = '/home/node/work/projects/ksz_v1/observed_map.npy'
    cmb_map_path = '/home/node/work/projects/ksz_v1/cmb_map_truth.npy'
    ksz_map_path = '/home/node/work/projects/ksz_v1/ksz_map_truth.npy'
    noise_map_path = '/home/node/work/projects/ksz_v1/noise_map.npy'
    catalog_path = '/home/node/work/projects/ksz_v1/halo_catalog.npz'
    observed_map = np.load(observed_map_path)
    cmb_map = np.load(cmb_map_path)
    ksz_map = np.load(ksz_map_path)
    noise_map = np.load(noise_map_path)
    catalog = np.load(catalog_path)
    print('--- Map Dimensions ---')
    print('Observed map shape: ' + str(observed_map.shape))
    print('CMB map shape: ' + str(cmb_map.shape))
    print('kSZ map shape: ' + str(ksz_map.shape))
    print('Noise map shape: ' + str(noise_map.shape))
    ix = catalog['ix_pixel']
    iy = catalog['iy_pixel']
    print('\n--- Pixel Indexing Verification ---')
    print('ix_pixel min/max: ' + str(ix.min()) + ' / ' + str(ix.max()))
    print('iy_pixel min/max: ' + str(iy.min()) + ' / ' + str(iy.max()))
    if ix.min() >= 0 and ix.max() < observed_map.shape[1] and iy.min() >= 0 and iy.max() < observed_map.shape[0]:
        print('Pixel indices are within map bounds.')
    else:
        print('WARNING: Pixel indices are out of bounds!')
    print('\n--- Map RMS Values ---')
    print('Observed map RMS: ' + str(np.round(np.std(observed_map), 2)) + ' uK (Expected: ~60 uK)')
    print('CMB map RMS: ' + str(np.round(np.std(cmb_map), 2)) + ' uK (Expected: ~57 uK)')
    print('kSZ map RMS: ' + str(np.round(np.std(ksz_map), 4)) + ' uK (Expected: ~0.14 uK)')
    print('Noise map RMS: ' + str(np.round(np.std(noise_map), 2)) + ' uK (Expected: 20 uK)')
    print('\n--- Halo Catalog Summary Statistics ---')
    keys_to_summarize = ['mass_Msun', 'tau', 'v_los_km_s', 'delta_T_ksz_uK']
    for key in keys_to_summarize:
        data = catalog[key]
        print(key + ':')
        print('  Min:  ' + str(np.min(data)))
        print('  Max:  ' + str(np.max(data)))
        print('  Mean: ' + str(np.mean(data)))
    N = 1200
    L_deg = 10.0
    L_rad = L_deg * np.pi / 180.0
    dx_rad = L_rad / N
    lx = np.fft.fftfreq(N, d=dx_rad) * 2 * np.pi
    ly = np.fft.fftfreq(N, d=dx_rad) * 2 * np.pi
    lx_2d, ly_2d = np.meshgrid(lx, ly)
    l_2d = np.sqrt(lx_2d**2 + ly_2d**2)
    l_2d_safe = np.copy(l_2d)
    l_2d_safe[0, 0] = 1e-5
    fwhm_arcmin = 1.4
    fwhm_rad = fwhm_arcmin * np.pi / (180.0 * 60.0)
    sigma_b = fwhm_rad / np.sqrt(8 * np.log(2))
    B_l_2d = np.exp(-l_2d_safe * (l_2d_safe + 1) * sigma_b**2 / 2.0)
    B_l_2d[0, 0] = 1.0
    D_l = 6000.0 * (l_2d_safe / 200.0)**2 / (1.0 + (l_2d_safe / 200.0)**1.2)**2.5 * np.exp(-(l_2d_safe / 4500.0)**2)
    Cl_CMB_2d = 2 * np.pi * D_l / (l_2d_safe * (l_2d_safe + 1))
    Cl_CMB_2d[0, 0] = 0.0
    sigma_noise_uK = 20.0
    omega_pix = dx_rad**2
    Cl_noise_2d = np.full_like(l_2d, sigma_noise_uK**2 * omega_pix)
    output_path = 'data/power_spectra_and_beam.npz'
    np.savez(output_path, l_2d=l_2d, B_l_2d=B_l_2d, Cl_CMB_2d=Cl_CMB_2d, Cl_noise_2d=Cl_noise_2d)
    print('\nSaved power spectra and beam function to ' + output_path)