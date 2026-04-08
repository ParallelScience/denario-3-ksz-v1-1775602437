# filename: codebase/step_3.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
import numpy as np
import scipy.stats
import os

def main():
    data_dir = 'data/'
    halo_catalog_path = '/home/node/work/projects/ksz_v1/halo_catalog.npz'
    catalog = np.load(halo_catalog_path)
    ra = catalog['ra']
    dec = catalog['dec']
    mass = catalog['mass_Msun']
    v_los_true = catalog['v_los_km_s']
    N = 1200
    density = np.zeros((N, N))
    x = ra * 120.0
    y = dec * 120.0
    x0 = np.floor(x).astype(int)
    y0 = np.floor(y).astype(int)
    dx = x - x0
    dy = y - y0
    for i in range(len(x)):
        xi = x0[i] % N
        yi = y0[i] % N
        xi1 = (xi + 1) % N
        yi1 = (yi + 1) % N
        m = mass[i]
        density[yi, xi] += m * (1 - dx[i]) * (1 - dy[i])
        density[yi, xi1] += m * dx[i] * (1 - dy[i])
        density[yi1, xi] += m * (1 - dx[i]) * dy[i]
        density[yi1, xi1] += m * dx[i] * dy[i]
    mean_dens = np.mean(density)
    delta = density / mean_dens - 1.0
    delta_k = np.fft.fft2(delta)
    pixel_size_rad = 0.5 * np.pi / (180 * 60)
    kx = np.fft.fftfreq(N, d=pixel_size_rad) * 2 * np.pi
    ky = np.fft.fftfreq(N, d=pixel_size_rad) * 2 * np.pi
    KX, KY = np.meshgrid(kx, ky)
    K2 = KX**2 + KY**2
    K2[0, 0] = 1e-10
    v_k = 1j * KX / K2 * delta_k
    v_k[0, 0] = 0
    v_rec_map = np.real(np.fft.ifft2(v_k))
    v_rec_map = (v_rec_map - np.mean(v_rec_map)) / np.std(v_rec_map) * 300.0
    v_rec = np.zeros(len(x))
    for i in range(len(x)):
        xi = x0[i] % N
        yi = y0[i] % N
        xi1 = (xi + 1) % N
        yi1 = (yi + 1) % N
        v_rec[i] = (v_rec_map[yi, xi] * (1 - dx[i]) * (1 - dy[i]) + v_rec_map[yi, xi1] * dx[i] * (1 - dy[i]) + v_rec_map[yi1, xi] * (1 - dx[i]) * dy[i] + v_rec_map[yi1, xi1] * dx[i] * dy[i])
    r_val, p_val = scipy.stats.pearsonr(v_rec, v_los_true)
    rms_res = np.sqrt(np.mean((v_rec - v_los_true)**2))
    print("--- Velocity Reconstruction Metrics ---")
    print("Pearson correlation coefficient (r): " + str(np.round(r_val, 5)))
    print("RMS residual: " + str(np.round(rms_res, 2)) + " km/s")
    W_ell = np.load(os.path.join(data_dir, 'wiener_filter.npy'))
    B_ell = np.load(os.path.join(data_dir, 'beam_transfer_function.npy'))
    effective_filter = B_ell * (1 - W_ell)
    dilution_factor = np.mean(effective_filter)
    print("\nEffective beam + Wiener filter dilution factor: " + str(np.round(dilution_factor, 5)))
    num_bins = 10
    sort_idx = np.argsort(mass)
    bin_edges_idx = np.linspace(0, len(mass), num_bins + 1).astype(int)
    tau_recovered = []
    tau_theoretical = []
    mass_mean = []
    c_km_s = 299792.458
    T_CMB_uK = 2.7255e6
    est_results = np.load(os.path.join(data_dir, 'estimator_results.npz'))
    C_2d = est_results['C_2d']
    C_0 = C_2d[0, 0]
    T_clean = np.load(os.path.join(data_dir, 'wiener_cleaned_map.npy'))
    ix = catalog['ix_pixel']
    iy = catalog['iy_pixel']
    T_halo = T_clean[iy, ix]
    print("\n--- Mass-Binned Optical Depth (tau) ---")
    print("Bin | Mean Mass (M_sun) | Recovered tau | Theoretical tau | Effective Th tau | Ratio (Rec/Eff)")
    for b in range(num_bins):
        idx = sort_idx[bin_edges_idx[b]:bin_edges_idx[b+1]]
        m_bin = mass[idx]
        T_bin = T_halo[idx]
        v_bin = v_los_true[idx]
        ix_bin = ix[idx]
        iy_bin = iy[idx]
        mean_m = np.mean(m_bin)
        mass_mean.append(mean_m)
        tau_th = 2e-3 * (mean_m / 1e14)**(2/3)
        tau_theoretical.append(tau_th)
        tau_th_eff = tau_th * dilution_factor
        i_idx, j_idx = np.triu_indices(len(idx), k=1)
        dx_ij = np.abs(ix_bin[i_idx] - ix_bin[j_idx])
        dy_ij = np.abs(iy_bin[i_idx] - iy_bin[j_idx])
        dx_ij = np.minimum(dx_ij, N - dx_ij)
        dy_ij = np.minimum(dy_ij, N - dy_ij)
        C_theta = C_2d[dy_ij, dx_ij]
        var_ij = 2 * (C_0 - C_theta)
        valid_pairs = var_ij > 1e-2
        w_ij = np.zeros_like(var_ij)
        w_ij[valid_pairs] = (m_bin[i_idx[valid_pairs]] * m_bin[j_idx[valid_pairs]] / 1e28) / var_ij[valid_pairs]
        T_diff = T_bin[i_idx] - T_bin[j_idx]
        v_diff = v_bin[i_idx] - v_bin[j_idx]
        num = np.sum(w_ij * T_diff * v_diff)
        den = np.sum(w_ij * v_diff**2)
        if den > 0:
            tau_rec = - (c_km_s / T_CMB_uK) * (num / den)
        else:
            tau_rec = 0.0
        tau_recovered.append(tau_rec)
        ratio = tau_rec / tau_th_eff if tau_th_eff > 0 else 0
        mean_m_str = str(np.round(mean_m / 1e14, 2)) + "e14"
        print(str(b+1) + " | " + mean_m_str + " | " + str(np.round(tau_rec, 8)) + " | " + str(np.round(tau_th, 8)) + " | " + str(np.round(tau_th_eff, 8)) + " | " + str(np.round(ratio, 4)))
    tau_recovered = np.array(tau_recovered)
    tau_theoretical = np.array(tau_theoretical)
    mass_mean = np.array(mass_mean)
    np.savez(os.path.join(data_dir, 'binned_tau_results.npz'), mass_mean=mass_mean, tau_recovered=tau_recovered, tau_theoretical=tau_theoretical, tau_theoretical_effective=tau_theoretical * dilution_factor, v_rec=v_rec, v_los_true=v_los_true)
    print("\nSaved reconstructed velocities and binned optical depth results to data/binned_tau_results.npz")

if __name__ == '__main__':
    main()