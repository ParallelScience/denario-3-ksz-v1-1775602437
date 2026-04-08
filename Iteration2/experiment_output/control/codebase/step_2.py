# filename: codebase/step_2.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
import numpy as np

def main():
    c_km_s = 299792.458
    T_CMB_uK = 2.7255e6
    obs_map_path = '/home/node/work/projects/ksz_v1/observed_map.npy'
    halo_catalog_path = '/home/node/work/projects/ksz_v1/halo_catalog.npz'
    data_dir = 'data/'
    cleaned_map_path = os.path.join(data_dir, 'wiener_cleaned_map.npy')
    obs_map = np.load(obs_map_path)
    cleaned_map = np.load(cleaned_map_path)
    catalog = np.load(halo_catalog_path)
    ix = catalog['ix_pixel']
    iy = catalog['iy_pixel']
    v_r = catalog['v_los_km_s']
    mass = catalog['mass_Msun']
    tau_true = catalog['tau']
    T_obs = obs_map[iy, ix]
    T_clean = cleaned_map[iy, ix]
    v_var = np.mean(v_r**2)
    tau_direct_unweighted_clean = - (c_km_s / (T_CMB_uK * v_var)) * np.mean(T_clean * v_r)
    tau_direct_unweighted_obs = - (c_km_s / (T_CMB_uK * v_var)) * np.mean(T_obs * v_r)
    M_weight = (mass / 1e14)**(2/3)
    num_direct_w = np.sum(M_weight * T_clean * v_r)
    den_direct_w = np.sum(M_weight**2 * v_r**2)
    A_hat_direct = - (c_km_s / T_CMB_uK) * (num_direct_w / den_direct_w)
    tau_direct_weighted_clean = A_hat_direct * np.mean(M_weight)
    num_direct_w_obs = np.sum(M_weight * T_obs * v_r)
    A_hat_direct_obs = - (c_km_s / T_CMB_uK) * (num_direct_w_obs / den_direct_w)
    tau_direct_weighted_obs = A_hat_direct_obs * np.mean(M_weight)
    N_pix = 1200
    C_2d = np.real(np.fft.ifft2(np.abs(np.fft.fft2(cleaned_map))**2)) / (N_pix * N_pix)
    C_0 = C_2d[0, 0]
    i_idx, j_idx = np.triu_indices(len(mass), k=1)
    dx = np.abs(ix[i_idx] - ix[j_idx])
    dy = np.abs(iy[i_idx] - iy[j_idx])
    dx = np.minimum(dx, N_pix - dx)
    dy = np.minimum(dy, N_pix - dy)
    C_theta = C_2d[dy, dx]
    var_ij = 2 * (C_0 - C_theta)
    var_ij = np.maximum(var_ij, 1e-6)
    mass_weight_i = M_weight[i_idx]
    mass_weight_j = M_weight[j_idx]
    w_ij = (mass[i_idx] * mass[j_idx] / 1e28) / var_ij
    T_diff_ij = T_clean[i_idx] - T_clean[j_idx]
    v_diff_ij = v_r[i_idx] - v_r[j_idx]
    num_pairwise = np.sum(w_ij * T_diff_ij * v_diff_ij)
    den_pairwise = np.sum(w_ij * v_diff_ij**2)
    F = np.sum(w_ij * (mass_weight_i + mass_weight_j)) / (2 * np.sum(w_ij))
    A_hat_pairwise = - (c_km_s / T_CMB_uK) * (num_pairwise / den_pairwise) / F
    tau_pairwise_clean = A_hat_pairwise * np.mean(M_weight)
    mean_tau_true = np.mean(tau_true)
    print("--- Recovered Mean Optical Depth (tau) ---")
    print("Catalog Mean tau: " + str(np.round(mean_tau_true, 8)))
    print("Direct Unweighted Estimator (Observed Map): " + str(np.round(tau_direct_unweighted_obs, 8)))
    print("Direct Unweighted Estimator (Cleaned Map):  " + str(np.round(tau_direct_unweighted_clean, 8)))
    print("Direct Mass-Weighted Estimator (Cleaned Map): " + str(np.round(tau_direct_weighted_clean, 8)))
    print("Pairwise Mass-Weighted Estimator (Cleaned Map): " + str(np.round(tau_pairwise_clean, 8)))
    np.savez(os.path.join(data_dir, 'estimator_results.npz'), T_obs=T_obs, T_clean=T_clean, tau_direct_unweighted_clean=tau_direct_unweighted_clean, tau_direct_unweighted_obs=tau_direct_unweighted_obs, tau_direct_weighted_clean=tau_direct_weighted_clean, tau_pairwise_clean=tau_pairwise_clean, mean_tau_true=mean_tau_true, C_2d=C_2d)
    np.savez(os.path.join(data_dir, 'pairwise_data.npz'), i_idx=i_idx, j_idx=j_idx, w_ij=w_ij, T_diff_ij=T_diff_ij, v_diff_ij=v_diff_ij)
    print("Saved per-halo temperatures and scalar results to data/estimator_results.npz")
    print("Saved computed pairwise statistics to data/pairwise_data.npz")

if __name__ == '__main__':
    main()