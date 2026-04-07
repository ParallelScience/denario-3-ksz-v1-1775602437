# filename: codebase/step_4.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
import numpy as np
import os

def aperture_sum(map_data, ix_arr, iy_arr, r_pix):
    ny, nx = map_data.shape
    r_int = int(np.ceil(r_pix))
    y, x = np.ogrid[-r_int:r_int+1, -r_int:r_int+1]
    mask_circle = x**2 + y**2 <= r_pix**2
    sums = np.zeros(len(ix_arr))
    for i in range(len(ix_arr)):
        cx, cy = ix_arr[i], iy_arr[i]
        y_indices = cy + np.arange(-r_int, r_int+1)
        x_indices = cx + np.arange(-r_int, r_int+1)
        valid_y = (y_indices >= 0) & (y_indices < ny)
        valid_x = (x_indices >= 0) & (x_indices < nx)
        if not (valid_y.any() and valid_x.any()):
            continue
        y_min_sub = np.argmax(valid_y)
        y_max_sub = len(valid_y) - np.argmax(valid_y[::-1])
        x_min_sub = np.argmax(valid_x)
        x_max_sub = len(valid_x) - np.argmax(valid_x[::-1])
        map_y_min = y_indices[y_min_sub]
        map_y_max = y_indices[y_max_sub-1] + 1
        map_x_min = x_indices[x_min_sub]
        map_x_max = x_indices[x_max_sub-1] + 1
        sub_map = np.zeros((2*r_int+1, 2*r_int+1))
        sub_map[y_min_sub:y_max_sub, x_min_sub:x_max_sub] = map_data[map_y_min:map_y_max, map_x_min:map_x_max]
        sums[i] = np.sum(sub_map[mask_circle])
    return sums

def fmt_sci(val):
    return np.format_float_scientific(val, precision=2)

def fmt_float(val):
    return str(np.round(val, 2))

if __name__ == '__main__':
    residual_map_path = 'data/residual_map.npy'
    catalog_path = '/home/node/work/projects/ksz_v1/halo_catalog.npz'
    transfer_function_path = 'data/transfer_function.npy'
    residual_map = np.load(residual_map_path)
    catalog = np.load(catalog_path)
    transfer_function = np.load(transfer_function_path)[0]
    if transfer_function == 0:
        transfer_function = 1e-10
    ix = catalog['ix_pixel']
    iy = catalog['iy_pixel']
    v_los = catalog['v_los_km_s']
    log_mass = catalog['log_mass']
    mass = catalog['mass_Msun']
    true_tau = catalog['tau']
    c_km_s = 299792.458
    T_CMB_uK = 2.7255e6
    radius_pix = 1.5 * 1.4 / 0.5
    T_res = aperture_sum(residual_map, ix, iy, radius_pix)
    n_bins = 10
    bins = np.linspace(13, 15, n_bins + 1)
    bin_indices = np.digitize(log_mass, bins) - 1
    bin_indices[bin_indices == n_bins] = n_bins - 1
    tau_hat_bins = np.zeros(n_bins)
    tau_hat_err_bins = np.zeros(n_bins)
    tau_true_mean_bins = np.zeros(n_bins)
    mass_mean_bins = np.zeros(n_bins)
    n_boot = 500
    bootstrap_distributions = np.zeros((n_bins, n_boot))
    print("--- Optical Depth Estimation ---")
    print("Transfer function applied: " + str(np.round(transfer_function, 4)))
    for b in range(n_bins):
        mask = bin_indices == b
        n_halos = np.sum(mask)
        if n_halos == 0:
            continue
        T_res_b = T_res[mask]
        v_los_b = v_los[mask]
        numerator = np.sum(T_res_b * v_los_b)
        denominator = np.sum(v_los_b**2)
        if denominator == 0:
            tau_raw = 0.0
        else:
            tau_raw = - (c_km_s / T_CMB_uK) * (numerator / denominator)
        tau_corrected = tau_raw / transfer_function
        tau_hat_bins[b] = tau_corrected
        tau_true_mean_bins[b] = np.mean(true_tau[mask])
        mass_mean_bins[b] = np.mean(mass[mask])
        boot_tau = np.zeros(n_boot)
        for i in range(n_boot):
            idx = np.random.randint(0, n_halos, n_halos)
            num_boot = np.sum(T_res_b[idx] * v_los_b[idx])
            den_boot = np.sum(v_los_b[idx]**2)
            if den_boot == 0:
                boot_tau[i] = 0.0
            else:
                boot_tau[i] = - (c_km_s / T_CMB_uK) * (num_boot / den_boot) / transfer_function
        tau_hat_err_bins[b] = np.std(boot_tau)
        bootstrap_distributions[b, :] = boot_tau
        print("Bin " + str(b+1) + " (logM " + str(np.round(bins[b], 1)) + "-" + str(np.round(bins[b+1], 1)) + "): N=" + str(n_halos) + ", True tau = " + fmt_sci(tau_true_mean_bins[b]) + ", Est tau = " + fmt_sci(tau_hat_bins[b]) + " +/- " + fmt_sci(tau_hat_err_bins[b]))
    print("\n--- Null Test (Shuffled Velocities) ---")
    null_tau_distributions = np.zeros((n_bins, n_boot))
    for i in range(n_boot):
        v_los_shuffled = np.random.permutation(v_los)
        for b in range(n_bins):
            mask = bin_indices == b
            if np.sum(mask) == 0:
                continue
            T_res_b = T_res[mask]
            v_los_b = v_los_shuffled[mask]
            numerator = np.sum(T_res_b * v_los_b)
            denominator = np.sum(v_los_b**2)
            if denominator == 0:
                tau_raw = 0.0
            else:
                tau_raw = - (c_km_s / T_CMB_uK) * (numerator / denominator)
            null_tau_distributions[b, i] = tau_raw / transfer_function
    null_means = np.mean(null_tau_distributions, axis=1)
    null_stds = np.std(null_tau_distributions, axis=1)
    for b in range(n_bins):
        if null_stds[b] > 0:
            sig = abs(tau_hat_bins[b] - null_means[b]) / null_stds[b]
        else:
            sig = 0.0
        print("Bin " + str(b+1) + " Null Test: Mean = " + fmt_sci(null_means[b]) + ", Std = " + fmt_sci(null_stds[b]) + ", Significance = " + fmt_float(sig) + " sigma")
    print("\n--- Pairwise kSZ Estimator (Cross-check) ---")
    pairwise_tau_bins = np.zeros(n_bins)
    for b in range(n_bins):
        mask = bin_indices == b
        n_h = np.sum(mask)
        if n_h < 2:
            continue
        T_res_b = T_res[mask]
        v_los_b = v_los[mask]
        sum_Tv = np.sum(T_res_b * v_los_b)
        sum_T = np.sum(T_res_b)
        sum_v = np.sum(v_los_b)
        sum_v2 = np.sum(v_los_b**2)
        numerator = n_h * sum_Tv - sum_T * sum_v
        denominator = n_h * sum_v2 - sum_v**2
        if denominator == 0:
            tau_raw_pw = 0.0
        else:
            tau_raw_pw = - (c_km_s / T_CMB_uK) * (numerator / denominator)
        tau_corrected_pw = tau_raw_pw / transfer_function
        pairwise_tau_bins[b] = tau_corrected_pw
        print("Bin " + str(b+1) + ": Pairwise Est tau = " + fmt_sci(tau_corrected_pw) + " (Direct Est = " + fmt_sci(tau_hat_bins[b]) + ")")
    print("\n--- Global Pairwise kSZ Estimator ---")
    N_total = len(T_res)
    sum_Tv_g = np.sum(T_res * v_los)
    sum_T_g = np.sum(T_res)
    sum_v_g = np.sum(v_los)
    sum_v2_g = np.sum(v_los**2)
    numerator_global = N_total * sum_Tv_g - sum_T_g * sum_v_g
    denominator_global = N_total * sum_v2_g - sum_v_g**2
    if denominator_global == 0:
        tau_raw_pw_global = 0.0
    else:
        tau_raw_pw_global = - (c_km_s / T_CMB_uK) * (numerator_global / denominator_global)
    tau_corrected_pw_global = tau_raw_pw_global / transfer_function
    print("Global Pairwise Est tau = " + fmt_sci(tau_corrected_pw_global))
    output_path = 'data/optical_depth_results.npz'
    np.savez(output_path, mass_mean_bins=mass_mean_bins, tau_true_mean_bins=tau_true_mean_bins, tau_hat_bins=tau_hat_bins, tau_hat_err_bins=tau_hat_err_bins, bootstrap_distributions=bootstrap_distributions, null_tau_distributions=null_tau_distributions, pairwise_tau_bins=pairwise_tau_bins, tau_corrected_pw_global=np.array([tau_corrected_pw_global]))
    print("\nSaved optical depth results to " + output_path)