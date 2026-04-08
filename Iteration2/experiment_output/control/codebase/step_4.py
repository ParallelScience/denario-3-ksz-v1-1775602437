# filename: codebase/step_4.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
from scipy.optimize import curve_fit
import time
from matplotlib.gridspec import GridSpecFromSubplotSpec
plt.rcParams['text.usetex'] = False
def get_binned_tau_jk(catalog, T_map, C_2d, num_bins=10, N_jk=100):
    ra = catalog['ra']
    dec = catalog['dec']
    mass = catalog['mass_Msun']
    v_los = catalog['v_los_km_s']
    ix = catalog['ix_pixel']
    iy = catalog['iy_pixel']
    N_pix = 1200
    C_0 = C_2d[0, 0]
    ra_bins = np.linspace(0, 10, 11)
    dec_bins = np.linspace(0, 10, 11)
    ra_idx = np.digitize(ra, ra_bins[1:-1])
    dec_idx = np.digitize(dec, dec_bins[1:-1])
    subregion_idx = ra_idx * 10 + dec_idx
    sort_idx = np.argsort(mass)
    bin_edges_idx = np.linspace(0, len(mass), num_bins + 1).astype(int)
    tau_bins = np.zeros(num_bins)
    tau_jk = np.zeros((N_jk, num_bins))
    mass_mean = np.zeros(num_bins)
    c_km_s = 299792.458
    T_CMB_uK = 2.7255e6
    for b in range(num_bins):
        idx = sort_idx[bin_edges_idx[b]:bin_edges_idx[b+1]]
        m_bin = mass[idx]
        T_bin = T_map[iy[idx], ix[idx]]
        v_bin = v_los[idx]
        ix_bin = ix[idx]
        iy_bin = iy[idx]
        sub_bin = subregion_idx[idx]
        mass_mean[b] = np.mean(m_bin)
        i_idx, j_idx = np.triu_indices(len(idx), k=1)
        dx_ij = np.abs(ix_bin[i_idx] - ix_bin[j_idx])
        dy_ij = np.abs(iy_bin[i_idx] - iy_bin[j_idx])
        dx_ij = np.minimum(dx_ij, N_pix - dx_ij)
        dy_ij = np.minimum(dy_ij, N_pix - dy_ij)
        C_theta = C_2d[dy_ij, dx_ij]
        var_ij = 2 * (C_0 - C_theta)
        valid_pairs = var_ij > 1e-2
        w_ij = np.zeros_like(var_ij)
        w_ij[valid_pairs] = (m_bin[i_idx[valid_pairs]] * m_bin[j_idx[valid_pairs]] / 1e28) / var_ij[valid_pairs]
        T_diff = T_bin[i_idx] - T_bin[j_idx]
        v_diff = v_bin[i_idx] - v_bin[j_idx]
        num = np.sum(w_ij * T_diff * v_diff)
        den = np.sum(w_ij * v_diff**2)
        tau_bins[b] = - (c_km_s / T_CMB_uK) * (num / den) if den > 0 else 0
        sub_i = sub_bin[i_idx]
        sub_j = sub_bin[j_idx]
        for k in range(N_jk):
            keep = (sub_i != k) & (sub_j != k)
            num_k = np.sum(w_ij[keep] * T_diff[keep] * v_diff[keep])
            den_k = np.sum(w_ij[keep] * v_diff[keep]**2)
            tau_jk[k, b] = - (c_km_s / T_CMB_uK) * (num_k / den_k) if den_k > 0 else 0
    err_bins = np.sqrt((N_jk - 1) / N_jk * np.sum((tau_jk - np.mean(tau_jk, axis=0))**2, axis=0))
    return mass_mean, tau_bins, err_bins
def compute_overall_tau_exact(T_map, v_array, catalog, N_jk=100):
    N_pix = 1200
    C_2d = np.real(np.fft.ifft2(np.abs(np.fft.fft2(T_map))**2)) / (N_pix * N_pix)
    C_0 = C_2d[0, 0]
    mass = catalog['mass_Msun']
    ix = catalog['ix_pixel']
    iy = catalog['iy_pixel']
    ra = catalog['ra']
    dec = catalog['dec']
    ra_bins = np.linspace(0, 10, 11)
    dec_bins = np.linspace(0, 10, 11)
    ra_idx = np.digitize(ra, ra_bins[1:-1])
    dec_idx = np.digitize(dec, dec_bins[1:-1])
    subregion_idx = ra_idx * 10 + dec_idx
    i_idx, j_idx = np.triu_indices(len(mass), k=1)
    dx_ij = np.abs(ix[i_idx] - ix[j_idx])
    dy_ij = np.abs(iy[i_idx] - iy[j_idx])
    dx_ij = np.minimum(dx_ij, N_pix - dx_ij)
    dy_ij = np.minimum(dy_ij, N_pix - dy_ij)
    C_theta = C_2d[dy_ij, dx_ij]
    var_ij = 2 * (C_0 - C_theta)
    var_ij = np.maximum(var_ij, 1e-6)
    M_weight = (mass / 1e14)**(2/3)
    mass_weight_i = M_weight[i_idx]
    mass_weight_j = M_weight[j_idx]
    w_ij = (mass[i_idx] * mass[j_idx] / 1e28) / var_ij
    T_halo = T_map[iy, ix]
    T_diff = T_halo[i_idx] - T_halo[j_idx]
    v_diff = v_array[i_idx] - v_array[j_idx]
    c_km_s = 299792.458
    T_CMB_uK = 2.7255e6
    sum_w = np.sum(w_ij)
    num = np.sum(w_ij * T_diff * v_diff)
    den = np.sum(w_ij * v_diff**2)
    F = np.sum(w_ij * (mass_weight_i + mass_weight_j)) / (2 * sum_w) if sum_w > 0 else 1
    A_hat = - (c_km_s / T_CMB_uK) * (num / den) / F if den > 0 else 0
    tau_all = A_hat * np.mean(M_weight)
    sub_i = subregion_idx[i_idx]
    sub_j = subregion_idx[j_idx]
    tau_jk = np.zeros(N_jk)
    for k in range(N_jk):
        keep = (sub_i != k) & (sub_j != k)
        sum_w_k = np.sum(w_ij[keep])
        num_k = np.sum(w_ij[keep] * T_diff[keep] * v_diff[keep])
        den_k = np.sum(w_ij[keep] * v_diff[keep]**2)
        F_k = np.sum(w_ij[keep] * (mass_weight_i[keep] + mass_weight_j[keep])) / (2 * sum_w_k) if sum_w_k > 0 else 1
        A_hat_k = - (c_km_s / T_CMB_uK) * (num_k / den_k) / F_k if den_k > 0 else 0
        tau_jk[k] = A_hat_k * np.mean(M_weight)
    err_all = np.sqrt((N_jk - 1) / N_jk * np.sum((tau_jk - np.mean(tau_jk))**2))
    return tau_all, err_all
def extract_stacked_cutouts(T_map, catalog, num_bins=3, cutout_size=21):
    mass = catalog['mass_Msun']
    v_los = catalog['v_los_km_s']
    ix = catalog['ix_pixel']
    iy = catalog['iy_pixel']
    sort_idx = np.argsort(mass)
    bin_edges_idx = np.linspace(0, len(mass), num_bins + 1).astype(int)
    half_size = cutout_size // 2
    N_pix = T_map.shape[0]
    stacks = []
    for b in range(num_bins):
        idx = sort_idx[bin_edges_idx[b]:bin_edges_idx[b+1]]
        stack = np.zeros((cutout_size, cutout_size))
        count = 0
        for i in idx:
            x = ix[i]
            y = iy[i]
            if x >= half_size and x < N_pix - half_size and y >= half_size and y < N_pix - half_size:
                cutout = T_map[y-half_size:y+half_size+1, x-half_size:x+half_size+1]
                weight = -v_los[i] / 300.0
                stack += cutout * weight
                count += 1
        if count > 0:
            stack /= count
        stacks.append(stack)
    return stacks
def main():
    data_dir = 'data/'
    halo_catalog_path = '/home/node/work/projects/ksz_v1/halo_catalog.npz'
    obs_map_path = '/home/node/work/projects/ksz_v1/observed_map.npy'
    catalog = np.load(halo_catalog_path)
    obs_map = np.load(obs_map_path)
    cleaned_map = np.load(os.path.join(data_dir, 'wiener_cleaned_map.npy'))
    perfect_cleaned_map = np.load(os.path.join(data_dir, 'perfect_cleaned_map.npy'))
    est_results = np.load(os.path.join(data_dir, 'estimator_results.npz'))
    C_2d_clean = est_results['C_2d']
    binned_res = np.load(os.path.join(data_dir, 'binned_tau_results.npz'))
    v_rec = binned_res['v_rec']
    v_los_true = binned_res['v_los_true']
    W_ell = np.load(os.path.join(data_dir, 'wiener_filter.npy'))
    B_ell = np.load(os.path.join(data_dir, 'beam_transfer_function.npy'))
    dilution_factor = np.mean(B_ell * (1 - W_ell))
    mass_mean, tau_bins, err_bins = get_binned_tau_jk(catalog, cleaned_map, C_2d_clean, num_bins=10, N_jk=100)
    def model(M, A, alpha):
        return A * (M / 1e14)**alpha
    try:
        popt, pcov = curve_fit(model, mass_mean, tau_bins, sigma=err_bins, absolute_sigma=True, p0=[2e-3 * dilution_factor, 0.66], bounds=([-np.inf, 0], [np.inf, 3.0]))
        A_fit, alpha_fit = popt
        err_A, err_alpha = np.sqrt(np.diag(pcov))
    except Exception:
        def model_fixed(M, A):
            return A * (M / 1e14)**(2/3)
        popt, pcov = curve_fit(model_fixed, mass_mean, tau_bins, sigma=err_bins, absolute_sigma=True, p0=[2e-3 * dilution_factor])
        A_fit = popt[0]
        alpha_fit = 2/3
        err_A = np.sqrt(np.diag(pcov))[0]
        err_alpha = 0.0
    configs = [
        ('Observed Map + True v', obs_map, v_los_true),
        ('Cleaned Map + True v', cleaned_map, v_los_true),
        ('Perfect Cleaned + True v', perfect_cleaned_map, v_los_true),
        ('Cleaned Map + Rec v', cleaned_map, v_rec)
    ]
    snr_labels = []
    snr_values = []
    for name, T_map, v_arr in configs:
        t_val, e_val = compute_overall_tau_exact(T_map, v_arr, catalog, N_jk=100)
        snr = t_val / e_val if e_val > 0 else 0
        snr_labels.append(name)
        snr_values.append(snr)
    stacks = extract_stacked_cutouts(cleaned_map, catalog, num_bins=3, cutout_size=21)
    fig = plt.figure(figsize=(16, 12))
    ax1 = plt.subplot(2, 2, 1)
    ax1.errorbar(mass_mean, tau_bins, yerr=err_bins, fmt='o', color='black', label='Recovered tau', capsize=4)
    M_smooth = np.logspace(np.log10(np.min(mass_mean)), np.log10(np.max(mass_mean)), 100)
    tau_th_eff = 2e-3 * (M_smooth / 1e14)**(2/3) * dilution_factor
    tau_fit = model(M_smooth, A_fit, alpha_fit)
    ax1.plot(M_smooth, tau_th_eff, 'b--', label='Theoretical')
    ax1.plot(M_smooth, tau_fit, 'r-', label='Best Fit')
    ax1.set_xscale('log')
    ax1.set_xlabel('Halo Mass M_500 [M_sun]')
    ax1.set_ylabel('Effective Optical Depth tau')
    ax1.set_title('(a) tau-M Scaling Relation')
    ax1.legend()
    ax1.grid(True, which='both', ls='--', alpha=0.5)
    ax2 = plt.subplot(2, 2, 2)
    r_val, _ = scipy.stats.pearsonr(v_rec, v_los_true)
    ax2.scatter(v_los_true, v_rec, alpha=0.3, s=10, color='purple')
    min_v = min(np.min(v_los_true), np.min(v_rec))
    max_v = max(np.max(v_los_true), np.max(v_rec))
    ax2.plot([min_v, max_v], [min_v, max_v], 'k--', alpha=0.5)
    ax2.set_xlabel('True Velocity [km/s]')
    ax2.set_ylabel('Reconstructed Velocity [km/s]')
    ax2.set_title('(b) Velocity Reconstruction (r = ' + str(np.round(r_val, 3)) + ')')
    ax2.grid(True, ls='--', alpha=0.5)
    gs = fig.add_gridspec(2, 2)
    gs_c = GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[1, 0], wspace=0.3)
    titles = ['Low Mass', 'Mid Mass', 'High Mass']
    vmin = min([np.min(s) for s in stacks])
    vmax = max([np.max(s) for s in stacks])
    vmax_sym = max(abs(vmin), abs(vmax))
    for i in range(3):
        ax_c = fig.add_subplot(gs_c[0, i])
        im = ax_c.imshow(stacks[i], cmap='coolwarm', origin='lower', extent=[-5.25, 5.25, -5.25, 5.25], vmin=-vmax_sym, vmax=vmax_sym)
        ax_c.set_title(titles[i])
        plt.colorbar(im, ax=ax_c, fraction=0.046, pad=0.04)
    fig.text(0.25, 0.48, '(c) Stacked kSZ Cutouts', ha='center', fontsize=12, fontweight='bold')
    ax4 = plt.subplot(2, 2, 4)
    x_pos = np.arange(len(snr_labels))
    ax4.bar(x_pos, snr_values, color=['gray', 'blue', 'green', 'orange'], alpha=0.7)
    ax4.set_xticks(x_pos)
    ax4.set_xticklabels(['Obs+TrueV', 'Clean+TrueV', 'Perf+TrueV', 'Clean+RecV'], rotation=15)
    ax4.set_ylabel('Detection SNR')
    ax4.set_title('(d) SNR Improvement Comparison')
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    timestamp = int(time.time())
    plot_filename = os.path.join(data_dir, 'ksz_analysis_summary_' + str(timestamp) + '.png')
    fig.savefig(plot_filename, dpi=300)
    print('Saved to ' + plot_filename)
if __name__ == '__main__':
    main()