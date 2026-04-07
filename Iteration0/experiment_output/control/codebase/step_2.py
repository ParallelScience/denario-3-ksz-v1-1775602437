# filename: codebase/step_2.py
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

if __name__ == '__main__':
    observed_map_path = '/home/node/work/projects/ksz_v1/observed_map.npy'
    ksz_map_path = '/home/node/work/projects/ksz_v1/ksz_map_truth.npy'
    catalog_path = '/home/node/work/projects/ksz_v1/halo_catalog.npz'
    ps_path = 'data/power_spectra_and_beam.npz'
    observed_map = np.load(observed_map_path)
    ksz_map = np.load(ksz_map_path)
    catalog = np.load(catalog_path)
    ix = catalog['ix_pixel']
    iy = catalog['iy_pixel']
    ps_data = np.load(ps_path)
    B_l_2d = ps_data['B_l_2d']
    Cl_CMB_2d = ps_data['Cl_CMB_2d']
    Cl_noise_2d = ps_data['Cl_noise_2d']
    numerator = Cl_CMB_2d * (B_l_2d**2)
    denominator = numerator + Cl_noise_2d
    W_l = np.zeros_like(numerator)
    mask = denominator > 0
    W_l[mask] = numerator[mask] / denominator[mask]
    observed_map_fft = np.fft.fft2(observed_map)
    filtered_cmb_fft = observed_map_fft * W_l
    filtered_cmb = np.real(np.fft.ifft2(filtered_cmb_fft))
    residual_map = observed_map - filtered_cmb
    ksz_map_fft = np.fft.fft2(ksz_map)
    filtered_ksz_cmb_fft = ksz_map_fft * W_l
    filtered_ksz_cmb = np.real(np.fft.ifft2(filtered_ksz_cmb_fft))
    residual_ksz = ksz_map - filtered_ksz_cmb
    radius_pix = 1.5 * 1.4 / 0.5
    true_ksz_sums = aperture_sum(ksz_map, ix, iy, radius_pix)
    residual_ksz_sums = aperture_sum(residual_ksz, ix, iy, radius_pix)
    transfer_function = np.sum(true_ksz_sums * residual_ksz_sums) / np.sum(true_ksz_sums**2)
    np.save('data/residual_map.npy', residual_map)
    np.save('data/transfer_function.npy', np.array([transfer_function]))
    print("--- Wiener Filter Construction ---")
    print("Wiener filter W(l) min/max: " + str(np.round(W_l.min(), 4)) + " / " + str(np.round(W_l.max(), 4)))
    print("Residual map RMS: " + str(np.round(np.std(residual_map), 2)) + " uK")
    print("Filtered CMB map RMS: " + str(np.round(np.std(filtered_cmb), 2)) + " uK")
    print("\n--- Transfer Function ---")
    print("Aperture radius: " + str(np.round(radius_pix, 2)) + " pixels")
    print("Mean absolute true kSZ aperture sum: " + str(np.round(np.mean(np.abs(true_ksz_sums)), 4)) + " uK")
    print("Mean absolute residual kSZ aperture sum: " + str(np.round(np.mean(np.abs(residual_ksz_sums)), 4)) + " uK")
    print("Signal Transfer Function (Attenuation Factor): " + str(np.round(transfer_function, 4)))
    print("\nSaved residual_map.npy and transfer_function.npy to data/")