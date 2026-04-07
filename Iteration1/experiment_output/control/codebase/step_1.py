# filename: codebase/step_1.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

def load_and_summarize_data(data_dir, observed_map_path, cmb_map_path, ksz_map_path, noise_map_path, catalog_path):
    observed_map = np.load(observed_map_path)
    cmb_map = np.load(cmb_map_path)
    ksz_map = np.load(ksz_map_path)
    noise_map = np.load(noise_map_path)
    catalog = np.load(catalog_path)
    print("--- Map Summary Statistics ---")
    maps = {"Observed Map": observed_map, "CMB Truth Map": cmb_map, "kSZ Truth Map": ksz_map, "Noise Map": noise_map}
    for name, m in maps.items():
        rms = np.sqrt(np.mean(m**2))
        print(name + ":")
        print("  RMS: " + str(rms) + " uK")
        print("  Min: " + str(np.min(m)) + " uK")
        print("  Max: " + str(np.max(m)) + " uK")
        print("")
    print("--- Catalog Summary Statistics ---")
    fields = ['mass_Msun', 'log_mass', 'v_los_km_s', 'tau', 'delta_T_ksz_uK']
    for field in fields:
        data = catalog[field]
        rms = np.sqrt(np.mean(data**2))
        print("Field: " + field)
        print("  RMS: " + str(rms))
        print("  Min: " + str(np.min(data)))
        print("  Max: " + str(np.max(data)))
        print("")
    pixel_size_deg = 10.0 / 1200.0
    df_catalog = pd.DataFrame({k: catalog[k] for k in catalog.keys()})
    df_catalog['x_deg_from_pix'] = df_catalog['ix_pixel'] * pixel_size_deg
    df_catalog['y_deg_from_pix'] = df_catalog['iy_pixel'] * pixel_size_deg
    processed_catalog_path = os.path.join(data_dir, "processed_halo_catalog.csv")
    df_catalog.to_csv(processed_catalog_path, index=False)
    print("Processed catalog saved to " + processed_catalog_path)
    return observed_map, cmb_map, ksz_map, noise_map, df_catalog

def plot_diagnostic_maps(data_dir, observed_map, cmb_map, ksz_map, noise_map):
    mpl.rcParams['text.usetex'] = False
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    map_plot_info = [(observed_map, "Observed Map", axes[0, 0]), (cmb_map, "CMB Truth Map", axes[0, 1]), (ksz_map, "kSZ Truth Map", axes[1, 0]), (noise_map, "Noise Map", axes[1, 1])]
    for m, title, ax in map_plot_info:
        rms = np.sqrt(np.mean(m**2))
        im = ax.imshow(m, cmap='viridis', origin='lower', extent=[0, 10, 0, 10])
        ax.set_title(title + " (RMS: " + str(round(rms, 2)) + " uK)")
        ax.set_xlabel("RA (deg)")
        ax.set_ylabel("Dec (deg)")
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label("Temperature (uK)")
    plt.tight_layout()
    timestamp = str(int(time.time()))
    plot_filename = "diagnostic_maps_1_" + timestamp + ".png"
    plot_filepath = os.path.join(data_dir, plot_filename)
    plt.savefig(plot_filepath, dpi=300)
    plt.close()
    print("Diagnostic maps plot saved to " + plot_filepath)

if __name__ == '__main__':
    data_directory = "data/"
    obs_map_path = "/home/node/work/projects/ksz_v1/observed_map.npy"
    cmb_truth_path = "/home/node/work/projects/ksz_v1/cmb_map_truth.npy"
    ksz_truth_path = "/home/node/work/projects/ksz_v1/ksz_map_truth.npy"
    noise_truth_path = "/home/node/work/projects/ksz_v1/noise_map.npy"
    halo_catalog_path = "/home/node/work/projects/ksz_v1/halo_catalog.npz"
    obs_map, cmb_truth, ksz_truth, noise_truth, catalog_df = load_and_summarize_data(data_directory, obs_map_path, cmb_truth_path, ksz_truth_path, noise_truth_path, halo_catalog_path)
    plot_diagnostic_maps(data_directory, obs_map, cmb_truth, ksz_truth, noise_truth)