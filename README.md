# Near-Moon_detectability
Data files to accompany 'Detectability of RSOs in the Near-Moon Environment: Optical and Short-wave Infrared' paper. CC-BY-4.0

[Detectability of RSOs in the Near-Moon Environment: Optical and Short-wave Infrared](https://arxiv.org) Cooke et al. 2026.

## head_data_filtered_clean.csv
Data file containing the per-pointing averaged frame statistics

### Column names
- date: Date of pointing (night of)
- date_idx: Date index
- time: Time (UTC)
- cam: Camera used (CMOS or SWIR)
- exp_time: Exposure time of pointing (sec)
- X: Average airmass of pointing (assuming $X = \sec(Z)$&#8203;)
- m_s: Moon separation (deg)
- m_i: Moon illumnation fraction (%)
- bg_mean: Background mean (ADU) - airmass and exposure time corrected
- bg_rms: Background rms (ADU) - airmass and exposure time corrected

## Analysis fits
Fits and parameters for background statistics as a function of $m_i$ and $m_s$
- CMOS_fit_params_update.npz: Fit parameters for CMOS $B_{\rm rms}$ as a function of $B_{\rm mean}$ (also includes CMOS G-band zero point)
- CMOS_bg_mean_rbf_spline_update.pkl: RBF interpolation fit for CMOS $B_{\rm mean}$ as a function of $m_i$ and $m_s$
- CMOS_bg_rms_rbf_spline_update.pkl: RBF interpolation fit for CMOS $B_{\rm rms}$ as a function of $m_i$ and $m_s$
- CMOS_bg_mean_scaling_params.npz: Normalisation constants for CMOS $B_{\rm mean}$ RBF interpolation
- CMOS_bg_rms_scaling_params.npz: Normalisation constants for CMOS $B_{\rm rms}$ RBF interpolation

- SWIR_fit_params_update.npz: Fit parameters for SWIR $B_{\rm rms}$ as a function of $B_{\rm mean}$ (also includes SWIR J-band zero point)
- SWIR_bg_mean_rbf_spline_update.pkl: RBF interpolation fit for SWIR $B_{\rm mean}$ as a function of $m_i$ and $m_s$
- SWIR_bg_rms_rbf_spline_update.pkl: RBF interpolation fit for SWIR $B_{\rm rms}$ as a function of $m_i$ and $m_s$
- SWIR_bg_mean_scaling_params.npz: Normalisation constants for SWIR $B_{\rm mean}$ RBF interpolation
- SWIR_bg_rms_scaling_params.npz: Normalisation constants for SWIR $B_{\rm rms}$ RBF interpolation

## near_moon_detectability_clean.py
Python code to convert analysis fits to target detectability as a function of lunar conditions. Produces data and plots providing limiting magnitude of detectable targets for both CMOS and SWIR in two regimes.

### Output files (compressed)
- CMOS_near_moon_detectability_out1_update.csv.zip: CMOS limiting G-band magnitude as a function of $m_i$, $m_s$ and $r$
- CMOS_near_moon_detectability_out2_update.csv.zip: CMOS limiting G-band magnitude as a function of $m_i$, $m_s$ and $n$
- SWIR_near_moon_detectability_out1_update.csv.zip: SWIR limiting J-band magnitude as a function of $m_i$, $m_s$ and $r$
- SWIR_near_moon_detectability_out2_update.csv.zip: SWIR limiting J-band magnitude as a function of $m_i$, $m_s$ and $n$

#### Column names
- moon_illum: Moon illumination (%)
- moon_sep: Moon separation (degrees)
- target_rate: Relative target rate ($''$/s)
- n_frames: Number of frames per stack
- limiting_mag: Limiting magnitude (G/J band for CMOS/SWIR)
