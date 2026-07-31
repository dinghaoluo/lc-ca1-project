# full pipeline run

This run sheet records the wrapper order I use to rebuild the repo-local outputs from preprocessing through figure generation, together with the commands called at each stage. Historical scripts under `_archive/` are outside this run sheet. Run commands from the repository root unless noted.

## before running

Start from the repo root. The path check prints the locations the scripts will use for generated data, figures, raw recordings, MATLAB spike-pipeline outputs, and Suite2p registrations. I check this before a long run because a wrong mount or stale path can waste hours before the first failure.

```powershell
# go to the repo root
cd Z:\Dinghao\code\lc-ca1-project

# print the main input and output roots used by the pipeline
python -u -c 'import sys; sys.path.insert(0, ''utils''); import project_paths as pp; print(pp.DATA_ROOT); print(pp.FIGURES_ROOT); print(pp.MICEEXP_ROOT); print(pp.MATLAB_SPIKE_PIPELINE_STEM); print(pp.SUITE2P_REGISTRATION_STEM)'
```

Commands below use `python -u` so progress prints immediately.

## environment and GPU

The Python environment needs `numpy`, `scipy`, `pandas`, `matplotlib`, `h5py`, `scikit-learn`, `tqdm`, and MATLAB-file readers such as `mat73`.

```powershell
# check the Python executable used by this terminal
where.exe python

# check CuPy import status from the same Python
python -u -c 'import sys, importlib.util; print(sys.executable); print(importlib.util.find_spec(''cupy''))'

# check NVIDIA driver and visible GPU status
nvidia-smi

# check CUDA device count and GPU name
python -u -c 'import cupy as cp; print(cp.cuda.runtime.getDeviceCount()); print(cp.cuda.runtime.getDeviceProperties(0)[''name''].decode())'
```

Successful GPU detection prints:

```text
GPU status: AVAILABLE (NVIDIA GeForce RTX 3090)
```

Without a visible GPU, the same check prints:

```text
GPU status: NOT AVAILABLE
```

CPU and GPU runs write the same set of outputs. CuPy shortens convolution-heavy steps, including dF/F computation and smoothing.

## MATLAB

The Python LC/HPC ephys extractors use existing MATLAB spike-pipeline `.mat` files. MATLAB is only needed when those files need to be regenerated.

```matlab
% refresh MATLAB spike-pipeline functions
cd('Z:\Dinghao\code\lc-ca1-project\preprocessing\matlab_spike_pipeline')
clear functions
rehash

% write MATLAB spike-pipeline outputs under data/matlab_spike_pipeline/
runspikepipeline
runspikepipeline_pix
```

## preprocessing data outputs

```powershell
python -u preprocessing\preprocessing_run_pipeline.py
```

The wrapper runs:

```powershell
# behaviour session pickles under data\behaviour\all_experiments\
python -u preprocessing\behaviour\process_behaviour.py --all

# aligned LC spike arrays under data\lc_ephys\all_sessions\
python -u preprocessing\ephys\lc\lc_all_extract.py

# aligned HPC spike arrays, all-trial trains, and run-trial trains under data\hpc_ephys\all_sessions\
python -u preprocessing\ephys\hpc\hpc_all_extract.py

# passive/Raphi HPC aligned spike arrays under data\hpc_ephys\all_sessions_raphi\
# A041 onward all_trains uses raw SciPy convolution for first-lick analyses
python -u preprocessing\ephys\hpc\hpc_all_extract_raphi.py

# suite2p registrations under data\suite2p_registration\, with raw-side mirroring where applicable
python -u preprocessing\imaging\suite2p_registration.py

# regular, immobile, and single-pixel axon-GCaMP session outputs
python -u preprocessing\imaging\axon_gcamp\lchpc_axon_all_extract.py
python -u preprocessing\imaging\axon_gcamp\lchpc_axon_all_extract_immobile.py
python -u preprocessing\imaging\axon_gcamp\lchpc_single_pixel_extract.py

# dLight, GRABNE, tone, and nLight processed traces and response maps
python -u preprocessing\imaging\dlight\hpc_dlight_lc_opto_extract.py
python -u preprocessing\imaging\grabne\hpc_grabne_lc_opto_extract.py
python -u preprocessing\imaging\grabne\hpc_grabne_tone_extract.py
python -u preprocessing\imaging\nlight\hpc_nlight_lc_opto_extract.py

# extraction figures for axon-GCaMP, dLight, GRABNE, tone, and nLight branches
python -u preprocessing\imaging\axon_gcamp\plot_lchpc_axon_all_extract.py
python -u preprocessing\imaging\axon_gcamp\plot_lchpc_axon_all_extract_immobile.py
python -u preprocessing\imaging\dlight\plot_hpc_dlight_lc_opto_extract.py
python -u preprocessing\imaging\grabne\plot_hpc_grabne_lc_opto_extract.py
python -u preprocessing\imaging\grabne\plot_hpc_grabne_tone_extract.py
python -u preprocessing\imaging\nlight\plot_hpc_nlight_lc_opto_extract.py
```

Single-session commands:

```powershell
# single-session axon-GCaMP pixel outputs with a smaller median-filter chunk
python -u preprocessing\imaging\axon_gcamp\lchpc_single_pixel_extract.py --recording A101i-20241031-01 --median-chunk-size 500

# single-session axon-GCaMP pixel outputs with the default median-filter chunk
python -u preprocessing\imaging\axon_gcamp\lchpc_single_pixel_extract.py --recording A101i-20241031-01

# single-session dLight LC-opto outputs plus inline extraction plots
python -u preprocessing\imaging\dlight\hpc_dlight_lc_opto_extract.py --recording A126i-20250605-02 --plots

# single-session nLight LC-opto outputs plus inline extraction plots
python -u preprocessing\imaging\nlight\hpc_nlight_lc_opto_extract.py --recording A171i-20260116-01 --plots
```

Interactive tools:

```powershell
# tiff-stack movie exports
python -u preprocessing\imaging\convert_movie_tif_gui.py

# FibreSight, including the former in-repository FibreSegger workflow, is now a standalone application; use the GUI from
# https://github.com/dinghaoluo/fibre-sight
# the older GUIs are archived under _archive/preprocessing/imaging/fibre_segger_GUI/ and _archive/preprocessing/imaging/fibre_sight/in_repo_workbench/
```

## core analysis data builders

```powershell
python -u analysis\core_analysis_run_pipeline.py
```

The wrapper runs:

```powershell
# behaviour speed/lick, pupil, reward-to-first-lick, and off-target run-bout outputs
python -u analysis\behaviour\analyse_speed_licks.py
python -u analysis\behaviour\analyse_pupil_size.py
python -u analysis\behaviour\first_lick_since_last_reward.py
python -u analysis\behaviour\off_target_run_bouts.py

# waveform, ACG, spike-time, ISI, identity UMAP, and LC cell-profile outputs
python -u analysis\lc\lc_all_waveforms_acgs.py
python -u analysis\lc\lc_all_spikes_isis.py
python -u analysis\lc\lc_all_identity_umap.py
python -u analysis\lc\lc_all_profiles.py

# profile pickles and waveform summaries for HPC cells
python -u analysis\hpc\hpc_all_profiles.py
python -u analysis\hpc\hpc_all_profiles_raphi.py
python -u analysis\hpc\hpc_all_waveforms.py

# axon-GCaMP profile table and plotting data
python -u analysis\imaging\lchpc_axon_all_profiles.py
```

The two HPC profile tables can also be rebuilt with alternative pre-run PyrUp/PyrDown classifications. These commands add columns such as `class_1s` and `class_half_s`; the original `class` columns remain unchanged.

```powershell
python -u analysis\hpc\hpc_all_profiles.py
python -u analysis\hpc\hpc_all_profiles_raphi.py
```

## controls

```powershell
python -u analysis\behaviour\behaviour_run_pipeline.py --branch controls
```

The wrapper runs:

```powershell
# behaviour-control figures for HPC and LC sessions
python -u analysis\behaviour_control\hpc_opto_speed_controls.py
python -u analysis\behaviour_control\lc_controls.py
python -u analysis\behaviour_control\lc_opto_speed_controls.py
```

## optogenetics summaries

```powershell
python -u analysis\optogenetics\optogenetics_run_pipeline.py
```

The wrapper runs:

```powershell
# ephys optogenetics behaviour summary figures
python -u analysis\optogenetics\summarise_opto.py

# imaging-animal optogenetics behaviour summary figures
python -u analysis\imaging\optogenetics\summarise_opto_imaging.py
```

## HPC analyses

```powershell
python -u analysis\hpc\hpc_run_pipeline.py --branch analyses
```

The wrapper runs:

```powershell
# first-lick distance and time for LC and LC-terminal stimulation
python -u analysis\hpc\behaviour\hpc_lc_stim_lick_comp.py --metric all --experiment all

# early/late first-lick HPC pyramidal-cell profile outputs
python -u analysis\hpc\first_lick_analysis\hpc_pyr_early_vs_late_lick_profiles_full.py

# behaviour-correlation data and figures for LC-stim CA1 PyrUp/PyrDown cells
python -u analysis\hpc\stim_ctrl\analyse_lc_stim_ca1_pyr_up_down_behaviour_corr.py
python -u analysis\hpc\stim_ctrl\plot_lc_stim_ca1_pyr_up_down_behaviour_corr.py

# CA1 place/time-cell metrics and sequence figures
python -u analysis\hpc\sequence_analysis\build_hpc_sequence_place_time_metrics.py
python -u plotting\hpc\plot_hpc_sequence_place_time.py

# effect-latency data and figures for LC-stim CA1 cells
python -u analysis\hpc\stim_ctrl\estimate_lc_stim_ca1_effect_latency.py
python -u analysis\hpc\stim_ctrl\plot_lc_stim_ca1_effect_latency.py
```

The alternative pre-run-window comparison is separate from the wrapper:

```powershell
python -u analysis\hpc\first_lick_analysis\hpc_pyr_early_vs_late_lick_profiles_alt_prewindows.py --repo-inputs
```

## LC analyses

```powershell
python -u analysis\lc\lc_analysis_run_pipeline.py
```

The wrapper runs:

```powershell
# run/cue/reward alignment and early/late first-lick response-profile outputs for LC cells
python -u analysis\lc\alignment_analysis\analyse_alignment_with_heatmap_run_cue_rew_aligned.py
python -u analysis\lc\first_lick_analysis\all_earlyvlate_ro_peak_fixed_threshold.py

# first-lick and control-vs-stim lick outputs for LC opto sessions
python -u analysis\lc\behaviour\lc_opto_first_lick_profile.py
python -u analysis\lc\behaviour\lc_opto_ctrl_vs_stim_lick_properties.py

# example-session lick and speed figures for LC and passive/Raphi sessions
python -u analysis\lc\behaviour\lc_example_session_licks.py
python -u analysis\lc\behaviour\lc_example_session_licks_passive_raphi.py
python -u analysis\lc\behaviour\lc_example_session_speed.py
python -u analysis\lc\behaviour\lc_example_session_speed_passive_raphi.py

# good-trial, lick-distance, lick-history, and lick-time outputs for LC opto sessions
python -u analysis\lc\behaviour\lc_opto_good_trial_percentage.py
python -u analysis\lc\behaviour\lc_opto_lick_distance.py
python -u analysis\lc\behaviour\hpc_lc_stim_lick_distance.py
python -u analysis\lc\behaviour\lc_opto_lick_history.py
python -u analysis\lc\behaviour\lc_opto_lick_history_comparison.py
python -u analysis\lc\behaviour\lc_opto_lick_time.py

# cue-start, run-bout, and single-trial example figures
python -u analysis\lc\behaviour\plot_cue_start_difference.py
python -u analysis\lc\behaviour\plot_run_bouts.py
python -u analysis\lc\behaviour\plot_run_bouts_examples.py
python -u analysis\lc\behaviour\plot_single_trial_example.py

# burst-amplitude behavioural-covariate outputs for LC cells
python -u analysis\lc\GLM\burst_amplitude_vs_baseline_rate.py
python -u analysis\lc\GLM\burst_amplitude_vs_reward_interval.py
python -u analysis\lc\GLM\burst_amplitude_vs_reward_interval_binned.py

# behaviour-encoding GLM permutation outputs for LC cells
python -u analysis\lc\GLM\glm_lc_beh_permutation.py
python -u analysis\lc\GLM\glm_lc_beh_permutation_full.py

# tonic firing-rate FFT output
python -u analysis\lc\GLM\tonic_fft_lc.py

# opto stimulation-response output for LC cells
python -u analysis\lc\ephys_opto\analyse_stim_response.py

# run-onset vs run-bout burst-amplitude output
python -u analysis\lc\run_onset_v_run_bout\lc_run_onset_vs_bout_burst_amplitude.py

# optogenetic tagging latency output
python -u analysis\lc\tagging_analysis\tagging_latency.py
```

## imaging analyses

```powershell
python -u analysis\imaging\imaging_analysis_run_pipeline.py
```

The wrapper runs:

```powershell
# axon-GCaMP run/cue/reward alignment heatmap outputs
python -u analysis\imaging\alignment_analysis\analyse_alignment_with_heatmap_run_cue_rew_aligned.py

# dLight expression, inhibition-control, spatial-dispersion, and LC-opto release-profile outputs
python -u analysis\imaging\controls\dlight_expression_control.py
python -u analysis\imaging\dLight_inhibition\hpc_dlight_lc_inh_stim_ctrl_run.py
python -u analysis\imaging\dLight_stim_dispersion\single_roi_binned_dilation_spatial_tau.py
python -u analysis\imaging\optogenetics\dlight_lc_opto_release_stim_ctrl.py

# dLight release-probability, significant-release proportion, and release-site heterogeneity outputs
python -u analysis\imaging\release_probability\prop_signif_release_dlight_stim.py
python -u analysis\imaging\release_probability\release_probability_dlight_stim.py
python -u analysis\imaging\release_probability\release_site_heterogeneity_dlight_stim.py

# mean and time-resolved ROI-vs-neuropil response-index outputs
python -u analysis\imaging\ROI_vs_neuropil\roi_vs_neuropil_ri_mean.py
python -u analysis\imaging\ROI_vs_neuropil\roi_vs_neuropil_ri_over_time.py

# nLight expression, spatial-dispersion, and LC-opto release-profile outputs
python -u analysis\imaging\controls\nlight_expression_control.py
python -u analysis\imaging\nLight_stim_dispersion\single_roi_binned_dilation_spatial_tau.py
python -u analysis\imaging\optogenetics\nlight_lc_opto_release_stim_ctrl.py

# nLight release-probability, significant-release proportion, and release-site heterogeneity outputs
python -u analysis\imaging\release_probability\prop_signif_release_nlight_stim.py
python -u analysis\imaging\release_probability\release_probability_nlight_stim.py
python -u analysis\imaging\release_probability\release_site_heterogeneity_nlight_stim.py

# mean and time-resolved nLight ROI-vs-neuropil response-index outputs
python -u analysis\imaging\ROI_vs_neuropil\roi_vs_neuropil_ri_mean_nlight.py
python -u analysis\imaging\ROI_vs_neuropil\roi_vs_neuropil_ri_over_time_nlight.py
```

## pharmacology

```powershell
python -u analysis\pharmacology\pharmacology_run_pipeline.py
```

The wrapper runs:

```powershell
# pharmacology summary figures under figures\pharmacology\
python -u analysis\pharmacology\summarise_all_drugs.py
```

## plotting scripts

```powershell
python -u plotting\plotting_run_pipeline.py
```

The wrapper runs:

```powershell
# behaviour speed/lick, pupil, reward-to-first-lick, example, immobile, and trial-by-trial figures
python -u plotting\behaviour\plot_speed_licks.py
python -u plotting\behaviour\plot_pupil_size.py
python -u plotting\behaviour\plot_first_lick_since_last_reward.py
python -u plotting\behaviour\plot_example_session.py
python -u plotting\behaviour\plot_example_trials.py
python -u plotting\behaviour\plot_immobile.py
python -u plotting\behaviour\plot_speeds.py
python -u plotting\behaviour\plot_trial_by_trial.py

# ctrl/stim, raster, heatmap, information, pre/post, and PyrUp/PyrDown figures for HPC cells
python -u plotting\hpc\plot_all_ctrl_stim_profiles.py
python -u plotting\hpc\plot_all_ctrl_stim_rasters.py
python -u plotting\hpc\plot_all_pyr_heatmap_dist.py
python -u plotting\hpc\plot_all_pyr_info_ctrl_stim.py
python -u plotting\hpc\plot_all_pyr_pre_post_ratio.py
python -u plotting\hpc\plot_all_pyr_pre_post_raw_change.py
python -u plotting\hpc\plot_run_onset_pyr_up_down_profiles.py
python -u plotting\hpc\plot_run_onset_pyr_up_down_profiles_raphi.py
python -u plotting\hpc\plot_hpc_sequence_place_time.py

# imaging reference, lick, heatmap, trace, whole-field, and dLight/nLight LC-opto figures
python -u plotting\imaging\example_sess_refs_release_tiff.py
python -u plotting\imaging\plot_ref_channel_16bit_maps.py
python -u plotting\imaging\plot_dlight_lc_opto_single_axon_stim_profiles.py
python -u plotting\imaging\plot_nlight_lc_opto_single_axon_stim_profiles.py
python -u plotting\imaging\plot_lick_profile.py
python -u plotting\imaging\plot_lick_profile_to_pumps.py
python -u plotting\imaging\plot_pooled_heatmap_axon_gcamp.py
python -u plotting\imaging\plot_raw_traces_axon_gcamp.py
python -u plotting\imaging\plot_raw_traces_axon_gcamp_example_trials.py
python -u plotting\imaging\plot_sorted_heatmaps_grids.py
python -u plotting\imaging\plot_sorted_heatmaps_rois.py
python -u plotting\imaging\plot_std_heatmap.py
python -u plotting\imaging\plot_whole_field.py
python -u plotting\imaging\summarise_dlight_lc_opto_all.py
python -u plotting\imaging\summarise_dlight_lc_opto_ctrl_inh.py
python -u plotting\imaging\summarise_nlight_lc_opto_all.py

# cell ACG, identity, profile, ISI, first-lick, raster, run-onset, waveform, and tagging figures for LC cells
python -u plotting\lc\plot_acgs_and_heatmap.py
python -u plotting\lc\plot_lc_tagged_vs_putative.py
python -u plotting\lc\plot_ctrl_stim_profiles.py
python -u plotting\lc\plot_isis.py
python -u plotting\lc\plot_lc_first_lick_sensitive_profiles.py
python -u plotting\lc\plot_rasters_1st_lick_ordered_early_late_only.py
python -u plotting\lc\plot_rasters_run_cue_rew_aligned.py
python -u plotting\lc\plot_runonset_burst_and_non_burst_profiles.py
python -u plotting\lc\plot_single_cell_acg.py
python -u plotting\lc\plot_single_cell_waveform.py
python -u plotting\lc\plot_tagged_example_good_bad_raster.py
python -u plotting\lc\plot_tagging_responses.py
python -u plotting\lc\plot_lc_trial_profiles.py

# rasters ordered by reward, lick, lick/reward sensitivity, and reward-to-run-onset context for LC cells
python -u plotting\lc\rasters\all_raster_last_reward_to_current_trial.py
python -u plotting\lc\rasters\all_raster_last_rew_ordered.py
python -u plotting\lc\rasters\all_raster_lick_ordered.py
python -u plotting\lc\rasters\all_raster_lick_ordered_raster_only.py
python -u plotting\lc\rasters\all_raster_lick_reward_sensitivity.py
python -u plotting\lc\rasters\all_raster_rew_ordered.py
python -u plotting\lc\rasters\all_raster_reward_to_run_onset_ordered.py
```
