# LC-CA1 project

This repository contains the working code for my projects on the locus coeruleus (LC) and hippocampal CA1: behaviour, electrophysiology, two-photon imaging, optogenetics, pharmacology, and modelling. It is a reorganised and cleaned-up version of [`code-mpfi-dinghao`](https://github.com/dinghaoluo/code-mpfi-dinghao); the older repository remains online in its original, unorganised form.

Raw recordings remain on MPFI storage. Several input paths therefore point to `Z:/Dinghao/...`, and some older scripts retain lab-specific conventions. The main paths are collected in `utils/project_paths.py`. Generated data, figures, and logs are kept under `data/`, `figures/`, and `outputs/`; all three directories are ignored by Git.

## repository layout

- `preprocessing/`: behaviour, electrophysiology, imaging, and MATLAB preprocessing
- `analysis/`: profile builders and branch-specific analyses
- `plotting/`: figure scripts that read saved analysis outputs
- `utils/`: shared path, IO, analysis, plotting, and imaging functions
- `modelling/`: the LC-DA-CA1 model and its support code
- `rec_list.py`: recording lists used across preprocessing and analysis
- `_archive/`: older and exploratory scripts outside the current run path

## running the pipeline

[`FULL_PIPELINE_RUN.md`](FULL_PIPELINE_RUN.md) records the wrapper order and the commands called at each stage, from preprocessing to figure generation. Most commands are run from the repository root. The broad sequence is:

1. Build the behaviour, electrophysiology, and imaging preprocessing outputs.
2. Build the shared LC, HPC, behaviour, and imaging profile tables.
3. Run the branch-specific analyses.
4. Regenerate figures from the saved outputs.

The run sheet also records the Python environment and GPU checks. MATLAB is only needed when the spike-pipeline `.mat` files have to be regenerated.

## branch notes

### behaviour

`preprocessing/behaviour/process_behaviour.py` writes one behaviour pickle per recording under `data/behaviour/all_experiments/<dataset>/`. The behaviour analyses build speed, lick, pupil, reward-history, and run-bout outputs from those session files.

### LC and HPC electrophysiology

The LC and HPC extractors write aligned spike arrays under `data/lc_ephys/` and `data/hpc_ephys/`. Their profile builders combine those session arrays with behaviour and waveform information; the Raphi/passive HPC recordings have a separate extractor and profile table.

### imaging

The imaging branches cover axon GCaMP, dLight, GRABNE, and nLight recordings. Extraction scripts save processed traces, response maps, and plotting inputs under the paths in `utils/project_paths.py`. Inline plotting is off by default in the newer extractors; `--plots` regenerates those figures during extraction.

### optogenetics, pharmacology, and modelling

Optogenetic summaries combine the behavioural, electrophysiological, and imaging branches. `analysis/pharmacology/summarise_all_drugs.py` builds the pharmacology summaries. The standalone LC-DA-CA1 model lives under `modelling/`.

## interactive tools

- `preprocessing/imaging/convert_movie_tif_gui.py`: TIFF-stack movie exports
- FibreSight is now a standalone application, and the older in-repository FibreSegger and FibreSight GUIs are archived. Use the current GUI at [dinghaoluo/fibre-sight](https://github.com/dinghaoluo/fibre-sight); the historical code is retained under `_archive/preprocessing/imaging/fibre_segger_GUI/` and `_archive/preprocessing/imaging/fibre_sight/in_repo_workbench/`.

## historical code

The current pipeline wrappers do not call scripts under `_archive/`. Those files remain for reference and may still contain old paths, duplicated helpers, unfinished branches, or imports from the earlier repository layout.
