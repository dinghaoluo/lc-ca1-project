# Fibre Sight workbench

Fibre Sight is a small workbench for labelling channel-2 axon ROIs, training the U-Net, running saved-model predictions, and correcting the resulting ROIs by hand. The older MSER route remains available before a trained model exists and on difficult images.

Both routes use the same ROI dictionary as the rest of the imaging pipeline: each ROI contains `xpix` and `ypix`. MSER proposals, saved labels, and model predictions open in the same editor. Fixed ROIs are written before a fresh MSER pass, so a decision made by hand survives another segmentation; predicted ROIs can be deleted, merged, or saved as the next set of training labels.

The defaults use the dLight HPC LC-opto data, training recipe, and model on the MPFI setup; several paths are local, and the data remain outside the repository.

## run

From the repository root:

```powershell
conda activate fibre-sight
python preprocessing\imaging\fibre_sight\fibre_sight_workbench_gui.py
```

## training data

Each usable session needs:

```text
<source root>/<session>/processed_data/ref_mat_ch2.npy
<source root>/<session>/processed_data/*_ROI_dict.npy
```

The channel-2 reference image and ROI dict must describe the same recording. `scan labelled sessions` writes the CSV used for training, including session paths, image and ROI summaries, inclusion state, and train, validation, or test split. Training settings are stored in YAML; trained checkpoints use `.pt`.

## the cycle

1. Load a `ref_mat_ch2.npy` image, or load an existing ROI dict to continue curation.
2. Before a model exists, run MSER, tune the visible segmentation controls, and fix ROIs that should survive another MSER pass.
3. Delete or merge poor proposals, then save the ROI dict beside the processed session.
4. Scan the labelled sessions and train the small U-Net. The diagnostics window runs held-out scoring and writes label previews and prediction overlays.
5. Load the saved `.pt` checkpoint and predict on a channel-2 image.
6. Adjust `strictness` or `min ROI size`; the saved confidence map is reused, so this step does not rerun the model.
7. Curate the prediction in the same editor, save it, and include the corrected labels in a later training round when needed.

## current model and diagnostics

The model registry points to `data/imaging/fibre_sight/runs/dlight_hpc_lc_opto_unet/best.pt`, with threshold `0.25`, minimum size `45`, and TTA enabled. Existing diagnostics live under:

```text
data/imaging/fibre_sight/model_comparison.txt
data/imaging/fibre_sight/runs/dlight_hpc_lc_opto_unet/history.csv
data/imaging/fibre_sight/runs/dlight_hpc_lc_opto_unet/history.png
data/imaging/fibre_sight/runs/dlight_hpc_lc_opto_unet/evaluate_test_t025_m45_tta.csv
data/imaging/fibre_sight/runs/dlight_hpc_lc_opto_unet/evaluation_examples_test.png
data/imaging/fibre_sight/diagnostics/final_model_prediction_overlays.png
data/imaging/fibre_sight/diagnostics/first_training_diagnostics/
```

## controls

The main and advanced MSER panels contain the values I changed while comparing proposals with hand-labelled ROIs.

- `strictness` sets the model-confidence threshold. Higher values usually keep fewer ROIs; lower values recover fainter ones.
- `min ROI size` removes small connected components after prediction.
- The prediction pass uses image flips and averaging by default.
- `show model confidence` leaves the confidence map under the ROI outlines while thresholds are adjusted.
- `ROI on / ROI off` hides the outlines when the raw channel-2 image needs a clean look.
