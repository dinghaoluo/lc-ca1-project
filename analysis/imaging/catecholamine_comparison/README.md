# Catecholamine comparison

This directory contains the imaging analyses that compare dLight3.6- and nLightG-reported LC-CA1 axonal signals. The scripts sit beside the data domain that they analyse, rather than in a thesis-specific wrapper; the chapter panel set is one possible composition of these measurements, while the underlying outputs remain reusable elsewhere in the repository.

## Assembly notes

When the upstream natural-event products have changed, rebuild them first:

```text
python analysis/imaging/nlight_vs_dlight/quantify_nlight_vs_dlight.py
```

Then run the current Chapter 3 package from the repository root:

```text
python analysis/imaging/catecholamine_comparison/build_catecholamine_session_manifest.py
python analysis/imaging/catecholamine_comparison/analyse_natural_threshold_sensitivity.py
python analysis/imaging/catecholamine_comparison/build_catecholamine_comparison_evidence.py
python Z:/Dinghao/code/render_profile_with_synapses.py
python Z:/Dinghao/code/assemble_chapter3_panels.py
```

The final source tables are written to `data/catecholamine_comparison/representative_inputs/`; the final PNG, PDF and SVG files are written to `Z:/Dinghao/code/lc_ca1_chapter3_assembly/outputs/`. `assemble_chapter3_panels.py` assembles the optical panel set, invokes the EM assembler for the EM panel set, and removes stale assembly exports, and refreshes `analysis_manifest.csv` after every successful run.

The assembled figures use Dinghao's opposing half-violin convention through `plot_violin_with_scatter`: every pale circular point is one session, all points have the same size and style, the maximum width of each half-violin is fixed across panels, and the black horizontal mark is the cohort median. Natural-event traces are summarised first within session and then as cohort mean +/- SEM, so a session with more candidate regions does not acquire greater weight; the representative stimulation traces align each trial to its first finite post-stimulation frame, the shutter-reopening point, before calculating mean +/- SEM.

## Figure logic

Figure 18 establishes the principal result under experimentally timed LC stimulation: fit-qualified nLightG traces decay more slowly than dLight3.6 traces, and the separation remains after each decay constant is divided by the published reporter off-kinetic. Figure 19 then asks where those evoked responses lie, using the Chapter 2 distance-ring and spatial-decay procedure in both cohorts; nLightG responses are detected in fewer LC-axon ROIs and have a shorter fitted pixel-domain spatial constant. Figure 20 tests whether the optogenetic kinetic separation is also detected at natural run onset, where it is not, while exposing the threshold dependence and green-channel measurement controls rather than hiding them. the EM panel set closes at a different scale: vesicles cluster within mitochondria-excluded cytoplasm in all four reconstructed profiles, while their proximity to nearby synapse centres and the unequal loading of those centres vary between profiles. A reserved panel remains for an author-supplied EM image showing LC axons near a local CA1 synapse.

This order is intentional. The longer optogenetic nLightG trace is the strongest functional result, the spatial analysis characterises the same stimulation cohort, natural events provide a less constrained comparison rather than a competing premise, and EM supplies anatomical context without being made to explain the optical measurements.

## Interpretation boundaries

The optogenetic decay result concerns reporter traces. It does not isolate release, clearance or reporter chemistry, because the measured decay also contains expression, acquisition and fit-selection effects; only 8 of 13 included nLightG time courses and 34 of 34 dLight3.6 time courses meet the decay-fit criterion. The sensor-normalised ratio shows that the observed separation is not explained by the published nLightG and dLight3.6 off-kinetics alone, but it is not an *in vivo* transmitter calibration.

The Figure 19 response criterion reproduces the implemented Chapter 2 rule: a one-sample greater-than-zero test with mean response index above 0.1. Circular shuffling belongs to the natural run-onset classifier in Figure 20 and is not substituted for the Chapter 2 stimulation criterion. Spatial distances remain in native pixels because the available field geometries do not support a single defensible physical conversion across every recording.

The EM null draws 500 deterministic patterns uniformly from the segmented cytoplasm, with mitochondria excluded, and retains one profile-level statistic from each pattern. The nearest-neighbour analysis tests organisation within the available axonal cytoplasm itself; the proximity and loading analyses then ask whether that organisation is related to nearby synapse centres. The four profiles are the units shown, not independent animals. The present bundle does not contain feasible-centre masks eroded by vesicle radius or an independent annotation set, so those controls remain explicit limitations rather than panels that merely report unavailable data. Nearby neuropil synapse centres are not proven contacts made by the reconstructed axon.

Earlier diagnostic scripts remain in this directory as method history, but their outputs are not part of the four-figure publication package. `figures/catecholamine_comparison/evidence/` contains diagnostic products only; the stale `figures/nlight_vs_dlight_complete/` material must not be used for Chapter 3 writing.

The standalone Results narrative, figure legends and exact interpretation limits are in `catecholamine_comparison_results_report.md`. `catecholamine_comparison_evidence_report.md` records the final evidence chain and the checks needed to reproduce it.


For Figure 18c, extract_representative_post_shutter_traces.py re-extracts whole-field traces from the registered movies. Each trial is baseline-subtracted over the preceding one second, then starts at its first observed frame after shutter reopening. Trials are aligned before the mean and SEM are calculated; samples overlapping the next stimulation or lacking a complete ten-second observation window are excluded. The compact output is data/catecholamine_comparison/representative_inputs/representative_post_shutter_traces.npz.
