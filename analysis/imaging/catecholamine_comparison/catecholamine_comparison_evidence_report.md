# Evidence record for the catecholamine-comparison publication package

Updated 31 August 2026.

## Scope

This record supports the four-figure Chapter 3 package and nothing broader: Figure 18, optogenetic persistence; Figure 19, spatial confinement; Figure 20, natural run-onset events; and the EM panel set, EM vesicle organisation. The analytical additions made for this package are the sensor-normalised decay ratio, the Chapter 2-style spatial analysis requested for the nLightG cohort, and the mitochondria-excluded cytoplasmic nearest-neighbour analysis for all four EM profiles. Earlier diagnostic outputs remain in the repository for provenance, but they do not enter the assembly output directory, the chapter narrative or the handoff manifest.

## Evidence chain

| Figure | Question | Principal source | Cohort and unit | Supported conclusion |
| --- | --- | --- | --- | --- |
| 18 | Does the stimulation-evoked reporter trace persist for different durations? | `roi_neuropil_timecourse_metrics.csv`, `sensor_normalised_session_metrics.csv`, `opto_trace_sensitivity.csv` | 34 dLight3.6 and 13 nLightG included sessions; 34 and 8 decay fits | Fit-qualified nLightG traces decay more slowly, and published reporter off-kinetics do not explain the separation |
| 19 | How does the evoked response change with distance from labelled LC axons? | `spatial_dilation_curves.csv`, `spatial_tau_metrics.csv` | 34 dLight3.6 and 13 nLightG sessions | nLightG responses meet the common criterion in fewer ROIs and have a shorter fitted pixel-domain reach |
| 20 | Is the same kinetic separation detected at natural run onset? | `nlight_vs_dlight_single_state/per_session_metrics.csv`, natural threshold tables | 57 dLight3.6 sessions from 4 animals and 26 nLightG sessions from 2 animals | Recruitment depends on the classifier setting; a between-cohort natural decay or rise-time difference is not detected |
| 21 | Are vesicles uniformly arranged within each reconstructed profile, and is that organisation related to nearby synapse centres? | `em_profile_null_statistics.csv`, `em_profile_null_distributions.npz` | Four profiles, each tested against 500 null patterns | Vesicles cluster within mitochondria-excluded cytoplasm in every profile; proximity and loading effects remain profile-specific, and no animal-level conclusion is available |

The logic is sequential. Experimentally timed stimulation supplies the clearest temporal result; the same stimulation cohort then supplies the spatial result; natural events test whether that temporal distinction generalises to behaviour; serial-section EM closes with an anatomical question that is deliberately not used to explain reporter kinetics.

## Verification of the two new analyses

### Sensor-normalised decay

For every fit-qualified session, `roi_tau_sensor_ratio` was checked against:

```text
roi_tau_s / sensor_tau_off_s
```

The maximum absolute discrepancy is constrained to be below `1e-12` in `build_catecholamine_comparison_evidence.py`. The published constants used here are 0.33 s for dLight3.6 and 0.19 s for nLightG. The cohort medians are 4.081 [3.855, 4.279] for dLight3.6 and 9.881 [8.734, 10.671] for nLightG; U = 272, session p = 1.69 x 10^-8, rank-biserial correlation = 1.0, and animal exact p = .0571. These values are stored as the `roi_tau_sensor_ratio` row of `data/catecholamine_comparison/evidence/statistical_results.csv`.

This ratio asks whether the fitted separation merely follows the two published off-kinetics. It does not convert fluorescence into transmitter concentration and does not remove cohort differences in expression, acquisition or fit selection.

### Chapter 2-style spatial analysis

The spatial analysis uses the 0-1 s response window, one-pixel distance rings around labelled LC-axon ROIs and the same exponential fitting procedure in both reporter cohorts. Session ring profiles enter the cohort mean and SEM with equal session weight. The implemented Chapter 2 response rule was recovered from the Chapter 2 code and reproduced directly: a one-sample greater-than-zero test together with mean response index above 0.1. It is not a shuffle-based stimulation classifier. Circular shuffling remains confined to the natural run-onset analysis.

The spatial decay constant is 5.856 [4.719, 8.545] pixels for dLight3.6 and 3.804 [3.339, 4.392] pixels for nLightG; U = 70, session p = 3.45 x 10^-4, rank-biserial correlation = -.683, and animal exact p = .200. The response-qualified fraction is .918 [.850, 1.000] and .0938 [.0476, .250], respectively; U = 4, session p = 2.19 x 10^-7, rank-biserial correlation = -.982, and animal exact p = .0286. Both rows are stored in `statistical_results.csv`.

Spatial distances remain in pixels. Of 47 available fields, 45 share the common 512-by-512, zoom-three geometry; `A152i-20251113-01` and `A152i-20251114-01` do not. A mixed or inferred micrometre conversion would appear more biological while being less defensible.

## Controls and inference boundaries

### Optogenetic traces

The merged representative stimulation panel begins at time zero defined separately for each trial as its first finite post-stimulation frame, the first moment at which the shutter has reopened and a reporter signal is observed, and continues to 10 s on that aligned display axis. Each trial is re-zeroed to its own prestimulation baseline before alignment; uncertainty is SEM across the 32 dLight3.6 or 18 nLightG aligned trials. Only the green dLight3.6 and blue nLightG reporter traces are shown. The separate neuropil analysis tests whether surrounding fluorescence accounts for persistence; the decay-fit cohort is smaller for nLightG, 8 fits from 13 included time courses, than for dLight3.6, 34 from 34.

### Natural events

Natural traces are likewise summarised within session before cohort mean and SEM. The primary response criterion uses the 95th percentile of a circular-shuffle distribution and an effect-size threshold of .05, with all nine percentile-by-threshold combinations displayed. At the primary setting the session comparison is p = .00293, but the exact animal comparison is p = .533 and only two nLightG animals are available. Natural decay gives session p = .698 and animal p = 1.000; rise time gives session p = .549 and animal p = .533. These results are absence of detected differences, not equivalence.

The publication controls remain in the reporter domain: baseline fluorescence, pre-onset noise, usable trial count and motion-regression fit. The red-channel response is retained in the compact source table for provenance but is not plotted. The half-violins use a common maximum width so differences in panel geometry do not reflect density or sample size.

### EM profiles

Each EM null pattern samples uniformly, with replacement, from the segmented axonal cytoplasm after mitochondrial exclusion. One median vesicle nearest-neighbour distance, one median distance to a nearby synapse centre and one loading coefficient are retained per pattern. Observed nearest-neighbour medians are 52.01, 39.08, 57.21 and 56.74 nm in P1-P4, against null medians of 73.14, 62.50, 105.12 and 149.12 nm; all four profiles give empirical p = q = .001996. P1 and P2 have proximity p = .00200, while P3 and P4 have p = .623 and .723. Loading p values are .00200, .0339, .631 and .0220 for P1-P4. The source table contains the corresponding BH-adjusted q values.

The current inputs do not contain vesicle-radius-eroded feasible-centre masks, so boundary sensitivity cannot be calculated, and they do not contain an independent annotation set, so annotation agreement cannot be calculated. These limitations remain documented but are not used as empty figure panels. the EM panel sete is instead reserved for an author-supplied EM image showing LC axons near a local CA1 synapse. The fixed 1,000-nm inclusion distance selects nearby neuropil synapse centres; it does not prove that the reconstructed axon forms those synapses.

## Inference hierarchy

Animals are the biological replicates for the imaging comparisons, sessions are nested observations, and ROIs or candidate regions are measurement units. The figures show session distributions because fit selection and field heterogeneity should remain visible, but each cross-cohort comparison also reports an exact test on animal medians. Session-level p values are not substituted for animal replication. In the EM analysis, the profile is the observed unit and every null pattern contributes one profile-level statistic; simulated coordinates are not treated as independent observations.

The strongest evidence is the directionally complete session separation in optogenetic decay, reinforced by the sensor-normalised ratio, and the response-qualified fraction difference in Figure 19, which is also supported by the exact animal-level test. The spatial decay and optogenetic decay effects remain limited by small animal counts, while the natural-event recruitment contrast remains threshold-dependent. The EM findings remain profile-specific.

## Reproducibility map

The final builders are:

```text
analysis/imaging/catecholamine_comparison/build_catecholamine_comparison_evidence.py
Z:/Dinghao/code/assemble_chapter3_panels.py
analysis/em/render_profile_with_synapses.py
Z:/Dinghao/code/assemble_em_panels.py
```

The panel-level source tables are:

```text
data/catecholamine_comparison/representative_inputs/panel_source_table.csv
data/catecholamine_comparison/representative_inputs/representative_fields.csv
data/catecholamine_comparison/representative_inputs/summary_statistics.csv
data/catecholamine_comparison/representative_inputs/persistence_statistics.csv
data/catecholamine_comparison/representative_inputs/em_profile_null_statistics.csv
data/catecholamine_comparison/representative_inputs/analysis_manifest.csv
```

The final figures are:

```text
Z:/Dinghao/code/lc_ca1_chapter3_assembly/outputs/figure_18_optogenetic_persistence.*
Z:/Dinghao/code/lc_ca1_chapter3_assembly/outputs/figure_19_spatial_confinement.*
Z:/Dinghao/code/lc_ca1_chapter3_assembly/outputs/figure_20_natural_events.*
Z:/Dinghao/code/lc_ca1_chapter3_assembly/outputs/figure_21_em_vesicles.*
```

Every figure is exported as a 300-dpi PNG, an editable PDF and an SVG. `analysis_manifest.csv` records size and SHA-256 for the builders, written reports, principal source tables and final figure files. The assembly script removes retired assembly outputs before it writes this manifest, so a successful build leaves exactly four figure stems in the assembly output directory.
