# circaTime 0.99.1

* `circularMAE()` and `circularError()` now accept a `period` argument (default
  24) and compute `pmin(abs(d), period - abs(d))` instead of hardcoding 24 hours;
  `benchmarkPhase()`, `transferPhase()`, and `circularAlign()` now pass their own
  `period` through to these metrics, so a non-24h cycle is scored on the requested
  cycle length. No behaviour change for the default 24h case.
* `benchmarkPhase()` and `transferPhase()` now also pass their `period` through to
  each engine's `fit` for engines that accept one (Molecular Timetable, ZeitZeiger,
  TimeSignatR), so the predictions themselves are produced on the requested cycle
  length rather than just the error metrics. tauFisher is inherently 24-hour and
  has no `period` argument; requesting a non-24h cycle with `"taufisher"` now
  emits a warning instead of silently returning 24h-scale predictions.

# circaTime 0.99.0

* Initial pre-release version.
* Four wrapped circadian phase-inference engines: Molecular Timetable
  (reimplemented from the original paper), tauFisher, ZeitZeiger, and
  TimeSignatR (the latter three wrap their authors' own R packages, installed
  separately since they are not on CRAN/Bioconductor).
* `benchmarkPhase()`: shared cross-validated benchmarking harness across all
  engines.
* `estimatePhase()`: apply a bundled reference model to new samples.
* Calibrated uncertainty layer: `bootstrapPhaseCI()`, `clockAmplitude()`,
  `calibrationCurve()`/`plotCalibration()`, `deriveAmplitudeThreshold()`.
* `standardizeGeneIds()`: translate Ensembl/Entrez/RefSeq IDs and known
  circadian gene aliases to gene symbols, with an AnnotationDbi/org.*.eg.db
  fallback for non-clock genes.
* `mapOrthologs()`: one-to-one mouse<->human ortholog mapping via `babelgene`;
  wired into `estimatePhase()` as `data_species`/`map_orthologs` for scoring
  samples against a reference model trained on the other species.
* `transferPhase()`: cross-dataset transfer benchmark - fit an engine on one
  `SummarizedExperiment` and score it on a different one (e.g. microarray-trained,
  RNA-seq-scored), reusing `benchmarkPhase()`'s engine registry and circular-error
  metrics.
* `estimatePhase(platform = ...)`: bundled reference models are now keyed by data
  type as well as species/tissue/method. `"microarray"` (GSE54650, the previous
  default and only option) and `"rnaseq"` (GSE54651) are both bundled for mouse.
* `estimatePhase(method = "consensus")`: combines every engine bundled for a
  given species/tissue/platform instead of making the user pick one. Two
  weighting schemes via `consensus_weights`: `"accuracy"` (default; each engine
  weighted by `1/cvMAE` using its own held-out error *on that tissue*) and
  `"equal"` (plain circular mean). The weights, each engine's own phase estimate,
  and a per-sample engine-agreement measure are all reported in
  `metadata(se)$circaTime_consensus` rather than hidden. Benchmarked against
  each single engine with real donor-split cross-validation on GSE54650
  (`benchmarks/run_consensus_benchmark.R`): pooled mean circular MAE 0.723h
  (accuracy-weighted consensus) vs. 0.726h (equal-weighted) vs. 0.738h
  (Molecular Timetable alone) vs. 0.804h (tauFisher alone) - both consensus
  schemes beat committing to a single engine everywhere, though the margin
  between the two weighting schemes themselves is narrow.
* `circularWeightedMean()`: weighted circular mean with a resultant-length
  ("concentration") agreement measure, the primitive behind consensus combining.
* `PhaseReferenceModel` gained a `cvMAE` slot (that model's own held-out
  cross-validated circular MAE) and `availableReferenceModels()` now reports it.
* `clockCoherence()`: dynamic replacement for the fixed 16-gene amplitude QC
  flag. Selects the informative gene panel *from the input dataset itself*
  (each gene's cosinor fit is tested against a phase-permuted null, panel =
  genes passing a BH-FDR cutoff) instead of a fixed core-clock panel, so
  dataset-specific rhythmic genes - not just the ones unaffected by aging per
  the legacy flag's own documented limitation - carry the QC signal.
* `estimatePhase(qc_method = "dynamic")`: uses `clockCoherence()` in place of
  `clockAmplitude()` for the `low_amplitude` flag, with a real-cohort
  requirement (`qc_min_samples`, default 10) since the panel is fit from the
  call's own samples; falls back to `qc_method = "legacy"` below that, recorded
  in `metadata(se)$circaTime_qc`. `qc_method = "legacy"` remains the default -
  no behaviour change for existing callers.
* Bug fix: `predictTauFisher()` (and therefore `estimatePhase(method =
  "taufisher")`/`"consensus"`) failed with a cryptic `no applicable method for
  'predict'` error in any session that predicts without having trained first -
  i.e. every real use of a bundled reference model. Caused by `nnet`'s S3
  method registration never being loaded (tauFisher's own NAMESPACE does not
  import `nnet`); fixed by requiring `nnet`'s namespace explicitly before
  predicting.
* `bayesianPhasePosterior()`: Phase 3 Step 3.3. A von Mises posterior predictive
  density and credible interval, built directly from the same residual pool
  [`bootstrapPhaseCI()`]/`estimatePhase()` already use rather than a new,
  independent inference method - concentration ([circularConcentration()], new)
  is measured from the pool, converted to the von Mises kappa parameter, and
  used to derive a parametric credible interval
  ([circularCredibleInterval()], new) and, on request, a full density grid.
* `estimatePhase()` now also returns `phase_credible_lo`/`phase_credible_hi`
  (the Bayesian interval above) alongside the existing empirical
  `phase_lo`/`phase_hi`; `metadata(se)$circaTime_bayesian` records the fitted
  `kappa`/`concentration`. Validated against the existing bootstrap interval on
  real held-out GSE54650 predictions (`benchmarks/run_bayesian_ci_benchmark.R`):
  pooled empirical coverage 88.9% (Bayesian) vs. 87.2% (bootstrap) at nominal
  80%, mean widths 2.72h vs. 2.78h.
