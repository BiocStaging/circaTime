# circaTime

A unified [Bioconductor](https://bioconductor.org) interface for inferring circadian
phase (time of sample collection relative to the internal clock) from bulk and
donor-level pseudobulk transcriptomes.

**Status: in development, not yet submitted to Bioconductor.** See
[Installation](#installation) below for how to install in the meantime.

## Why

Tens of thousands of public RNA-seq samples have no recorded collection time, but
circadian phase can often be recovered computationally from expression data alone.
Several published methods do this, but they live in separate repositories, in
different languages, take incompatible inputs, produce incompatible outputs, and have
never been benchmarked against each other on shared, labelled ground truth.
`circaTime` does not introduce a new phase-inference algorithm — it wraps existing
published methods behind one `SummarizedExperiment`-based interface, adds a shared
cross-validation benchmarking harness, and adds a calibrated uncertainty layer:
bootstrap and Bayesian confidence intervals validated for real coverage on labelled
data, plus a quality flag (fixed core-clock amplitude, or a dynamic per-dataset
coherence panel) for samples whose clock is too damped for a phase estimate to be
trustworthy. A consensus mode combines every engine bundled for a sample instead of
requiring you to pick one, and gene-ID/cross-species handling lets a reference model
trained on one species and platform score data from another.

## Installation

`circaTime` is still in development and has not yet been submitted to Bioconductor.
Until it is accepted, install the development version directly from GitHub:

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("dralperenuysal/circaTime")
```

Three of the wrapped engines ([`tauFisher`](https://github.com/micnngo/tauFisher),
[`TimeSignatR`](https://github.com/braunr/TimeSignatR),
[`zeitzeiger`](https://github.com/hugheylab/zeitzeiger)) live only on GitHub, so they
are `Suggests`, not hard dependencies — `circaTime` installs and works (Molecular
Timetable engine included) without them. Install whichever ones you need separately:

```r
remotes::install_github("micnngo/tauFisher")
remotes::install_github("braunr/TimeSignatR")
remotes::install_github("hugheylab/zeitzeiger")
```

Non-symbol gene identifiers (Ensembl/Entrez/RefSeq) and cross-species ortholog
mapping (`standardizeGeneIds()`, `mapOrthologs()`) are also optional: install
`AnnotationDbi`, `org.Mm.eg.db`/`org.Hs.eg.db`, and `babelgene` if you need them.

**Once `circaTime` is accepted into Bioconductor**, this section will be replaced with
the standard:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("circaTime")
```

## Usage

### Estimate phase for new samples

```r
library(circaTime)
library(SummarizedExperiment)

# See what reference models are currently bundled (species/tissue/method)
availableReferenceModels()

se <- estimatePhase(se, method = "timetable", tissue = "Liv")
colData(se)[, c("phase", "phase_lo", "phase_hi", "clock_amp", "phase_flag")]
```

`estimatePhase()` never hard-errors on a difficult individual sample — it sets
`phase = NA` and an informative `phase_flag` (`"low_amplitude"`, `"gene_coverage"`, or
`"failed"`) instead. Only malformed input (an unsupported tissue/species/method
combination, or missing gene identifiers) is a hard error.

Beyond the bootstrap `phase_lo`/`phase_hi`, every call also returns a Bayesian
`phase_credible_lo`/`phase_credible_hi` (von Mises posterior, same residual pool). Set
`method = "consensus"` to combine every engine bundled for a sample's
(species, tissue, platform) via a weighted circular mean instead of committing to one,
and `qc_method = "dynamic"` to flag low-quality samples with `clockCoherence()` (a
per-dataset gene panel) instead of the fixed 16-gene amplitude flag — useful when the
core clock genes themselves are not the ones affected (e.g. aging cohorts).

### Benchmark a method against known ground truth

```r
bm <- benchmarkPhase(se_truth,
  truth_col = "ZT",
  methods = c("timetable", "taufisher")
)
bm$summary   # one row per method: circular MAE, circular correlation, n
bm$residuals # one row per held-out sample
```

### Check whether circadian time is confounding a design

You don't have to be a circadian researcher to care about this one. If your samples'
collection times went unrecorded and they happen to be unbalanced across your
comparison — cases collected in the morning, controls in the afternoon — every
rhythmic gene will show an apparent condition effect that is really a clock effect.

```r
phaseConfounding(se, ~ condition)
#>        term        type  n n_groups statistic effect_size p_value
#> 1 condition categorical 20        2  339.6259    8.132343   0.001
```

`effect_size` is in hours: here the two groups' mean collection times sit about eight
hours apart. Phase that varies randomly *within* each group is not a confounder, and
the test says so. All comparisons are circular (ZT23 and ZT1 are two hours apart, not
twenty-two) and p-values come from permutation rather than parametric approximations.

## What's wrapped so far

| Method | Reference | Status |
|---|---|---|
| Molecular Timetable | Ueda et al. 2004, *PNAS* | Reimplemented from the method description (no maintained software release exists) |
| tauFisher | Duan, Ngo et al. 2024, *Nat Commun* | Wraps the authors' own R package (optional install, see above) |
| ZeitZeiger | Hughey, Hastie & Butte 2016, *Nucleic Acids Res* | Wraps the authors' own R package (optional install, see above) |
| TimeSignatR | Braun et al. 2018, *PNAS* | Wraps the authors' own R package (optional install, see above) |

Reference models currently bundled for `estimatePhase()` are trained on GSE54650
(Zhang et al. 2014, *PNAS*; mouse microarray) and GSE54651 (same study, RNA-seq).
The microarray set covers 8 of GSE54650's 12 tissues (Cer, Hrt, Hyp, Kid, Liv,
Lun, Mus, WFat) and two engines — Molecular Timetable and tauFisher, plus their
`"consensus"` combination — while the RNA-seq set bundles Molecular Timetable only
(select with `estimatePhase(platform = "rnaseq")`). ZeitZeiger and TimeSignatR are
wrapped and usable through their `fit*()`/`predict*()` functions and
`benchmarkPhase()`/`transferPhase()`, but are deliberately not bundled as
pre-trained `estimatePhase()` models (their fitted objects are much larger for no
measured accuracy benefit). Call `availableReferenceModels()` for the exact current
list. Additional labelled cohorts (human whole blood, mouse liver single-cell
pseudobulk) used for benchmarking and transfer testing, but not yet bundled as
trained reference models, are distributed via the companion
[`circaTimeData`](https://github.com/dralperenuysal/circaTimeData) package.

## License

[Artistic-2.0](https://opensource.org/licenses/Artistic-2.0), matching Bioconductor's
standard license requirements.

## Citation

A manuscript describing `circaTime`'s benchmarking results and uncertainty layer is in
preparation. Citation details will be added here once available.
