# Reference identity selection in OpenASChem 1.3.0

## What the default mode does

1. Extract the target monoisotopic XIC and align the three standard replicates
   using the existing Python/Numba backend. Enumerate local maxima with relative
   prominence >=0.10; the large enumeration cap is 10000, not a top3 strategy.
2. Apply existing standard filters to each candidate: each standard quantity
   above the configured minimum, each half-prominence width <=10 seconds, and
   mean standard boundaries strictly within 50–300 seconds.
3. Sum the individual smoothed apex heights across the three standards. Define
   relative height against the strongest accepted candidate of that compound.
   Retain candidates with relative height >=0.5 by default. This denominator is
   calculated before isotope vetoes; weak alternatives are not promoted merely
   because a dominant interference was removed.
4. Collect raw intensity-weighted mass centroids within each candidate's standard
   boundaries and quantify carbon-13, paired carbon-13, chlorine-37, bromine-81
   and sulfur-34 channels. Extraction uses the configured m/z unit and tolerance;
   the tested defaults are positive mode, singly charged [M+H]+ and ±3 ppm.
5. Calculate isotope mismatch from absolute log-ratio deviations with a 0.005
   stabilizer. Carbon-13 is always scored. Cl/Br/S channels are scored when the
   theoretical ratio is >=0.025; queried channel overlap is accounted for.
   Summarize over the three standards by the median. Paired carbon-13 and the
   channel one carbon-isotope increment below the target also support vetoes.
6. Estimate the mass-error center as the median of candidate centroid errors
   among candidates with relative height >=0.8 in that pool. Score a candidate as
   `-isotope_loss -abs(mass_error_ppm-center_ppm) +ln(relative_height)`.
   All three coefficients are 1. There is no additional hard mass-error cutoff.
7. Check M+1 explanations using the parent library, independently of its pool
   labels. The parent must have the expected mass within the extraction window;
   parent/target median chromatographic correlation must be >=0.98, target/parent
   and next-isotope/target ratios must agree within factors 1.3 and 1.4, and the
   parent-isotope mass explanation must have residual <=1 ppm and improve on the
   target formula's residual by >=1 ppm. Such a candidate is vetoed.
8. Choose the highest-scoring eligible candidate. If a competing peak is >5 s
   away, has >=half the chosen peak's height and a score gap <0.15, defer the
   compound as ambiguous instead of making an automatic assignment.
9. Match the chosen reference to experimental runs and apply the existing
   quantification, background, FC, replicate and statistical filters. The score
   and veto do not inspect experimental FC, target names or pool membership.

## Parent libraries and formula checks

By default, parent candidates come from the supplied analysis library. Supply
`--identity_parent_library` to search a larger library without adding its
compounds to targeted hit screening. Required columns are `formula` and either
`SGC ID for Component` or the historical `SGC ID for Compotent` spelling.
`MONOISOTOPIC_MASS`, if supplied, must be the neutral monoisotopic mass. Otherwise
neutral parent masses are calculated from formulas. Other columns are optional.
The analysis library retains the normal OpenASChem input schema and mass
convention; source masses are never automatically corrected.

The library check compares the analysis library's neutral mass with its formula.
A discrepancy >0.001 Da is marked for review, not automatically treated as an
error or used to overwrite the input; salts and charged representations need
chemical interpretation. Duplicate or invalid parent formulas are handled before
selection; invalid formulas produce a clear input error.

## Outputs and interpretation

The main report and result CSVs contain accepted references and final hits.
The companion review HTML contains all qualifying candidate peaks with their
selection status, and only final accepted hits carry hit flags. The decisions
CSV also lists compounds without an accepted standard peak. No-standard-peak,
isotope-excluded and ambiguous decisions do not establish chemical absence.

Selection annotations and `ReferencePeakID` preserve the original multi-digit
candidate identity through CSV export and HTML. Reference rank is an identifier
within the enumerated candidate set, not an identity-score rank. Scores are
heuristics rather than calibrated probabilities, posterior odds or FDR estimates.
Missing isotope signal can penalize a genuine weak compound, while a mass-matched
fragment or isomer may still be indistinguishable from the intended compound.

The defaults were optimized on one mixed-standard benchmark. Changing the
ambiguity margin or standard evidence can change the accepted list. A benchmark
with zero B-pool calls partly reflects deferral of uncertain assignments, not
proof of zero identification error. Transfer to the tested real AS-MS dataset
preserved the previous hit list, but does not establish improved specificity
without independent truth labels. Negative acquisition and Da windows are
supported by the extraction code; the published pilot counts refer to positive
mode at ±3 ppm. Area quantification retains the corrected original-RT integration.

Use `--peak_selection top1` for the previous standard peak selection behavior.
This changes reference assignment only; it does not restore removed CV filtering,
the former index-based width cutoff, R dependencies, or automatic static PDFs.
