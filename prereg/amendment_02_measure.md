# Amendment 2 to the Part II pre-registration (25 Nov 2026)

Frozen text: prereg/partII_prereg.md (commit da8bdae); Amendment 1 (commit 26af051). The original text is unchanged.

## Why
- D4 measured the self-agreement of the published-style table (reconstructed method, 1,000 half-size splits by trait x LD block; decomposition/results/d4/). The labelled-tissue overlap is unstable: for 8 of 15 traits its 95% interval spans 0 to 1 or 0 to 0.5 (Age at Menarche, BMD, Diastolic BP, FEV:FVC, HLSRC, Hair pigment, Chronotype, FVC). A measure that varies this much between halves of the authors' own data cannot separate the readings for those traits.
- The correlation between the two halves' 31-tissue profiles is stable: Spearman-Brown corrected median 0.80 to 0.98 for 13 of 15 traits (Chronotype 0.52 and FVC 0.51, with 20 and 41 units).
- This decision uses only the published-style table. No reference-only table has been built.

## Change
- Main measure (section 7): per trait, the Pearson correlation between the two tables' 31-tissue profiles. Everything else in section 7 is unchanged: the same 1,000 splits (decomposition/results/d4/splits.parquet); main = reference-only half 1 vs published-style half 2; self-agreement = published-style half 1 vs half 2; luck = reference-only half 1, built after shuffling units across traits, vs published-style half 2; 95% percentile intervals; and the three readings.
- The labelled-tissue overlap is kept as a secondary measure and reported in full for every trait.
- The Spearman-Brown corrected correlation is reported as a full-size description.
- Predictions 2 and 3 are checked against the readings of the main measure.
