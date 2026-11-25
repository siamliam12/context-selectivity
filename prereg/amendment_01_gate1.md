# Amendment 1 to the Part II pre-registration (25 Nov 2026)

Frozen text: prereg/partII_prereg.md (commit da8bdae). This amendment adds to it; the original text is unchanged.

## What happened
- Gate 1 was run as pre-registered on 25 Nov 2026 (sets A, B and C; SD by n and n-1). No table met the frozen criteria, so Gate 1 = Unresolved (decomposition/logs/gate1.log).
- Following the pre-registered diagnosis, we tested candidate reconstructions of the authors' method: score file (standardized, raw, re-standardized), trait-level summary (median, mean, weighted mean, mean per credible set), track handling, SD convention, and variant sets (cut-offs 0.1 to 0.5; promoted, demoted, prioritized; credible-set membership; SNPs only). Every try is logged in decomposition/results/d3/attempts_detail.csv.
- We searched every released folder (borzoi_veps, processed_borzoi_parquets, _20260810_tmp_files) for the figure's inputs. None were found.

## What we found
1. The figure follows the mean of the variants' Z-scores, not the median stated in the caption. With the median, at most 9 of 15 traits' labelled tissues match; with the mean, up to 15 of 15.
2. The closest variant set is the differential variants (|dPIP| > 0.25) that lie in a Sniff 95% credible set. With it, the four printed values that agree with the figure are reproduced within 0.03 (2.41, 1.06, 1.50, 0.56 against 2.42, 1.08, 1.48, 0.59), every trait has the same top tissue, and 12 of 15 traits' labelled tissues match.
3. The remaining gaps are trait-specific (largest: MPV, 0.2 to 0.3 too low) and move between traits as settings change. The likeliest cause is that the figure was made from a fine-mapping run that differs from the released files; the authors' temporary folder shows runs with a fixed and with a trait-specific GWAS sample size.
4. Two values printed in the text do not match the figure (section 5).

## Decision
- The figure is reported as not reproduced under the frozen criteria, with the reconstructed method and its remaining gaps reported in full.
- Part E proceeds with the reconstructed method as the published-style table: mean across a trait's variants; set B restricted to variants in a Sniff credible set (fig4c_reconstructed_set.csv); SD dividing by n; everything else as pre-registered.
- Robustness: every Part E result is also computed with two nearby reconstructions: the mean over set A (all 15 traits' labels match), and set B in Sniff credible sets averaged per credible set (closest overall to the figure's dots). A conclusion is drawn only if it holds across all three.
- All other rules, readings and predictions are unchanged. Prediction 1 (set A) is recorded as not confirmed.
