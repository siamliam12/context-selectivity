# Part II Pre-registration

## 1. About this pre-registration

| Item | Entry |
| --- | --- |
| Project | Context selectivity, Part II: decomposition of a published tissue-attribution figure |
| Authors | MD Siam Ahmed(me), Dr. Ziliang Zong (PI) |
| Draft written | 16 Nov 2026 |
| Frozen on | To be completed at the freeze commit |
| Commit hash | To be completed |
| External timestamp | OSF registration link, to be completed |

**Question.** How much of the tissue pattern in Sniff Fig. 4C can be rebuilt from where genes are normally expressed at each locus, with no information about the variant's alleles?

**What has been seen before freezing.** The variant lists, positions and letters; the score file's track names; GTEx and ENCODE sample labels and gene expression; and, to check the raw-to-standardized score conversion, the scores of 20,000 random variants that are not in the figure. No trait-by-tissue table has been built, and the scores of the 5,704 figure variants have not been read.

**Scope.** This document covers Gate 1 (reproducing the figure), the reference-only rebuild (Part E: E1 measured expression, E2 the model's own reference prediction) and the genome-wide check (D5). The allele-substitution test has its own pre-registration. Part I (the measured eQTL reference) is reported separately.

## 2. Target figure and data

- **Figure.** Fig. 4C of Srivastava et al., Sniff, bioRxiv v1, 14 July 2025, doi 10.1101/2025.07.09.663936. Checked on 16 Nov 2026: only v1 is public, no revised version was found, and papers posted in 2026 still cite it as a preprint. The target is therefore preprint Fig. 4C. If a revised figure becomes public before any table is built, this document is amended first.
- **Scores.** The released Borzoi scores, scores\_indel\_snp\_distmatched.h5 (standardized), from one replicate (f0), L2 score (logD2) over the 196,608-base window around the variant, indels scored with the stitch option and then quantile-normalized to the SNP distribution. The file's positions are hg19, and the scores were computed on hg19 (every REF letter fits hg19; 25% fit hg38 at the same positions, which is chance).
- **Tracks.** 89 RNA tracks: 87 GTEx tracks plus 2 GTEx-made K562 tracks labelled bone\_marrow. They are score columns 6296 to 6384 and form 31 tissue groups (track\_map.csv).
- **Fine-mapping.** 660 Sniff and baseline PolyFun-SuSiE files for 15 UK Biobank traits (gs://seqnn-share/sniff/sniff\_results\_ukb\_15/), checksummed.
- **Trait labels.** Taken from the file names (trait\_map.csv), not from supplementary table S2, whose trait column is shifted against its rows.

## 3. Variants

- **Differential PIP** = Sniff PIP minus baseline PIP, per variant and trait.
- **Set A (promoted):** difference above 0.25. 3,230 variant-trait rows, 3,160 unique variants; identical to supplementary table S3.
- **Set B (promoted plus demoted):** set A plus differences below -0.25. 5,830 rows, 5,704 unique variants.
- **Set C (Sniff-prioritized):** the set named in the panel's own title and defined in the paper's Methods: PIP below 0.1 under both SuSiE and PolyFun-Baseline, and Sniff raises it by at least 0.25. Built from S3 (set\_c\_variants.csv), a subset of set A: 1,172 variant-trait rows and 1,158 unique variants, indels included, exactly the paper's count.
- A variant found in two traits is kept under each trait.
- **Indels** are kept in Gate 1 and Part E. The allele test uses SNPs only (decision P11).
- **Positions.** Each variant's scores are read from its row in the score file (score\_row in variants.csv), found by hg19 position plus letters. Every script that uses positions runs the permanent position check (coord\_check.py) on every variant.

## 4. Building the tables (the table function)

One function builds every trait-by-tissue table in Parts D and E. For a given set of variants:

1. **Per variant:** average the scores of the tracks in each of the 31 groups (mean).
2. **Per variant:** Z-score the 31 group values (subtract their mean, divide by their SD).
3. **Per trait and group:** take the median Z across the trait's variants.

The result is a 15 x 31 table. A trait's **labelled tissues** are the groups with median Z above 0.5.

- **SD convention.** For Gate 1, tables are built both ways (dividing by n, and by n - 1), because the authors do not say which they used. Every new table divides by n.
- **Zero spread.** A variant whose 31 values are all equal is left out of that table, and the count is reported.

## 5. Gate 1: reproducing the figure

**The printed values.** The paper's text gives six values. Four agree with the figure:

| Trait | Tissue group | Printed value |
| --- | --- | --- |
| Cholesterol | liver | 2.42 |
| Haemoglobin (MCH) | bone\_marrow | 1.08 |
| Monocyte count (MC) | blood | 1.48 |
| Diastolic BP | heart | 0.59 |

Two disagree with the figure:

- Text: reticulocyte count (HLSRC) with bone marrow, 0.78. The figure shows HLSRC's bone-marrow dot at about 0.82; the dot at 0.78 is MPV's bone marrow.
- Text: monocyte count (MC) with spleen, 0.50. The figure shows MC's spleen dot at about 0.83; MC's dot at 0.50 is lung, and the spleen dot at 0.50 is MPV's.

The figure is the target, so these two cells are checked against the figure readings. The mismatch between text and figure is reported as a finding about the paper, whatever Gate 1 shows.

**The figure readings.** Read by eye from the published image on 16 Nov 2026, before any table was built, and frozen (fig4c\_labels.csv). Coloured dots are the tissues the authors labelled (median Z above 0.5): 17 dots in 11 traits, each with its approximate height. Height, FEV:FVC, FVC and BMD have none. Two dots sit on the 0.5 line (MPV spleen, MC lung); they are ignored in the label match, as are groups within 0.02 of the line in the rebuilt table.

**What is built.** Six tables: sets A, B and C, each with both SD conventions.

| Label | Meaning |
| --- | --- |
| Exact | One table hits the four agreeing values within ±0.02, every labelled dot's height within ±0.05 of the figure reading, and the labelled tissues for all 15 traits |
| Partial | The four values match, but labelled tissues differ in 1 to 3 traits or 1 to 2 dot heights miss; or tissues and heights all match, but 1 of the four values misses |
| Unresolved | Anything worse |

- If more than one table qualifies: prefer set A, then set B, then set C; then n over n - 1. The chosen set and SD convention are used for every later table.
- **Exact or Partial:** Part E runs. After Partial, the gaps are stated with every later result.
- **Unresolved = non-reproduction.** It means the authors' released scores and stated method do not reproduce the figure. It is not a claim that the authors made an error. The cause is looked for in a fixed order: variant set; then tissue grouping and the Z steps; then positions. This includes a check using S2's shifted trait names. At least one hour is spent looking for a bug of our own; then the authors are emailed; then it is written up as a finding. Part E does not run. Part I and the genome-wide check still do.
- If the authors say which set they used before Gate 1 is run, their answer decides.

## 6. The reference-only rebuild (Part E)

Two versions are built. In both, the variant's alleles play no part; only its location is used. Each runs only after Gate 1 is Exact or Partial, using the set chosen there.

**E1, measured expression**

- **Genes.** All genes overlapping 98,304 bases either side of the variant's hg38 position (the score's 196,608-base window), from GTEx's v8 gene model (GENCODE v26), any gene type. Check version: the nearest gene start only.
- **Expression per group.** GTEx v8 median TPM, averaged over each group's sub-tissues with weights matching its tracks (tissue\_map.csv). bone\_marrow uses the K562 column from 8 ENCODE polyA-plus experiments x 2 replicates (GENCODE V29). Check version: GTEx v10.
- **Per locus and group:** sum TPM over the window's genes, then log(x + 1).
- **Rescaling:** each group is standardized (mean 0, SD 1) across the genome-wide sample loci (section 10), mirroring how the authors standardized each track across all variants.
- Then table-function steps 2 and 3: Z across the 31 groups per variant, median per trait.

**E2, the model's own reference prediction**

- Borzoi replicate f0, run on the hg19 reference sequence only, with the release's options (reverse complement, no shifts), for the 89 tracks (model outputs 7522 onwards).
- **Per locus and track:** the predicted coverage summed over the output window, then log(x + 1).
- **Rescaling:** each track is standardized across the Borzoi reference predictions for a seeded random 10,000 of the genome-wide sample loci (section 10).
- Then table-function steps 1 to 3.

## 7. Reading the result

Each version (E1, E2) is read separately against the published-style table.

- **Measure, per trait:** overlap = groups labelled in both tables divided by groups labelled in either. One table names none: 0. Neither names any: 1.
- **Resampling unit:** trait x LD block. The genome is divided into about 1,700 approximately independent LD blocks (Berisa and Pickrell 2016, European, hg19), and each variant belongs to the block containing its hg19 position.
- **Half-size reading (main).** 1,000 random splits of each trait's units into two halves. For each split:
  - **main** = overlap of the reference-only table from half 1 with the published-style table from half 2;
  - **self-agreement** = overlap of the published-style tables from half 1 and half 2;
  - **luck** = overlap of a reference-only table from half 1, built after shuffling units across traits (keeping each trait's number of units), with the published-style table from half 2.
- **Intervals:** 95% percentile intervals over the 1,000 values.
- **Readings, per trait:**

| Main interval | Reading |
| --- | --- |
| Overlaps the luck interval | Explains little |
| Above luck, and overlaps the self-agreement interval | Explains essentially all |
| In between | Explains part, reported as a position between the two lines |

- **Overall:** the number of traits in each reading is reported; the most common reading is the overall reading used to check the predictions.
- **Also reported:** the full-size overlap, as a description only.
- **Secondary:** per trait, the correlation between the two halves' 31-group profiles, corrected with the Spearman-Brown formula; and the same for reference-only versus published.
- Every trait is reported. Traits with few units are marked, not dropped, and each trait's number of units is reported.

## 8. The authors' six pairs

The relevance standard is the authors' own claims: the six trait-tissue pairs printed in the paper (author\_pairs.csv, the same six values as Gate 1). No independent human labelling is done.

- For each of E1 and E2, each pair is reported as yes or no: does that trait's group have median Z above 0.5 in the full-size reference-only table?
- The count out of 6 is reported beside the count expected by chance: for each pair, the share of that trait's 31 groups that are labelled in the reference-only table, summed over the six pairs.
- Two of the six values in the paper's text disagree with the figure (section 5). The pairs themselves (HLSRC with bone marrow, MC with spleen) are labelled in the figure, so the pair check is unaffected.
- When writing up, each pair is backed by standard physiology citations. bone\_marrow is described as the K562 cell line standing in for bone marrow, not as bone marrow tissue.

## 9. Cell-level test (Part E)

For each of the 465 cells (15 traits x 31 groups) of each reference-only table:

- Units are shuffled across traits, keeping each trait's number of units, 10,000 times.
- p = (1 + number of shuffles with median Z at least as high as observed) / (1 + 10,000).
- Benjamini-Hochberg correction at 5% across the 465 cells.

## 10. Genome-wide check (D5)

This runs whatever the Gate 1 outcome, because it uses only the released scores and expression.

- **Sample.** 50,000 variants drawn at random from the score file's 19,534,182 rows, leaving out the 5,704 figure variants; SNPs and indels. Each is converted to hg38 with the C2 method; failures are counted and dropped. The position check runs on every variant.
- **Score profile.** The standardized scores for the 89 tracks, averaged into the 31 groups, then Z-scored across the groups.
- **Expression profile.** The window genes and steps of E1 (section 6): sum TPM, log(x + 1), each group standardized across this sample, then Z-scored across the groups. Variants with zero spread (for example, no gene in the window) are left out and counted.
- **Per variant:** the Spearman correlation between the two 31-group profiles. **Main statistic:** the median across variants.
- **Null.** Each variant is paired with the expression profile of another variant on a different chromosome and in the same tenth of overall expression level (mean log expression across groups), repeated 1,000 times; each round gives one median.
- **Main result:** the real median minus the mean of the null medians, with a 95% interval from redrawing LD blocks 1,000 times.
- **Reading:** an interval above 0 means score profiles follow their own locus's expression beyond the general tissue pattern. The size is reported. There is no pass or fail.
- **Splits fixed in advance:** the top 10% of variants by their largest absolute group score, versus the rest; variants inside a gene body versus outside; SNPs versus indels.
- **Extra:** the share of variants whose top-scoring group is also their top-expression group, against the same share under the null.
- **Check versions:** nearest gene only; GTEx v10.

## 11. Predictions

Written before any result. They do not change any method or reading; each is checked and reported afterwards, right or wrong.

| # | Question | Prediction | Main reason |
| --- | --- | --- | --- |
| 1 | Which of sets A, B and C reproduces Fig. 4C? | Set A (promoted only) | S3, titled "differentially prioritized", holds promoted variants only; the abstract and the 4D example describe promoted variants |
| 2 | E1 (measured expression) reading | Explains part | Genes near trait loci are preferentially expressed in trait-relevant tissues (Finucane et al. 2018); GTEx medians and the ENCODE K562 source differ from the model's training samples |
| 3 | E2 (model reference) reading | Explains essentially all | Same model and tracks; models track gene-level expression well but variant effects between individuals poorly (Huang et al. 2023); the log(x + 1) score gives larger values where coverage is high |
| 4 | Author pairs above 0.5 in E1 | 4 of 6 | Likely: cholesterol-liver, both red-cell traits with bone\_marrow (K562), monocytes-blood; spleen and heart sit nearer the line |
| 5 | Genome-wide: real minus chance above 0 | Yes, clearly | The log(x + 1) score mechanism; a variant changes expression only where its gene is expressed |
| 6 | Genome-wide: top 10% of scores vs rest | Stronger in the top 10% | Small scores are mostly noise |
| 7 | Top-scoring group equals top-expression group more than chance | Yes | Follows from 5 |

## 12. What each outcome permits

| Outcome | Allowed statement | Not allowed |
| --- | --- | --- |
| Reference-only recovers most of the pattern | The published tissue assignment largely reflects where the locus is expressed, not an effect of the variant | Any claim that the model ignores tissue, or about the paper's fine-mapping, laboratory validation or replication |
| Reference-only recovers little | The pattern carries information beyond where the locus is expressed; the allele test is needed to say whether it depends on the alleles | Any claim about which tissue is biologically correct |
| Published values not reproduced (Gate 1 Unresolved) | A computational reproducibility result, reported as such | Any decomposition verdict |

The project succeeds if every analysis ends in a defensible, pre-registered outcome, including non-reproduction. It does not succeed or fail on whether a published claim is overturned.

## 13. Amendments, seeds and software

- **Amendments.** Any change after the freeze is appended as a dated amendment with its reason. The original text is never edited.
- **Seeds.** Master seed 20261116 for every random draw (splits, shuffles, samples, nulls). Each analysis derives its own seed from it in a fixed way, written in its script.
- **Software.** The Python environment "decomposition" (versions recorded in its environment file); Borzoi through baskerville 1.0.0; all scripts in decomposition/scripts/, each writing its output to decomposition/logs/.
- **Groups for later analyses.** Mutation class, CpG status and GC content (101-base window), computed in hg19, are frozen in variants.csv for the allele test.

## 14. Frozen files and checksums

All in decomposition/data/frozen/, read-only, with fingerprints in frozen.sha256:

| File | Contents |
| --- | --- |
| trait\_map.csv | 15 traits: file code, figure label, full name |
| variants.csv | 5,704 variants: hg19 and hg38 positions, the authors' letters, score\_row, groups |
| variant\_traits.csv | 5,830 variant-trait rows (sets A and B) |
| set\_c\_variants.csv | Set C, the Sniff-prioritized variants, built from S3 |
| track\_map.csv | 89 tracks: score column, group, model output |
| tissue\_map.csv | 42 sub-tissue rows with weights for the 31 groups |
| author\_pairs.csv | The authors' six trait-tissue pairs and the values printed in the text |
| fig4c\_labels.csv | The figure's 17 labelled dots with approximate heights, read by eye on 16 Nov 2026 |

Raw inputs (fine-mapping files, genomes, GTEx and ENCODE files) are fingerprinted in finemap.sha256, genome.sha256, gtex.sha256 and encode.sha256.

**Checksums at the freeze:** the full contents of frozen.sha256 are pasted here at the freeze commit.

```
103fa6d089423fb6c50047cfba5c636cbd76501bddd71fbe3c06931b973059e1  trait_map.csv
dc0a0936e51dbd7fb89e38ae51f1cd0c694de864c067b2317949d5612941c073  variants.csv
815301389141404f0af9ffb5bbde514b1dedd34fbf6512e8e7d6b673425872e5  variant_traits.csv
87f87b4d8febce15435a696de1b59080cf1836362f228ee6cd0af0d7763c60a6  track_map.csv
b5f407b18e3ae16690d8f0e6a5c01deb2e40b0b9eb9d8d47c840313621deb513  tissue_map.csv
b96c8d3b2fe89d9d492542dc5bf63448787cfd31443b6ff091f1f1ff7c97e9d6  author_pairs.csv
289a5d1bb5bd4c8d3681611f7b8d28e2600bf6e8d414621231e2c22b520aa180  fig4c_labels.csv
3e9d3712bc227b9f214924bec7c7c8b6715fa36705f712a5a768a9db75aa4090  set_c_variants.csv
```
