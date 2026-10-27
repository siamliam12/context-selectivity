README

Files

1. scores.h5 
This file contains Borzoi variant effect predictions across 6,385 functional genomic datasets for 19,534,182 autosomal SNPs and indels. 
Variants were scored with the L2 score (specified as the logD2 score in the Borzoi codebase). 
Indels were scored using the "--stitch" option, which implements a stitched scoring strategy. 
The command used to score variants is: 

# python hound_snp.py \
         --float16 \
         -u \
         --stats logD2 \
         --shifts 0 \
         --rc \
         --indel stitch \
         --require gpu \
         <params_file> \
         <model_file> \
         <vcf_file>

We used a single Borzoi replicate to score all variants, the model was accessed via: https://storage.googleapis.com/seqnn-share/borzoi/f0/model0_best.h5
The parameter file required to score variants was accessed via: https://storage.googleapis.com/seqnn-share/borzoi/params.json

# The keys in the h5 file include: 
# 'chr'
# 'pos'
# 'snp'
# 'ref_allele'
# 'alt_allele'
# 'logD2'
# 'target_ids'
# 'target_labels'

# The 'logD2' key maps to a 19,534,182 x 6,385 array containing Borzoi L2 scores, computed as described in the manuscript. Each of the 6,385 columns represent scores for a single dataset, corresponding to datasets listed in arrays "target_ids" and "target_labels". No additional normalizations have been performed on this data.


2. scores_indel_snp_distmatched.h5 
This file contains Borzoi variant effect predictions in which the distributions of indels and SNPs have been aligned (see manuscript for details, we recommend using this score set for tasks such as fine mapping which jointly consider indels and SNPs). The scores have further been standardized (mean=0, variance=1). The metadata for this file is the same as the metadata for scores.h5 (described above).


For processed data (PCs used in Sniff finemapping), please refer to https://storage.googleapis.com/seqnn-share/sniff/borzoi_102_annotation_set
