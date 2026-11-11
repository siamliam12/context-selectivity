# Tests the permanent position check: it must pass on correct data and fail on broken data.
import pandas as pd
from coord_check import check_table

v = pd.read_csv("decomposition/data/interim/variant_context_hg38.csv")

# 1. Correct positions: must pass.
check_table(v, label="real table")

# 2. Every position moved by one base: must fail.
shifted = v.copy()
shifted["BP_hg38"] = shifted["BP_hg38"] + 1
try:
    check_table(shifted, label="shifted by +1")
    print("PROBLEM: the shifted table passed; the check is not working")
except AssertionError as e:
    print("shifted by +1 correctly failed:", str(e)[:120])

# 3. REF and ALT swapped (SNPs only): must fail.
swapped = v[~v["is_indel"]].copy()
swapped[["REF", "ALT"]] = swapped[["ALT", "REF"]].values
try:
    check_table(swapped, label="REF/ALT swapped")
    print("PROBLEM: the swapped table passed; the check is not working")
except AssertionError as e:
    print("REF/ALT swapped correctly failed:", str(e)[:120])