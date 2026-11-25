# The table function (D1; mean added by Amendment 1): builds every trait-by-tissue table in Parts D and E.
import pandas as pd

def group_means(long):
    """Rows of (variant, group, score) -> one row per variant, one column per group: the mean score."""
    return long.groupby(["variant", "group"])["score"].mean().unstack("group")

def variant_z(gm, ddof=0):
    """Z-score each variant across its groups. ddof=0 divides by n, ddof=1 by n-1.
    Variants whose values are all equal are dropped; returns how many."""
    mu, sd = gm.mean(axis=1), gm.std(axis=1, ddof=ddof)
    keep = sd > 0
    z = gm[keep].sub(mu[keep], axis=0).div(sd[keep], axis=0)
    return z, int((~keep).sum())

def trait_table(z, membership, agg="mean"):
    """membership = rows of (variant, label). agg = 'mean' (Amendment 1) or 'median' (the caption's wording).
    Returns label, group, trait_z, n_variants."""
    zl = z.stack().rename("z").reset_index()
    m = membership[["variant", "label"]].merge(zl, on="variant")
    return (m.groupby(["label", "group"])["z"]
             .agg(trait_z=agg, n_variants="size").reset_index())

def build_table(long, membership, ddof=0, agg="mean"):
    z, n_dropped = variant_z(group_means(long), ddof)
    return trait_table(z, membership, agg), n_dropped