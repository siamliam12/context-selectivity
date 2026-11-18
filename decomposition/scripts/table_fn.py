# The table function (D1): builds every trait-by-tissue table in Parts D and E.
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

def trait_table(z, membership):
    """membership = rows of (variant, label); a variant may sit under several traits.
    Returns label, group, median_z, n_variants."""
    zl = z.stack().rename("z").reset_index()
    m = membership[["variant", "label"]].merge(zl, on="variant")
    return (m.groupby(["label", "group"])["z"]
             .agg(median_z="median", n_variants="size").reset_index())

def build_table(long, membership, ddof=0):
    z, n_dropped = variant_z(group_means(long), ddof)
    return trait_table(z, membership), n_dropped