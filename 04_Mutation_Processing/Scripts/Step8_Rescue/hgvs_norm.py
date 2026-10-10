"""
hgvs_norm.py — shared HGVS query normalisation
==============================================
Prepares a legacy IDbases HGVS description for submission to an external
validator (VariantValidator, Mutalyzer).

Two corrections are applied, in order:

1. Reference-wrapper removal
   Legacy descriptions carry a genomic reference with a transcript selector,
   e.g. NG_007385.1(NM_000022.4):c.646G>A. VariantValidator cannot align an
   NM_ transcript against an NG_ RefSeqGene record and returns
   "TranscriptMapper(...): No transcript exons" for every such query. The
   transcript-only form NM_000022.4:c.646G>A validates normally. Both NC_ and
   NG_ wrappers are stripped.

2. Redundant deleted/duplicated bases
   HGVS recommendations from 2016 onward deprecate spelling out the deleted or
   duplicated sequence: c.2643delA is written c.2643del. VariantValidator
   rejects the legacy form with "VariantSyntaxError: Removing redundant
   reference bases from variant description" rather than validating it.
   delins is left intact, since the inserted sequence is required there; only
   the redundant deleted portion of a delins is removed.

This module deliberately performs no I/O and reads no environment variables,
so it can be imported from any pipeline script without side effects.
"""

import re

__all__ = [
    "normalize_hgvs_query",
    "strip_reference_wrapper",
    "normalize_deldup",
    "describe_normalization",
]

_WRAPPER_RE = re.compile(
    r"^N[CG]_[0-9]+\.[0-9]+\((N[MGR]_[0-9]+\.[0-9]+)\):([cnr]\..+)$"
)

_DELINS_RE = re.compile(r"del[ACGTUacgtu]+(ins[ACGTUacgtu]+)$")

_DELDUP_RE = re.compile(r"(del|dup)[ACGTUacgtu]+$")


def strip_reference_wrapper(hgvs):
    s = (hgvs or "").strip()
    m = _WRAPPER_RE.match(s)
    if m:
        return f"{m.group(1)}:{m.group(2)}"
    return s


def normalize_deldup(hgvs):
    s = (hgvs or "").strip()
    m = _DELINS_RE.search(s)
    if m:
        return s[: m.start()] + "del" + m.group(1)
    m = _DELDUP_RE.search(s)
    if m:
        return s[: m.start()] + m.group(1).lower()
    return s


def normalize_hgvs_query(hgvs):
    return normalize_deldup(strip_reference_wrapper(hgvs))


def describe_normalization(hgvs):
    """
    Return a semicolon-separated audit label for the corrections applied.

    'deldup_shortened' is scientifically load-bearing and must not be treated as
    a cosmetic fix. Dropping the deleted sequence from c.1293delC discards the
    source record's assertion about *which* base was deleted while retaining its
    assertion about the position. Where that base disagreed with the reference,
    a validator will now accept the description and return the reference base
    instead. The position is being trusted over the allele, which is a curation
    decision that belongs in manual review rather than in a silent merge.
    """
    raw = (hgvs or "").strip()
    stripped = strip_reference_wrapper(raw)
    final = normalize_deldup(stripped)
    applied = []
    if stripped != raw:
        applied.append("wrapper_stripped")
    if final != stripped:
        applied.append("deldup_shortened")
    return ";".join(applied)
