"""
source_change.py
================
Maps pipeline rows back to the source change they describe.

A source change is one position in an allele's notation lists: the i-th genomic
token together with the i-th coding token of the same allele in all_mutations.tsv.
Its id is GENE|accession|allele|genomic token|coding token ("-" for a missing
level). The same token pair on both alleles of one record (homozygous) gets the
first allele, as the NM tracks already do when they keep the first ok row.

Track A rows carry the BED-sanitised genomic notation as sysname, so they are
looked up by that; every other track carries the coding notation and its allele.
"""

import csv
import re

TRACK_A = "NG_IDRefseq"


def sanitise(notation):
    return re.sub(r"[^A-Za-z0-9_.>\-]", "", re.sub(r"[,\s;]+", "_", notation))


def allele_text(value):
    v = str(value or "").strip()
    return v[:-2] if v.endswith(".0") else v


class SourceChanges:
    def __init__(self, all_mutations_tsv):
        records = {}
        with open(all_mutations_tsv, encoding="utf-8") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                level = {"genomic": "g", "coding": "c"}.get(row["variant_type"].strip())
                if level is None:
                    continue
                key = (row["gene"].strip().upper(), row["accession"].strip())
                alleles = records.setdefault(key, {})
                alleles.setdefault(allele_text(row["allele_num"]), {"g": [], "c": []})[level].append(row["notation"].strip())
        self.by_genomic = {}
        self.by_coding = {}
        for (gene, acc), alleles in records.items():
            canonical = {}
            for allele, levels in alleles.items():
                g_tokens, c_tokens = levels["g"], levels["c"]
                for i in range(max(len(g_tokens), len(c_tokens))):
                    g_tok = g_tokens[i] if i < len(g_tokens) else "-"
                    c_tok = c_tokens[i] if i < len(c_tokens) else "-"
                    first = canonical.setdefault((g_tok, c_tok), allele)
                    entry = (f"{gene}|{acc}|{first}|{g_tok}|{c_tok}", first)
                    if g_tok != "-":
                        self.by_genomic.setdefault((gene, acc, sanitise(g_tok)), entry)
                    if c_tok != "-":
                        self.by_coding.setdefault((gene, acc, allele, c_tok), entry)
                        self.by_coding.setdefault((gene, acc, "", c_tok), entry)

    def lookup(self, gene, accession, allele, sysname, track):
        gene, accession, sysname = str(gene).strip().upper(), str(accession).strip(), str(sysname).strip()
        if track == TRACK_A:
            return self.by_genomic.get((gene, accession, sysname), ("", ""))
        return (self.by_coding.get((gene, accession, allele_text(allele), sysname))
                or self.by_coding.get((gene, accession, "", sysname), ("", "")))


def split_id(source_change_id):
    parts = source_change_id.split("|", 4)
    return parts if len(parts) == 5 else ["", "", "", "", ""]
