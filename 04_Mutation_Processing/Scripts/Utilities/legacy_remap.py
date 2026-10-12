"""
legacy_remap.py
===============
Legacy c. numbering: test whether a source change reproduces its source protein on a
transcript at a shifted position, and remap source c. notations with a per-change table.

Test
----
A source c. change (substitution, deletion, duplication, insertion, delins) placed at
c.(N + k) on the transcript CDS reconciles with the source protein when

    substitution   the reference base matches, the codon translates to the source
                   reference and alternative residues, and the codon number equals the
                   source protein position                         (evidence: substitution)
    frameshift     the first changed residue is the source position and residue
                                                                   (evidence: frameshift_start)
    in-frame       the source position lies within the changed span and carries the
                   source residue                                  (evidence: inframe_indel)

Remap table
-----------
legacy_remap.tsv, one row per (gene, accession, source c.): the remapped c. on the
MANE transcript. Built by Utilities/build_legacy_remap.py from data; read here with
load_remap and applied with remap.
"""

import csv
import re
from pathlib import Path

_BASES = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: _AA[16 * i + 4 * j + k]
         for i, a in enumerate(_BASES) for j, b in enumerate(_BASES) for k, c in enumerate(_BASES)}
THREE = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
         "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
         "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "*"}
RE_C_SUB = re.compile(r"^c\.(\d+)([ACGT])>([ACGT])$", re.I)
RE_C_INDEL = re.compile(r"^c\.(\d+)(?:_(\d+))?(delins|del|dup|ins)([ACGT]*)$", re.I)
RE_P = re.compile(r"^p\.\(?([A-Z][a-z]{2}|[A-Z*])(\d+)(.*?)\)?$")
RE_POS = re.compile(r"(?P<sign>[-*]?)(?P<num>\d+)(?P<intron>[+-]\d+)?")
REMAP_COLS = ["gene", "accession", "source_c", "remapped_c", "shift", "evidence", "mane_nm"]


def one_letter(res):
    if res in ("X", "*", "Ter"):
        return "*"
    return THREE.get(res) or (res if len(res) == 1 else None)


def translate(seq):
    out = []
    for i in range(0, len(seq) - 2, 3):
        aa = CODON.get(seq[i:i + 3], "X")
        out.append(aa)
        if aa == "*":
            break
    return "".join(out)


def parse_c(token):
    t = token.strip().replace(" ", "")
    m = RE_C_SUB.match(t)
    if m:
        return "sub", int(m.group(1)), int(m.group(1)), m.group(2).upper(), m.group(3).upper()
    m = RE_C_INDEL.match(t)
    if m:
        a = int(m.group(1))
        return m.group(3).lower(), a, int(m.group(2) or a), m.group(4).upper(), ""
    return None


def parse_p(protein):
    m = RE_P.match(protein.strip())
    if not m:
        return None
    ref, rest = one_letter(m.group(1)), m.group(3)
    if ref is None:
        return None
    if "fs" in rest:
        return "fs", int(m.group(2)), ref, None
    if re.search(r"del|dup|ins", rest):
        return "inframe", int(m.group(2)), ref, None
    alt = ref if rest in ("=", "") else one_letter(rest)
    if alt is None:
        return None
    return "sub", int(m.group(2)), ref, alt


def mutate(seq, cds0, change, k):
    kind, a, b, s1, s2 = change
    i, j = cds0 + a - 1 + k, cds0 + b - 1 + k
    if i < cds0 or j >= len(seq):
        return None
    if kind == "sub":
        return seq[:i] + s2 + seq[i + 1:] if seq[i] == s1 else None
    if kind == "del":
        if s1 and a == b and len(s1) > 1:
            j = i + len(s1) - 1
        if s1 and seq[i:j + 1] != s1:
            return None
        return seq[:i] + seq[j + 1:]
    if kind == "dup":
        if s1 and seq[i:j + 1] != s1:
            return None
        return seq[:j + 1] + seq[i:j + 1] + seq[j + 1:]
    if kind == "ins":
        return seq[:i + 1] + s1 + seq[i + 1:] if s1 else None
    if kind == "delins":
        return seq[:i] + s1 + seq[j + 1:] if s1 else None
    return None


def reconcile(seq, cds0, change, protein, k):
    mut = mutate(seq, cds0, change, k)
    if mut is None:
        return None
    wt, mt = translate(seq[cds0:]), translate(mut[cds0:])
    kind, pos, ref, alt = protein
    if not 0 < pos <= len(wt) or wt[pos - 1] != ref:
        return None
    n = next((x for x in range(min(len(wt), len(mt))) if wt[x] != mt[x]), None)
    if kind == "sub":
        if change[0] != "sub" or (change[1] + k - 1) // 3 + 1 != pos:
            return None
        if n is None:
            return "substitution" if alt == ref else None
        return "substitution" if n + 1 == pos and mt[n] == alt else None
    if n is None:
        return None
    if kind == "fs":
        return "frameshift_start" if n + 1 == pos else None
    span = (len(change[3]) if change[3] else change[2] - change[1] + 1) // 3
    return "inframe_indel" if n + 1 <= pos <= n + 1 + span else None


def shift_notation(c_notation, k):
    if not c_notation.startswith("c."):
        return None
    bad = []

    def repl(m):
        if m.group("sign"):
            bad.append(m.group(0))
            return m.group(0)
        new = int(m.group("num")) + k
        if new <= 0:
            bad.append(m.group(0))
        return f"{new}{m.group('intron') or ''}"

    body = re.split(r"(del|dup|ins|>)", c_notation[2:], maxsplit=1)
    shifted = RE_POS.sub(repl, body[0])
    return None if bad else "c." + shifted + "".join(body[1:])


def load_remap(path):
    table = {}
    if not path or not Path(path).exists():
        return table
    with open(path, encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            table[(r["gene"].strip().upper(), r["accession"].strip(), r["source_c"].strip())] = r["remapped_c"].strip()
    return table


def remap(table, gene, accession, c_notation):
    return table.get((gene.strip().upper(), accession.strip(), c_notation.strip()))
