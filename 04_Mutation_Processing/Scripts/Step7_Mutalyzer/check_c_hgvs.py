"""
check_c_hgvs.py
===============
Cross-checks the c_hgvs column in mutalyzer_results.tsv against
the c. notations in all_mutations.tsv.

Logic:
  - For each row in mutalyzer_results.tsv with a non-empty c_hgvs,
    find the matching accession in all_mutations.tsv (same gene + accession).
  - Extract the c. notation from all_mutations.tsv (variant_type=coding rows).
  - Compare the c. position and change from both sources.
  - Classify the result as:
      match         — c. number and change agree exactly
      position_only — c. position agrees but change differs (transcript version mismatch candidate)
      mismatch      — c. position disagrees (coordinate error candidate)
      no_idbase_c   — no c. notation found in all_mutations.tsv for this accession
      intronic_skip — skipped because c_hgvs is intronic (not meaningful to compare)

Output: cross_check_c_hgvs.csv with columns:
  gene, accession, hgvs_input, c_hgvs_mutalyzer, c_notation_idbase,
  c_pos_mutalyzer, c_pos_idbase, result, note

Saved to: Output/Step7_Mutalyzer/cross_check_c_hgvs.csv
"""

import os
import re
import csv

MUTALYZER_TSV = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Output\Step7_Mutalyzer\mutalyzer_results.tsv"
ALL_MUT_TSV   = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Output\Step1_Extraction\all_mutations.tsv"
OUT_CSV       = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Output\Step7_Mutalyzer\cross_check_c_hgvs.csv"

RE_INTRONIC = re.compile(r"c\.\*?\-?\d+[+\-]\d+")
RE_C_FULL   = re.compile(r"c\.(\*?\-?\d+(?:[+\-]\d+)?)(.*)")


def extract_c_core(c_hgvs):
    """
    Extract just the c. notation part from a full HGVS string.
    NC_000021.9(NM_000383.4):c.22C>T  ->  c.22C>T
    NM_000383.4:c.22C>T               ->  c.22C>T
    c.22C>T                           ->  c.22C>T
    """
    if not c_hgvs:
        return ""
    m = re.search(r"(c\.[^\s,;]+)", c_hgvs)
    return m.group(1) if m else ""


def parse_c_position(c_notation):
    """
    Extract the numeric position from a c. notation for comparison.
    c.22C>T       -> '22'
    c.55+552G>A   -> '55+552'
    c.1119del     -> '1119'
    c.*3A>G       -> '*3'
    Returns empty string if unparseable.
    """
    m = RE_C_FULL.match(c_notation)
    return m.group(1) if m else ""


def main():
    print("Loading all_mutations.tsv...")
    idbase_c = {}
    with open(ALL_MUT_TSV, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("variant_type", "").strip() != "coding":
                continue
            notation = row.get("notation", "").strip()
            if not notation.startswith("c."):
                continue
            key = (row["gene"].strip(), row["accession"].strip())
            if key not in idbase_c:
                idbase_c[key] = []
            idbase_c[key].append(notation)

    print(f"  {len(idbase_c)} gene/accession pairs with c. notations")

    print("Loading mutalyzer_results.tsv...")
    with open(MUTALYZER_TSV, encoding="utf-8") as f:
        mut_rows = list(csv.DictReader(f, delimiter="\t"))
    print(f"  {len(mut_rows)} rows")

    results = []
    n_match = n_pos_only = n_mismatch = n_no_idbase = n_intronic = n_no_c = 0

    for row in mut_rows:
        if row.get("status") != "ok":
            continue

        c_hgvs = row.get("c_hgvs", "").strip()
        if not c_hgvs:
            n_no_c += 1
            continue

        if RE_INTRONIC.search(c_hgvs):
            n_intronic += 1
            continue

        gene      = row["gene"].strip()
        accession = row["accession"].strip()
        key       = (gene, accession)

        c_core_mut = extract_c_core(c_hgvs)
        pos_mut    = parse_c_position(c_core_mut)

        idbase_notations = idbase_c.get(key, [])
        if not idbase_notations:
            results.append({
                "gene":               gene,
                "accession":          accession,
                "hgvs_input":         row.get("hgvs_input",""),
                "c_hgvs_mutalyzer":   c_hgvs,
                "c_notation_idbase":  "",
                "c_pos_mutalyzer":    pos_mut,
                "c_pos_idbase":       "",
                "result":             "no_idbase_c",
                "note":               "No coding c. notation found in all_mutations.tsv for this accession",
            })
            n_no_idbase += 1
            continue

        best_result   = None
        best_priority = 99

        for idb_c in idbase_notations:
            pos_idb = parse_c_position(idb_c)

            if c_core_mut.lower() == idb_c.lower():
                r = "match"
                priority = 0
                note = "Exact match — coordinate chain verified"
            elif pos_mut and pos_idb and pos_mut == pos_idb:
                r = "position_only"
                priority = 1
                note = "Position matches but change differs — possible transcript version mismatch"
            else:
                r = "mismatch"
                priority = 2
                note = f"Position mismatch — Mutalyzer: {pos_mut}, IDbases: {pos_idb} — possible coordinate error"

            if priority < best_priority:
                best_priority = priority
                best_result = {
                    "gene":               gene,
                    "accession":          accession,
                    "hgvs_input":         row.get("hgvs_input",""),
                    "c_hgvs_mutalyzer":   c_hgvs,
                    "c_notation_idbase":  idb_c,
                    "c_pos_mutalyzer":    pos_mut,
                    "c_pos_idbase":       pos_idb,
                    "result":             r,
                    "note":               note,
                }

        results.append(best_result)
        if best_result["result"] == "match":
            n_match += 1
        elif best_result["result"] == "position_only":
            n_pos_only += 1
        else:
            n_mismatch += 1

    fieldnames = ["gene","accession","hgvs_input","c_hgvs_mutalyzer","c_notation_idbase",
                  "c_pos_mutalyzer","c_pos_idbase","result","note"]

    with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    total = n_match + n_pos_only + n_mismatch + n_no_idbase
    print(f"\nResults written to {OUT_CSV}")
    print(f"  match          : {n_match}  ({100*n_match/total:.1f}%)" if total else "")
    print(f"  position_only  : {n_pos_only}  ({100*n_pos_only/total:.1f}%)" if total else "")
    print(f"  mismatch       : {n_mismatch}  ({100*n_mismatch/total:.1f}%)" if total else "")
    print(f"  no_idbase_c    : {n_no_idbase}")
    print(f"  intronic_skip  : {n_intronic}")
    print(f"  no_c_hgvs      : {n_no_c}")
    print(f"  Total compared : {total}")


if __name__ == "__main__":
    main()
