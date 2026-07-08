import csv
import os
import re

THESIS_DIR = r"C:\Users\BornLoser\Desktop\Assignment\Thesis"
STEP8_OUT  = os.path.join(THESIS_DIR, "04_Mutation_Processing", "Output", "Step8_Merging")

IN_PATH  = os.path.join(STEP8_OUT, "unresolved_variants.tsv")
OUT_PATH = os.path.join(STEP8_OUT, "unresolved_disposition.tsv")
SUMMARY  = os.path.join(STEP8_OUT, "unresolved_disposition_summary.tsv")

WATSON = {"A": "T", "T": "A", "C": "G", "G": "C"}


def rc(s):
    return "".join(WATSON.get(c, "N") for c in reversed(s))


def strip_desc(s):
    s = (s or "").strip()
    return s.rsplit(":", 1)[1] if ":" in s else s


def classify(row):
    err = row.get("errors") or ""
    mut = row.get("mut_type", "")

    if "ENOSELECTORFOUND" in err:
        return {
            "category": "ENOSELECTORFOUND",
            "subcategory": "nm_not_in_ng",
            "root_cause": "Specified NM transcript not annotated in the chosen NG_ reference (NM version obsolescence or wrong NG/NM pairing)",
            "rescue_path": "resubmit against MANE Select NM directly or via VariantValidator; or use the NM version annotated in the NG",
            "disposition": "RESCUABLE_AUTO",
        }

    if "EINTRONIC" in err:
        return {
            "category": "EINTRONIC",
            "subcategory": "intronic_c_notation",
            "root_cause": "Intronic position not resolvable in c. context",
            "rescue_path": "VariantValidator API with MANE transcript",
            "disposition": "RESCUABLE_AUTO",
        }

    if "ESYNTAXUEOF" in err:
        return {
            "category": "ESYNTAXUEOF",
            "subcategory": "missing_inserted_sequence",
            "root_cause": "delins with no inserted sequence (data loss in original curation)",
            "rescue_path": "none",
            "disposition": "UNRESCUABLE",
        }

    if "EINSERTIONRANGE" in err:
        return {
            "category": "EINSERTIONRANGE",
            "subcategory": "non_consecutive_ins_range",
            "root_cause": "Insertion range positions not consecutive in transcript",
            "rescue_path": "manual HGVS rewrite (likely off-by-one) then resubmit",
            "disposition": "RESCUABLE_MANUAL",
        }

    if "ERANGEREVERSED" in err:
        return {
            "category": "ERANGEREVERSED",
            "subcategory": "start_greater_than_end",
            "root_cause": "Intronic offset range written in wrong order",
            "rescue_path": "manual HGVS rewrite to canonical form then resubmit",
            "disposition": "RESCUABLE_MANUAL",
        }

    if "ESEQUENCEMISMATCH" in err:
        m = re.search(
            r"`([ACGT]+)`\s*was not found in the reference sequence;\s*`([ACGT]+)`\s*was found",
            err,
        )
        if not m:
            return {
                "category": "ESEQUENCEMISMATCH",
                "subcategory": "unparseable",
                "root_cause": "Reference mismatch with unparseable error details",
                "rescue_path": "manual inspection",
                "disposition": "UNKNOWN",
            }
        submitted, actual = m.group(1), m.group(2)
        if len(submitted) == 1 and len(actual) == 1:
            if WATSON[submitted] == actual:
                return {
                    "category": "ESEQUENCEMISMATCH",
                    "subcategory": "single_base_complement",
                    "root_cause": f"Submitted ref {submitted} is complement of actual {actual}; curator may have recorded genomic + strand base on minus-strand gene",
                    "rescue_path": "none (would require fabricating correct base)",
                    "disposition": "UNRESCUABLE_DATA_ERROR",
                }
            return {
                "category": "ESEQUENCEMISMATCH",
                "subcategory": "single_base_wrong",
                "root_cause": f"Submitted ref {submitted} does not match GRCh38 ({actual}); curation error or obsolete reference",
                "rescue_path": "none",
                "disposition": "UNRESCUABLE_DATA_ERROR",
            }
        if submitted == rc(actual):
            return {
                "category": "ESEQUENCEMISMATCH",
                "subcategory": "multibase_rc_match",
                "root_cause": f"Submitted {submitted} is reverse complement of actual {actual}",
                "rescue_path": "manual review; likely recorded in wrong orientation",
                "disposition": "RESCUABLE_MANUAL",
            }
        if len(submitted) == len(actual) and (
            submitted.startswith(actual[:3]) or submitted.endswith(actual[-3:])
        ):
            return {
                "category": "ESEQUENCEMISMATCH",
                "subcategory": "positional_offset",
                "root_cause": f"Submitted {submitted} shares prefix/suffix with actual {actual}; probable off-by-one coordinate",
                "rescue_path": "manual: strip base from HGVS or shift position then resubmit",
                "disposition": "RESCUABLE_MANUAL",
            }
        if len(submitted) != len(actual):
            return {
                "category": "ESEQUENCEMISMATCH",
                "subcategory": "length_mismatch",
                "root_cause": f"Length of submitted ({len(submitted)}) differs from actual ({len(actual)})",
                "rescue_path": "none",
                "disposition": "UNRESCUABLE_DATA_ERROR",
            }
        return {
            "category": "ESEQUENCEMISMATCH",
            "subcategory": "completely_different",
            "root_cause": f"Submitted {submitted} completely different from actual {actual}",
            "rescue_path": "none",
            "disposition": "UNRESCUABLE_DATA_ERROR",
        }

    return {
        "category": "OTHER",
        "subcategory": "uncategorised",
        "root_cause": err[:200],
        "rescue_path": "manual review",
        "disposition": "UNKNOWN",
    }


def main():
    with open(IN_PATH, encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    groups = {}
    for r in rows:
        key = strip_desc(r.get("hgvs_input", "")) or r.get("sysname", "")
        gid = (r["gene"], key)
        g = groups.setdefault(gid, {"rep": r, "acc": set()})
        a = r.get("accession")
        if a:
            g["acc"].add(a)

    distinct = []
    for (gene, key), g in groups.items():
        rep = dict(g["rep"])
        rep["dedup_key"] = key
        rep["patient_count"] = len(g["acc"])
        rep["accession_list"] = ";".join(sorted(g["acc"]))
        distinct.append(rep)

    base_fields = list(rows[0].keys())
    out_fields = ["dedup_key", "patient_count", "accession_list"] + \
        [c for c in base_fields if c not in ("dedup_key", "patient_count", "accession_list")] + \
        ["category", "subcategory", "root_cause", "rescue_path", "disposition"]

    counts = {}
    disp_counts = {}
    with open(OUT_PATH, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=out_fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in distinct:
            c = classify(r)
            r.update(c)
            w.writerow(r)
            key = (c["category"], c["subcategory"])
            counts[key] = counts.get(key, 0) + 1
            disp_counts[c["disposition"]] = disp_counts.get(c["disposition"], 0) + 1

    total = len(distinct)
    with open(SUMMARY, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["category", "subcategory", "count", "pct_of_unresolved"])
        for (cat, sub), n in sorted(counts.items(), key=lambda x: -x[1]):
            w.writerow([cat, sub, n, f"{100*n/total:.1f}%"])
        w.writerow([])
        w.writerow(["disposition", "", "count", "pct_of_unresolved"])
        for disp, n in sorted(disp_counts.items(), key=lambda x: -x[1]):
            w.writerow([disp, "", n, f"{100*n/total:.1f}%"])

    print(f"Input rows           : {len(rows)}")
    print(f"Distinct unresolved  : {total}")
    print(f"Wrote                -> {OUT_PATH}")
    print(f"Summary              -> {SUMMARY}")
    print()
    print("Disposition summary (distinct variants):")
    for disp, n in sorted(disp_counts.items(), key=lambda x: -x[1]):
        print(f"  {disp:24s} {n:4d}  ({100*n/total:.1f}%)")


if __name__ == "__main__":
    main()
