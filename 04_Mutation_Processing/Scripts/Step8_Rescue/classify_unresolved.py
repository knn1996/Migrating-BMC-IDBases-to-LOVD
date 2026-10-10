import csv
import json
import os
import re

IN_PATH  = os.environ["IN_PATH"]
OUT_PATH = os.environ["OUT_PATH"]
SUMMARY  = os.environ["SUMMARY"]

WATSON = {"A": "T", "T": "A", "C": "G", "G": "C"}
RE_IVS = re.compile(r"\bIVS\d+", re.IGNORECASE)


def _err_codes(err_str):
    if not err_str:
        return set()
    try:
        items = json.loads(err_str)
        if isinstance(items, list):
            return {d.get("code", "") for d in items if isinstance(d, dict)}
    except (json.JSONDecodeError, ValueError):
        pass
    return set(re.findall(r"\b(E[A-Z]+)\b", err_str))


def _has_ivs(row):
    return bool(
        RE_IVS.search(row.get("sysname") or "")
        or RE_IVS.search(row.get("hgvs_input") or "")
    )


def rc(s):
    return "".join(WATSON.get(c, "N") for c in reversed(s))


def strip_desc(s):
    s = (s or "").strip()
    return s.rsplit(":", 1)[1] if ":" in s else s


def classify(row):
    err    = row.get("errors") or ""
    codes  = _err_codes(err)
    status = row.get("status", "").strip()

    if status == "no_result":
        return {
            "category": "NO_REF",
            "subcategory": "no_usable_reference",
            "root_cause": "No NG_ or NM_ track produced any result for this variant; "
                          "gene likely absent from all reference FASTA directories",
            "rescue_path": "none",
            "disposition": "UNRESCUABLE",
        }

    if "ENOSELECTORFOUND" in err or "ENOSELECTORFOUND" in codes:
        return {
            "category": "ENOSELECTORFOUND",
            "subcategory": "nm_not_in_ng",
            "root_cause": "Specified NM transcript not annotated in the chosen NG_ reference "
                          "(NM version obsolescence or wrong NG/NM pairing)",
            "rescue_path": "resubmit against MANE Select NM directly or via VariantValidator; "
                           "or use the NM version annotated in the NG",
            "disposition": "RESCUABLE_AUTO",
        }

    if "EINTRONIC" in err or "EINTRONIC" in codes:
        if _has_ivs(row):
            return {
                "category": "EINTRONIC",
                "subcategory": "IVS-syntax",
                "root_cause": "Legacy IVS notation (e.g. IVS3+1G>T) is not valid HGVS; "
                              "Mutalyzer cannot resolve intronic position from NG_ context",
                "rescue_path": "VariantValidator API with bare NM_:c. after IVS normalisation",
                "disposition": "RESCUABLE_AUTO",
            }
        return {
            "category": "EINTRONIC",
            "subcategory": "intronic_c_notation",
            "root_cause": "Intronic position not resolvable in c. context via Mutalyzer NG_ route",
            "rescue_path": "VariantValidator API with MANE transcript",
            "disposition": "RESCUABLE_AUTO",
        }

    if "ESYNTAXUEOF" in err or "ESYNTAXUEOF" in codes:
        return {
            "category": "ESYNTAXUEOF",
            "subcategory": "missing_inserted_sequence",
            "root_cause": "delins with no inserted sequence (data loss in original curation)",
            "rescue_path": "none",
            "disposition": "UNRESCUABLE",
        }

    if "EINSERTIONRANGE" in err or "EINSERTIONRANGE" in codes:
        return {
            "category": "EINSERTIONRANGE",
            "subcategory": "non_consecutive_ins_range",
            "root_cause": "Insertion range positions not consecutive in transcript",
            "rescue_path": "manual HGVS rewrite (likely off-by-one) then resubmit",
            "disposition": "RESCUABLE_MANUAL",
        }

    if "ERANGEREVERSED" in err or "ERANGEREVERSED" in codes:
        return {
            "category": "ERANGEREVERSED",
            "subcategory": "start_greater_than_end",
            "root_cause": "Intronic offset range written in wrong order",
            "rescue_path": "manual HGVS rewrite to canonical form then resubmit",
            "disposition": "RESCUABLE_MANUAL",
        }

    if "ESEQUENCEMISMATCH" in err or "ESEQUENCEMISMATCH" in codes:
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
            if WATSON.get(submitted) == actual:
                return {
                    "category": "ESEQUENCEMISMATCH",
                    "subcategory": "single_base_complement",
                    "root_cause": f"Submitted ref {submitted} is complement of actual {actual}; "
                                  "curator may have recorded genomic + strand base on minus-strand gene",
                    "rescue_path": "none (would require fabricating correct base)",
                    "disposition": "UNRESCUABLE_DATA_ERROR",
                }
            return {
                "category": "ESEQUENCEMISMATCH",
                "subcategory": "single_base_wrong",
                "root_cause": f"Submitted ref {submitted} does not match GRCh38 ({actual}); "
                              "curation error or obsolete reference",
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
                "root_cause": f"Submitted {submitted} shares prefix/suffix with actual {actual}; "
                              "probable off-by-one coordinate from legacy reference",
                "rescue_path": "apply empirical sstart offset from lrg_offset_results.csv then resubmit",
                "disposition": "RESCUABLE_AUTO",
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
