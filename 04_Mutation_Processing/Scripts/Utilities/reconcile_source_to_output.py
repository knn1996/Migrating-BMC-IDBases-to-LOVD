import pandas as pd, re, os

D = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Output"
S8 = os.path.join(D, "Step8_Merging")

src = pd.read_csv(os.path.join(D, "all_mutations.tsv"), sep="\t", dtype=str).fillna("")
mer = pd.read_csv(os.path.join(S8, "merged_variants.tsv"), sep="\t", dtype=str).fillna("")
unr = pd.read_csv(os.path.join(S8, "unresolved_disposition.tsv"), sep="\t", dtype=str).fillna("")

def gpos(s):
    return set(int(x) for x in re.findall(r"g\.(\d+)", s))

mer_gkeys = set()
mer_cov = set()
for _, r in mer.iterrows():
    mer_cov.add((r["gene"], r["accession"]))
    for p in gpos(r["sysname"]) | gpos(r["g_hgvs"]):
        mer_gkeys.add((r["gene"], r["accession"], p))

unr_cov = set(zip(unr["gene"], unr["accession"]))

buckets = {"resolved_exact": 0, "resolved_patient_present": 0,
           "unresolved": 0, "no_notation": 0, "unaccounted": 0}
unacc = []
for _, r in src.iterrows():
    ga = (r["gene"], r["accession"])
    if r["sysname"].strip() == "":
        buckets["no_notation"] += 1
        continue
    ps = gpos(r["sysname"])
    if any((r["gene"], r["accession"], p) in mer_gkeys for p in ps):
        buckets["resolved_exact"] += 1
    elif ga in mer_cov:
        buckets["resolved_patient_present"] += 1
    elif ga in unr_cov:
        buckets["unresolved"] += 1
    else:
        buckets["unaccounted"] += 1
        unacc.append(r)

print("source rows:", len(src))
for k, v in buckets.items():
    print(f"  {k:26} {v:5}  ({100*v/len(src):4.1f}%)")
print("  " + "-" * 38)
print(f"  {'sum':26} {sum(buckets.values()):5}")
print()
print("distinct source patient-genes:", src.groupby(['gene','accession']).ngroups)
print("  present in merged:   ", len(mer_cov & set(zip(src['gene'], src['accession']))))
print("  present in unresolved:", len(unr_cov & set(zip(src['gene'], src['accession']))))

if unacc:
    out = os.path.join(S8, "reconcile_unaccounted.tsv")
    pd.DataFrame(unacc).to_csv(out, sep="\t", index=False)
    print("\nwrote", out, "rows:", len(unacc))
