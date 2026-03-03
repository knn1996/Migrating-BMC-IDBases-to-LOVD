"""
generate_bed.py
===============
Generates LiftOver-ready BED files from idbase *pub.html mutation databases.

For each *base folder in idbase/:
  - Parses DNA mutation positions from Feature dna /loc lines
  - Maps IDRefSeq positions to hg16/17/18 chromosomal coordinates via
    100%-identity BLAST hits in First hits.csv
  - Writes GENE_hgXX.BED files ready for UCSC LiftOver

Usage:  python generate_bed.py
"""

import os, re, csv
from collections import defaultdict

# ── Configuration ──────────────────────────────────────────────────────────────
THESIS_DIR = r"C:\Users\BornLoser\Desktop\Assignment\Thesis"
IDBASE_DIR = os.path.join(THESIS_DIR, "02_Source_Database", "idbase")
BED_DIR    = os.path.join(THESIS_DIR, "03_BED_Files", "BED")
CSV_PATH   = os.path.join(THESIS_DIR, "04_Mutation_Processing", "Output", "First hits.csv")

# Special gene-name aliases in First hits.csv
ALIASES = {
    "LRRC8A":  ["LRRC8_DNA", "LRRC8A_DNA"],
    "IL12RB1": ["IL12RB_DNA"],
}

RE_STRIP = re.compile(r"<[^>]+>")
RE_LOC   = re.compile(r"/loc:\s+IDRefSeq:\s+\w+:\s*(\d+)(?:\.\.(\d+))?", re.IGNORECASE)

# ── Helpers ────────────────────────────────────────────────────────────────────
def load_blast_map(path):
    """Return {QSEQID_UPPER: [{"hg","chrom","qstart","qend","sstart","send"}]}
    for rows with pident == 100."""
    blast_map = defaultdict(list)
    skipped = 0
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            try:
                if float(row["pident"]) != 100.0:
                    skipped += 1; continue
            except (ValueError, KeyError):
                continue
            fname = row["file"].upper()
            hg = "hg16" if "HG16" in fname else "hg17" if "HG17" in fname else "hg18"
            blast_map[row["qseqid"].strip().upper()].append({
                "hg": hg, "chrom": row["sseqid"].strip(),
                "qstart": int(row["qstart"]), "qend": int(row["qend"]),
                "sstart": int(row["sstart"]), "send":  int(row["send"]),
            })
    print(f"[BLAST] {len(blast_map)} genes at 100% identity ({skipped} partial skipped)")
    return blast_map


def blast_key(gene, blast_map):
    """Try standard naming conventions; return matching key or None."""
    for c in [f"{gene}_DNA", gene] + ALIASES.get(gene, []):
        if c.upper() in blast_map:
            return c.upper()
    return None


def idref_to_genomic(pos, entry):
    """Convert 1-based IDRefSeq position to chromosomal coordinate."""
    offset = pos - entry["qstart"]
    return entry["sstart"] + offset if entry["sstart"] <= entry["send"] \
           else entry["sstart"] - offset


def parse_mutations(pub_html_path):
    """Parse all Feature dna /loc positions from a pub.html file.
    Returns [{"accession","sysname","pos_start","pos_end"}].
    - HTML tags stripped before parsing.
    - All alleles per entry captured (in_dna not reset after /loc).
    - Fallback sysname from /name + coords when no Systematic name exists.
    - Duplicate (accession, pos_start, pos_end) rows deduplicated.
    """
    lines = RE_STRIP.sub("", open(pub_html_path, encoding="utf-8",
                                  errors="replace").read()).splitlines()
    mutations, seen = [], set()
    acc = sys_name = feat_name = ""
    in_dna = False

    for line in lines:
        s = line.strip()
        if s.startswith("Accession"):
            m = re.match(r"Accession\s+(\S+)", s)
            if m:
                acc, sys_name, feat_name, in_dna = m.group(1), "", "", False
        elif s.startswith("Systematic name"):
            m = re.match(r"Systematic name\s+(.+)", s)
            if m:
                sys_name = m.group(1).strip()
        elif re.match(r"Feature\s+dna;", s):
            in_dna, feat_name = True, ""
        elif re.match(r"Feature\s+(rna|aa);", s):
            in_dna = False
        elif in_dna:
            m = re.match(r"/name:\s+(\S+)", s)
            if m:
                feat_name = m.group(1)
            m = RE_LOC.search(s)
            if m:
                p1 = int(m.group(1))
                p2 = int(m.group(2)) if m.group(2) else p1
                if sys_name:
                    name = sys_name
                elif feat_name:
                    coord = f"{p1}..{p2}" if p2 != p1 else str(p1)
                    name = f"g.{coord}{feat_name}"
                else:
                    name = ""
                key = (acc, p1, p2)
                if key not in seen:
                    seen.add(key)
                    mutations.append({"accession": acc, "sysname": name,
                                      "pos_start": p1, "pos_end": p2})
                # in_dna intentionally NOT reset -- captures all alleles
    return mutations


def make_bed_name(accession, sysname):
    safe = re.sub(r"[,\s;]+", "_", sysname)
    safe = re.sub(r"[^A-Za-z0-9_.>\-]", "", safe)
    return f"{accession}_{safe}" if safe else accession


def write_bed(gene, hg, entry, mutations, pub_path):
    """Write BED file; return number of lines written."""
    strand = "+" if entry["sstart"] <= entry["send"] else "-"
    header = [
        f"# BED file for {gene} ({hg}) - LiftOver format",
        f"# Source: {pub_path}",
        f"# Chrom: {entry['chrom']}  BLAST: qstart={entry['qstart']} "
        f"qend={entry['qend']} sstart={entry['sstart']} send={entry['send']}",
    ]
    bed_lines = []
    for mut in mutations:
        # Linear extrapolation beyond qend is valid -- no OOB filter applied.
        g1 = idref_to_genomic(mut["pos_start"], entry)
        g2 = idref_to_genomic(mut["pos_end"],   entry)
        cs = min(g1, g2) - 1          # 0-based start
        ce = max(g1, g2)              # half-open end
        if mut["pos_start"] == mut["pos_end"]:
            ce = cs + 1               # SNV: exactly 1 bp
        bed_lines.append(
            f"{entry['chrom']}\t{cs}\t{ce}\t"
            f"{make_bed_name(mut['accession'], mut['sysname'])}\t0\t{strand}"
        )
    with open(os.path.join(BED_DIR, f"{gene}_{hg}.BED"), "w", encoding="utf-8") as f:
        f.write("\n".join(header + bed_lines) + "\n")
    return len(bed_lines)


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    os.makedirs(BED_DIR, exist_ok=True)
    blast_map = load_blast_map(CSV_PATH)

    folders = sorted(
        d for d in os.listdir(IDBASE_DIR)
        if os.path.isdir(os.path.join(IDBASE_DIR, d))
        and d.lower().endswith("base")
        and d.lower() != "immunomebase"
    )

    total_files = total_muts = 0
    skipped = []

    for folder in folders:
        gene = re.sub(r"base$", "", folder, flags=re.IGNORECASE)
        folder_path = os.path.join(IDBASE_DIR, folder)

        pub = next((f for f in os.listdir(folder_path)
                    if f.lower().endswith("pub.html")), None)
        if not pub:
            skipped.append((gene, "no pub.html")); continue

        key = blast_key(gene, blast_map)
        if not key:
            skipped.append((gene, "no 100% BLAST match")); continue

        mutations = parse_mutations(os.path.join(folder_path, pub))
        if not mutations:
            skipped.append((gene, "no DNA mutations")); continue

        for entry in blast_map[key]:
            n = write_bed(gene, entry["hg"], entry, mutations,
                          os.path.join(folder_path, pub))
            print(f"  WROTE {gene}_{entry['hg']}.BED  ({n} mutations)")
            total_files += 1
            total_muts  += n

    print(f"\nBED files: {total_files}  |  Total entries: {total_muts}"
          f"  |  Skipped genes: {len(skipped)}")
    for g, reason in skipped:
        print(f"  SKIP {g}: {reason}")


if __name__ == "__main__":
    main()