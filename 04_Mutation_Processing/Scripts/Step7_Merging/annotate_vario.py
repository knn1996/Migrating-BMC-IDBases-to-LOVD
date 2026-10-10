import os
import re
import xml.etree.ElementTree as ET
import pandas as pd

IN_PATH   = os.environ["IN_PATH"]
OUT_PATH  = os.environ["OUT_PATH"]
VARIO_OWX = os.environ["VARIO_OWX"]

DNA_COL_IN   = "VariantOnTranscript/DNA"
DNA_COL_ALT  = "VariantOnGenome/DNA"
PROT_COL_IN  = "VariantOnTranscript/Protein"
DNA_COL_OUT  = "VariO/DNA"
PROT_COL_OUT = "VariO/Protein"

SUB      = "VariO_0136"
TS_PUR   = "VariO_0315"
TS_PYR   = "VariO_0314"
TV       = "VariO_0316"
DEL      = "VariO_0141"
INS      = "VariO_0142"
INDEL    = "VariO_0143"
INV      = "VariO_0145"
P_MIS    = "VariO_0308"
P_NON    = "VariO_0310"
P_SYN    = "VariO_0343"
P_TERM   = "VariO_0309"
P_AADEL  = "VariO_0014"
P_AAINS  = "VariO_0018"
P_AAINDEL= "VariO_0022"

FALLBACK = {
    "VariO_0136": "DNA substitution",
    "VariO_0315": "purine transition",
    "VariO_0314": "pyrimidine transition",
    "VariO_0316": "transversion",
    "VariO_0141": "DNA deletion",
    "VariO_0142": "DNA insertion",
    "VariO_0143": "DNA indel",
    "VariO_0145": "DNA inversion",
    "VariO_0308": "missense variation",
    "VariO_0310": "nonsense variation",
    "VariO_0343": "synonymous variation",
    "VariO_0309": "termination codon change",
    "VariO_0014": "amino acid deletion",
    "VariO_0018": "amino acid insertion",
    "VariO_0022": "amino acid indel",
}

TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
PURINE_TS   = {("A", "G"), ("G", "A")}
SUB_RE      = re.compile(r"([ACGT]+)>([ACGT]+)")
AA_RE       = re.compile(r"^([A-Za-z]{3})(\d+)([A-Za-z]{3})$")


def load_labels(path):
    labels = {}
    try:
        root = ET.parse(path).getroot()
    except Exception:
        return labels
    for aa in root.iter():
        if aa.tag.split("}")[-1] != "AnnotationAssertion":
            continue
        prop = subj = lit = None
        for e in aa:
            t = e.tag.split("}")[-1]
            if t == "AnnotationProperty":
                prop = e.get("IRI") or e.get("abbreviatedIRI") or ""
            elif t == "IRI":
                subj = e.text
            elif t == "Literal":
                lit = e.text
        if prop and lit and subj and prop.endswith("label"):
            m = re.search(r"(VariO_\d+)", subj or "")
            if m:
                labels[m.group(1)] = lit
    return labels


LABELS = load_labels(VARIO_OWX)


def term(tid):
    name = LABELS.get(tid) or FALLBACK.get(tid, tid)
    return f"{name} ({tid.replace('_', ':')})"


def join(ids):
    return ";".join(term(i) for i in ids if i)


def dna_terms(dna, genomic):
    s = (dna or "").strip() or (genomic or "").strip()
    if not s:
        return ""
    low = s.lower()
    if "delins" in low:
        return join([INDEL])
    if "inv" in low:
        return join([INV])
    if ">" in s:
        ids = [SUB]
        m = SUB_RE.search(s.upper())
        if m:
            r, a = m.group(1), m.group(2)
            if len(r) == 1 and len(a) == 1:
                if (r, a) in TRANSITIONS:
                    ids.append(TS_PUR if (r, a) in PURINE_TS else TS_PYR)
                else:
                    ids.append(TV)
        return join(ids)
    if "dup" in low or "ins" in low:
        return join([INS])
    if "del" in low:
        return join([DEL])
    return ""


def protein_terms(p):
    body = (p or "").strip()
    if not body:
        return ""
    if body.startswith("p."):
        body = body[2:]
    body = body.strip().lstrip("(").rstrip(")").strip()
    low = body.lower()
    if body == "=":
        return join([P_SYN])
    if not body or body == "?":
        return ""
    if "fs" in low:
        return ""
    if body.endswith("="):
        return join([P_SYN])
    if "ext" in low:
        return join([P_TERM])
    if "met1" in low or body.endswith("?"):
        return ""
    if "delins" in low or ("del" in low and "ins" in low):
        return join([P_AAINDEL])
    if "del" in low:
        return join([P_AADEL])
    if "dup" in low or "ins" in low:
        return join([P_AAINS])
    if body.endswith("Ter") or body.endswith("*"):
        return join([P_NON])
    m = AA_RE.match(body)
    if m:
        if m.group(1).lower() == m.group(3).lower():
            return join([P_SYN])
        return join([P_MIS])
    return ""


def main():
    df = pd.read_csv(IN_PATH, sep="\t", dtype=str).fillna("")
    dna = df[DNA_COL_IN] if DNA_COL_IN in df.columns else pd.Series([""] * len(df))
    gen = df[DNA_COL_ALT] if DNA_COL_ALT in df.columns else pd.Series([""] * len(df))
    prot = df[PROT_COL_IN] if PROT_COL_IN in df.columns else pd.Series([""] * len(df))
    df[DNA_COL_OUT] = [dna_terms(d, g) for d, g in zip(dna, gen)]
    df[PROT_COL_OUT] = [protein_terms(x) for x in prot]
    df.to_csv(OUT_PATH, sep="\t", index=False)
    n = len(df)
    print(f"rows                 : {n}")
    print(f"VariO/DNA filled     : {(df[DNA_COL_OUT] != '').sum()}")
    print(f"VariO/Protein filled : {(df[PROT_COL_OUT] != '').sum()}")
    print(f"ontology labels used : {len(LABELS)}")


if __name__ == "__main__":
    main()
