#!/usr/bin/env bash
shopt -s nullglob

IN_DIR="${IN_DIR}"
OUT_DIR="${OUT_DIR}"
CHAIN_DIR="${CHAIN_DIR}"
UNMAP_DIR="${OUT_DIR}/unmapped"

LIFTOVER="${CHAIN_DIR}/liftOver"

mkdir -p "$OUT_DIR" "$UNMAP_DIR"

for f in "$IN_DIR"/*_hg16.BED; do
    gene=$(basename "$f" _hg16.BED)
    "$LIFTOVER" "$f" "$CHAIN_DIR/hg16ToHg38.over.chain" "$OUT_DIR/${gene}_hg38_16.BED" "$UNMAP_DIR/${gene}_unmap_16.txt"
done
for f in "$IN_DIR"/*_hg17.BED; do
    gene=$(basename "$f" _hg17.BED)
    "$LIFTOVER" "$f" "$CHAIN_DIR/hg17ToHg38.over.chain" "$OUT_DIR/${gene}_hg38_17.BED" "$UNMAP_DIR/${gene}_unmap_17.txt"
done
for f in "$IN_DIR"/*_hg18.BED; do
    gene=$(basename "$f" _hg18.BED)
    "$LIFTOVER" "$f" "$CHAIN_DIR/hg18ToHg38.over.chain" "$OUT_DIR/${gene}_hg38_18.BED" "$UNMAP_DIR/${gene}_unmap_18.txt"
done
