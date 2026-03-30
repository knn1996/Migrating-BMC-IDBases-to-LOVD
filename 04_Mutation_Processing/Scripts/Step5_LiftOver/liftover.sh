shopt -s nullglob

HG38_DIR="/home/khoi1996/Documents/BED/BED_IDRefseq/hg38"
UNMAP_DIR="/home/khoi1996/Documents/BED/BED_IDRefseq/unmapped"

for f in *_hg16.BED; do
    gene=$(basename $f _hg16.BED)
    ./liftOver "$f" hg16ToHg38.over.chain "$HG38_DIR/${gene}_hg38_16.BED" "$UNMAP_DIR/${gene}_unmap_16.txt"
done
for f in *_hg17.BED; do
    gene=$(basename $f _hg17.BED)
    ./liftOver "$f" hg17ToHg38.over.chain "$HG38_DIR/${gene}_hg38_17.BED" "$UNMAP_DIR/${gene}_unmap_17.txt"
done
for f in *_hg18.BED; do
    gene=$(basename $f _hg18.BED)
    ./liftOver "$f" hg18ToHg38.over.chain "$HG38_DIR/${gene}_hg38_18.BED" "$UNMAP_DIR/${gene}_unmap_18.txt"
done
