#!/bin/bash
# crosscheck.sh
# Compare hg38 coordinates lifted from hg16, hg17, hg18
# Reports genes where coordinates DIFFER between versions

MISMATCHES=0

for f in *_hg38_18.BED; do
    gene=$(basename $f _hg38_18.BED)

    f16="${gene}_hg38_16.BED"
    f17="${gene}_hg38_17.BED"
    f18="${gene}_hg38_18.BED"

    # Skip if any file is missing
    if [[ ! -f "$f16" || ! -f "$f17" || ! -f "$f18" ]]; then
        echo "MISSING  $gene — one or more hg38 files not found"
        continue
    fi

    # Strip comment lines, sort, compare
    clean16=$(grep -v "^#" "$f16" | sort)
    clean17=$(grep -v "^#" "$f17" | sort)
    clean18=$(grep -v "^#" "$f18" | sort)

    if [[ "$clean16" == "$clean17" && "$clean17" == "$clean18" ]]; then
        echo "OK       $gene — all three assemblies agree"
    else
        echo "MISMATCH $gene — coordinates differ between assemblies!"
        MISMATCHES=$((MISMATCHES + 1))

        # Show the diff between hg16 and hg18 as example
        echo "  --- hg16 vs hg18 ---"
        diff <(echo "$clean16") <(echo "$clean18") | head -10
        echo ""
    fi
done

echo ""
echo "=== DONE: $MISMATCHES gene(s) with mismatches ==="
