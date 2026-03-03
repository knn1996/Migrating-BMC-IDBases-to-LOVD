# ============================================================
# Generate LiftOver BED files from idbase *pub.html files
# Using First hits.csv for chromosome coordinate mapping
# ============================================================

$thesisDir = "C:\Users\BornLoser\Desktop\Assignment\Thesis"
$idbaseDir = "$thesisDir\02_Source_Database\idbase"
$bedDir    = "$thesisDir\03_BED_Files\BED"
$csvPath   = "$thesisDir\04_Mutation_Processing\Output\First hits.csv"

# ---- Load BLAST coordinate map ----
# Key: GENENAME_DNA (uppercase), Value: list of {hg, chrom, sstart, send, qstart, qend}
$blastMap = @{}
$csvLines = Get-Content $csvPath | Select-Object -Skip 1
foreach ($line in $csvLines) {
    $f = $line -split ","
    if ($f.Count -lt 11) { continue }
    $file    = $f[0].Trim()   # e.g. ADA_vs_HG16.txt
    $qseqid  = $f[1].Trim()   # e.g. ADA_DNA
    $sseqid  = $f[2].Trim()   # e.g. chr20
    $qstart  = [int]$f[7].Trim()
    $qend    = [int]$f[8].Trim()
    $sstart  = [long]$f[9].Trim()
    $send    = [long]$f[10].Trim()

    # Determine HG version from filename
    $hg = "UNKNOWN"
    if ($file -match "HG16") { $hg = "hg16" }
    elseif ($file -match "HG17") { $hg = "hg17" }
    elseif ($file -match "HG18") { $hg = "hg18" }

    $key = $qseqid.ToUpper()
    if (-not $blastMap.ContainsKey($key)) { $blastMap[$key] = @() }
    $blastMap[$key] += [PSCustomObject]@{
        hg     = $hg
        chrom  = $sseqid
        sstart = $sstart
        send   = $send
        qstart = $qstart
        qend   = $qend
    }
}

Write-Host "Loaded $($blastMap.Count) gene entries from BLAST map"

# ---- Helper: convert IDRefSeq position to genomic coordinate ----
function Convert-ToGenomic {
    param($idrefPos, $blastEntry)
    # offset within the query alignment
    $offset = $idrefPos - $blastEntry.qstart  # 0-based offset
    if ($blastEntry.sstart -le $blastEntry.send) {
        # Forward strand
        return $blastEntry.sstart + $offset
    } else {
        # Reverse strand
        return $blastEntry.sstart - $offset
    }
}

# ---- Process each *base folder ----
$baseFolders = Get-ChildItem $idbaseDir -Directory | Where-Object { $_.Name -like "*base" -and $_.Name -ne "ImmunomeBase" }

$totalBED = 0
$skipped  = 0

foreach ($folder in $baseFolders) {
    $geneName = $folder.Name -replace "base$", ""
    
    # Find *pub.html
    $pubFile = Get-ChildItem $folder.FullName -Filter "*pub.html" | Select-Object -First 1
    if (-not $pubFile) {
        Write-Host "  SKIP $geneName - no pub.html found"
        $skipped++
        continue
    }

    # Read content, strip HTML tags for parsing
    $rawContent = [System.IO.File]::ReadAllText($pubFile.FullName)
    $cleanContent = $rawContent -replace '<[^>]+>', ''
    $lines = $cleanContent.Split("`n")

    # Get the DNA sequence ID referenced in the header
    # Look for "Sequence" line to find the D#### IDRefSeq accession
    $seqLine = ($lines | Where-Object { $_ -match "^Sequence\s" } | Select-Object -First 1)
    
    # Extract IDRefSeq D#### accession (e.g., D0001, D0002, NCF1_DNA)
    $dnaAccession = $null
    if ($seqLine -match "(?i)IDRefSeq:\s*([A-Za-z0-9_]+)") {
        $dnaAccession = $matches[1].ToUpper()
    }

    # Build lookup key for BLAST map
    # Most genes: GENENAME_DNA (e.g., ADA_DNA)
    # Special: NCF1_DNA uses "NCF1_DNA" as qseqid directly
    $lookupKeys = @()
    
    # Try standard: GENENAME_DNA
    $stdKey = ($geneName + "_DNA").ToUpper()
    $lookupKeys += $stdKey
    
    # Also try the DNA accession if it's NCF1_DNA style
    if ($dnaAccession -and $dnaAccession -match "_DNA$") {
        $lookupKeys += $dnaAccession
    }
    
    # Special cases
    if ($geneName -eq "CD40L") { $lookupKeys += "CD40L" }
    if ($geneName -eq "LRRC8A") { $lookupKeys += "LRRC8_DNA"; $lookupKeys += "LRRC8A_DNA" }
    if ($geneName -eq "IL12RB1") { $lookupKeys += "IL12RB_DNA" }
    if ($geneName -eq "FUCT1") { $lookupKeys += "FUCT1_DNA" }

    # Find matching BLAST entries
    $blastEntries = $null
    foreach ($key in $lookupKeys) {
        if ($blastMap.ContainsKey($key)) {
            $blastEntries = $blastMap[$key]
            break
        }
    }

    if (-not $blastEntries) {
        Write-Host "  SKIP $geneName - no BLAST coordinate mapping (uses EMBL/GenBank not in First hits.csv)"
        $skipped++
        continue
    }

    # ---- Parse mutations from pub.html ----
    # We look for blocks: Feature dna; N followed by /loc: IDRefSeq: D####: POSITION
    # We also capture the accession ID (for naming) and the change type
    
    $mutations = [System.Collections.Generic.List[PSCustomObject]]::new()
    $currentAccession = $null
    $currentSysName   = $null
    $inDnaFeature     = $false
    
    for ($i = 0; $i -lt $lines.Count; $i++) {
        $line = $lines[$i].Trim()
        
        # Capture mutation accession (e.g., "Accession       A0052")
        if ($line -match "^Accession\s+(\S+)") {
            $currentAccession = $matches[1]
            $currentSysName   = ""
            $inDnaFeature     = $false
        }
        
        # Capture systematic name
        if ($line -match "^Systematic name\s+(.+)") {
            $currentSysName = $matches[1].Trim()
        }
        
        # Detect start of DNA feature block
        if ($line -match "^Feature\s+dna;") {
            $inDnaFeature = $true
        }
        
        # If in dna feature, look for /loc: IDRefSeq: D####: POSITION
        if ($inDnaFeature -and $line -match "/loc:\s+(?:IDRefSeq|IdRefSeq):\s+[A-Za-z0-9_]+:\s*([\d]+)(?:\.\.([\d]+))?") {
            $posStart = [int]$matches[1]
            $posEnd   = if ($matches[2]) { [int]$matches[2] } else { $posStart }
            
            # Only add if position is within query alignment range of at least one blast entry
            $mutations.Add([PSCustomObject]@{
                accession = $currentAccession
                sysname   = $currentSysName
                posStart  = $posStart
                posEnd    = $posEnd
            })
            $inDnaFeature = $false  # Reset after capturing position
        }
        
        # Exit DNA feature if we hit a non-feature line (new section)
        if ($line -match "^Feature\s+rna;" -or $line -match "^Feature\s+aa;") {
            $inDnaFeature = $false
        }
    }

    if ($mutations.Count -eq 0) {
        Write-Host "  SKIP $geneName - no DNA mutations parsed from $($pubFile.Name)"
        $skipped++
        continue
    }

    # ---- Write BED files for each HG version ----
    $hgVersions = $blastEntries | ForEach-Object { $_.hg } | Sort-Object -Unique
    
    foreach ($hg in $hgVersions) {
        $hgEntry = $blastEntries | Where-Object { $_.hg -eq $hg } | Select-Object -First 1
        $bedFile = Join-Path $bedDir "$geneName`_$hg.BED"
        
        $bedLines = [System.Collections.Generic.List[string]]::new()
        $bedLines.Add("# BED file for $geneName ($hg) - LiftOver format")
        $bedLines.Add("# Source: $($pubFile.FullName)")
        $bedLines.Add("# Chrom: $($hgEntry.chrom), BLAST: qstart=$($hgEntry.qstart) qend=$($hgEntry.qend) sstart=$($hgEntry.sstart) send=$($hgEntry.send)")
        
        $strand = if ($hgEntry.sstart -le $hgEntry.send) { "+" } else { "-" }
        
        foreach ($mut in $mutations) {
            # Check position is within BLAST alignment range
            if ($mut.posStart -lt $hgEntry.qstart -or $mut.posEnd -gt $hgEntry.qend) {
                continue  # Position outside mapped region
            }
            
            $gStart = Convert-ToGenomic $mut.posStart $hgEntry
            $gEnd   = Convert-ToGenomic $mut.posEnd $hgEntry
            
            # BED is 0-based, half-open [start, end)
            # For reverse strand, gStart > gEnd after conversion; normalize
            $chromStart = [Math]::Min($gStart, $gEnd) - 1  # convert to 0-based
            $chromEnd   = [Math]::Max($gStart, $gEnd)       # end is exclusive, so +0 for SNV
            
            # For SNVs (single position), end = start + 1
            if ($mut.posStart -eq $mut.posEnd) {
                $chromEnd = $chromStart + 1
            }
            
            # Name: accession_systematicname (sanitize)
            $name = "$($mut.accession)"
            if ($mut.sysname) {
                $safeName = $mut.sysname -replace '[,\s;]+', '_' -replace '[^A-Za-z0-9_.>-]', ''
                $name = "$($mut.accession)_$safeName"
            }
            
            $bedLines.Add("$($hgEntry.chrom)`t$chromStart`t$chromEnd`t$name`t0`t$strand")
        }
        
        # Write BED file
        [System.IO.File]::WriteAllLines($bedFile, $bedLines)
        $lineCount = $bedLines.Count - 3  # subtract header lines
        Write-Host "  WROTE $geneName $hg`: $lineCount mutations -> $bedFile"
        $totalBED++
    }
}

Write-Host ""
Write-Host "=== DONE ==="
Write-Host "BED files written: $totalBED"
Write-Host "Genes skipped (no BLAST mapping): $skipped"