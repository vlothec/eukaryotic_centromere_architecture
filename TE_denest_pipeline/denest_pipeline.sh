#!/usr/bin/env bash
# =============================================================================
# denest_pipeline.sh -- Iteratively peel nested intact transposable elements.
#
# Input coordinates are BED (0-based, half-open): [start, end).
# A "layer" is one excision event.  Layer 0 removes TEs already annotated in
# the supplied BED/GFF.  Later layers remove full-length BLAST matches exposed
# after the previous layer was excised.
# =============================================================================

set -euo pipefail

log() {
    printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2
}

usage() {
    cat <<'USAGE'
Usage:
  denest_pipeline.sh -g genome.fa (-r chr:start-end | -R rois.bed) \
    (-q intact_tes.fa | -b intact_tes.bed | -a edta.gff3) [options]

Required inputs:
  -g FILE  Genome assembly FASTA.
  -r TEXT  One ROI as chr:start-end (BED coordinates, end excluded).
  -R FILE  Multiple ROIs in BED format.  Columns 1--3 are required; column 4
           is an optional unique ROI name.
  -q FILE  Intact TE sequences (one FASTA record per TE), used as BLAST query.
  -b FILE  Intact TE coordinates in BED.  They are removed as layer 0 and
           extracted from the genome to make BLAST queries.
  -a FILE  EDTA GFF3. Intact LTR-RTs are identified with
           extract_intact_ltr_edta.sh and removed as layer 0.

Options:
  -o DIR   Output directory [./te_denest_out].
  -c NUM   Minimum fraction of a query aligned (0--100) [80].
  -p NUM   Starting nucleotide identity percentage [95].
  -s NUM   Stopping nucleotide identity percentage [80].
  -d NUM   Identity decrease after a threshold is exhausted [5].
  -m NUM   Maximum number of BLAST-derived layers per ROI [100].
  -t NUM   BLAST threads [8].
  -C/-P/-S/-D  Legacy aliases for -c/-p/-s/-d, respectively.
  -h       Show this help.

Dependencies: samtools, bedtools, BLAST+ (makeblastdb and blastn), Python 3.
USAGE
}

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

require_command() {
    command -v "$1" >/dev/null 2>&1 || die "Required program not found: $1"
}

fasta_length() {
    awk '!/^>/ { length_sum += length($0) } END { print length_sum + 0 }' "$1"
}

write_empty_bed() {
    : > "$1"
}

# Convert a single chr:start-end ROI to a four-column BED row.
make_single_roi() {
    local roi="$1"
    local output="$2"
    if [[ ! "$roi" =~ ^([^:]+):([0-9]+)-([0-9]+)$ ]]; then
        die "ROI must be chr:start-end with BED coordinates: $roi"
    fi
    local chromosome="${BASH_REMATCH[1]}"
    local start="${BASH_REMATCH[2]}"
    local end="${BASH_REMATCH[3]}"
    (( start < end )) || die "ROI start must be smaller than end: $roi"
    printf '%s\t%s\t%s\t%s\n' "$chromosome" "$start" "$end" \
        "${chromosome}_${start}_${end}" > "$output"
}

# Ensure every ROI has a safe, unique name for directory and FASTA identifiers.
normalise_rois() {
    awk 'BEGIN { OFS="\t" }
        /^#/ || NF < 3 { next }
        $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $2 >= $3 {
            print "Invalid ROI BED row " NR > "/dev/stderr"; exit 1
        }
        {
            name = (NF >= 4 && $4 != ".") ? $4 : $1 "_" $2 "_" $3
            gsub(/[^A-Za-z0-9_.-]/, "_", name)
            if (++seen[name] > 1) {
                print "ROI names must be unique: " name > "/dev/stderr"; exit 1
            }
            print $1, $2, $3, name
        }' "$1" > "$2"
}

make_layer() {
    local roi_dir="$1"
    local layer_number="$2"
    local source_fasta="$3"
    local local_bed="$4"
    local roi_name="$5"
    local identity="$6"
    local hit_count="$7"

    local layer_dir
    layer_dir=$(printf '%s/layer_%03d' "$roi_dir" "$layer_number")
    mkdir -p "$layer_dir"

    # Store raw hits separately from merged intervals: raw hits count candidate
    # TEs, whereas merged intervals define the bases actually excised.
    cp "$local_bed" "$layer_dir/candidates.bed"
    bedtools sort -i "$local_bed" | bedtools merge -i - > "$layer_dir/remove_merged.bed"
    local removed_bases
    removed_bases=$(awk '{ sum += $3 - $2 } END { print sum + 0 }' "$layer_dir/remove_merged.bed")
    (( removed_bases > 0 )) || return 1

    python3 "$HELPER" excise \
        --fasta "$source_fasta" \
        --intervals "$layer_dir/remove_merged.bed" \
        --output "$layer_dir/denested.fasta" \
        --name "${roi_name}_layer_${layer_number}"

    python3 "$HELPER" project \
        --mapping "$roi_dir/current_to_original.bed" \
        --removals "$layer_dir/remove_merged.bed" \
        --next-map "$layer_dir/current_to_original.bed" \
        --removed-original "$layer_dir/removed_original_roi.bed" \
        --name "$roi_name"
    cp "$layer_dir/current_to_original.bed" "$roi_dir/current_to_original.bed"

    local remaining_bases
    remaining_bases=$(fasta_length "$layer_dir/denested.fasta")
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$layer_number" "$identity" "$hit_count" "$removed_bases" "$remaining_bases" \
        "$(basename "$layer_dir")" >> "$roi_dir/layer_summary.tsv"
    cat "$layer_dir/removed_original_roi.bed" >> "$roi_dir/all_removed_original_roi.bed"
    printf '%s\n' "$layer_dir"
}

GENOME=""
SINGLE_ROI=""
ROI_BED=""
QUERY_FASTA=""
TE_BED=""
EDTA_GFF=""
OUTDIR="./te_denest_out"
MIN_COVERAGE=80
IDENTITY_START=95
IDENTITY_STOP=80
IDENTITY_STEP=5
MAX_LAYERS=100
THREADS=8

while getopts ':g:r:R:q:b:a:o:c:p:s:d:m:t:C:P:S:D:h' option; do
    case "$option" in
        g) GENOME="$OPTARG" ;;
        r) SINGLE_ROI="$OPTARG" ;;
        R) ROI_BED="$OPTARG" ;;
        q) QUERY_FASTA="$OPTARG" ;;
        b) TE_BED="$OPTARG" ;;
        a) EDTA_GFF="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        c) MIN_COVERAGE="$OPTARG" ;;
        C) MIN_COVERAGE="$OPTARG" ;;
        p) IDENTITY_START="$OPTARG" ;;
        P) IDENTITY_START="$OPTARG" ;;
        s) IDENTITY_STOP="$OPTARG" ;;
        S) IDENTITY_STOP="$OPTARG" ;;
        d) IDENTITY_STEP="$OPTARG" ;;
        D) IDENTITY_STEP="$OPTARG" ;;
        m) MAX_LAYERS="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        h) usage; exit 0 ;;
        :) die "Option -$OPTARG requires an argument" ;;
        *) usage; exit 1 ;;
    esac
done

[[ -n "$GENOME" ]] || die "Specify a genome with -g"
[[ -n "$SINGLE_ROI" || -n "$ROI_BED" ]] || die "Specify -r or -R"
[[ -z "$SINGLE_ROI" || -z "$ROI_BED" ]] || die "Use only one of -r and -R"
input_count=0
[[ -n "$QUERY_FASTA" ]] && ((input_count+=1))
[[ -n "$TE_BED" ]] && ((input_count+=1))
[[ -n "$EDTA_GFF" ]] && ((input_count+=1))
(( input_count == 1 )) || die "Specify exactly one of -q, -b, or -a"
[[ -f "$GENOME" ]] || die "Genome FASTA not found: $GENOME"
[[ -z "$ROI_BED" || -f "$ROI_BED" ]] || die "ROI BED not found: $ROI_BED"
[[ -z "$QUERY_FASTA" || -f "$QUERY_FASTA" ]] || die "TE FASTA not found: $QUERY_FASTA"
[[ -z "$TE_BED" || -f "$TE_BED" ]] || die "TE BED not found: $TE_BED"
[[ -z "$EDTA_GFF" || -f "$EDTA_GFF" ]] || die "EDTA GFF3 not found: $EDTA_GFF"
[[ "$MIN_COVERAGE" =~ ^[0-9]+$ ]] && (( MIN_COVERAGE > 0 && MIN_COVERAGE <= 100 )) || die "-c must be 1--100"
[[ "$IDENTITY_START" =~ ^[0-9]+$ && "$IDENTITY_STOP" =~ ^[0-9]+$ && "$IDENTITY_STEP" =~ ^[0-9]+$ ]] || die "Identity settings must be integers"
(( IDENTITY_START >= IDENTITY_STOP && IDENTITY_STOP > 0 && IDENTITY_STEP > 0 )) || die "Invalid identity settings"
[[ "$MAX_LAYERS" =~ ^[0-9]+$ ]] && (( MAX_LAYERS > 0 )) || die "-m must be positive"

for program in samtools bedtools makeblastdb blastn python3; do
    require_command "$program"
done

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
HELPER="$SCRIPT_DIR/denest_helpers.py"
[[ -f "$HELPER" ]] || die "Missing helper script: $HELPER"
mkdir -p "$OUTDIR"

ROI_FILE="$OUTDIR/rois.normalised.bed"
if [[ -n "$SINGLE_ROI" ]]; then
    make_single_roi "$SINGLE_ROI" "$ROI_FILE"
else
    normalise_rois "$ROI_BED" "$ROI_FILE"
fi
[[ -s "$ROI_FILE" ]] || die "No valid ROIs were provided"

samtools faidx "$GENOME"

# Build a single query FASTA and a global layer-0 TE BED when annotations are
# supplied. A FASTA-only run has no layer 0; it begins with BLAST discovery.
SEED_BED="$OUTDIR/intact_tes_layer0.bed"
if [[ -n "$QUERY_FASTA" ]]; then
    cp "$QUERY_FASTA" "$OUTDIR/intact_te_queries.fasta"
    write_empty_bed "$SEED_BED"
elif [[ -n "$TE_BED" ]]; then
    awk 'BEGIN { OFS="\t" } !/^#/ && NF >= 3 { print $1, $2, $3, (NF >= 4 ? $4 : "TE"), ".", (NF >= 6 ? $6 : "+") }' \
        "$TE_BED" | bedtools sort -i - > "$SEED_BED"
    bedtools getfasta -fi "$GENOME" -bed "$SEED_BED" -s -name > "$OUTDIR/intact_te_queries.fasta"
else
    bash "$SCRIPT_DIR/extract_intact_ltr_edta.sh" "$EDTA_GFF" ltr \
        | awk 'BEGIN { OFS="\t" } { print $1, $4 - 1, $5, $9, ".", $7 }' \
        | bedtools sort -i - > "$OUTDIR/intact_te_queries.bed"
    bedtools getfasta -fi "$GENOME" -bed "$OUTDIR/intact_te_queries.bed" -s -name > "$OUTDIR/intact_te_queries.fasta"
    bash "$SCRIPT_DIR/extract_intact_ltr_edta.sh" "$EDTA_GFF" cut \
        | bedtools sort -i - > "$SEED_BED"
fi
[[ -s "$OUTDIR/intact_te_queries.fasta" ]] || die "No intact TE query sequences were produced"

printf 'roi\tlayer\tidentity\tcandidate_TE_hits\tremoved_bases\tremaining_bases\tdirectory\n' > "$OUTDIR/all_roi_layer_summary.tsv"
log "Starting de-nesting for $(wc -l < "$ROI_FILE" | tr -d ' ') ROI(s)."

while IFS=$'\t' read -r chromosome roi_start roi_end roi_name; do
    roi_dir="$OUTDIR/$roi_name"
    mkdir -p "$roi_dir"
    log "ROI $roi_name: ${chromosome}:${roi_start}-${roi_end}"

    samtools faidx "$GENOME" "${chromosome}:$((roi_start + 1))-${roi_end}" \
        | awk -v name="$roi_name" '/^>/ { print ">" name; next } { print }' > "$roi_dir/original.fasta"
    printf '%s\t0\t%s\n' "$roi_name" "$((roi_end - roi_start))" > "$roi_dir/current_to_original.bed"
    write_empty_bed "$roi_dir/all_removed_original_roi.bed"
    printf 'layer\tidentity\tcandidate_TE_hits\tremoved_bases\tremaining_bases\tdirectory\n' > "$roi_dir/layer_summary.tsv"

    current_fasta="$roi_dir/original.fasta"
    layer=0

    # Layer 0 is clipped to this ROI and converted to ROI-local coordinates.
    awk -v chr="$chromosome" -v roi_start="$roi_start" -v roi_end="$roi_end" -v name="$roi_name" '
        $1 == chr && $2 < roi_end && $3 > roi_start {
            start = ($2 > roi_start) ? $2 - roi_start : 0
            end = ($3 < roi_end) ? $3 - roi_start : roi_end - roi_start
            print name "\t" start "\t" end
        }' "$SEED_BED" > "$roi_dir/layer_000_candidates.bed"
    seed_hits=$(wc -l < "$roi_dir/layer_000_candidates.bed" | tr -d ' ')
    if (( seed_hits > 0 )); then
        layer_dir=$(make_layer "$roi_dir" "$layer" "$current_fasta" "$roi_dir/layer_000_candidates.bed" "$roi_name" "annotation" "$seed_hits")
        current_fasta="$layer_dir/denested.fasta"
        layer=$((layer + 1))
    fi

    identity="$IDENTITY_START"
    blast_layers=0
    # If layer 0 removed the complete ROI, skip BLAST database creation on an
    # empty sequence and still write a valid final FASTA and summary.
    while (( identity >= IDENTITY_STOP && blast_layers < MAX_LAYERS )) && \
          (( $(fasta_length "$current_fasta") > 0 )); do
        work_dir="$roi_dir/blast_identity_${identity}"
        mkdir -p "$work_dir"
        makeblastdb -in "$current_fasta" -dbtype nucl -out "$work_dir/current" >/dev/null

        # qlen permits an explicit post-filter on one HSP's query coverage.
        blastn -query "$OUTDIR/intact_te_queries.fasta" -db "$work_dir/current" \
            -perc_identity "$identity" -num_threads "$THREADS" \
            -outfmt '6 qseqid qlen sseqid pident length qstart qend sstart send evalue bitscore' \
            > "$work_dir/hits.tsv"

        awk -v min_coverage="$MIN_COVERAGE" -v name="$roi_name" 'BEGIN { OFS="\t" }
            $5 / $2 * 100 >= min_coverage {
                start = ($8 < $9) ? $8 - 1 : $9 - 1
                end = ($8 < $9) ? $9 : $8
                print name, start, end, $1, $4, $5 / $2 * 100
            }' "$work_dir/hits.tsv" > "$work_dir/candidates.bed"

        hit_count=$(wc -l < "$work_dir/candidates.bed" | tr -d ' ')
        if (( hit_count == 0 )); then
            log "  ${identity}% identity exhausted."
            # Do not reset the final cutoff to itself; that would loop forever.
            if (( identity == IDENTITY_STOP )); then
                break
            fi
            identity=$((identity - IDENTITY_STEP))
            (( identity < IDENTITY_STOP )) && identity="$IDENTITY_STOP"
            continue
        fi

        cp "$work_dir/candidates.bed" "$roi_dir/layer_$(printf '%03d' "$layer")_candidates.bed"
        layer_dir=$(make_layer "$roi_dir" "$layer" "$current_fasta" "$work_dir/candidates.bed" "$roi_name" "$identity" "$hit_count") || break
        current_fasta="$layer_dir/denested.fasta"
        log "  Layer $layer: removed $hit_count qualifying TE hit(s) at ${identity}% identity."
        layer=$((layer + 1))
        blast_layers=$((blast_layers + 1))
        # Keep the same cutoff: excision can reveal another nested layer.
    done

    if (( blast_layers == MAX_LAYERS )); then
        log "  Reached maximum BLAST-derived layers ($MAX_LAYERS)."
    fi
    cp "$current_fasta" "$roi_dir/final_denested.fasta"
    awk -v roi="$roi_name" 'NR > 1 { print roi "\t" $0 }' "$roi_dir/layer_summary.tsv" >> "$OUTDIR/all_roi_layer_summary.tsv"
done < "$ROI_FILE"

log "Finished. Per-layer sequences and coordinates are in: $OUTDIR"
