#!/bin/bash

# LDLR MPRA Full Dataset Analysis
# This script runs the complete fiberseq_MPRA pipeline with all LDLR variant files

set -e  # Exit on any error

echo "=========================================="
echo "LDLR MPRA Full Dataset Analysis"
echo "=========================================="

# Define paths
ANALYSIS_DIR="/gscratch/stergachislab/bmallo/large_home/plasmid_fiberseq/fiberseq_MPRA_test/fiberseq_MPRA_full_analysis"
WT_BED="/gscratch/stergachislab/bmallo/large_home/FiberHMM/FiberHMM_samples/LDLR_MPRA/SNV_sorted_reads/fp_WT_bed/LDLR_WT_fp.bed"
VARIANT_DIR="/gscratch/stergachislab/bmallo/large_home/FiberHMM/FiberHMM_samples/LDLR_MPRA/LDLR_pure_SNVs/fp_beds"
TSV_FILE="/gscratch/stergachislab/bmallo/large_home/plasmid_fiberseq/figure_scripts/figure_7/GRCh38_LDLR.tsv"

# Create analysis directories
mkdir -p "$ANALYSIS_DIR"/{control,null,results_standard,results_fdr}
cd "$ANALYSIS_DIR"

echo "Analysis directory: $ANALYSIS_DIR"
echo "WT BED file: $WT_BED"
echo "Variant directory: $VARIANT_DIR"
echo "TSV file: $TSV_FILE"
echo ""

# Check variant files
echo "Checking variant files..."
VARIANT_COUNT=$(find "$VARIANT_DIR" -name "*_fp.bed" | wc -l)
echo "Found $VARIANT_COUNT variant BED files"

if [[ $VARIANT_COUNT -eq 0 ]]; then
    echo "ERROR: No variant BED files found in $VARIANT_DIR"
    exit 1
fi

echo "Sample variant files:"
find "$VARIANT_DIR" -name "*_fp.bed" | head -5
echo ""

# Check WT file size
WT_SIZE=$(stat -c%s "$WT_BED" 2>/dev/null || stat -f%z "$WT_BED" 2>/dev/null || echo "unknown")
WT_LINES=$(wc -l < "$WT_BED" 2>/dev/null || echo "unknown")
echo "WT BED file: $WT_SIZE bytes, $WT_LINES lines"
echo ""

# Step 1: Generate control data with more iterations for production
echo "=========================================="
echo "Step 1: Generating control data (production)"
echo "=========================================="

# Check if control already exists
if [[ -f "control/control.pkl" ]]; then
    echo "Control data already exists. Skipping generation."
    echo "To regenerate, delete control/control.pkl"
else
    python -m fiberseq_MPRA generate-control \
        --input-bed "$WT_BED" \
        --output control/control.pkl \
        --ref-seq-length 4718 \
        --subsample-size 5000 \
        --num-iterations 10000 \
        --row-min 1 \
        --row-max 400 \
        --col-min 3000 \
        --col-max 3600 \
        --bin-size 10 \
        --processes 128 \
        --parallel
    
    echo "✓ Control data generation complete"
fi
echo ""

# Step 2: Generate null distributions with more iterations for production
echo "=========================================="
echo "Step 2: Generating null distributions (production)"
echo "=========================================="

# Check if null distributions already exist
if [[ -f "null/null_means.pkl" && -f "null/null_stds.pkl" ]]; then
    echo "Null distributions already exist. Skipping generation."
    echo "To regenerate, delete null/*.pkl files"
else
    python -m fiberseq_MPRA generate-null \
        --input-bed "$WT_BED" \
        --output-means null/null_means.pkl \
        --output-stds null/null_stds.pkl \
        --ref-seq-length 4718 \
        --subsample-size 5000 \
        --num-iterations 10000 \
        --row-min 1 \
        --row-max 400 \
        --col-min 3000 \
        --col-max 3600 \
        --bin-size 10 \
        --processes 128 \
        --parallel
    
    echo "✓ Null distribution generation complete"
fi
echo ""

# Check generated files
echo "Generated preprocessing files:"
ls -la control/ null/
echo ""

# Step 3: Run main analysis (standard)
echo "=========================================="
echo "Step 3: Running main analysis (standard p-values)"
echo "=========================================="

python -m fiberseq_MPRA main \
    --input-dir "$VARIANT_DIR" \
    --control-pkl control/control.pkl \
    --null-means-pkl null/null_means.pkl \
    --null-stds-pkl null/null_stds.pkl \
    --output-dir results_standard/ \
    --ref-seq-length 4718 \
    --row-min 1 \
    --row-max 400 \
    --col-min 3000 \
    --col-max 3600 \
    --bin-size 10 \
    --vmin 2 \
    --vmax 10 \
    --xticks 50 \
    --cmap magma \
    --plot-both \
    --processes 128 \
    --tsv-file "$TSV_FILE" \
    --min-coverage 2000

echo "✓ Standard analysis complete"
echo ""

# Step 4: Run main analysis with FDR correction
echo "=========================================="
echo "Step 4: Running analysis with FDR correction"
echo "=========================================="

python -m fiberseq_MPRA main \
    --input-dir "$VARIANT_DIR" \
    --control-pkl control/control.pkl \
    --null-means-pkl null/null_means.pkl \
    --null-stds-pkl null/null_stds.pkl \
    --output-dir results_fdr/ \
    --ref-seq-length 4718 \
    --row-min 1 \
    --row-max 400 \
    --col-min 3000 \
    --col-max 3600 \
    --bin-size 10 \
    --vmin 2 \
    --vmax 10 \
    --xticks 50 \
    --cmap magma \
    --plot-both \
    --processes 128 \
    --tsv-file "$TSV_FILE" \
    --min-coverage 2000 \
    --apply-fdr \
    --fdr-threshold 0.05

echo "✓ FDR-corrected analysis complete"
echo ""

# Generate summary report
echo "=========================================="
echo "Analysis Summary Report"
echo "=========================================="

echo "Input Data:"
echo "- WT BED file: $WT_SIZE bytes, $WT_LINES reads"
echo "- Variant files processed: $VARIANT_COUNT files"
echo "- TSV reference data: $(wc -l < "$TSV_FILE") entries"
echo ""

echo "Results Generated:"
STANDARD_PLOTS=$(find results_standard/ -name "*.png" 2>/dev/null | wc -l)
FDR_PLOTS=$(find results_fdr/ -name "*.png" 2>/dev/null | wc -l)
echo "- Standard p-value plots: $STANDARD_PLOTS"
echo "- FDR-corrected plots: $FDR_PLOTS"
echo ""

echo "Output Directories:"
echo "- Standard results: results_standard/"
echo "- FDR-corrected results: results_fdr/"
echo "- Control data: control/"
echo "- Null distributions: null/"
echo ""

echo "File sizes:"
du -sh control/ null/ results_standard/ results_fdr/ 2>/dev/null || echo "Could not get directory sizes"
echo ""

# Create analysis log
cat > analysis_log.txt << EOF
LDLR MPRA Full Dataset Analysis Log
Generated: $(date)

Input Files:
- WT BED: $WT_BED ($WT_SIZE bytes, $WT_LINES reads)
- Variant directory: $VARIANT_DIR ($VARIANT_COUNT files)
- TSV reference: $TSV_FILE

Parameters Used:
- Reference sequence length: 4718
- Subsample size: 5000
- Iterations: 1000
- Row range: 1-400 (footprint sizes)
- Column range: 3000-3600 (genomic positions)
- Bin size: 10
- Minimum coverage: 1000 reads
- Processes: 8

Results:
- Standard p-value plots: $STANDARD_PLOTS
- FDR-corrected plots: $FDR_PLOTS
- FDR threshold: 0.05

Directories:
- Analysis root: $ANALYSIS_DIR
- Standard results: results_standard/
- FDR results: results_fdr/
- Control data: control/
- Null distributions: null/
EOF

echo "Analysis log saved to: analysis_log.txt"
echo ""

echo "=========================================="
echo "Full dataset analysis completed successfully!"
echo "=========================================="
echo ""
echo "Next steps:"
echo "1. Review plots in results_standard/ and results_fdr/"
echo "2. Compare results with and without FDR correction"
echo "3. Examine specific variants of interest"
echo "4. Consider adjusting parameters if needed"
echo "5. Generate publication-quality figures"
