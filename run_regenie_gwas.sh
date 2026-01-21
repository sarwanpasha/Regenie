#!/bin/bash

###############################################################################
# REGENIE GWAS Pipeline
# A portable pipeline for genome-wide association studies using REGENIE
###############################################################################

set -e  # Exit on error

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored messages
print_info() { echo -e "${BLUE}[INFO]${NC} $1"; }
print_success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
print_warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }
print_error() { echo -e "${RED}[ERROR]${NC} $1"; }

# Function to display usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Required Arguments:
  -b, --bfile PREFIX       Path to PLINK binary file prefix (without .bed/.bim/.fam)
  -p, --pheno FILE        Path to phenotype file
  -o, --outdir DIR        Output directory for results

Optional Arguments:
  -c, --covar COLS        Covariate column names (comma-separated, e.g., "COV1,COV2")
                          Default: COV1,COV2
  -y, --pheno-col NAME    Phenotype column name (default: PHENO)
  --maf FLOAT             Minor allele frequency threshold (default: 0.01)
  --geno FLOAT            Genotype missingness threshold (default: 0.05)
  --mind FLOAT            Sample missingness threshold (default: 0.05)
  --hwe FLOAT             Hardy-Weinberg p-value threshold (default: 1e-6)
  --bsize INT             Block size for REGENIE (default: 1000 for step1, 400 for step2)
  --pthresh FLOAT         P-value threshold for Firth regression (default: 0.01)
  --threads INT           Number of threads (default: 4)
  --no-plots              Skip plot generation
  -h, --help              Display this help message

Example:
  $0 -b /path/to/plink_data -p /path/to/pheno.txt -o ./results -c "AGE,SEX,PC1,PC2"

Phenotype File Format:
  Tab or space-delimited file with header containing:
  - FID: Family ID
  - IID: Individual ID
  - PHENO: Phenotype (or custom name specified with --pheno-col)
  - Covariates: One or more covariate columns

EOF
    exit 1
}

# Default parameters
PLINK_PREFIX=""
PHENO_FILE=""
OUT_DIR=""
COVAR_COLS="COV1,COV2"
PHENO_COL="PHENO"
MAF=0.01
GENO=0.05
MIND=0.05
HWE=1e-6
BSIZE_STEP1=1000
BSIZE_STEP2=400
PTHRESH=0.01
THREADS=4
SKIP_PLOTS=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -b|--bfile)
            PLINK_PREFIX="$2"
            shift 2
            ;;
        -p|--pheno)
            PHENO_FILE="$2"
            shift 2
            ;;
        -o|--outdir)
            OUT_DIR="$2"
            shift 2
            ;;
        -c|--covar)
            COVAR_COLS="$2"
            shift 2
            ;;
        -y|--pheno-col)
            PHENO_COL="$2"
            shift 2
            ;;
        --maf)
            MAF="$2"
            shift 2
            ;;
        --geno)
            GENO="$2"
            shift 2
            ;;
        --mind)
            MIND="$2"
            shift 2
            ;;
        --hwe)
            HWE="$2"
            shift 2
            ;;
        --bsize)
            BSIZE_STEP1="$2"
            BSIZE_STEP2="$2"
            shift 2
            ;;
        --pthresh)
            PTHRESH="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --no-plots)
            SKIP_PLOTS=true
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            print_error "Unknown option: $1"
            usage
            ;;
    esac
done

# Validate required arguments
if [[ -z "$PLINK_PREFIX" ]] || [[ -z "$PHENO_FILE" ]] || [[ -z "$OUT_DIR" ]]; then
    print_error "Missing required arguments!"
    usage
fi

# Check if input files exist
if [[ ! -f "${PLINK_PREFIX}.bed" ]] || [[ ! -f "${PLINK_PREFIX}.bim" ]] || [[ ! -f "${PLINK_PREFIX}.fam" ]]; then
    print_error "PLINK binary files not found at: ${PLINK_PREFIX}"
    exit 1
fi

if [[ ! -f "$PHENO_FILE" ]]; then
    print_error "Phenotype file not found: ${PHENO_FILE}"
    exit 1
fi

# Create output directory
mkdir -p "$OUT_DIR"

# Print configuration
print_info "======================================"
print_info "REGENIE GWAS Pipeline Configuration"
print_info "======================================"
print_info "PLINK prefix: $PLINK_PREFIX"
print_info "Phenotype file: $PHENO_FILE"
print_info "Output directory: $OUT_DIR"
print_info "Phenotype column: $PHENO_COL"
print_info "Covariate columns: $COVAR_COLS"
print_info "MAF threshold: $MAF"
print_info "Genotype missingness: $GENO"
print_info "Sample missingness: $MIND"
print_info "HWE threshold: $HWE"
print_info "Threads: $THREADS"
print_info "======================================"

# Load modules (if needed - modify based on your HPC environment)
if command -v module &> /dev/null; then
    print_info "Loading required modules..."
    module load Plink/1.9.10 2>/dev/null || print_warning "Plink module not found, using system plink"
    module load R/4.5.2 2>/dev/null || print_warning "R module not found, using system R"
fi

# Step 1: Quality Control
print_info "Step 1: Running quality control..."
plink \
  --bfile "$PLINK_PREFIX" \
  --maf "$MAF" \
  --geno "$GENO" \
  --hwe "$HWE" \
  --mind "$MIND" \
  --make-bed \
  --threads "$THREADS" \
  --out "${OUT_DIR}/qc_data"

if [[ $? -eq 0 ]]; then
    print_success "Quality control completed"
else
    print_error "Quality control failed"
    exit 1
fi

# Step 2: Create high-quality SNP list for REGENIE Step 1
print_info "Step 2: Creating high-quality SNP list..."
plink \
  --bfile "${OUT_DIR}/qc_data" \
  --maf 0.05 \
  --geno 0.01 \
  --hwe 1e-15 \
  --indep-pairwise 1000 100 0.9 \
  --threads "$THREADS" \
  --out "${OUT_DIR}/qc_pass"

if [[ $? -eq 0 ]]; then
    print_success "SNP list created"
else
    print_error "SNP pruning failed"
    exit 1
fi

# Step 3: Prepare phenotype and covariate files
print_info "Step 3: Preparing phenotype and covariate files..."

R --vanilla <<EOF
# Read phenotype file
pheno <- read.table("${PHENO_FILE}", header=TRUE, stringsAsFactors=FALSE)

# Check if required columns exist
required_cols <- c("FID", "IID", "${PHENO_COL}")
if (!all(required_cols %in% colnames(pheno))) {
    missing <- required_cols[!required_cols %in% colnames(pheno)]
    stop("Missing required columns in phenotype file: ", paste(missing, collapse=", "))
}

# Prepare phenotype file for REGENIE
pheno_out <- pheno[, c("FID", "IID", "${PHENO_COL}")]
write.table(pheno_out, 
            "${OUT_DIR}/phenotype_regenie.txt", 
            row.names=FALSE, 
            col.names=TRUE, 
            quote=FALSE, 
            sep="\t")

# Prepare covariate file
covar_cols <- unlist(strsplit("${COVAR_COLS}", ","))
covar_cols <- trimws(covar_cols)  # Remove whitespace

# Check if covariate columns exist
missing_covars <- covar_cols[!covar_cols %in% colnames(pheno)]
if (length(missing_covars) > 0) {
    warning("Missing covariate columns: ", paste(missing_covars, collapse=", "))
    covar_cols <- covar_cols[covar_cols %in% colnames(pheno)]
}

if (length(covar_cols) > 0) {
    covar_out <- pheno[, c("FID", "IID", covar_cols)]
    write.table(covar_out, 
                "${OUT_DIR}/covariates_regenie.txt", 
                row.names=FALSE, 
                col.names=TRUE, 
                quote=FALSE, 
                sep="\t")
    cat("Covariates included:", paste(covar_cols, collapse=", "), "\n")
} else {
    # Create empty covariate file if none specified
    covar_out <- pheno[, c("FID", "IID")]
    write.table(covar_out, 
                "${OUT_DIR}/covariates_regenie.txt", 
                row.names=FALSE, 
                col.names=TRUE, 
                quote=FALSE, 
                sep="\t")
    cat("No covariates included\n")
}

cat("Phenotype and covariate files prepared successfully\n")
EOF

if [[ $? -eq 0 ]]; then
    print_success "Phenotype and covariate files prepared"
else
    print_error "Failed to prepare phenotype/covariate files"
    exit 1
fi

# Step 4: REGENIE Step 1 (fit null model)
print_info "Step 4: Running REGENIE Step 1 (null model fitting)..."
regenie \
  --step 1 \
  --bed "${OUT_DIR}/qc_data" \
  --extract "${OUT_DIR}/qc_pass.prune.in" \
  --phenoFile "${OUT_DIR}/phenotype_regenie.txt" \
  --covarFile "${OUT_DIR}/covariates_regenie.txt" \
  --bsize "$BSIZE_STEP1" \
  --bt \
  --lowmem \
  --lowmem-prefix "${OUT_DIR}/tmp_rg" \
  --threads "$THREADS" \
  --force-step1 \
  --out "${OUT_DIR}/fit_out"

if [[ $? -eq 0 ]]; then
    print_success "REGENIE Step 1 completed"
else
    print_error "REGENIE Step 1 failed"
    exit 1
fi

# Step 5: REGENIE Step 2 (association testing)
print_info "Step 5: Running REGENIE Step 2 (association testing)..."
regenie \
  --step 2 \
  --bed "${OUT_DIR}/qc_data" \
  --phenoFile "${OUT_DIR}/phenotype_regenie.txt" \
  --covarFile "${OUT_DIR}/covariates_regenie.txt" \
  --bt \
  --firth --approx \
  --pThresh "$PTHRESH" \
  --bsize "$BSIZE_STEP2" \
  --threads "$THREADS" \
  --pred "${OUT_DIR}/fit_out_pred.list" \
  --out "${OUT_DIR}/test_out"

if [[ $? -eq 0 ]]; then
    print_success "REGENIE Step 2 completed"
else
    print_error "REGENIE Step 2 failed"
    exit 1
fi

# Step 6: Generate plots
if [[ "$SKIP_PLOTS" == false ]]; then
    print_info "Step 6: Generating plots..."
    
    Rscript - <<'EOF'
    library(qqman)
    
    # Get arguments from environment
    out_dir <- Sys.getenv("OUT_DIR")
    pheno_col <- Sys.getenv("PHENO_COL")
    
    # Find REGENIE output file
    result_file <- list.files(out_dir, pattern="test_out_.*\\.regenie$", full.names=TRUE)[1]
    
    if (is.na(result_file) || !file.exists(result_file)) {
        stop("REGENIE output file not found")
    }
    
    cat("Reading results from:", result_file, "\n")
    
    # Read REGENIE output
    results <- read.table(result_file, header=TRUE)
    
    cat("Total variants:", nrow(results), "\n")
    cat("Columns:", colnames(results), "\n")
    
    # Prepare for qqman
    results$CHR <- as.numeric(gsub("chr", "", results$CHROM))
    results$P <- 10^(-results$LOG10P)
    
    # Remove non-autosomal and invalid values
    results <- results[results$CHR %in% 1:22, ]
    results <- results[!is.na(results$P) & is.finite(results$LOG10P), ]
    
    cat("Valid autosomal variants:", nrow(results), "\n")
    
    # Calculate lambda
    chisq <- qchisq(1 - results$P, 1)
    lambda <- median(chisq, na.rm=TRUE) / qchisq(0.5, 1)
    
    # Calculate statistics
    n_gws <- sum(results$P < 5e-8, na.rm=TRUE)
    n_suggestive <- sum(results$P < 1e-5, na.rm=TRUE)
    
    cat("Lambda GC:", lambda, "\n")
    cat("Genome-wide significant variants:", n_gws, "\n")
    cat("Suggestive variants:", n_suggestive, "\n")
    
    # Manhattan Plot
    png(file.path(out_dir, "manhattan_plot.png"), 
        width=1800, height=900, res=150, type="cairo")
    
    par(mar=c(5.5, 5.5, 4.5, 2), mgp=c(3.5, 1, 0))
    
    manhattan(results, 
              chr="CHR", bp="GENPOS", snp="ID", p="P",
              col=c("#3B82F6", "#F59E0B", "#10B981", "#EC4899", 
                    "#8B5CF6", "#EF4444", "#06B6D4", "#F97316",
                    "#14B8A6", "#A855F7", "#84CC16"),
              suggestiveline=-log10(1e-5),
              genomewideline=-log10(5e-8),
              cex=0.8, cex.axis=1.4, cex.lab=1.6,
              main=paste0("Genome-Wide Association Study\nλ = ", round(lambda, 3), 
                         " | Genome-wide: ", n_gws, " | Suggestive: ", n_suggestive),
              cex.main=1.8,
              xlab="Chromosome",
              ylab=expression(bold(-log[10](italic(P)))),
              ylim=c(0, max(-log10(results$P), na.rm=TRUE) * 1.1),
              lwd=2.5)
    
    dev.off()
    
    # QQ Plot
    png(file.path(out_dir, "qq_plot.png"), 
        width=1000, height=1000, res=150, type="cairo")
    
    par(mar=c(5.5, 5.5, 4.5, 2), mgp=c(3.5, 1, 0))
    
    qq(results$P, col="#3B82F6", pch=19, cex=1,
       cex.axis=1.4, cex.lab=1.6,
       main=paste0("Quantile-Quantile Plot\nλ = ", round(lambda, 3)),
       cex.main=1.8,
       xlab=expression(bold(Expected~~-log[10](italic(P)))),
       ylab=expression(bold(Observed~~-log[10](italic(P)))),
       xlim=c(0, max(-log10(ppoints(length(results$P))))),
       ylim=c(0, max(-log10(results$P), na.rm=TRUE)))
    
    text(x=max(-log10(ppoints(length(results$P))))*0.15,
         y=max(-log10(results$P), na.rm=TRUE)*0.95,
         labels=bquote(bold(lambda == .(round(lambda, 3)))),
         cex=2.2, col="red3", font=2)
    
    text(x=max(-log10(ppoints(length(results$P))))*0.15,
         y=max(-log10(results$P), na.rm=TRUE)*0.88,
         labels=paste0("N = ", format(nrow(results), big.mark=",")),
         cex=1.6, col="gray30", font=2)
    
    dev.off()
    
    # Combined Plot
    png(file.path(out_dir, "combined_plot.png"), 
        width=2200, height=1000, res=150, type="cairo")
    
    layout(matrix(c(1,1,1,2,2), nrow=1))
    
    par(mar=c(5, 5, 4, 1))
    manhattan(results, chr="CHR", bp="GENPOS", snp="ID", p="P",
              col=c("#3B82F6", "#F59E0B", "#10B981"),
              suggestiveline=-log10(1e-5), genomewideline=-log10(5e-8),
              cex=0.7, cex.axis=1.3, cex.lab=1.5,
              main="Manhattan Plot", cex.main=1.7, lwd=2.5,
              xlab="Chromosome",
              ylab=expression(bold(-log[10](italic(P)))))
    
    par(mar=c(5, 5, 4, 2))
    qq(results$P, col="#3B82F6", pch=19, cex=0.9,
       cex.axis=1.3, cex.lab=1.5, main="Q-Q Plot", cex.main=1.7,
       xlab=expression(bold(Expected~~-log[10](italic(P)))),
       ylab=expression(bold(Observed~~-log[10](italic(P)))))
    
    mtext(paste0("GWAS Results Summary (λ = ", round(lambda, 3), 
                 " | Significant: ", n_gws, ")"),
          side=3, line=-2, outer=TRUE, cex=1.5, font=2)
    
    dev.off()
    
    # Summary Report
    cat("\n===========================================\n")
    cat("           GWAS RESULTS SUMMARY\n")
    cat("===========================================\n")
    cat("Total variants analyzed:", format(nrow(results), big.mark=","), "\n")
    cat("Genomic inflation (λ):", round(lambda, 4), "\n")
    cat("Minimum p-value:", format(min(results$P), scientific=TRUE), "\n")
    cat("Maximum -log10(P):", round(max(-log10(results$P)), 2), "\n")
    cat("Genome-wide significant (P<5e-8):", n_gws, "\n")
    cat("Suggestive (P<1e-5):", n_suggestive, "\n")
    cat("\n===========================================\n")
    cat("           TOP 20 VARIANTS\n")
    cat("===========================================\n")
    
    top_variants <- results[order(results$P), ][1:min(20, nrow(results)), 
                                                 c("ID", "CHR", "GENPOS", "P", "LOG10P")]
    print(top_variants, row.names=FALSE)
    
    # Save summary to file
    sink(file.path(out_dir, "gwas_summary.txt"))
    cat("GWAS RESULTS SUMMARY\n")
    cat("===========================================\n")
    cat("Analysis date:", format(Sys.Date()), "\n")
    cat("Total variants:", format(nrow(results), big.mark=","), "\n")
    cat("Genomic inflation (λ):", round(lambda, 4), "\n")
    cat("Genome-wide significant (P<5e-8):", n_gws, "\n")
    cat("Suggestive (P<1e-5):", n_suggestive, "\n\n")
    cat("TOP 20 VARIANTS\n")
    cat("===========================================\n")
    print(top_variants, row.names=FALSE)
    sink()
    
    cat("\n===========================================\n")
    cat("           PLOTS GENERATED\n")
    cat("===========================================\n")
    cat("  ✓ manhattan_plot.png\n")
    cat("  ✓ qq_plot.png\n")
    cat("  ✓ combined_plot.png\n")
    cat("  ✓ gwas_summary.txt\n")
    cat("===========================================\n")
EOF

    export OUT_DIR
    export PHENO_COL
    
    if [[ $? -eq 0 ]]; then
        print_success "Plots generated successfully"
    else
        print_warning "Plot generation encountered issues"
    fi
else
    print_info "Skipping plot generation (--no-plots flag set)"
fi

# Final summary
print_success "======================================"
print_success "REGENIE GWAS Pipeline Completed!"
print_success "======================================"
print_info "Results saved to: $OUT_DIR"
print_info "Key output files:"
print_info "  - ${OUT_DIR}/test_out_*.regenie (association results)"
print_info "  - ${OUT_DIR}/manhattan_plot.png"
print_info "  - ${OUT_DIR}/qq_plot.png"
print_info "  - ${OUT_DIR}/combined_plot.png"
print_info "  - ${OUT_DIR}/gwas_summary.txt"
print_success "======================================"