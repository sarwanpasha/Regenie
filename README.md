# REGENIE GWAS Pipeline

A portable, user-friendly pipeline for running genome-wide association studies (GWAS) using REGENIE with automatic quality control, association testing, and visualization.

## Features

- ✅ Automated quality control of genotype data
- ✅ REGENIE two-step association analysis
- ✅ Automatic generation of Manhattan and QQ plots
- ✅ Comprehensive summary statistics
- ✅ Flexible parameter configuration
- ✅ Portable and easy to use

## Requirements

### Software Dependencies
- [PLINK 1.9+](https://www.cog-genomics.org/plink/)
- [REGENIE v3.0+](https://rgcgithub.github.io/regenie/)
- [R 4.0+](https://www.r-project.org/)
- R packages: `qqman`

### System Requirements
- Linux/Unix environment
- Sufficient disk space for temporary files and results
- At least 8GB RAM (depends on dataset size)

## Installation

### 1. Clone the repository
```bash
git clone https://github.com/yourusername/regenie-gwas-pipeline.git
cd regenie-gwas-pipeline
```

### 2. Make the script executable
```bash
chmod +x run_regenie_gwas.sh
```

### 3. Install R package (if not already installed)
```R
install.packages("qqman")
```

### 4. Verify dependencies
```bash
# Check PLINK
plink --version

# Check REGENIE
regenie --version

# Check R
R --version
```

## Quick Start

### Basic Usage
```bash
./run_regenie_gwas.sh \
  -b /path/to/your/plink_data \
  -p /path/to/your/phenotype.txt \
  -o ./results
```

### With Custom Covariates
```bash
./run_regenie_gwas.sh \
  -b /path/to/your/plink_data \
  -p /path/to/your/phenotype.txt \
  -o ./results \
  -c "AGE,SEX,PC1,PC2,PC3"
```

### Advanced Usage
```bash
./run_regenie_gwas.sh \
  -b /path/to/your/plink_data \
  -p /path/to/your/phenotype.txt \
  -o ./results \
  -c "AGE,SEX,PC1,PC2,PC3,PC4,PC5" \
  -y "DISEASE_STATUS" \
  --maf 0.01 \
  --geno 0.05 \
  --hwe 1e-6 \
  --threads 8
```

## Input File Formats

### PLINK Binary Files
Standard PLINK binary format (`.bed`, `.bim`, `.fam`). Specify the prefix without extensions.

Example:
```
mydata.bed
mydata.bim
mydata.fam
```
Use: `-b mydata`

### Phenotype File
Tab or space-delimited text file with header. Required columns:
- `FID`: Family ID
- `IID`: Individual ID
- `PHENO`: Phenotype (or custom column name)
- Additional covariate columns

Example (`phenotype.txt`):
```
FID    IID    PHENO    AGE    SEX    PC1        PC2
FAM1   ID1    1        45     1      0.023      -0.012
FAM2   ID2    0        52     2      -0.015     0.031
FAM3   ID3    1        38     1      0.008      0.019
```

**Phenotype Encoding:**
- Binary traits: 0 = control, 1 = case
- Missing: NA or -9

## Command-Line Options

### Required Arguments
| Option | Description |
|--------|-------------|
| `-b, --bfile PREFIX` | Path to PLINK binary file prefix |
| `-p, --pheno FILE` | Path to phenotype file |
| `-o, --outdir DIR` | Output directory |

### Optional Arguments
| Option | Default | Description |
|--------|---------|-------------|
| `-c, --covar COLS` | COV1,COV2 | Comma-separated covariate column names |
| `-y, --pheno-col NAME` | PHENO | Phenotype column name |
| `--maf FLOAT` | 0.01 | Minor allele frequency threshold |
| `--geno FLOAT` | 0.05 | Genotype missingness threshold |
| `--mind FLOAT` | 0.05 | Sample missingness threshold |
| `--hwe FLOAT` | 1e-6 | Hardy-Weinberg p-value threshold |
| `--bsize INT` | 1000/400 | Block size for REGENIE |
| `--pthresh FLOAT` | 0.01 | P-value threshold for Firth regression |
| `--threads INT` | 4 | Number of threads |
| `--no-plots` | - | Skip plot generation |
| `-h, --help` | - | Display help message |

## Output Files

The pipeline generates the following outputs in the specified output directory:

### Association Results
- `test_out_*.regenie`: Main REGENIE association results
- `gwas_summary.txt`: Text summary of results

### Quality Control Files
- `qc_data.{bed,bim,fam}`: QC-filtered genotype data
- `qc_pass.prune.in`: List of pruned SNPs for Step 1

### Intermediate Files
- `phenotype_regenie.txt`: Formatted phenotype file
- `covariates_regenie.txt`: Formatted covariate file
- `fit_out*`: REGENIE Step 1 outputs

### Visualizations
- `manhattan_plot.png`: Manhattan plot with genome-wide significant variants
- `qq_plot.png`: Q-Q plot showing genomic inflation
- `combined_plot.png`: Side-by-side Manhattan and Q-Q plots

## Example Workflows

### 1. Basic Case-Control Study
```bash
./run_regenie_gwas.sh \
  -b ./data/mydata \
  -p ./data/pheno.txt \
  -o ./results/basic_gwas
```

### 2. Study with Principal Components
```bash
./run_regenie_gwas.sh \
  -b ./data/mydata \
  -p ./data/pheno.txt \
  -o ./results/pc_adjusted \
  -c "PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
```

### 3. Stringent QC Parameters
```bash
./run_regenie_gwas.sh \
  -b ./data/mydata \
  -p ./data/pheno.txt \
  -o ./results/stringent_qc \
  -c "AGE,SEX,PC1,PC2,PC3" \
  --maf 0.05 \
  --geno 0.02 \
  --mind 0.02 \
  --hwe 1e-10
```

### 4. High-Performance Computing
```bash
./run_regenie_gwas.sh \
  -b ./data/mydata \
  -p ./data/pheno.txt \
  -o ./results/hpc_run \
  -c "AGE,SEX,BMI,PC1,PC2,PC3,PC4,PC5" \
  --threads 16 \
  --bsize 2000
```

## Interpreting Results

### Manhattan Plot
- **Y-axis**: -log10(P-value)
- **Red line**: Genome-wide significance (P < 5×10⁻⁸)
- **Blue line**: Suggestive significance (P < 1×10⁻⁵)
- Points above the red line indicate genome-wide significant associations

### QQ Plot
- **Lambda (λ)**: Genomic inflation factor
  - λ ≈ 1.0: No inflation (ideal)
  - λ > 1.05: Possible population stratification or polygenic signal
  - λ < 0.95: Possible deflation
- Deviation from diagonal indicates true associations or systematic bias

### Summary Statistics
Check `gwas_summary.txt` for:
- Total variants analyzed
- Genomic inflation factor (λ)
- Number of genome-wide significant variants
- Top 20 variants by P-value

## Troubleshooting

### Common Issues

**1. "PLINK binary files not found"**
- Ensure you specify the prefix without file extensions
- Check that `.bed`, `.bim`, and `.fam` files all exist

**2. "Missing required columns in phenotype file"**
- Verify that `FID`, `IID`, and phenotype column exist
- Check that column names match exactly (case-sensitive)

**3. "Covariates not found"**
- Ensure covariate column names match those in the phenotype file
- Use exact column names (case-sensitive)

**4. Memory issues**
- Reduce `--bsize` parameter
- Use `--lowmem` flag (already included in pipeline)
- Increase available RAM

**5. Module load errors**
- Comment out the `module load` lines if not using HPC
- Ensure software is in your PATH

## HPC/Cluster Usage

For HPC environments, you can submit the pipeline as a job:

```bash
#!/bin/bash
#SBATCH --job-name=regenie_gwas
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8

module load Plink/1.9.10
module load R/4.5.2

./run_regenie_gwas.sh \
  -b /path/to/data \
  -p /path/to/pheno.txt \
  -o /path/to/results \
  -c "COV1,COV2,COV3" \
  --threads 8
```

## Citation

If you use this pipeline, please cite:

**REGENIE:**
- Mbatchou et al. (2021) Nature Genetics. "Computationally efficient whole-genome regression for quantitative and binary traits"

**PLINK:**
- Chang et al. (2015) GigaScience. "Second-generation PLINK: rising to the challenge of larger and richer datasets"

**qqman:**
- Turner (2014) "qqman: an R package for visualizing GWAS results using Q-Q and manhattan plots"

## License

MIT License - see LICENSE file for details

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request

## Support

For issues and questions:
- Open an issue on GitHub
- Check existing issues for solutions
- Provide example data and error messages

## Authors

[Your Name] - Initial work

## Acknowledgments

- REGENIE development team
- PLINK development team
- qqman package authors
