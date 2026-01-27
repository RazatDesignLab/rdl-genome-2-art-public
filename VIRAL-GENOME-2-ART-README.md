# Viral Genome-2-Art: SARS-CoV-2 Quantum Wave Visualization

Transform SARS-CoV-2 viral genome variants into stunning quantum wave interference patterns and artistic visualizations.

## Overview

This project extends the genome-2-art concept to viral genomes, specifically SARS-CoV-2 data from NCBI's public datasets. It converts genetic variants in viral samples into quantum wave patterns, creating beautiful scientific art while revealing the underlying structure and diversity of viral genomes.

## Features

- **Circular Genome Visualization**: Display viral variants as radial spikes around a circular genome backbone
- **Quantum Wave Fields**: Generate interference patterns from multiple viral samples using wave physics
- **Sample Comparison**: Compare variant positions and patterns across different viral samples
- **Multi-format Output**: Export visualizations as HTML, PNG, and SVG formats
- **Statistics Analysis**: Detailed breakdown of variant counts and heterozygosity rates

## Public Dataset: NCBI SARS-CoV-2 Coronaviridae Collection

### Dataset Information

This tool uses publicly available SARS-CoV-2 genomic data from the **NCBI SRA (Sequence Read Archive) AWS Registry of Open Data**. The dataset is part of the **ACTIV TRACE initiative** and contains a comprehensive collection of SARS-CoV-2 variant data.

**Source**: NCBI SRA Coronaviridae Datasets on AWS
**URL**: https://registry.opendata.aws/ncbi-covid-19/
**AWS Bucket**: `s3://sra-pub-sars-cov2/`
**Documentation**: https://www.ncbi.nlm.nih.gov/sra/docs/sra-aws-download/

### Dataset Details

The SARS-CoV-2 dataset includes:

- **VCF Directory** (`vcf/`): Variant Call Format files generated from SRA data, organized by run accession ID
- **Original Data** (`sra-src/`): Submitted sequence files in their original format
- **Processed Data** (`run/`): Normalized sequence data accessible via SRA Toolkit
- **Metadata**: Available in parquet format for analysis with AWS Athena Service

### Data Processing Pipeline

The dataset uses the **SARS-CoV-2 Variant Calling Pipeline** which:

- Processes Sequence Read Archive (SRA) files optimized for viral variant identification
- Identifies genetic variations in viral samples compared to the reference genome (NC_045512.2)
- Presents variations in standard VCF (Variant Call Format) files
- Updates data every 6 hours with daily AWS bucket updates
- Provides detailed variant annotations using SnpEff

### Access Information

**No AWS Account Required** - Data can be accessed freely using:

```bash
# List available data
aws s3 ls s3://sra-pub-sars-cov2/ --no-sign-request

# List VCF files
aws s3 ls s3://sra-pub-sars-cov2/vcf/ --no-sign-request

# Download specific sample VCF
aws s3 cp s3://sra-pub-sars-cov2/vcf/DRR259112/DRR259112.ref.snpeff.vcf . --no-sign-request
```

**Web Console Access**: https://s3.console.aws.amazon.com/s3/buckets/sra-pub-sars-cov2/

## Installation & Requirements

### Prerequisites

- Python 3.8+
- Virtual environment (recommended)
- AWS CLI (for data download)

### Python Dependencies

```bash
pip install -r requirements.txt
```

Required packages:

- `pandas` - Data manipulation and analysis
- `numpy` - Numerical computations
- `plotly` - Interactive visualizations
- `kaleido` - Static image export

### Setup

1. Clone or download the repository
2. Create and activate virtual environment:
   ```bash
   python -m venv .venv
   source .venv/bin/activate  # On Windows: .venv\Scripts\activate
   ```
3. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```
4. Create public genome data directory:
   ```bash
   mkdir -p public_genomes
   ```

## Usage

### 1. Download SARS-CoV-2 Data

Download sample VCF files to the `public_genomes/` directory:

```bash
# Download multiple samples
aws s3 cp s3://sra-pub-sars-cov2/vcf/DRR259112/DRR259112.ref.snpeff.vcf public_genomes/ --no-sign-request
aws s3 cp s3://sra-pub-sars-cov2/vcf/DRR259113/DRR259113.ref.snpeff.vcf public_genomes/ --no-sign-request
aws s3 cp s3://sra-pub-sars-cov2/vcf/DRR272391/DRR272391.ref.snpeff.vcf public_genomes/ --no-sign-request
```

### 2. Run the Visualization

```bash
# Activate virtual environment
source .venv/bin/activate

# Run the viral genome art generator
python viral-genome-2-art.py
```

### 3. View Results

The script generates:

- `viral-genome-2-art.html` - Interactive web visualization
- `art/viral-genome-wave.png` - High-resolution wave field image
- `art/viral-genome-wave.svg` - Vector format wave field

## Understanding the Visualizations

### Circular Genome Plot

- **Circle**: Represents the ~30kb SARS-CoV-2 genome
- **Position markers**: Every 5kb around the circle
- **Radial spikes**: Individual variants, length indicates impact
- **Colors**: Encode variant frequency and amplitude properties

### Quantum Wave Field

- **Wave interference**: Combined patterns from all samples
- **Colors**: Scientific colorscale (cyan-blue to magenta)
- **Intensity**: Represents quantum wave amplitude
- **Patterns**: Show genetic relationships between samples

### Sample Comparison

- **Scatter plot**: Variant positions across samples
- **Colors**: Distinguish between samples and variant types
- **Y-axis**: Each sample on separate row
- **X-axis**: Genomic position (base pairs)

### Statistics Table

- **Sample names**: Run accession identifiers
- **Variant counts**: Total number of variants per sample
- **Heterozygous percentage**: Proportion of mixed genotypes

## Scientific Interpretation

### Wave Properties

Each nucleotide is assigned wave characteristics:

- **A (Adenine)**: Low frequency (1.5), high amplitude (1.0)
- **T (Thymine)**: Medium frequency (2.8), high amplitude (0.9)
- **C (Cytosine)**: High frequency (4.2), medium amplitude (0.7)
- **G (Guanine)**: Very high frequency (5.5), low amplitude (0.5)

### Variant Types

- **SNVs (Single Nucleotide Variants)**: Point mutations
- **Indels**: Insertions and deletions
- **Homozygous**: Same allele (pure waves)
- **Heterozygous**: Mixed alleles (interference patterns)

### Biological Significance

- **High variant density**: Indicates rapidly evolving genome regions
- **Wave interference**: Shows genetic relationships between samples
- **Heterozygosity patterns**: Reveal mixed infections or sequencing artifacts

## Data Format Details

### VCF File Structure

The input VCF files contain:

- **CHROM**: Chromosome (NC_045512.2 for SARS-CoV-2)
- **POS**: Position in genome (1-29,903)
- **REF**: Reference allele
- **ALT**: Alternative allele
- **QUAL**: Variant quality score
- **FILTER**: Quality filters applied
- **INFO**: Variant annotations (SnpEff)
- **FORMAT**: Sample format specification
- **Sample data**: Genotype and quality metrics

### Reference Genome

- **NC_045512.2**: SARS-CoV-2 reference genome
- **Length**: 29,903 base pairs
- **Genes**: ORF1ab, S, ORF3a, E, M, ORF6, ORF7a, ORF8, N, ORF10

## Example Analysis

Using the included sample data:

```
🧬 Analysis Summary:
   • DRR272391: 10 variants, 10.0% heterozygous, span 28,640 bp
   • DRR259112: 36 variants, 86.1% heterozygous, span 28,854 bp
   • DRR259113: 20 variants, 70.0% heterozygous, span 28,908 bp
```

This reveals:

- **DRR259112**: Most variants (36), high heterozygosity - possibly mixed sample
- **DRR272391**: Fewest variants (10), low heterozygosity - clean sample
- **DRR259113**: Intermediate complexity (20 variants, 70% heterozygous)

## Customization

### Resolution Settings

Modify wave field resolution in the code:

```python
wave_field, wx, wy = create_viral_wave_field(df_list, sample_names, resolution=600)
```

### Color Schemes

Customize the scientific colorscale:

```python
scientific_colorscale = [
    [0.0, '#000011'],   [0.1, '#001133'],   [0.2, '#003366'],
    [0.3, '#0066AA'],   [0.4, '#00BFFF'],   [0.5, '#7F5FFF'],
    [0.6, '#BF3FFF'],   [0.7, '#FF00FF'],   [0.8, '#FF66FF'],
    [0.9, '#FFAAFF'],   [1.0, '#FFFFFF']
]
```

### Wave Properties

Adjust nucleotide wave mapping:

```python
nucleotide_wave_map = {
    'A': {'freq': 1.5, 'amp': 1.0},
    'T': {'freq': 2.8, 'amp': 0.9},
    'C': {'freq': 4.2, 'amp': 0.7},
    'G': {'freq': 5.5, 'amp': 0.5}
}
```

## Troubleshooting

### Common Issues

1. **No VCF files found**: Ensure files are in `public_genomes/` directory
2. **Color errors**: Check that RGB values are within 0-255 range
3. **Memory issues**: Reduce resolution parameter for large datasets
4. **Export failures**: Ensure kaleido is properly installed

### Data Quality

- **Low variant counts**: May indicate high-quality reference strain
- **High heterozygosity**: Could suggest mixed infections or artifacts
- **Missing positions**: Normal for high-quality viral samples

## Contributing

Contributions welcome! Areas for enhancement:

- Additional viral genome support (influenza, HIV, etc.)
- Interactive parameter tuning
- Phylogenetic tree integration
- Real-time data streaming from NCBI

## Citation

If using this tool for research, please cite:

- The NCBI SRA Coronaviridae dataset
- The ACTIV TRACE initiative
- Original genome-2-art methodology

## License

This project is open source. The underlying SARS-CoV-2 data is freely available through NCBI's Registry of Open Data on AWS.

## Acknowledgments

- **NCBI SRA team** for maintaining the public SARS-CoV-2 dataset
- **ACTIV TRACE initiative** for coordinated variant surveillance
- **AWS Registry of Open Data** for free data access
- **Plotly team** for excellent visualization tools

---

_Transform viral genomics into quantum art - where science meets creativity._
