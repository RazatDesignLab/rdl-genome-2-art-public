import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio
import kaleido
import argparse
import sys
import os

# ─── Color Palette Definitions ───────────────────────────────────────────────
COLOR_PALETTES = {
    'scientific': [
        [0.0, '#000000'],   # Black for background
        [0.1, '#001133'],   # Very dark blue
        [0.2, '#003366'],   # Dark blue-cyan
        [0.3, '#0066AA'],   # Medium blue-cyan
        [0.4, '#00BFFF'],   # Pure cyan blue
        [0.5, '#7F5FFF'],   # Cyan-magenta blend
        [0.6, '#BF3FFF'],   # Blue-magenta blend
        [0.7, '#FF00FF'],   # Pure magenta
        [0.8, '#FF66FF'],   # Light magenta
        [0.9, '#FFAAFF'],   # Very light magenta
        [1.0, '#010000'],   # Near-black for highlights
    ],
    'ocean': [
        [0.0, '#000000'],
        [0.1, '#000A1A'],
        [0.2, '#001A3A'],
        [0.3, '#003366'],
        [0.4, '#005599'],
        [0.5, '#0077BB'],
        [0.6, '#0099DD'],
        [0.7, '#33BBEE'],
        [0.8, '#66DDFF'],
        [0.9, '#AAEEFF'],
        [1.0, '#DDFAFF'],
    ],
    'fire': [
        [0.0, '#000000'],
        [0.1, '#1A0000'],
        [0.2, '#4D0000'],
        [0.3, '#800000'],
        [0.4, '#B30000'],
        [0.5, '#E60000'],
        [0.6, '#FF3300'],
        [0.7, '#FF6600'],
        [0.8, '#FF9933'],
        [0.9, '#FFCC66'],
        [1.0, '#FFE599'],
    ],
    'nebula': [
        [0.0, '#000000'],
        [0.1, '#0D0015'],
        [0.2, '#1A002A'],
        [0.3, '#330055'],
        [0.4, '#550088'],
        [0.5, '#7700AA'],
        [0.6, '#9933CC'],
        [0.7, '#BB66DD'],
        [0.8, '#DD99EE'],
        [0.9, '#EEBBF5'],
        [1.0, '#F8DDFA'],
    ],
    'earth': [
        [0.0, '#000000'],
        [0.1, '#1A0F00'],
        [0.2, '#332200'],
        [0.3, '#554400'],
        [0.4, '#776622'],
        [0.5, '#998844'],
        [0.6, '#AA9955'],
        [0.7, '#BBAA66'],
        [0.8, '#CCBB88'],
        [0.9, '#DDCCAA'],
        [1.0, '#EEDDCC'],
    ],
    'aurora': [
        [0.0, '#000000'],
        [0.1, '#001111'],
        [0.2, '#002222'],
        [0.3, '#004433'],
        [0.4, '#006644'],
        [0.5, '#008855'],
        [0.6, '#00AA66'],
        [0.7, '#33CC77'],
        [0.8, '#66EE88'],
        [0.9, '#99FF99'],
        [1.0, '#CCFFCC'],
    ],
    'infrared': [
        [0.0, '#000000'],
        [0.1, '#0A0010'],
        [0.2, '#1A0033'],
        [0.3, '#2D0055'],
        [0.4, '#440077'],
        [0.5, '#5500AA'],
        [0.6, '#7733BB'],
        [0.7, '#9966CC'],
        [0.8, '#BB99DD'],
        [0.9, '#CCAAEE'],
        [1.0, '#E0CCF5'],
    ],
}


def get_palette_legend_colors(palette_name):
    """Extract representative low and high colors from a palette for legends."""
    palette = COLOR_PALETTES[palette_name]
    low_color = palette[4][1]   # Position 0.4 — low intensity
    high_color = palette[7][1]  # Position 0.7 — high intensity
    return low_color, high_color


def load_ancestry_data(filepath):
    """Load AncestryDNA raw data file.

    Handles both tab-delimited (real AncestryDNA exports) and
    whitespace-delimited (mock/test) files robustly.
    """
    df = pd.read_csv(filepath, sep=r'\s+', comment='#',
                     names=['rsid', 'chromosome', 'position', 'allele1', 'allele2'],
                     engine='python', dtype={'chromosome': str})

    # Drop header rows read as data
    df = df[df['rsid'] != 'rsid']

    # Clean and convert data
    df['chromosome'] = df['chromosome'].astype(str).str.replace('chr', '')
    df['position'] = pd.to_numeric(df['position'], errors='coerce')
    df = df.dropna(subset=['position'])
    df['position'] = df['position'].astype(int)

    # Filter to standard chromosomes
    valid_chroms = [str(i) for i in range(1, 23)] + ['X', 'Y']
    df = df[df['chromosome'].isin(valid_chroms)]

    return df

def genotype_to_wave_properties(allele1, allele2):
    """Convert genotype to wave properties (amplitude and frequency)"""
    # Map nucleotides to wave characteristics
    nucleotide_wave_map = {
        'A': {'freq': 1.0, 'amp': 1.0},      # Low frequency, high amplitude
        'T': {'freq': 2.0, 'amp': 0.8},      # Medium frequency, medium amplitude
        'C': {'freq': 3.0, 'amp': 0.6},      # High frequency, low amplitude
        'G': {'freq': 4.0, 'amp': 0.4},      # Very high frequency, very low amplitude
        '--': {'freq': 0.0, 'amp': 0.0}      # No wave
    }

    wave1 = nucleotide_wave_map.get(str(allele1), {'freq': 0, 'amp': 0})
    wave2 = nucleotide_wave_map.get(str(allele2), {'freq': 0, 'amp': 0})

    # Combine alleles to create unique wave signature
    # Heterozygous creates beat frequencies, homozygous creates pure tones
    if allele1 == allele2:  # Homozygous - pure wave
        frequency = wave1['freq']
        amplitude = wave1['amp']
        phase_shift = 0
    else:  # Heterozygous - creates interference patterns
        frequency = (wave1['freq'] + wave2['freq']) / 2
        amplitude = (wave1['amp'] + wave2['amp']) / 2
        # Create phase difference for interference
        phase_shift = abs(wave1['freq'] - wave2['freq']) * np.pi / 4

    return frequency, amplitude, phase_shift

def create_circos_plot(df, max_snps=2000):
    """Create a CIRCOS-style plot showing genome organization"""

    print(f"Creating CIRCOS plot from genome data...")

    if len(df) > max_snps:
        df_sample = pd.concat([
            group.sample(n=max(1, min(len(group), max_snps // 24)))
            for _, group in df.groupby('chromosome')
        ]).reset_index(drop=True)
    else:
        df_sample = df

    # Calculate chromosome statistics
    chrom_stats = []
    for chrom in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12',
                  '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']:
        chrom_data = df_sample[df_sample['chromosome'] == chrom]
        if len(chrom_data) > 0:
            # Calculate heterozygosity rate for this chromosome
            het_rate = sum(1 for _, row in chrom_data.iterrows()
                          if str(row['allele1']) != str(row['allele2']) and
                          str(row['allele1']) != '--' and str(row['allele2']) != '--') / len(chrom_data)

            # Calculate average position (for relative sizing)
            avg_pos = chrom_data['position'].mean() / 1e6  # Convert to Mb
            max_pos = chrom_data['position'].max() / 1e6

            chrom_stats.append({
                'chromosome': chrom,
                'het_rate': het_rate,
                'snp_count': len(chrom_data),
                'avg_pos_mb': avg_pos,
                'max_pos_mb': max_pos
            })

    # Create circular layout
    n_chroms = len(chrom_stats)
    angles = np.linspace(0, 2*np.pi, n_chroms, endpoint=False)

    # Create the CIRCOS plot traces
    traces = []

    # 1. Outer ring - Chromosome boundaries
    outer_r = 10
    inner_r = 8

    for i, (angle, stats) in enumerate(zip(angles, chrom_stats)):
        # Create chromosome arc
        arc_angles = np.linspace(angle - 0.12, angle + 0.12, 50)

        # Outer boundary
        x_outer = outer_r * np.cos(arc_angles)
        y_outer = outer_r * np.sin(arc_angles)

        # Inner boundary
        x_inner = inner_r * np.cos(arc_angles)
        y_inner = inner_r * np.sin(arc_angles)

        # Chromosome arc
        x_arc = np.concatenate([x_outer, x_inner[::-1], [x_outer[0]]])
        y_arc = np.concatenate([y_outer, y_inner[::-1], [y_outer[0]]])

        traces.append(go.Scatter(
            x=x_arc, y=y_arc,
            fill='toself',
            fillcolor=f'rgba({50 + i*8}, {100 + i*6}, {150 + i*4}, 0.7)',
            line=dict(color='white', width=1),
            name=f'Chr {stats["chromosome"]}',
            hovertemplate=f'Chromosome {stats["chromosome"]}<br>SNPs: {stats["snp_count"]}<br>Heterozygosity: {stats["het_rate"]:.1%}<extra></extra>',
            showlegend=False
        ))

        # Add chromosome labels with better positioning for dark background
        label_r = outer_r + 1.5
        label_x = label_r * np.cos(angle)
        label_y = label_r * np.sin(angle)

        traces.append(go.Scatter(
            x=[label_x], y=[label_y],
            mode='text',
            text=[stats['chromosome']],
            textfont=dict(size=12, color='white'),  # Changed to white for dark background
            showlegend=False,
            hoverinfo='skip'
        ))

    # 2. Inner ring - Heterozygosity visualization
    het_r = 6
    for i, (angle, stats) in enumerate(zip(angles, chrom_stats)):
        # Heterozygosity bar height
        bar_height = stats['het_rate'] * 2  # Scale for visibility

        bar_angles = np.linspace(angle - 0.08, angle + 0.08, 20)
        x_bar_outer = (het_r + bar_height) * np.cos(bar_angles)
        y_bar_outer = (het_r + bar_height) * np.sin(bar_angles)
        x_bar_inner = het_r * np.cos(bar_angles)
        y_bar_inner = het_r * np.sin(bar_angles)

        x_bar = np.concatenate([x_bar_outer, x_bar_inner[::-1], [x_bar_outer[0]]])
        y_bar = np.concatenate([y_bar_outer, y_bar_inner[::-1], [y_bar_outer[0]]])

        # Color based on heterozygosity rate
        color_intensity = stats['het_rate']
        color = f'rgba({int(255*color_intensity)}, {int(150*(1-color_intensity))}, {int(100*color_intensity)}, 0.8)'

        traces.append(go.Scatter(
            x=x_bar, y=y_bar,
            fill='toself',
            fillcolor=color,
            line=dict(color='white', width=0.5),
            showlegend=False,
            hovertemplate=f'Chr {stats["chromosome"]} Heterozygosity: {stats["het_rate"]:.1%}<extra></extra>'
        ))

    return traces

def create_snp_distribution_plot(df, max_snps=3000):
    """Create circular SNP distribution plot showing individual SNPs radiating from chromosomes"""

    print(f"Creating SNP distribution plot...")

    if len(df) > max_snps:
        df_sample = pd.concat([
            group.sample(n=max(1, min(len(group), max_snps // 24)))
            for _, group in df.groupby('chromosome')
        ]).reset_index(drop=True)
    else:
        df_sample = df

    traces = []

    # Create chromosome positions (same as CIRCOS for consistency)
    chrom_list = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12',
                  '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']

    chrom_angles = {}
    for i, chrom in enumerate(chrom_list):
        chrom_angles[chrom] = i * 2 * np.pi / len(chrom_list)

    # Add chromosome markers
    for chrom in chrom_list:
        angle = chrom_angles[chrom]
        label_r = 10
        label_x = label_r * np.cos(angle)
        label_y = label_r * np.sin(angle)

        traces.append(go.Scatter(
            x=[label_x], y=[label_y],
            mode='markers+text',
            marker=dict(size=15, color='#2C3E50', symbol='circle'),
            text=[chrom],
            textposition='middle center',
            textfont=dict(size=10, color='white', family="Arial Black"),
            name=f'Chr {chrom}',
            showlegend=False,
            hovertemplate=f'Chromosome {chrom}<extra></extra>'
        ))

    # Create SNP rays for each chromosome
    het_snps_x, het_snps_y = [], []
    hom_snps_x, hom_snps_y = [], []

    for _, row in df_sample.iterrows():
        chrom = str(row['chromosome'])
        if chrom in chrom_angles:
            # Get chromosome angle
            base_angle = chrom_angles[chrom]

            # Create slight angular spread based on position within chromosome
            position_spread = (row['position'] % 1000000) / 1000000 * 0.3 - 0.15  # ±0.15 radians
            snp_angle = base_angle + position_spread

            # Distance based on position magnitude
            distance = 2 + (row['position'] % 50000000) / 50000000 * 6  # Between 2-8 units

            snp_x = distance * np.cos(snp_angle)
            snp_y = distance * np.sin(snp_angle)

            # Classify as heterozygous or homozygous
            is_heterozygous = (str(row['allele1']) != str(row['allele2']) and
                             str(row['allele1']) != '--' and str(row['allele2']) != '--')

            if is_heterozygous:
                het_snps_x.append(snp_x)
                het_snps_y.append(snp_y)
            else:
                hom_snps_x.append(snp_x)
                hom_snps_y.append(snp_y)

    # Add heterozygous SNPs (red/orange dots)
    if het_snps_x:
        traces.append(go.Scatter(
            x=het_snps_x, y=het_snps_y,
            mode='markers',
            marker=dict(size=3, color='#E74C3C', opacity=0.7),
            name='Heterozygous',
            showlegend=False,  # Will add legend below plot
            hovertemplate='Heterozygous SNP<extra></extra>'
        ))

    # Add homozygous SNPs (green dots)
    if hom_snps_x:
        traces.append(go.Scatter(
            x=hom_snps_x, y=hom_snps_y,
            mode='markers',
            marker=dict(size=2, color='#27AE60', opacity=0.6),
            name='Homozygous',
            showlegend=False,  # Will add legend below plot
            hovertemplate='Homozygous SNP<extra></extra>'
        ))

    # Add connecting lines from center to chromosomes
    for chrom in chrom_list:
        if any(df_sample['chromosome'] == chrom):
            angle = chrom_angles[chrom]

            # Line from center to chromosome
            line_x = [0, 10 * np.cos(angle)]
            line_y = [0, 10 * np.sin(angle)]

            traces.append(go.Scatter(
                x=line_x, y=line_y,
                mode='lines',
                line=dict(color='#BDC3C7', width=1, dash='dot'),
                showlegend=False,
                hoverinfo='skip'
            ))

    return traces

def create_enhanced_quantum_iris(df, max_snps=3000, resolution=500):
    """Create enhanced iris-like visualization with rich detail - HIGH RESOLUTION"""

    print(f"Creating Enhanced Quantum Genome Iris from {min(len(df), max_snps)} SNPs...")

    if len(df) > max_snps:
        # Sample across all chromosomes proportionally
        df_sample = pd.concat([
            group.sample(n=max(1, min(len(group), max_snps // 24)))
            for _, group in df.groupby('chromosome')
        ]).reset_index(drop=True)
    else:
        df_sample = df

    # HIGH RESOLUTION coordinate grid
    x = np.linspace(-15, 15, resolution)
    y = np.linspace(-15, 15, resolution)
    X, Y = np.meshgrid(x, y)

    # Convert to polar coordinates for iris features
    R = np.sqrt(X**2 + Y**2)
    THETA = np.arctan2(Y, X)
    THETA[THETA < 0] += 2*np.pi  # Ensure theta is 0 to 2π

    # Initialize multiple wave fields for complexity
    primary_wave = np.zeros_like(X, dtype=complex)
    harmonic_wave = np.zeros_like(X, dtype=complex)
    radial_wave = np.zeros_like(X, dtype=complex)

    print(f"Generating high-resolution wave interference from {len(df_sample)} genetic variants...")

    # Create primary wave interference
    for i, (_, row) in enumerate(df_sample.iterrows()):
        frequency, amplitude, phase_shift = genotype_to_wave_properties(row['allele1'], row['allele2'])

        if amplitude > 0:  # Only process valid genotypes
            # Position wave source based on chromosome (circular arrangement)
            chrom_num = int(row['chromosome']) if row['chromosome'].isdigit() else (23 if row['chromosome'] == 'X' else 24)
            chrom_angle = (chrom_num - 1) * 2 * np.pi / 24

            # Radial position based on genomic position within chromosome
            normalized_pos = (row['position'] % 10000000) / 10000000  # Cycle every 10M bases
            radius = 8 + 4 * normalized_pos  # Between radius 8-12

            source_x = radius * np.cos(chrom_angle)
            source_y = radius * np.sin(chrom_angle)

            # Calculate distance from each grid point to wave source
            r = np.sqrt((X - source_x)**2 + (Y - source_y)**2)
            r = np.maximum(r, 0.1)  # Avoid division by zero

            # PRIMARY WAVE
            wave_contribution = amplitude * np.exp(1j * (frequency * r + phase_shift)) / r
            primary_wave += wave_contribution

            # HARMONIC WAVES: Add overtones for richer patterns
            harmonic_freq = frequency * 1.5  # Harmonic frequency
            harmonic_contribution = amplitude * 0.4 * np.exp(1j * (harmonic_freq * r + phase_shift * 1.5)) / (r * 1.2)
            harmonic_wave += harmonic_contribution

            # RADIAL IRIS FIBERS: Create radial patterns from each chromosome sector
            if str(row['allele1']) != str(row['allele2']):  # Heterozygous creates radial fibers
                # Calculate angle from this chromosome's sector
                angle_diff = THETA - chrom_angle
                angle_diff = np.where(angle_diff > np.pi, angle_diff - 2*np.pi, angle_diff)
                angle_diff = np.where(angle_diff < -np.pi, angle_diff + 2*np.pi, angle_diff)

                # Create radial fiber pattern
                radial_strength = amplitude * 0.6 * np.exp(-(angle_diff**2) / 0.3)  # Angular focus
                radial_modulation = np.cos(frequency * R * 2 + phase_shift)  # Radial waves
                radial_fiber = radial_strength * radial_modulation
                radial_wave += radial_fiber

    # Combine all wave components
    total_wave = primary_wave + harmonic_wave + radial_wave

    # Calculate enhanced wave properties
    intensity = np.abs(total_wave)**2
    phase_field = np.angle(total_wave)

    # Create iris anatomical features
    # 1. Pupil region (R < 2)
    pupil_mask = R < 2.0

    # 2. Iris region (2 <= R <= 10)
    iris_mask = (R >= 2.0) & (R <= 10.0)

    # 3. Limbal ring (R around 10)
    limbal_mask = (R >= 9.5) & (R <= 10.5)

    # Enhanced iris pattern creation
    iris_pattern = intensity.copy()

    # Apply pupil (black center)
    iris_pattern[pupil_mask] = 0

    # Apply iris region enhancements
    iris_pattern = iris_pattern * iris_mask  # Only show in iris region

    # Add collarette enhancement (textural boundary)
    collarette_enhancement = 0.5 * np.exp(-((R - 4.0)**2) / 0.5) * iris_mask
    iris_pattern += collarette_enhancement

    # Add radial fiber details using phase information
    radial_fibers = 0.3 * (1 + np.cos(8 * THETA + phase_field)) * iris_mask
    iris_pattern += radial_fibers

    # Add concentric rings (iris muscle structure)
    concentric_rings = 0.2 * np.sin(R * 2) * iris_mask
    iris_pattern += concentric_rings

    # Add limbal darkening
    limbal_darkening = 0.4 * np.exp(-((R - 10.0)**2) / 0.3) * limbal_mask
    iris_pattern += limbal_darkening

    # Add fine texture details
    fine_texture = 0.15 * np.sin(THETA * 20) * np.cos(R * 6) * iris_mask
    iris_pattern += fine_texture

    # Add crypts and furrows (heterozygous creates more texture)
    heterozygous_count = sum(1 for _, row in df_sample.iterrows()
                           if str(row['allele1']) != str(row['allele2']) and
                           str(row['allele1']) != '--' and str(row['allele2']) != '--')
    texture_complexity = heterozygous_count / len(df_sample)

    crypts_furrows = texture_complexity * 0.25 * np.random.normal(0, 1, iris_pattern.shape) * iris_mask
    iris_pattern += crypts_furrows

    # Final pupil masking
    iris_pattern[pupil_mask] = 0
    iris_pattern[R > 10.5] = 0  # Outside eye boundary

    # Normalize patterns
    if np.max(intensity) > 0:
        intensity = intensity / np.max(intensity)
    if np.max(iris_pattern) > 0:
        iris_pattern = iris_pattern / np.max(iris_pattern)

    print(f"   • High-res wave intensity range: {np.min(intensity):.4f} to {np.max(intensity):.4f}")
    print(f"   • Enhanced iris pattern range: {np.min(iris_pattern):.4f} to {np.max(iris_pattern):.4f}")
    print(f"   • Genetic texture complexity: {texture_complexity:.1%}")

    # Phase field for wave plot
    phase_normalized = (phase_field + np.pi) / (2 * np.pi)  # Normalize to 0-1
    phase_normalized[pupil_mask] = 0
    phase_normalized[R > 10.5] = 0

    return iris_pattern, phase_normalized, x, y, len(df_sample), intensity, texture_complexity

def analyze_genome_uniqueness(df):
    """Analyze what makes this genome's wave pattern unique"""

    # Count genotype patterns
    genotype_counts = {}
    wave_signatures = []

    for _, row in df.iterrows():
        genotype = f"{row['allele1']}{row['allele2']}"
        genotype_counts[genotype] = genotype_counts.get(genotype, 0) + 1

        freq, amp, phase = genotype_to_wave_properties(row['allele1'], row['allele2'])
        wave_signatures.append((freq, amp, phase))

    # Calculate uniqueness metrics
    total_snps = len(df)
    heterozygous_ratio = sum(1 for _, row in df.iterrows()
                           if str(row['allele1']) != str(row['allele2']) and
                           str(row['allele1']) != '--' and str(row['allele2']) != '--') / total_snps

    # Wave complexity (variety of frequencies)
    unique_frequencies = len(set(sig[0] for sig in wave_signatures if sig[1] > 0))

    # Chromosome distribution
    chrom_distribution = df['chromosome'].value_counts().to_dict()

    return {
        'total_snps': total_snps,
        'heterozygous_ratio': heterozygous_ratio,
        'genotype_counts': genotype_counts,
        'unique_frequencies': unique_frequencies,
        'chromosome_distribution': chrom_distribution,
        'wave_complexity_score': heterozygous_ratio * unique_frequencies
    }

# Usage example
if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Genome-2-Art: Transform genetic data into quantum wave art'
    )
    parser.add_argument('-i', '--input', default='AncestryDNA_mock.txt',
                        help='Input file path (default: AncestryDNA_mock.txt)')
    parser.add_argument('-p', '--palette', default='scientific',
                        choices=COLOR_PALETTES.keys(),
                        help='Color palette name (default: scientific)')
    parser.add_argument('--list-palettes', action='store_true',
                        help='Display available palettes and exit')
    args = parser.parse_args()

    # Handle --list-palettes
    if args.list_palettes:
        print("Available color palettes:")
        for name in COLOR_PALETTES:
            low, high = get_palette_legend_colors(name)
            print(f"  {name:12s}  low={low}  high={high}")
        sys.exit(0)

    # Load genome data
    df = load_ancestry_data(args.input)

    print(f"Loaded {len(df):,} SNPs from your genome")

    # Analyze uniqueness
    uniqueness = analyze_genome_uniqueness(df)
    print(f"\n🧬 Your Genome's Enhanced Signature:")
    print(f"   • Total genetic variants: {uniqueness['total_snps']:,}")
    print(f"   • Heterozygous ratio: {uniqueness['heterozygous_ratio']:.1%}")
    print(f"   • Wave complexity score: {uniqueness['wave_complexity_score']:.2f}")
    print(f"   • Unique wave frequencies: {uniqueness['unique_frequencies']}")

    # Create all four visualizations
    print(f"\n🔬 Creating comprehensive genome visualization suite...")

    # Get all plot data
    circos_traces = create_circos_plot(df)
    snp_traces = create_snp_distribution_plot(df)
    iris_pattern, phase_pattern, x, y, snps_used, intensity, texture_complexity = create_enhanced_quantum_iris(df, resolution=600)

    # Get selected colorscale from palette
    selected_colorscale = COLOR_PALETTES[args.palette]

    # Create the 2x2 figure with adjusted spacing - bottom plots moved down 3%
    fig = make_subplots(
        rows=2, cols=2,
        row_heights=[1, 1.2],
        subplot_titles=(
            "Genome CIRCOS Plot",
            "SNP Distribution",
            "Quantum Genome Iris",
            "Wave Phase Field"
        ),
        specs=[
            [{"type": "scatter"}, {"type": "scatter"}],
            [{"type": "heatmap"}, {"type": "heatmap"}]
        ],
        horizontal_spacing=0.02,  # Much tighter horizontal spacing
        vertical_spacing=0.11   # Increased by 3% to move bottom plots down
    )

    # Add CIRCOS plot traces (top left)
    for trace in circos_traces:
        fig.add_trace(trace, row=1, col=1)

    # Add SNP distribution plot traces (top right)
    for trace in snp_traces:
        fig.add_trace(trace, row=1, col=2)

    # Add iris pattern (bottom left) with scientific colorscale
    fig.add_trace(
        go.Heatmap(
            z=iris_pattern,
            x=x,
            y=y,
            colorscale=selected_colorscale,
            name="Quantum Iris",
            hovertemplate="Position: (%{x:.1f}, %{y:.1f})<br>Iris Pattern: %{z:.4f}<extra></extra>",
            showscale=False
        ),
        row=2, col=1
    )

    # Invert phase pattern but keep background areas black (mapped to 0)
    wave_display = 1.0 - phase_pattern
    wave_display[phase_pattern == 0] = 0

    # Add wave phase field (bottom right) with inverted colorscale for more contrast
    fig.add_trace(
        go.Heatmap(
            z=wave_display,
            x=x,
            y=y,
            colorscale=selected_colorscale,
            name="Wave Phase",
            hovertemplate="Position: (%{x:.1f}, %{y:.1f})<br>Wave Phase: %{z:.4f}<extra></extra>",
            showscale=False
        ),
        row=2, col=2
    )

    # Update layout with black background and improved text visibility
    fig.update_layout(
        title=dict(
            text="<b>Genome-2-Art: From Scientific Graphing to Artistic Quantum Waves</b>",
            x=0.5,
            font=dict(size=20, color='white'),  # White title for dark background
            pad=dict(t=20, b=20)
        ),
        paper_bgcolor='#000000',  # Black background
        plot_bgcolor='#000000',   # Black plot background
        font=dict(color='white', size=12),  # White text
        height=900,   # Slightly smaller for tighter layout
        width=1200,   # Slightly smaller for tighter layout
        showlegend=False,  # Remove side legends
        margin=dict(t=100, b=80, l=50, r=50)  # Adjusted margins
    )

    # Configure axes for each subplot with standard zoom
    axis_range = [-12, 12]  # Standard range for top plots
    bottom_range = [-15, 15]  # Standard range for bottom plots

    # CIRCOS plot (row 1, col 1)
    fig.update_xaxes(
        title_text="", showticklabels=False, showgrid=False, zeroline=False,
        showline=True, linecolor='lightgrey', linewidth=1,
        scaleanchor="y", scaleratio=1, range=axis_range,
        row=1, col=1
    )
    fig.update_yaxes(
        title_text="", showticklabels=False, showgrid=False, zeroline=False,
        showline=True, linecolor='lightgrey', linewidth=1,
        range=axis_range,
        row=1, col=1
    )

    # SNP distribution plot (row 1, col 2)
    fig.update_xaxes(
        title_text="", showticklabels=False, showgrid=False, zeroline=False,
        showline=True, linecolor='lightgrey', linewidth=1,
        scaleanchor="y2", scaleratio=1, range=axis_range,
        row=1, col=2
    )
    fig.update_yaxes(
        title_text="", showticklabels=False, showgrid=False, zeroline=False,
        showline=True, linecolor='lightgrey', linewidth=1,
        range=axis_range,
        row=1, col=2
    )

    # Iris plot (row 2, col 1) - 5% zoom increase
    fig.update_xaxes(
        title_text="", showticklabels=False, showgrid=False, zeroline=False,
        showline=True, linecolor='lightgrey', linewidth=1,
        scaleanchor="y3", scaleratio=1, range=bottom_range,
        row=2, col=1
    )
    fig.update_yaxes(
        title_text="", showticklabels=False, showgrid=False, zeroline=False,
        showline=True, linecolor='lightgrey', linewidth=1,
        range=bottom_range,
        row=2, col=1
    )

    # Wave plot (row 2, col 2) - 5% zoom increase
    fig.update_xaxes(
        title_text="", showticklabels=False, showgrid=False, zeroline=False,
        showline=True, linecolor='lightgrey', linewidth=1,
        scaleanchor="y4", scaleratio=1, range=bottom_range,
        row=2, col=2
    )
    fig.update_yaxes(
        title_text="", showticklabels=False, showgrid=False, zeroline=False,
        showline=True, linecolor='lightgrey', linewidth=1,
        range=bottom_range,
        row=2, col=2
    )

    # Add black background for Wave Phase Field heatmap - must be after axes configuration
    fig.add_shape(
        type="rect",
        xref="x domain", yref="y domain",
        x0=0, x1=1, y0=0, y1=1,
        fillcolor="black", line_width=0,
        layer="below",
        row=2, col=2
    )

    # Update subplot titles with better styling for dark background
    for annotation in fig['layout']['annotations']:
        annotation['font'] = dict(size=14, color='white')
        annotation['y'] = annotation['y'] + 0.02  # Move titles slightly up to prevent overlap

    # Add custom legends with improved formatting and positioning
    # Top legends appear above the bottom plot titles

    # Legend positions
    top_legend_above_bottom_titles = 0.568  # Moved up by 2% for optimal positioning
    bottom_legend_y = -0.08  # Below bottom plots

    # Legend for CIRCOS plot (top left) - positioned above Iris plot title with padding
    fig.add_annotation(
        text="<b>●</b> Chromosomes &nbsp;&nbsp;&nbsp; <b>●</b> Heterozygosity Rate",
        x=0.24, y=top_legend_above_bottom_titles,
        xref="paper", yref="paper",
        showarrow=False,
        font=dict(size=12, color='#E8E8E8', family="Arial"),
        xanchor="center",
        bgcolor="rgba(0,0,0,0.8)",
        bordercolor="#333333",
        borderwidth=1,
        borderpad=8  # Added padding for better spacing from plot titles
    )

    # Legend for SNP Distribution (top right) - positioned above Wave Phase plot title with padding
    fig.add_annotation(
        text="<b style='color:#E74C3C'>●</b> Heterozygous SNPs &nbsp;&nbsp;&nbsp; <b style='color:#27AE60'>●</b> Homozygous SNPs",
        x=0.76, y=top_legend_above_bottom_titles,
        xref="paper", yref="paper",
        showarrow=False,
        font=dict(size=12, color='#E8E8E8', family="Arial"),
        xanchor="center",
        bgcolor="rgba(0,0,0,0.8)",
        bordercolor="#333333",
        borderwidth=1,
        borderpad=8  # Added padding for better spacing from plot titles
    )

    # Get dynamic legend colors from selected palette
    legend_low_color, legend_high_color = get_palette_legend_colors(args.palette)

    # Legend for Iris (bottom left) - below the plot
    fig.add_annotation(
        text=f"<b style='color:{legend_low_color}'>Low Intensity</b> ← → <b style='color:{legend_high_color}'>High Intensity</b><br><i>{args.palette.capitalize()} Genome Pattern Visualization</i>",
        x=0.24, y=bottom_legend_y,
        xref="paper", yref="paper",
        showarrow=False,
        font=dict(size=11, color='#E8E8E8', family="Arial"),
        xanchor="center",
        bgcolor="rgba(0,0,0,0.8)",
        bordercolor="#333333",
        borderwidth=1
    )

    # Legend for Wave Phase (bottom right) - below the plot
    fig.add_annotation(
        text=f"<b style='color:{legend_low_color}'>Wave Minimum</b> ← → <b style='color:{legend_high_color}'>Wave Maximum</b><br><i>Quantum Wave Interference Phases</i>",
        x=0.76, y=bottom_legend_y,
        xref="paper", yref="paper",
        showarrow=False,
        font=dict(size=11, color='#E8E8E8', family="Arial"),
        xanchor="center",
        bgcolor="rgba(0,0,0,0.8)",
        bordercolor="#333333",
        borderwidth=1
    )

    print(f"   • CIRCOS: Chromosome organization with heterozygosity visualization")
    print(f"   • SNP Distribution: Individual variants positioned by chromosome")
    print(f"   • Iris: Enhanced anatomical pattern from {snps_used:,} SNPs")
    print(f"   • Wave: High-resolution quantum interference (600x600 grid)")

    # Save comprehensive suite first
    fig.write_html("genome-2-art.html", include_plotlyjs="cdn",
                   config={
        'displayModeBar': True,
        'displaylogo': False,
        'responsive': True
    })

    print(f"   • Clean HTML output without embedded documentation")

    # Export the full suite as high-res images
    os.makedirs("art", exist_ok=True)

    # Create standalone Wave Phase Field figure for export
    wave_fig = go.Figure()
    wave_fig.add_trace(
        go.Heatmap(
            z=wave_display,
            x=x,
            y=y,
            colorscale=selected_colorscale,
            showscale=False
        )
    )

    wave_fig.update_layout(
        title="Wave Phase Field",
        paper_bgcolor='#000000',
        plot_bgcolor='#000000',
        font=dict(color='white', size=12),
        height=800,
        width=800,
        margin=dict(t=50, b=50, l=50, r=50),
        showlegend=False
    )

    wave_fig.update_xaxes(
        showticklabels=False, showgrid=False, zeroline=False,
        showline=False, range=[-15, 15]
    )
    wave_fig.update_yaxes(
        showticklabels=False, showgrid=False, zeroline=False,
        showline=False, range=[-15, 15]
    )

    try:
        wave_fig.write_image("art/genome-2-art.svg", width=1200, height=1200)
        print(f"   • Wave Phase Field SVG: art/genome-2-art.svg (vector format)")
    except Exception as e:
        print(f"   • SVG export failed: {e}")

    try:
        wave_fig.write_image("art/genome-2-art.png", width=1200, height=1200)
        print(f"   • Wave Phase Field PNG: art/genome-2-art.png (high-res raster)")
    except Exception as e:
        print(f"   • PNG export failed: {e}")

    print(f"\n💾 Export Complete:")
    print(f"   • Full suite HTML: genome-2-art.html")
    print(f"   • Color palette: {args.palette}")
    print(f"\n📐 Improvements Made:")
    print(f"   • Black background with white text for better contrast")
    print(f"   • Legends moved below each plot to prevent overlap")
    print(f"   • Tighter spacing with uniform plot sizing")
    print(f"   • {args.palette.capitalize()} colorscale for genome visualization")
    print(f"   • Wave Phase plot inverted to show more cyan blue")
    print(f"   • Consistent axis ranges for better visual balance")
    print(f"   • All quantum interference patterns preserved at full resolution")
    print(f"   • Color scheme matches scientific standards for genomic data")
