import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio
import kaleido
import os
import glob

def load_vcf_data(filepath):
    """Load SARS-CoV-2 VCF data file"""

    # Read VCF file, skipping header lines that start with ##
    with open(filepath, 'r') as file:
        lines = file.readlines()

    # Find the header line (starts with #CHROM)
    header_line = None
    data_start = 0
    for i, line in enumerate(lines):
        if line.startswith('#CHROM'):
            header_line = line.strip().split('\t')
            data_start = i + 1
            break

    if header_line is None:
        raise ValueError("Could not find VCF header line")

    # Parse data lines
    data = []
    for line in lines[data_start:]:
        if line.strip():
            data.append(line.strip().split('\t'))

    df = pd.DataFrame(data, columns=header_line)

    # Clean and convert data
    df['POS'] = pd.to_numeric(df['POS'])

    # Extract genotype information (GT field)
    # Split the last column (sample data) to get genotype
    sample_col = df.columns[-1]  # Last column is sample name
    df['GT'] = df[sample_col].str.split(':').str[0]  # First field is genotype

    # Parse REF and ALT alleles
    df['REF_base'] = df['REF']
    df['ALT_base'] = df['ALT']

    # Convert genotype to alleles
    df['allele1'] = df.apply(lambda row: row['REF_base'] if row['GT'].split('/')[0] == '0' else row['ALT_base'], axis=1)
    df['allele2'] = df.apply(lambda row: row['REF_base'] if row['GT'].split('/')[1] == '0' else row['ALT_base'], axis=1)

    # Add chromosome info (all SARS-CoV-2 is one chromosome)
    df['chromosome'] = '1'
    df['position'] = df['POS']

    print(f"Loaded {len(df)} variants from {filepath}")
    return df[['chromosome', 'position', 'allele1', 'allele2', 'REF_base', 'ALT_base', 'GT']].copy()

def genotype_to_wave_properties_viral(allele1, allele2):
    """Convert viral genotype to wave properties - optimized for viral mutations"""

    # Enhanced nucleotide mapping for viral genomes with more distinct patterns
    nucleotide_wave_map = {
        'A': {'freq': 1.5, 'amp': 1.0},      # Adenine - Low frequency, high amplitude
        'T': {'freq': 2.8, 'amp': 0.9},      # Thymine - Medium frequency, high amplitude
        'C': {'freq': 4.2, 'amp': 0.7},      # Cytosine - High frequency, medium amplitude
        'G': {'freq': 5.5, 'amp': 0.5},      # Guanine - Very high frequency, lower amplitude
        'N': {'freq': 0.1, 'amp': 0.1},      # Unknown - very low
        '-': {'freq': 0.0, 'amp': 0.0}       # Deletion
    }

    # Handle multi-character alleles (indels)
    def get_wave_for_allele(allele):
        if len(str(allele)) == 1:
            return nucleotide_wave_map.get(str(allele).upper(), {'freq': 0.5, 'amp': 0.5})
        else:
            # For indels, use length-based properties
            length = len(str(allele))
            return {'freq': min(6.0, 1.0 + length * 0.5), 'amp': min(1.0, 0.3 + length * 0.1)}

    wave1 = get_wave_for_allele(allele1)
    wave2 = get_wave_for_allele(allele2)

    # Create wave interference patterns
    if allele1 == allele2:  # Homozygous - pure wave
        frequency = wave1['freq']
        amplitude = wave1['amp']
        phase_shift = 0
    else:  # Heterozygous - creates interference patterns
        frequency = (wave1['freq'] + wave2['freq']) / 2
        amplitude = (wave1['amp'] + wave2['amp']) / 2 * 1.3  # Boost heterozygous visibility
        phase_shift = abs(wave1['freq'] - wave2['freq']) * np.pi / 3

    return frequency, amplitude, phase_shift

def create_viral_genome_circle(df, sample_name):
    """Create circular visualization of viral genome variants"""

    print(f"Creating viral genome circle for {sample_name}...")

    # SARS-CoV-2 genome is ~29,903 bases
    genome_length = df['position'].max()

    # Create circular positions
    angles = 2 * np.pi * df['position'] / genome_length

    traces = []

    # Base circle (genome backbone)
    circle_angles = np.linspace(0, 2*np.pi, 500)
    circle_r = 10
    circle_x = circle_r * np.cos(circle_angles)
    circle_y = circle_r * np.sin(circle_angles)

    traces.append(go.Scatter(
        x=circle_x, y=circle_y,
        mode='lines',
        line=dict(color='#34495E', width=3),
        name='Viral Genome',
        showlegend=False,
        hoverinfo='skip'
    ))

    # Add genome position markers every 5kb
    for pos in range(0, int(genome_length), 5000):
        angle = 2 * np.pi * pos / genome_length
        label_r = 12
        marker_x = label_r * np.cos(angle)
        marker_y = label_r * np.sin(angle)

        traces.append(go.Scatter(
            x=[marker_x], y=[marker_y],
            mode='markers+text',
            marker=dict(size=8, color='#7F8C8D'),
            text=[f'{pos//1000}kb'],
            textposition='middle center',
            textfont=dict(size=10, color='white'),
            showlegend=False,
            hovertemplate=f'Position: {pos:,}<extra></extra>'
        ))

    # Plot variants as radial spikes
    variant_x = []
    variant_y = []
    colors = []
    hover_texts = []

    for _, row in df.iterrows():
        angle = 2 * np.pi * row['position'] / genome_length

        # Get wave properties for visual encoding
        freq, amp, phase = genotype_to_wave_properties_viral(row['allele1'], row['allele2'])

        # Variant spike length based on impact/amplitude
        spike_length = 3 + amp * 4  # 3-7 units

        # Position on circle
        base_x = circle_r * np.cos(angle)
        base_y = circle_r * np.sin(angle)

        # Spike end position
        spike_x = (circle_r + spike_length) * np.cos(angle)
        spike_y = (circle_r + spike_length) * np.sin(angle)

        # Add spike line with proper color clamping
        r_val = max(0, min(255, int(255*freq/6)))
        g_val = max(0, min(255, int(255*(1-amp))))
        b_val = max(0, min(255, int(255*amp)))

        traces.append(go.Scatter(
            x=[base_x, spike_x], y=[base_y, spike_y],
            mode='lines',
            line=dict(
                color=f'rgba({r_val}, {g_val}, {b_val}, 0.8)',
                width=2
            ),
            showlegend=False,
            hovertemplate=f'Position: {row["position"]:,}<br>' +
                         f'Variant: {row["REF_base"]}→{row["ALT_base"]}<br>' +
                         f'Genotype: {row["GT"]}<extra></extra>'
        ))

    # Add sample name in center
    traces.append(go.Scatter(
        x=[0], y=[0],
        mode='text',
        text=[sample_name],
        textfont=dict(size=16, color='white', family="Arial Black"),
        showlegend=False,
        hoverinfo='skip'
    ))

    return traces

def create_viral_wave_field(df_list, sample_names, resolution=400):
    """Create wave interference field from multiple viral samples"""

    print(f"Creating viral wave interference field from {len(df_list)} samples...")

    # Create coordinate grid
    x = np.linspace(-12, 12, resolution)
    y = np.linspace(-12, 12, resolution)
    X, Y = np.meshgrid(x, y)

    # Convert to polar coordinates
    R = np.sqrt(X**2 + Y**2)
    THETA = np.arctan2(Y, X)
    THETA[THETA < 0] += 2*np.pi

    # Initialize wave field
    wave_field = np.zeros_like(X, dtype=complex)

    # Combine waves from all samples
    for sample_idx, (df, sample_name) in enumerate(zip(df_list, sample_names)):
        print(f"  Processing {sample_name}...")

        genome_length = df['position'].max()

        for _, row in df.iterrows():
            freq, amp, phase = genotype_to_wave_properties_viral(row['allele1'], row['allele2'])

            if amp > 0:
                # Position wave source on circle based on genomic position
                angle = 2 * np.pi * row['position'] / genome_length
                source_radius = 8 + sample_idx * 1.5  # Different radii for different samples

                source_x = source_radius * np.cos(angle)
                source_y = source_radius * np.sin(angle)

                # Calculate wave contribution
                r = np.sqrt((X - source_x)**2 + (Y - source_y)**2)
                r = np.maximum(r, 0.1)

                # Add wave with sample-specific phase offset
                sample_phase = sample_idx * np.pi / 4
                wave_contribution = amp * np.exp(1j * (freq * r + phase + sample_phase)) / r
                wave_field += wave_contribution

    # Calculate intensity
    intensity = np.abs(wave_field)**2

    # Normalize
    if np.max(intensity) > 0:
        intensity = intensity / np.max(intensity)

    # Create mask for visualization
    mask = R < 11
    intensity = intensity * mask

    return intensity, x, y

def create_variant_comparison_plot(df_list, sample_names):
    """Create comparison plot of variant positions across samples"""

    print(f"Creating variant comparison plot...")

    fig = go.Figure()

    colors = ['#E74C3C', '#3498DB', '#2ECC71', '#F39C12', '#9B59B6', '#1ABC9C']

    for i, (df, sample_name) in enumerate(zip(df_list, sample_names)):

        # Create y-values (one row per sample)
        y_values = [i] * len(df)

        # Color code by variant type
        variant_colors = []
        for _, row in df.iterrows():
            if len(row['REF_base']) == 1 and len(row['ALT_base']) == 1:
                variant_colors.append(colors[i % len(colors)])  # SNV
            else:
                variant_colors.append('#E67E22')  # Indel

        fig.add_trace(go.Scatter(
            x=df['position'],
            y=y_values,
            mode='markers',
            marker=dict(
                size=8,
                color=variant_colors,
                opacity=0.8,
                line=dict(width=1, color='white')
            ),
            name=sample_name,
            hovertemplate=f'{sample_name}<br>' +
                         'Position: %{x:,}<br>' +
                         '<extra></extra>'
        ))

    fig.update_layout(
        title="Viral Genome Variant Positions Across Samples",
        xaxis_title="Genomic Position (bp)",
        yaxis_title="Sample",
        yaxis=dict(
            tickmode='array',
            tickvals=list(range(len(sample_names))),
            ticktext=sample_names
        ),
        plot_bgcolor='black',
        paper_bgcolor='black',
        font=dict(color='white'),
        height=400
    )

    return fig

# Main execution
if __name__ == "__main__":

    # Find all VCF files in public_genomes directory
    vcf_files = glob.glob("public_genomes/*.vcf")

    if not vcf_files:
        print("No VCF files found in public_genomes directory!")
        exit(1)

    print(f"Found {len(vcf_files)} VCF files to process...")

    # Load all VCF data
    df_list = []
    sample_names = []

    for vcf_file in vcf_files:
        try:
            df = load_vcf_data(vcf_file)
            df_list.append(df)
            sample_name = os.path.basename(vcf_file).split('.')[0]  # Extract sample name
            sample_names.append(sample_name)

            print(f"   • {sample_name}: {len(df)} variants")

        except Exception as e:
            print(f"   ✗ Failed to load {vcf_file}: {e}")

    if not df_list:
        print("No valid VCF files could be loaded!")
        exit(1)

    print(f"\n🦠 Creating SARS-CoV-2 Genome Art from {len(df_list)} samples...")

    # Create visualizations
    print(f"\n🎨 Generating visualizations...")

    # Create the main subplot figure
    if len(df_list) == 1:
        # Single sample layout
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=(
                f"Viral Genome Circle - {sample_names[0]}",
                "Variant Positions",
                "Quantum Wave Field",
                "Variant Statistics"
            ),
            specs=[
                [{"type": "scatter"}, {"type": "scatter"}],
                [{"type": "heatmap"}, {"type": "scatter"}]
            ]
        )

        # Single sample circle
        circle_traces = create_viral_genome_circle(df_list[0], sample_names[0])
        for trace in circle_traces:
            fig.add_trace(trace, row=1, col=1)

    else:
        # Multi-sample layout
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=(
                "Sample Comparison",
                "Variant Position Comparison",
                "Combined Quantum Wave Field",
                "Statistics Summary"
            ),
            specs=[
                [{"type": "scatter"}, {"type": "scatter"}],
                [{"type": "heatmap"}, {"type": "table"}]
            ]
        )

        # Multi-sample circles (first 3 samples)
        display_samples = min(3, len(df_list))
        for i in range(display_samples):
            circle_traces = create_viral_genome_circle(df_list[i], sample_names[i])
            for trace in circle_traces:
                # Offset circles horizontally
                for j, point in enumerate(trace['x']):
                    if trace['x'] is not None:
                        trace['x'] = [x + (i-1)*25 for x in trace['x']]
                fig.add_trace(trace, row=1, col=1)

    # Variant comparison plot (row 1, col 2)
    comparison_fig = create_variant_comparison_plot(df_list, sample_names)
    for trace in comparison_fig.data:
        fig.add_trace(trace, row=1, col=2)

    # Wave field (row 2, col 1)
    wave_field, wx, wy = create_viral_wave_field(df_list, sample_names)

    scientific_colorscale = [
        [0.0, '#000011'],   [0.1, '#001133'],   [0.2, '#003366'],
        [0.3, '#0066AA'],   [0.4, '#00BFFF'],   [0.5, '#7F5FFF'],
        [0.6, '#BF3FFF'],   [0.7, '#FF00FF'],   [0.8, '#FF66FF'],
        [0.9, '#FFAAFF'],   [1.0, '#FFFFFF']
    ]

    fig.add_trace(
        go.Heatmap(
            z=wave_field,
            x=wx, y=wy,
            colorscale=scientific_colorscale,
            showscale=False,
            hovertemplate="Wave Intensity: %{z:.3f}<extra></extra>"
        ),
        row=2, col=1
    )

    # Statistics summary (row 2, col 2)
    stats_data = []
    for df, name in zip(df_list, sample_names):
        hetero_count = sum(1 for _, row in df.iterrows() if row['allele1'] != row['allele2'])
        stats_data.append([name, len(df), hetero_count, f"{100*hetero_count/len(df):.1f}%"])

    stats_df = pd.DataFrame(stats_data, columns=['Sample', 'Variants', 'Heterozygous', 'Het %'])

    fig.add_trace(
        go.Table(
            header=dict(values=list(stats_df.columns),
                       fill_color='#34495E',
                       font=dict(color='white', size=12)),
            cells=dict(values=[stats_df[col] for col in stats_df.columns],
                      fill_color='#2C3E50',
                      font=dict(color='white', size=11)),
        ),
        row=2, col=2
    )

    # Update layout
    fig.update_layout(
        title=dict(
            text="<b>Viral Genome-2-Art: SARS-CoV-2 Quantum Wave Visualization</b>",
            x=0.5,
            font=dict(size=18, color='white')
        ),
        paper_bgcolor='#000000',
        plot_bgcolor='#000000',
        font=dict(color='white', size=12),
        height=800,
        width=1200,
        showlegend=False
    )

    # Configure axes
    for row in [1, 2]:
        for col in [1, 2]:
            if row == 1:  # Top row scatter plots
                fig.update_xaxes(
                    showgrid=False, zeroline=False, showticklabels=False,
                    scaleanchor=f"y{'' if row==1 and col==1 else row*2+col-2}", scaleratio=1,
                    row=row, col=col
                )
                fig.update_yaxes(
                    showgrid=False, zeroline=False, showticklabels=False,
                    row=row, col=col
                )
            elif row == 2 and col == 1:  # Wave field
                fig.update_xaxes(showticklabels=False, showgrid=False, row=row, col=col)
                fig.update_yaxes(showticklabels=False, showgrid=False, row=row, col=col)

    # Save the visualization
    output_file = "viral-genome-2-art.html"
    fig.write_html(output_file, include_plotlyjs="cdn")

    print(f"\n💾 Viral Genome Art Complete!")
    print(f"   • HTML output: {output_file}")
    print(f"   • Samples analyzed: {len(df_list)}")
    print(f"   • Total variants: {sum(len(df) for df in df_list)}")

    # Create art directory and export wave field
    os.makedirs("art", exist_ok=True)

    # Export wave field as standalone art
    wave_fig = go.Figure()
    wave_fig.add_trace(
        go.Heatmap(
            z=wave_field,
            x=wx, y=wy,
            colorscale=scientific_colorscale,
            showscale=False
        )
    )

    wave_fig.update_layout(
        title="SARS-CoV-2 Quantum Wave Field",
        paper_bgcolor='black',
        plot_bgcolor='black',
        height=800, width=800,
        margin=dict(t=50, b=50, l=50, r=50),
        font=dict(color='white')
    )

    wave_fig.update_xaxes(showticklabels=False, showgrid=False)
    wave_fig.update_yaxes(showticklabels=False, showgrid=False)

    try:
        wave_fig.write_image("art/viral-genome-wave.png", width=1200, height=1200)
        print(f"   • Wave field art: art/viral-genome-wave.png")
    except Exception as e:
        print(f"   • PNG export failed: {e}")

    try:
        wave_fig.write_image("art/viral-genome-wave.svg", width=1200, height=1200)
        print(f"   • Wave field vector: art/viral-genome-wave.svg")
    except Exception as e:
        print(f"   • SVG export failed: {e}")

    print(f"\n🧬 Analysis Summary:")
    for df, name in zip(df_list, sample_names):
        hetero_ratio = sum(1 for _, row in df.iterrows() if row['allele1'] != row['allele2']) / len(df)
        genome_span = df['position'].max() - df['position'].min()
        print(f"   • {name}: {len(df)} variants, {hetero_ratio:.1%} heterozygous, span {genome_span:,} bp")
