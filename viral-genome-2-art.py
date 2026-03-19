import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio
import kaleido
import os
import glob

# ── Color palette matching genome-2-art ─────────────────────────────────────
SCIENTIFIC_COLORSCALE = [
    [0.0, '#000011'], [0.1, '#001133'], [0.2, '#003366'],
    [0.3, '#0066AA'], [0.4, '#00BFFF'], [0.5, '#7F5FFF'],
    [0.6, '#BF3FFF'], [0.7, '#FF00FF'], [0.8, '#FF66FF'],
    [0.9, '#FFAAFF'], [1.0, '#FFFFFF']
]

SAMPLE_COLORS = ['#E74C3C', '#3498DB', '#2ECC71', '#F39C12', '#9B59B6']

# Warm amber/gold — complementary to the cool blue-cyan-purple iris colorscale
PHASE_COLORSCALE = [
    [0.0,  '#000000'],
    [0.15, '#1A0800'],
    [0.35, '#7A2000'],
    [0.55, '#D45000'],
    [0.75, '#FF8C00'],
    [0.90, '#FFB800'],
    [1.0,  '#FFE066']
]


# ── Data loading ─────────────────────────────────────────────────────────────

def load_vcf_data(filepath):
    """Load SARS-CoV-2 VCF data file"""

    with open(filepath, 'r') as file:
        lines = file.readlines()

    header_line = None
    data_start = 0
    for i, line in enumerate(lines):
        if line.startswith('#CHROM'):
            header_line = line.strip().split('\t')
            data_start = i + 1
            break

    if header_line is None:
        raise ValueError("Could not find VCF header line")

    data = []
    for line in lines[data_start:]:
        if line.strip():
            data.append(line.strip().split('\t'))

    df = pd.DataFrame(data, columns=header_line)
    df['POS'] = pd.to_numeric(df['POS'])

    sample_col = df.columns[-1]
    df['GT'] = df[sample_col].str.split(':').str[0]

    df['REF_base'] = df['REF']
    df['ALT_base'] = df['ALT']

    df['allele1'] = df.apply(
        lambda row: row['REF_base'] if row['GT'].split('/')[0] == '0' else row['ALT_base'], axis=1
    )
    df['allele2'] = df.apply(
        lambda row: row['REF_base'] if row['GT'].split('/')[1] == '0' else row['ALT_base'], axis=1
    )

    df['chromosome'] = '1'
    df['position'] = df['POS']

    print(f"Loaded {len(df)} variants from {filepath}")
    return df[['chromosome', 'position', 'allele1', 'allele2', 'REF_base', 'ALT_base', 'GT']].copy()


def genotype_to_wave_properties_viral(allele1, allele2):
    """Convert viral genotype to wave properties"""

    nucleotide_wave_map = {
        'A': {'freq': 1.5, 'amp': 1.0},
        'T': {'freq': 2.8, 'amp': 0.9},
        'C': {'freq': 4.2, 'amp': 0.7},
        'G': {'freq': 5.5, 'amp': 0.5},
        'N': {'freq': 0.1, 'amp': 0.1},
        '-': {'freq': 0.0, 'amp': 0.0}
    }

    def get_wave(allele):
        if len(str(allele)) == 1:
            return nucleotide_wave_map.get(str(allele).upper(), {'freq': 0.5, 'amp': 0.5})
        else:
            length = len(str(allele))
            return {'freq': min(6.0, 1.0 + length * 0.5), 'amp': min(1.0, 0.3 + length * 0.1)}

    wave1 = get_wave(allele1)
    wave2 = get_wave(allele2)

    if allele1 == allele2:
        frequency = wave1['freq']
        amplitude = wave1['amp']
        phase_shift = 0
    else:
        frequency = (wave1['freq'] + wave2['freq']) / 2
        amplitude = (wave1['amp'] + wave2['amp']) / 2 * 1.3
        phase_shift = abs(wave1['freq'] - wave2['freq']) * np.pi / 3

    return frequency, amplitude, phase_shift


# ── Visualization functions (mirroring genome-2-art structure) ───────────────

def create_viral_circos_plot(df_list, sample_names, n_sectors=12):
    """CIRCOS-style plot: viral genome divided into position sectors, one arc per sector.
    Inner bars show per-sample heterozygosity — mirrors genome create_circos_plot."""

    genome_length = max(df['position'].max() for df in df_list)
    sector_size = genome_length / n_sectors
    arc_width = 2 * np.pi / n_sectors

    traces = []
    outer_r = 10
    inner_r = 8

    # Outer sector arcs (like chromosome arcs in genome version)
    for i in range(n_sectors):
        angle = i * arc_width
        arc_angles = np.linspace(angle - arc_width / 2 + 0.03, angle + arc_width / 2 - 0.03, 30)

        x_outer = outer_r * np.cos(arc_angles)
        y_outer = outer_r * np.sin(arc_angles)
        x_inner = inner_r * np.cos(arc_angles)
        y_inner = inner_r * np.sin(arc_angles)

        x_arc = np.concatenate([x_outer, x_inner[::-1], [x_outer[0]]])
        y_arc = np.concatenate([y_outer, y_inner[::-1], [y_outer[0]]])

        sector_start_kb = int(i * sector_size / 1000)
        sector_end_kb = int((i + 1) * sector_size / 1000)

        traces.append(go.Scatter(
            x=x_arc, y=y_arc,
            fill='toself',
            fillcolor=f'rgba({50 + i * 15}, {80 + i * 10}, {150 + i * 8}, 0.65)',
            line=dict(color='white', width=1),
            showlegend=False,
            hovertemplate=f'{sector_start_kb}kb–{sector_end_kb}kb<extra></extra>'
        ))

        # Sector position labels
        label_r = outer_r + 1.5
        traces.append(go.Scatter(
            x=[label_r * np.cos(angle)], y=[label_r * np.sin(angle)],
            mode='text',
            text=[f'{sector_start_kb}kb'],
            textfont=dict(size=10, color='white'),
            showlegend=False,
            hoverinfo='skip'
        ))

    # Inner bars: het rate per sample per sector (like genome het bars, stacked by sample)
    n_samples = len(df_list)
    bar_base_radii = [5.0 + s * 0.9 for s in range(n_samples)]
    bar_max_height = 0.75

    for s_idx, (df, sample_name) in enumerate(zip(df_list, sample_names)):
        bar_r = bar_base_radii[s_idx]
        color = SAMPLE_COLORS[s_idx % len(SAMPLE_COLORS)]

        for i in range(n_sectors):
            angle = i * arc_width
            sector_start = i * sector_size
            sector_end = (i + 1) * sector_size

            sector_data = df[(df['position'] >= sector_start) & (df['position'] < sector_end)]
            if len(sector_data) == 0:
                continue

            het_rate = sum(
                1 for _, row in sector_data.iterrows()
                if str(row['allele1']) != str(row['allele2'])
            ) / len(sector_data)

            bar_height = het_rate * bar_max_height
            bar_angles = np.linspace(angle - arc_width / 2 + 0.06, angle + arc_width / 2 - 0.06, 10)

            x_bar_outer = (bar_r + bar_height) * np.cos(bar_angles)
            y_bar_outer = (bar_r + bar_height) * np.sin(bar_angles)
            x_bar_inner = bar_r * np.cos(bar_angles)
            y_bar_inner = bar_r * np.sin(bar_angles)

            x_bar = np.concatenate([x_bar_outer, x_bar_inner[::-1], [x_bar_outer[0]]])
            y_bar = np.concatenate([y_bar_outer, y_bar_inner[::-1], [y_bar_outer[0]]])

            # Parse hex color to rgba
            r_c = int(color[1:3], 16)
            g_c = int(color[3:5], 16)
            b_c = int(color[5:7], 16)

            traces.append(go.Scatter(
                x=x_bar, y=y_bar,
                fill='toself',
                fillcolor=f'rgba({r_c},{g_c},{b_c},0.85)',
                line=dict(color='white', width=0.5),
                showlegend=False,
                hovertemplate=f'{sample_name}<br>{int(sector_start//1000)}–{int(sector_end//1000)}kb<br>Het: {het_rate:.0%}<extra></extra>'
            ))

    return traces


def create_viral_snp_distribution(df_list, sample_names):
    """Radial scatter of variants by genomic position — mirrors genome create_snp_distribution_plot."""

    genome_length = max(df['position'].max() for df in df_list)

    traces = []

    # Genome position markers every 5kb
    for pos in range(0, int(genome_length) + 1, 5000):
        angle = 2 * np.pi * pos / genome_length
        label_r = 10
        traces.append(go.Scatter(
            x=[label_r * np.cos(angle)], y=[label_r * np.sin(angle)],
            mode='markers+text',
            marker=dict(size=14, color='#2C3E50', symbol='circle'),
            text=[f'{pos // 1000}kb'],
            textposition='middle center',
            textfont=dict(size=9, color='white', family="Arial Black"),
            showlegend=False,
            hovertemplate=f'Position: {pos:,}<extra></extra>'
        ))
        # Spoke line
        traces.append(go.Scatter(
            x=[0, label_r * np.cos(angle)], y=[0, label_r * np.sin(angle)],
            mode='lines',
            line=dict(color='#BDC3C7', width=1, dash='dot'),
            showlegend=False, hoverinfo='skip'
        ))

    # Variant dots per sample
    het_colors = ['#E74C3C', '#E67E22', '#E91E63']
    hom_colors = ['#27AE60', '#2ECC71', '#00BCD4']

    for s_idx, (df, sample_name) in enumerate(zip(df_list, sample_names)):
        het_x, het_y, hom_x, hom_y = [], [], [], []

        for _, row in df.iterrows():
            base_angle = 2 * np.pi * row['position'] / genome_length
            spread = (row['position'] % 1000) / 1000 * 0.2 - 0.1
            snp_angle = base_angle + spread + s_idx * 0.06
            distance = 2 + (row['position'] % 25000) / 25000 * 6

            x = distance * np.cos(snp_angle)
            y = distance * np.sin(snp_angle)

            if str(row['allele1']) != str(row['allele2']):
                het_x.append(x); het_y.append(y)
            else:
                hom_x.append(x); hom_y.append(y)

        if het_x:
            traces.append(go.Scatter(
                x=het_x, y=het_y, mode='markers',
                marker=dict(size=6, color=het_colors[s_idx % 3], opacity=0.85),
                name=f'{sample_name} Het', showlegend=False,
                hovertemplate=f'{sample_name} Heterozygous<extra></extra>'
            ))
        if hom_x:
            traces.append(go.Scatter(
                x=hom_x, y=hom_y, mode='markers',
                marker=dict(size=5, color=hom_colors[s_idx % 3], opacity=0.75),
                name=f'{sample_name} Hom', showlegend=False,
                hovertemplate=f'{sample_name} Homozygous<extra></extra>'
            ))

    return traces


def create_viral_quantum_iris(df_list, sample_names, resolution=400):
    """Create quantum iris + wave phase field from viral genome wave interference.
    Uses the same iris anatomy as genome-2-art create_quantum_genome_iris,
    adapted for viral: genome position → angle, sample index → source radius."""

    n_total = sum(len(df) for df in df_list)
    print(f"Creating Viral Quantum Iris from {n_total} variants across {len(df_list)} samples...")

    genome_length = max(df['position'].max() for df in df_list)

    # HIGH RESOLUTION coordinate grid (same as genome version)
    x = np.linspace(-15, 15, resolution)
    y = np.linspace(-15, 15, resolution)
    X, Y = np.meshgrid(x, y)

    R = np.sqrt(X**2 + Y**2)
    THETA = np.arctan2(Y, X)
    THETA[THETA < 0] += 2 * np.pi

    # Three wave components (same as genome version)
    primary_wave  = np.zeros_like(X, dtype=complex)
    harmonic_wave = np.zeros_like(X, dtype=complex)
    radial_wave   = np.zeros_like(X, dtype=complex)

    for sample_idx, (df, sample_name) in enumerate(zip(df_list, sample_names)):
        print(f"  Processing {sample_name}...")
        sample_phase = sample_idx * np.pi / 4

        for _, row in df.iterrows():
            frequency, amplitude, phase_shift = genotype_to_wave_properties_viral(
                row['allele1'], row['allele2']
            )

            if amplitude > 0:
                # Genome position → angular placement (like chromosome angle in genome version)
                genome_angle = 2 * np.pi * row['position'] / genome_length
                # Radial position: sample index sets base radius, position modulates within it
                normalized_pos = (row['position'] % 5000000) / 5000000
                radius = (8 + sample_idx * 1.5) + 2 * normalized_pos

                source_x = radius * np.cos(genome_angle)
                source_y = radius * np.sin(genome_angle)

                r = np.sqrt((X - source_x)**2 + (Y - source_y)**2)
                r = np.maximum(r, 0.1)

                # PRIMARY WAVE
                primary_wave += amplitude * np.exp(1j * (frequency * r + phase_shift + sample_phase)) / r

                # HARMONIC WAVE
                harmonic_freq = frequency * 1.5
                harmonic_wave += (amplitude * 0.4
                                  * np.exp(1j * (harmonic_freq * r + phase_shift * 1.5))
                                  / (r * 1.2))

                # RADIAL IRIS FIBERS (heterozygous variants only)
                if str(row['allele1']) != str(row['allele2']):
                    angle_diff = THETA - genome_angle
                    angle_diff = np.where(angle_diff >  np.pi, angle_diff - 2 * np.pi, angle_diff)
                    angle_diff = np.where(angle_diff < -np.pi, angle_diff + 2 * np.pi, angle_diff)

                    radial_strength    = amplitude * 0.6 * np.exp(-(angle_diff**2) / 0.3)
                    radial_modulation  = np.cos(frequency * R * 2 + phase_shift)
                    radial_wave       += radial_strength * radial_modulation

    total_wave  = primary_wave + harmonic_wave + radial_wave
    intensity   = np.abs(total_wave)**2
    phase_field = np.angle(total_wave)

    # ── Iris anatomical features (identical to genome version) ───────────────
    pupil_mask  = R < 2.0
    iris_mask   = (R >= 2.0) & (R <= 10.0)
    limbal_mask = (R >= 9.5) & (R <= 10.5)

    iris_pattern = intensity.copy()
    iris_pattern[pupil_mask] = 0
    iris_pattern = iris_pattern * iris_mask

    # Collarette ring
    iris_pattern += 0.5 * np.exp(-((R - 4.0)**2) / 0.5) * iris_mask

    # Radial fibers
    iris_pattern += 0.3 * (1 + np.cos(8 * THETA + phase_field)) * iris_mask

    # Concentric rings (iris muscle structure)
    iris_pattern += 0.2 * np.sin(R * 2) * iris_mask

    # Limbal darkening
    iris_pattern += 0.4 * np.exp(-((R - 10.0)**2) / 0.3) * limbal_mask

    # Fine angular texture
    iris_pattern += 0.15 * np.sin(THETA * 20) * np.cos(R * 6) * iris_mask

    # Crypts and furrows driven by het rate
    all_rows = pd.concat(df_list)
    texture_complexity = sum(
        1 for _, row in all_rows.iterrows()
        if str(row['allele1']) != str(row['allele2'])
    ) / len(all_rows)

    iris_pattern += texture_complexity * 0.25 * np.random.normal(0, 1, iris_pattern.shape) * iris_mask

    # Final masking
    iris_pattern[pupil_mask] = 0
    iris_pattern[R > 10.5] = 0

    # Normalize
    if np.max(intensity) > 0:
        intensity = intensity / np.max(intensity)
    if np.max(iris_pattern) > 0:
        iris_pattern = iris_pattern / np.max(iris_pattern)

    print(f"   • Iris intensity range: {np.min(intensity):.4f} to {np.max(intensity):.4f}")
    print(f"   • Genetic texture complexity: {texture_complexity:.1%}")

    # Phase field: weight by amplitude so zero-amplitude regions stay at exactly 0
    # (maps to black on PHASE_COLORSCALE). All non-iris pixels are explicitly zeroed.
    phase_normalized = (phase_field + np.pi) / (2 * np.pi) * intensity
    phase_normalized[~iris_mask] = 0
    phase_normalized[pupil_mask] = 0
    phase_normalized[R > 10.5] = 0
    if np.max(phase_normalized) > 0:
        phase_normalized = phase_normalized / np.max(phase_normalized)

    return iris_pattern, phase_normalized, x, y


# ── Main ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":

    vcf_files = glob.glob("public_genomes/*.vcf")

    if not vcf_files:
        print("No VCF files found in public_genomes directory!")
        exit(1)

    print(f"Found {len(vcf_files)} VCF files to process...")

    df_list, sample_names = [], []

    for vcf_file in vcf_files:
        try:
            df = load_vcf_data(vcf_file)
            df_list.append(df)
            sample_name = os.path.basename(vcf_file).split('.')[0]
            sample_names.append(sample_name)
            print(f"   • {sample_name}: {len(df)} variants")
        except Exception as e:
            print(f"   ✗ Failed to load {vcf_file}: {e}")

    if not df_list:
        print("No valid VCF files could be loaded!")
        exit(1)

    print(f"\n🦠 Creating SARS-CoV-2 Genome Art from {len(df_list)} samples...")
    print(f"\n🎨 Generating visualizations...")

    # ── Build the four panels ─────────────────────────────────────────────────
    circos_traces = create_viral_circos_plot(df_list, sample_names)
    snp_traces    = create_viral_snp_distribution(df_list, sample_names)
    iris_pattern, phase_normalized, ix, iy = create_viral_quantum_iris(df_list, sample_names)

    # ── Assemble 2×2 subplot figure (identical structure to genome-2-art) ─────
    fig = make_subplots(
        rows=2, cols=2,
        row_heights=[1, 1.2],
        subplot_titles=(
            "Viral Genome CIRCOS Plot",
            "Variant Distribution",
            "Quantum Viral Iris",
            "Wave Phase Field"
        ),
        specs=[
            [{"type": "scatter"}, {"type": "scatter"}],
            [{"type": "heatmap"}, {"type": "heatmap"}]
        ],
        horizontal_spacing=0.02,
        vertical_spacing=0.11
    )

    for trace in circos_traces:
        fig.add_trace(trace, row=1, col=1)

    for trace in snp_traces:
        fig.add_trace(trace, row=1, col=2)

    fig.add_trace(
        go.Heatmap(
            z=iris_pattern, x=ix, y=iy,
            colorscale=SCIENTIFIC_COLORSCALE,
            showscale=False,
            name="Quantum Viral Iris",
            hovertemplate="Position: (%{x:.1f}, %{y:.1f})<br>Iris: %{z:.4f}<extra></extra>"
        ),
        row=2, col=1
    )

    # Show phase directly (no inversion): zero amplitude stays dark,
    # wave peaks glow warm amber — complementary to the cool blue iris.
    fig.add_trace(
        go.Heatmap(
            z=phase_normalized, x=ix, y=iy,
            colorscale=PHASE_COLORSCALE,
            zmin=0, zmax=1,
            showscale=False,
            name="Wave Phase Field",
            hovertemplate="Position: (%{x:.1f}, %{y:.1f})<br>Phase: %{z:.4f}<extra></extra>"
        ),
        row=2, col=2
    )

    # ── Layout (matching genome-2-art) ────────────────────────────────────────
    fig.update_layout(
        title=dict(
            text="<b>Viral Genome-2-Art: SARS-CoV-2 Quantum Wave Visualization</b>",
            x=0.5,
            font=dict(size=20, color='white'),
            pad=dict(t=20, b=20)
        ),
        paper_bgcolor='#000000',
        plot_bgcolor='#000000',
        font=dict(color='white', size=12),
        height=900,
        width=1200,
        showlegend=False,
        margin=dict(t=100, b=80, l=50, r=50)
    )

    # ── Axis configuration (mirrors genome-2-art exactly) ─────────────────────
    axis_range  = [-12, 12]
    bottom_range = [-15, 15]

    # CIRCOS (row 1, col 1)
    fig.update_xaxes(
        title_text="", showticklabels=False, showgrid=False, zeroline=False,
        showline=True, linecolor='lightgrey', linewidth=1,
        scaleanchor="y", scaleratio=1, range=axis_range,
        row=1, col=1
    )
    fig.update_yaxes(
        title_text="", showticklabels=False, showgrid=False, zeroline=False,
        showline=True, linecolor='lightgrey', linewidth=1,
        range=axis_range, row=1, col=1
    )

    # Variant distribution (row 1, col 2)
    fig.update_xaxes(
        title_text="", showticklabels=False, showgrid=False, zeroline=False,
        showline=True, linecolor='lightgrey', linewidth=1,
        scaleanchor="y2", scaleratio=1, range=axis_range,
        row=1, col=2
    )
    fig.update_yaxes(
        title_text="", showticklabels=False, showgrid=False, zeroline=False,
        showline=True, linecolor='lightgrey', linewidth=1,
        range=axis_range, row=1, col=2
    )

    # Iris heatmap (row 2, col 1)
    fig.update_xaxes(
        title_text="", showticklabels=False, showgrid=False, zeroline=False,
        showline=True, linecolor='lightgrey', linewidth=1,
        scaleanchor="y3", scaleratio=1, range=bottom_range,
        row=2, col=1
    )
    fig.update_yaxes(
        title_text="", showticklabels=False, showgrid=False, zeroline=False,
        showline=True, linecolor='lightgrey', linewidth=1,
        range=bottom_range, row=2, col=1
    )

    # Phase field heatmap (row 2, col 2)
    fig.update_xaxes(
        title_text="", showticklabels=False, showgrid=False, zeroline=False,
        showline=True, linecolor='lightgrey', linewidth=1,
        scaleanchor="y4", scaleratio=1, range=bottom_range,
        row=2, col=2
    )
    fig.update_yaxes(
        title_text="", showticklabels=False, showgrid=False, zeroline=False,
        showline=True, linecolor='lightgrey', linewidth=1,
        range=bottom_range, row=2, col=2
    )

    # ── Save ──────────────────────────────────────────────────────────────────
    output_file = "viral-genome-2-art.html"
    fig.write_html(output_file, include_plotlyjs="cdn")

    print(f"\n💾 Viral Genome Art Complete!")
    print(f"   • HTML output: {output_file}")
    print(f"   • Samples analyzed: {len(df_list)}")
    print(f"   • Total variants: {sum(len(df) for df in df_list)}")

    # Export standalone art images
    os.makedirs("art", exist_ok=True)

    iris_fig = go.Figure()
    iris_fig.add_trace(go.Heatmap(z=iris_pattern, x=ix, y=iy,
                                   colorscale=SCIENTIFIC_COLORSCALE, showscale=False))
    iris_fig.update_layout(
        title="SARS-CoV-2 Quantum Viral Iris",
        paper_bgcolor='black', plot_bgcolor='black',
        height=800, width=800,
        margin=dict(t=50, b=50, l=50, r=50),
        font=dict(color='white')
    )
    iris_fig.update_xaxes(showticklabels=False, showgrid=False)
    iris_fig.update_yaxes(showticklabels=False, showgrid=False)

    try:
        iris_fig.write_image("art/viral-genome-iris.png", width=1200, height=1200)
        print(f"   • Iris art: art/viral-genome-iris.png")
    except Exception as e:
        print(f"   • PNG export failed: {e}")

    try:
        iris_fig.write_image("art/viral-genome-iris.svg", width=1200, height=1200)
        print(f"   • Iris vector: art/viral-genome-iris.svg")
    except Exception as e:
        print(f"   • SVG export failed: {e}")

    print(f"\n🧬 Analysis Summary:")
    for df, name in zip(df_list, sample_names):
        het_ratio = sum(1 for _, row in df.iterrows()
                        if str(row['allele1']) != str(row['allele2'])) / len(df)
        genome_span = df['position'].max() - df['position'].min()
        print(f"   • {name}: {len(df)} variants, {het_ratio:.1%} heterozygous, span {genome_span:,} bp")
