# Changelog

All notable changes to the Genome-2-Art project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.1.0] - 2025-02-07

### Added

- **Enhanced Color Palette System**: 7 scientifically-inspired color palettes for quantum wave visualization
  - `scientific`: Cyan-magenta gradient (default)
  - `ocean`: Deep ocean blues to white
  - `fire`: Black to red to golden yellow
  - `nebula`: Deep space purples to light pink
  - `earth`: Dark earth tones to cream
  - `aurora`: Night sky teals to bright greens
  - `infrared`: Deep purples to pale lavender
- **Command-line Interface** via argparse:
  - `-i, --input` flag to specify input file (default: `AncestryDNA_mock.txt`)
  - `-p, --palette` flag to select color schemes (default: `scientific`)
  - `--list-palettes` option to display available palettes with representative colors
- **Dynamic Legend System**: Iris and Wave Phase Field legends automatically adapt colors and label text to match the selected palette
- **High-Resolution Export**: Black background support for both PNG and SVG exports across all palettes

### Fixed

- **Wave Phase Field Background Bug**: Inverted phase pattern (`1.0 - phase`) caused background areas to map to the top of the colorscale instead of black; now correctly masks background to 0 before display
- **Data Loading Compatibility**: Switched to `sep=r'\s+'` with `engine='python'` and `dtype={'chromosome': str}` to handle both tab-delimited (real AncestryDNA) and whitespace-delimited (mock) files; filters out header rows read as data
- **Pandas 3.0 Sampling Compatibility**: Replaced `groupby().apply()` sampling with `pd.concat` pattern to preserve chromosome column in pandas 3.0+
- **Legend Color Accuracy**: Legends now pull dynamic low/high colors from the active palette instead of hardcoded `#00BFFF`/`#FF00FF`
- **Export Background Consistency**: Standalone PNG/SVG `paper_bgcolor` and `plot_bgcolor` changed from `#111111` to `#000000` to match the main visualization
- **Export Block Scope**: Moved image export code inside `if __name__ == "__main__"` block (was previously at module level, causing `NameError` on import)
- **Deprecated Kaleido Engine Argument**: Removed `engine="kaleido"` from `write_image()` calls to silence deprecation warnings

### Improved

- **Genetic Data Processing**: Successfully handles 650,230+ SNPs from real AncestryDNA files
- **Privacy Protection**: `.gitignore` updated to exclude generated genome art output (`art/genome-2-art.png`, `art/genome-2-art.svg`) and git backups
- **Code Organization**: Centralized `COLOR_PALETTES` dictionary at module level with `get_palette_legend_colors()` helper

### Usage Examples

```bash
# Use default scientific palette with mock data
python genome-2-art_v2.py

# Use with real DNA data
python genome-2-art_v2.py -i AncestryDNA.txt

# Use ocean color palette
python genome-2-art_v2.py -i AncestryDNA.txt -p ocean

# List all available palettes
python genome-2-art_v2.py --list-palettes
```

### Files Changed

- `genome-2-art_v2.py`: Color palette system, CLI arguments, data loading fixes, export fixes
- `.gitignore`: Added generated art output and git backup exclusions

---

## [2.0.0] - Previous Version

### Features

- Initial quantum wave genome visualization
- CIRCOS plot generation
- SNP distribution visualization
- Enhanced quantum iris patterns
- Basic color scheme support

---

_For older versions and detailed commit history, see the Git log._
