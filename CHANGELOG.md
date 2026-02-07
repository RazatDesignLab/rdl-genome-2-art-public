# Changelog

All notable changes to the Genome-2-Art project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.1.0] - 2025-02-07

### Added

- **Enhanced Color Palette System**: 7 scientifically-inspired color palettes for quantum wave visualization
  - `scientific`: Cyan-magenta gradient for scientific accuracy
  - `ocean`: Deep ocean blues to white
  - `fire`: Black to red to golden yellow
  - `nebula`: Deep space purples to light pink
  - `earth`: Dark earth tones to cream
  - `aurora`: Night sky teals to bright greens
  - `infrared`: Deep purples to pale lavender
- **Command-line Color Palette Support**:
  - `-p, --palette` flag to select color schemes
  - `--list-palettes` option to display available palettes
- **Dynamic Legend System**: Legends now automatically adapt to show actual colors from selected palette
- **Enhanced Virtual Environment Support**: Complete setup with Python 3.13.2 compatibility
- **High-Resolution Export**: Black background support for both PNG and SVG exports

### Fixed

- **Critical Background Rendering Bug**: Wave Phase Field now displays proper black backgrounds instead of unwanted white backgrounds
- **Data Loading Compatibility**: Improved AncestryDNA file parsing for both mock and real genetic data files
- **Visualization Sampling Issues**: Fixed chromosome column preservation during data sampling for large datasets
- **Legend Color Accuracy**: Legends now reflect actual palette colors instead of hardcoded cyan/magenta
- **Export Background Consistency**: Standalone PNG/SVG files now maintain black backgrounds matching the main visualization

### Improved

- **Genetic Data Processing**: Successfully handles 650,230+ SNPs from real AncestryDNA files
- **Privacy Protection**: Virtual environment (`venv/`) properly excluded from Git commits via .gitignore
- **Code Organization**: Enhanced color palette management with centralized COLOR_PALETTES dictionary
- **User Experience**: Better error handling and debug information for data loading
- **Visual Quality**: Consistent black backgrounds across all plot types and export formats

### Technical Details

- **Python Version**: Compatible with Python 3.13.2
- **Dependencies**: All packages from requirements.txt installed successfully
  - kaleido (visualization export)
  - plotly (interactive plotting)
  - pandas (data manipulation)
  - numpy (numerical computing)
  - matplotlib (plotting backend)
  - seaborn (statistical visualization)
- **Data Security**: Real genetic data (AncestryDNA.txt) protected by .gitignore patterns
- **Performance**: Optimized sampling for large datasets while preserving chromosome information

### Usage Examples

```bash
# Use default scientific palette
python genome-2-art_v2.py -i AncestryDNA.txt

# Use ocean color palette
python genome-2-art_v2.py -i AncestryDNA.txt -p ocean

# List all available palettes
python genome-2-art_v2.py --list-palettes

# Use with virtual environment
source venv/bin/activate && python genome-2-art_v2.py -i AncestryDNA.txt -p fire
```

### Files Changed

- `genome-2-art_v2.py`: Major enhancements to color system and bug fixes
- `.gitignore`: Already properly configured for virtual environments
- `requirements.txt`: All dependencies verified and working
- `art/genome-2-art.png`: Export functionality improved with black backgrounds
- `art/genome-2-art.svg`: Vector export with proper background rendering

### Development Environment

- **Virtual Environment**: `venv/` directory created with Python 3.13.2
- **Package Management**: All dependencies installed via pip
- **Git Integration**: Virtual environment properly excluded from version control
- **Cross-platform**: Compatible with macOS, tested on Apple Silicon

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
