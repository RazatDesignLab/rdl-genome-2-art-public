# Security Policy

## Our Commitment

Genome-2-Art processes sensitive genetic information. While all processing happens locally on your machine and no data is transmitted externally, we take security seriously and appreciate responsible disclosure of any vulnerabilities.

## Supported Versions

| Version | Supported          |
| ------- | ------------------ |
| Latest  | :white_check_mark: |

## Reporting a Vulnerability

If you discover a security vulnerability, please report it responsibly:

### For General Security Issues

1. **Do NOT open a public issue** for security vulnerabilities
2. Email the maintainers directly or use GitHub's private vulnerability reporting feature
3. Include:
   - Description of the vulnerability
   - Steps to reproduce
   - Potential impact
   - Suggested fix (if any)

### Response Timeline

- **Acknowledgment**: Within 48 hours
- **Initial assessment**: Within 1 week
- **Resolution target**: Depends on severity, typically 2-4 weeks

## Security Considerations

### Genetic Data Privacy

This project is designed with privacy as a core principle:

- **Local processing only** — Your genetic data never leaves your computer
- **No network calls** — The visualization scripts make no external API requests
- **No telemetry** — We collect no usage data or analytics
- **No cloud dependencies** — All processing uses local Python libraries

### Best Practices for Users

1. **Never commit real genetic data** to any repository
2. **Use the mock data** (`AncestryDNA_mock.txt`) for testing and development
3. **Keep your raw DNA files** in a secure, non-synced location
4. **Review the .gitignore** to ensure your data patterns are excluded

### What We Consider Security Issues

- Code that could inadvertently transmit genetic data
- Vulnerabilities in dependencies that affect data handling
- Issues that could expose user data through generated outputs
- Flaws in the local processing pipeline

### What We Don't Consider Security Issues

- Theoretical attacks requiring physical access to the user's machine
- Issues in upstream dependencies with no practical exploit path
- Feature requests disguised as security concerns

## Dependencies

We use well-established scientific Python libraries:

- `pandas` — Data manipulation
- `numpy` — Numerical computing
- `plotly` — Visualization
- `kaleido` — Image export

We recommend keeping dependencies updated:

```bash
pip install --upgrade -r requirements.txt
```

## Acknowledgments

We appreciate security researchers who help keep this project safe. Contributors who responsibly disclose vulnerabilities will be acknowledged (with permission) in our release notes.

---

*Your genetic data is uniquely yours. We're committed to keeping it that way.*
