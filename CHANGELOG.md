# Changelog

## 0.1.2 - 2026-10-31

### Fixed

- `half_bit_threshold` and `frc_res(..., threshold="half_bit")` no longer fail with `module 'numpy' has no attribute 'int'` on NumPy 1.24 and newer ([#3](https://github.com/tiptenbrink/frc/issues/3)).
- On Python 3.8-3.10, DIPlib is capped at the last release with wheels for that Python, so installers no longer pick a DIPlib they cannot install.
- The license metadata says MIT, as `LICENSE` does.

### Changed

- NumPy, SciPy, DIPlib and rustfrc no longer have an upper version bound, so NumPy 2 is allowed.
- Releases are approved with a security key through [touchgate](https://github.com/borink-org/touchgate) and published with Trusted Publishing and PEP 740 attestations.

## 0.1.1 - 2024-04-03

### Added

- Python 3.12 support, by updating to rustfrc 1.1.5.

### Removed

- Python 3.7 support.

## 0.1.0 - 2022-04-14

- GitHub release for citation.
