# Changelog

## Unreleased

frc 0.3 computes the FRC in Rust, through rustfrc 2, and depends only on NumPy and rustfrc.

### Added

- A `seed` for `one_frc`. The same image and seed give the same curve on every CPU; without a seed, every call splits with a fresh one.
- `utility.tukey_window(length, alpha)`, the one-dimensional window that `tukey_square` uses.

### Changed

- `one_frc` and `two_frc` run in rustfrc 2: the split, the Fourier transforms and the ring sums. `one_frc` on a 4096 x 4096 image takes about 30 ms instead of about 1 s, on 8 threads of an AMD EPYC 9V45. Rings hold the same frequencies as DIPlib's `RadialSum`, which frc used before.
- `one_frc` truncates float images to int32 counts before the split, for both halves; 0.2 kept the fractional part in the second half.
- `gaussf` filters with NumPy. For sigmas of 0.8 and more it gives what DIPlib gave, to within 1e-14; for smaller sigmas, DIPlib used a Fourier method instead.
- `tukey_square` computes its window with NumPy, as SciPy did.
- frc requires rustfrc 2.0 or later, and no longer caps the Python version.

### Removed

- DIPlib and SciPy as dependencies, and with them DIPlib images as input: pass NumPy arrays, for example `np.asarray(dip_image)`.
- Python 3.10 support. It has reached its end of life, and rustfrc 2 requires Python 3.11.

### Added

- Support for Python 3.13 and 3.14, with DIPlib 3.5.2 and 3.6.

### Changed

- frc requires rustfrc 1.2 or later, which brings wheels for every CPython from 3.9 on.
- The oldest NumPy, SciPy and DIPlib allowed on each Python is the oldest release with wheels for it, so installers no longer try to build one from source.

### Removed

- Python 3.8 and 3.9 support. Both have reached their end of life, and no pytest that runs on them fixes GHSA-6w46-j5rx-g56g.

## 0.1.3 - 2026-11-01

### Changed

- Releases upload PEP 740 attestations, so PyPI shows the workflow that built and published each file. 0.1.2 was published without them.

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
