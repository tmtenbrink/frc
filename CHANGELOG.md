# Changelog

## 0.1.2 2026-10-31

* **Fix the half-bit threshold on NumPy 1.24+**: `half_bit_threshold` and `frc_res(..., threshold="half_bit")` no longer fail with `module 'numpy' has no attribute 'int'` ([#3](https://github.com/tmtenbrink/frc/issues/3)).
* NumPy, SciPy, DIPlib and rustfrc no longer have an upper version bound, so NumPy 2 is allowed.
* On Python 3.8-3.10, DIPlib is capped at the last release with wheels for that version, so installers no longer pick a DIPlib they cannot install.

### Internal

* Replace Poetry with uv and update GitHub Actions, pinned to commit SHAs
* Publish to PyPI with Trusted Publishing from a protected branch, see `docs/RELEASING.md`
* Fix the license metadata to MIT

## 0.1.1 2024-04-03

* **Python 3.12 support** (by updating to `rustfrc 1.1.5`)
* **Drop Python 3.7 support**

## 0.1.0 2022-04-14

* GitHub release for citation
