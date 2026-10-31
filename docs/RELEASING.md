# Releasing

frc is released to PyPI with [touchgate](https://github.com/borink-org/touchgate): a release publishes only after it was approved with the release key, a FIDO2 security key, by touch and PIN. Its README describes the process and what it guarantees; this file records how this repository is set up for it. touchgate releases frc to PyPI because `pyproject.toml` has a `[project]` table.

## A release

1. Start the `touchgate` workflow on `main` with the version, such as `0.1.2`:

   ```nu
   gh workflow run touchgate.yml --repo tmtenbrink/frc -f version=0.1.2
   ```

   It dates the `## Unreleased` section of `CHANGELOG.md`, sets the version in `pyproject.toml` and `uv.lock`, and opens a pull request from `release/<version>`.

2. Check out that branch, look it over, and approve it with the release key, which asks for its PIN and a touch. An approval lasts an hour, so merge once the checks pass, or approve again.

   ```nu
   git fetch origin
   git switch release/0.1.2
   touchgate approve
   ```

   Install touchgate with `cargo install --locked --git https://github.com/borink-org/touchgate`.

3. Merge the pull request with a merge commit. The `touchgate` workflow then starts `release.yml` on the `publish` branch. It checks the approval, waits for every check on the merge, builds the sdist and the wheel, signs their PEP 740 attestations and publishes them to PyPI, and only then tags `v<version>` and makes the GitHub release with the built files attached.

A release that failed halfway can run again: the upload skips what PyPI has already.

Every pull request adds a line under `## Unreleased` in `CHANGELOG.md`, or carries the `no-changelog` label. The `touchgate` workflow checks both, and that the newest released section of the changelog is the version in `pyproject.toml`.

## Settings

The release is only as safe as these settings. They are made by `tmtenbrink`, the repository's owner, in the browser; no tool or agent holds its credentials.

- `tiptenbrink`, whose `gh` token agents use, is a collaborator with write access, never admin.
- The `publish` branch holds the copy of `release.yml` that runs. A ruleset restricts its creation, updates and deletion and blocks force pushes, with only the repository admin allowed to bypass it. `main` needs a `release.yml` too, because GitHub starts a workflow by hand only if the default branch has a file at that path. Only the copy on `publish` ever releases: run from `main`, the publish job is refused by the `release` environment. Keep the `main` copy the same as the `publish` copy, so that `main` shows what runs: when `tmtenbrink` changes `publish`, make the same change on `main`.
- The `release` environment admits the `publish` branch alone, by its exact name.
- Immutable releases are on, GitHub Actions may create pull requests, merge commits are allowed, and a `no-changelog` label exists.
- On PyPI, `frc` trusts only the workflow `release.yml` of `tmtenbrink/frc` in the environment `release`. No PyPI API token for the project exists, and the repository has no `PYPI_API_TOKEN` secret.

`release.yml` and `touchgate.yml` run touchgate at a pinned commit, and `release.yml` holds the public release key and the uv version, `UV_VERSION`, that builds and uploads the release. Changing any of them in `release.yml` is a change to the `publish` branch, so `tmtenbrink` makes it, after reading what a new touchgate commit or uv release changes.

frc builds with `uv_build`, and uv builds with the copy built into itself only when the `uv_build` range in `pyproject.toml` holds its version; otherwise it downloads one from PyPI. `touchgate check --build` refuses to build in that case, so when `UV_VERSION` moves past the range, widen the range in the same change.
