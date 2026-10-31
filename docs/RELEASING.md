# Releasing

frc is released to PyPI by `.github/workflows/release.yml`, with Trusted Publishing and no API token. This file records the process and the settings it relies on.

## A release

1. On `main`, set the version in `pyproject.toml`, add a section `## <version> <date>` at the top of `CHANGELOG.md`, and push. Wait for `CI` to pass.
2. Start the `release` workflow on the `publish` branch with the full SHA of that commit:

   ```nu
   gh workflow run release.yml --ref publish -f commit=(git rev-parse HEAD)
   ```

3. The workflow checks that the commit is on `main`, that the version and the newest changelog section agree and that the version is not tagged yet. It waits for every check on the commit to pass, builds the sdist and the wheel, signs PEP 740 attestations for them, publishes both to PyPI, and only then tags the commit `v<version>` and makes the GitHub release with the changelog section and the built files.

A release that failed halfway can run again: the upload skips what PyPI has already.

## Settings

The release is only as safe as these settings. They are made by `tmtenbrink`, the repository's owner, in the browser; no tool or agent holds its credentials.

- `tiptenbrink`, whose `gh` token agents use, is a collaborator with write access, never admin.
- The `publish` branch holds the copy of `release.yml` that runs. A ruleset restricts its creation, updates and deletion and blocks force pushes, with only the repository admin allowed to bypass it. `main` needs a `release.yml` too, because GitHub starts a workflow by hand only if the default branch has a file at that path. Only the copy on `publish` ever releases: run from `main`, the publish job is refused by the `release` environment. Keep the `main` copy the same as the `publish` copy, so that `main` shows what runs: when `tmtenbrink` changes `publish`, make the same change on `main`.
- The `release` environment admits the `publish` branch alone, by its exact name.
- Immutable releases are on.
- On PyPI, `frc` trusts only the workflow `release.yml` of `tmtenbrink/frc` in the environment `release`. No PyPI API token for the project exists, and the repository has no `PYPI_API_TOKEN` secret.
