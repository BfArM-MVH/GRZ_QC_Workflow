# Contributing

## Pull requests

PR titles must follow [Conventional Commits](https://www.conventionalcommits.org/). The CI workflow `PR Title Lint` enforces this.

Allowed types: `feat`, `fix`, `chore`, `docs`, `refactor`, `test`, `ci`, `build`, `perf`, `style`, `revert`.

The PR title becomes the squash-merge commit subject and is parsed by release-please to determine the next version:

| Title prefix                             | Version bump | Example            |
| ---------------------------------------- | ------------ | ------------------ |
| `fix:`                                   | patch        | `1.3.1` -> `1.3.2` |
| `feat:`                                  | minor        | `1.3.1` -> `1.4.0` |
| `feat!:` (or `BREAKING CHANGE:` in body) | major        | `1.3.1` -> `2.0.0` |
| others (`chore`, `docs`, ...)            | no release   |                    |

Subject must start with a lowercase letter and not end with a period. WIP PRs are exempt while marked WIP.

## Releases

Releases are automated by [release-please](https://github.com/googleapis/release-please).

On every push to `main`, release-please opens (or updates) a release PR titled `chore(main): release <version>`. The PR aggregates conventional commits since the last release and updates:

- `CHANGELOG.md`
- `manifest.version` in `nextflow.config`
- `.release-please-manifest.json`

Merging the release PR creates a Git tag and GitHub Release. Do not manually edit `CHANGELOG.md` or bump `manifest.version` in PRs.

## Formatting/linting files

You can easily create a development environment with all the following software on Linux with:

```shell
conda create --name grzqc --file environment-dev.conda.linux-64.lock
```

### Prettier

```shell
prettier --check --write .
```

### Nextflow

```shell
find main.nf workflows/ subworkflows/local modules/local -exec nextflow lint -format -sort-declarations -harshil-alignment {} +
```

### Python

```shell
ruff format bin/
ruff check --fix --extend-select I bin/
```

## Maintenance

### Update Conda development environment lock file

```shell
conda lock \
  --kind explicit \
  --platform linux-64 \
  --filename-template 'environment-dev.conda.{platform}.lock' \
  --file environment-dev.yaml
```

### Update local module Docker/Singularity environments

Use the [Seqera Containers](https://seqera.io/containers) service to generate containers that match the Conda environment specification for the module.

Add all of the same packages at the same versions and then use "Get Container".
You can copy the new URL into the `main.nf` file.
The Docker URL is the second URL, the Singularity URL is the first (on top).
To generate the Singularity container, choose "Singularity" instead of "Docker" next to "Container settings".
You must then view the build logs and choose the "HTTPS" checkbox when it is ready, then copy that URL into the `main.nf` file for the local module.
