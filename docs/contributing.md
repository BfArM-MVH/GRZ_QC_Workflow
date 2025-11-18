# Contributing

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
