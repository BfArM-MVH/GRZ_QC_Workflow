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
