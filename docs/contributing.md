# Contributing

## Maintenance

### Update Conda development environment lock file

```shell
conda lock \
  --kind explicit \
  --platform linux-64 \
  --filename-template 'environment-dev.conda.{platform}.lock' \
  --file environment-dev.yaml
```
