## Summary

<!-- Brief description of the change and why -->

## Tasks

- [ ] PR title follows [Conventional Commits](https://www.conventionalcommits.org/) (`feat:`, `fix:`, `chore:`, `docs:`, `refactor:`, `test:`, `ci:`, `build:`, `perf:`, `style:`, `revert:`)
  - Use `feat!:` or include `BREAKING CHANGE:` in the body for breaking changes
  - Version bump and `CHANGELOG.md` are handled automatically by release-please on merge to `main`
- [ ] Request reviews from the relevant people

<!--
Version impact (decided by PR title prefix on squash merge):
  fix:    -> patch bump   (1.3.0 -> 1.3.1)
  feat:   -> minor bump   (1.3.0 -> 1.4.0)
  feat!:  -> major bump   (1.3.0 -> 2.0.0)
  others  -> no release
-->
