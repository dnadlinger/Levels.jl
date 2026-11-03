# [Developer documentation](@id dev_docs)

See the [Contributing guidelines](90-contributing.md) for proposing changes.
Run the commands below from the repository root with Julia 1.12 or newer.
The root project, tests and documentation share a Julia workspace and manifest.

## Local setup

Fork and clone the repository, then add the original as `upstream`:

```bash
git remote add upstream https://github.com/dnadlinger/Levels.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Maintainers working directly in the repository can use `origin` wherever this
guide says `upstream`.

## Tests

Run the suite used by CI:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Test items can be selected by name or tag. Pass the repository root to
`run_tests` so TestItemRunner finds the package project:

```bash
julia --project=test -e 'using TestItemRunner; TestItemRunner.run_tests("."; filter=ti->:fast in ti.tags, verbose=true)'
julia --project=test -e 'using TestItemRunner; TestItemRunner.run_tests("."; filter=ti->occursin("Spectro", ti.name), verbose=true)'
```

## Formatting and linting

Use an [EditorConfig](https://editorconfig.org) plugin in your editor.
Install [pre-commit](https://pre-commit.com), then enable the hook and run the
checks:

```bash
pre-commit install
pre-commit run -a
```

The hook runs the configured checks before each commit. JuliaFormatter's
version is pinned in `.pre-commit-config.yaml`; use that version for manual
formatting too, since line-breaking rules can differ between versions.

CI also checks links with [lychee](https://github.com/lycheeverse/lychee).
With lychee installed, run:

```bash
lychee --no-progress --config lychee.toml .
```

## Branches and pull requests

Start a branch from the current upstream `master`:

```bash
git fetch upstream
git switch master
git merge --ff-only upstream/master
git switch -c 42-add-answer-universe
```

Use the issue number and a short imperative description for the branch name.
For small changes without an issue, use a descriptive prefix such as `typo`
or `hotfix`; discuss larger changes in an issue first.

Keep commits focused, with informative titles in imperative or present tense.
Explain breaking changes in the commit body. Before opening a pull request,
run tests and pre-commit, and rebase onto any new upstream changes:

```bash
git fetch upstream
git rebase upstream/master
git push -u origin HEAD
```

## Building the documentation

Build once, or run a live preview:

```bash
julia --project=docs docs/make.jl
julia --project=docs docs/serve.jl
```

The workspace already points at the local package; no separate `Pkg.develop`
step is needed. The pages under "Examples" are rendered from the Literate.jl
scripts in `examples/`. Repository scaffolding comes from BestieTemplate/Copier;
`copier update` may overwrite changes to generated files.

## Releases

1. Create a `release-x.y.z` branch and update `version` in `Project.toml`.
2. In `CHANGELOG.md`, rename `Unreleased` to `[x.y.z] - yyyy-mm-dd`, add a new
   `Unreleased` section, and update the comparison links at the bottom.
3. Commit as `Release vx.y.z`, open a pull request, and merge after checks pass.
4. Comment `@JuliaRegistrator register` on the merged commit and follow the
   registry pull request returned by the bot.
5. After registration, check that TagBot creates the tag and release, then
   verify the documentation workflow and the version shown in the
   [stable docs](https://dnadlinger.github.io/Levels.jl/stable).
