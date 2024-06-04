# notes to maintain documentation

- built with [Documenter](https://documenter.juliadocs.org/stable/).
- deployed [here](https://crsl4.github.io/PhyloNetworks.jl/)
  (go to `dev/` or `stable/`)
  using GitHub and files committed to the `gh-pages` branch.

## how it works: overview

- `.github/workflows/documentation.yml` asks to start the doc project
  (installs R and dependencies like `PhyloPlots` & `Documenter`) and
  run `./docs/make.jl` after a successful test & build.
- the julia script `docs/make.jl` has these steps:
  1. check out the master version of PhyloPlots
  2. run `makedocs()` from `Documenter`: make the documentation.
     also runs all `jldoctest` blocks in the source files, to check that
     the output in the blocks matches the actual output.
     This steps also translate the "Documenter md" documentation files
     into html files.
  3. run `deploydocs(...)` also from Documenter:
     to push the files on github, gh-pages branch.

for now, docstrings are automatically used to build an entry for
- each internal thing that has a docstring (e.g. not exported in `src/PhyloNetworks.jl`)
- each public *type*
Therefore: any public *function* needs to be manually listed in `docs/src/lib/public.md`,
in a section to get a nice organization of all these manual entries.

## The "Documenter md" format

### Note on the format

The documentation pages are all written in this format. It is a regular md, but
with extra blocks of codes (as `@example` and `@setup`) that contain Julia
commands. These lines will be executed during the `makedoc()` process. See the
`Documenter` [doc](https://juliadocs.github.io/Documenter.jl/stable/man/syntax/)
for more details on the syntax. For instance, @example blocks with the same "name"
are run in the same session. Otherwise, an @example blocks with no name
is run in its own anonymous Modules.

### Setting up the plot environment

Some of these blocs may contain plots, which are going to be drawn during the
process, requiring the use of `PhyloPlots` along with `RCall`. Hence,
before the doc is built, `.github/workflows/ci.yml` installs `R` on the server,
sets up the julia environment with dependencies like `PhyloPlots` before
starting the build in itself.

### Directory of the plots

We chose to group all the output plots in the directory `assets/figures`.
Hence, the typical setup in a documentation page containing plots is:

    ```@setup name
    using PhyloPlots, RCall
    mkpath("../assets/figures")
    R"name <- function(x) file.path('..', 'assets', 'figures', x)" # hide
    ```

The `mkpath` command is there to ensure that the target directory does indeed
exist. In theory, it only needs to be called once (by the first documentation
page being built). However, as this order might be subject to change over time,
it could be safer to include it on every such page.

### Format of the plots

After trial and error and discussions, we chose to use only the *SVG* format
for the plot. This format should ensure that when a plot is drawn again,
identical, in a new build, then Git will recognize that it has not change, and
hence not save a new version of it. This should ensure that the repository does
not grow unreasonably at each build of the doc, i.e. after each push to
master. The typical commands to save and display a plot should hence be:

    ```@example name
    R"svg(name('my_useful_name.svg'), width=4, height=4)" # hide
    plot(net);
    R"dev.off()" # hide
    nothing # hide
    ```
    ![my_useful_name](../assets/figures/my_useful_name.svg)

**Warning**: this is not like an interactive session. If the same file name
is re-used by some other documentation page for some other plot, only the
final version of the plot will be committed by git, with possible unintended
consequences. Make sure to use different file names for plots that are supposed
to look different (across the whole site).

## to make a local version of the website

generally: to see deprecation warnings, add the option `--depwarn=yes`.

```shell
julia --project=docs/ -e 'using Pkg; Pkg.instantiate(); Pkg.develop(PackageSpec(path=pwd()))'
julia --project=docs/ --color=yes docs/make.jl
```

or interactively in `docs/`:

```shell
pkg> activate .
pkg> status # just to check
pkg> status --manifest
pkg> instantiate # after deleting Manifest.toml and undoing changes to Project.toml
pkg> # add RCall#master # in case some dependency causes an issue
pkg> dev PhyloNetworks
pkg> add PhyloPlots#master # to get the master branch: done by make.jl
julia> include("make.jl")
```

or, after project & manifest setup:
```shell
julia --project --color=yes -e 'include("make.jl")'
```

it will:
- test the `jldoctest` blocks of examples in the docstrings
- create or update a `build/` directory with html files.

To see the result, just open `docs/build/index.html` in a browser and follow the links.

## references

big difference:

    [blabla](@ref)
    [`blabla`](@ref)

The first version will look for a *section* header "blabla", to link to that section.
The secon version will look for a *function* named "blabla",
to link to the documentation for that function.
