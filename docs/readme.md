# notes to maintain documentation

- built with [Documenter](https://juliadocs.github.io/Documenter.jl).
- deployed [here](https://crsl4.github.io/PhyloNetworks.jl/)
  (go to `dev/` or `stable/`)
  using GitHub and files committed to the `gh-pages` branch.

## how it works: overview

- `.travis.yml` asks to start the doc project
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
Therefore: any public *function* needs to be manually listed in `docs/src/man/public.md`,
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
before the doc is built, the script `.travis.yml` installs `R` on the server,
sets up the julia environment with dependencies like `PhyloPlots` before
starting the build in itself.
Note that, for an unknown reason, `R` must be installed *outside* of (the former) `make.sh`,
in the main body of `.travis.sh`.

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
    plot(net, :R);
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

```shell
julia --project=docs/ -e 'using Pkg; Pkg.instantiate(); Pkg.develop(PackageSpec(path=pwd()))'
julia --project=docs/ --color=yes docs/make.jl
```

or interactively in `docs/`:

```shell
pkg> activate .
pkg> status # just to check
pkg> status --manifest
pkg> instantiate # after deleting Manifest.toml and undo changes to Project.toml
pkg> # add RCall#master # in case some dependency causes an issue
pkg> dev PhyloNetworks
pkg> add PhyloPlots#master # to get the master branch: done by make.jl
julia> include("make.jl")
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

## earlier versions

### weave

[`Weave`](https://github.com/mpastell/Weave.jl) was used to format the
documentation pages until PhyloNetworks v0.7.0
(see [v0.7.0 doc](http://crsl4.github.io/PhyloNetworks.jl/v0.7.0/)),
and the saving of the plots on the Git repository was handled with an
extra Travis environment variable DRAW_FIG.
Instructions about this previous previous setup can be found in this very
[docs/readme file, in v0.7.0](https://github.com/crsl4/PhyloNetworks.jl/blob/v0.7.0/docs/readme.md).
An extra step used file [make_weave.jl](https://github.com/crsl4/PhyloNetworks.jl/blob/v0.7.0/docs/src/man/src/make_weave.jl).

### mkdocs

[MkDocs](http://www.mkdocs.org) was used to format the documentation pages until
PhyloNetworks [v0.8.0](http://crsl4.github.io/PhyloNetworks.jl/v0.8.0/).
For this:
- `makedocs( ..., format = Markdown(), ...)`
  created vanilla "GitHub md" files only (no html conversion),
- `deploydocs` called `mkdocs` to turn the markdown files in `docs/build/*.md` into html files:

        deploydocs(
            repo = ...,
            deps= Deps.pip("pygments", "mkdocs==0.17.5", "mkdocs-material==2.9.4", "python-markdown-math"),
            make = () -> run(`mkdocs build`),
            target = "site" # which files get copied to gh-pages
        )

    problem: Travis uses python 2, but mkdocs v1.0 needs python 3.
  Versions of `mkdocs` and `mkdocs-material` specified manually to avoid conflicts,
  See https://discourse.julialang.org/t/mkdocs-material-in-documenter/13764/3
- In this conversion, the `mkdocs-material` package was used, for its "material" theme,
  via configuration in the `mkdocs.yml` file, in `docs/`:


```yml
# this is for MkDocs, to turn the .md files produced by Documenter into .html files
site_name:        PhyloNetworks.jl
repo_url:         https://github.com/crsl4/PhyloNetworks.jl
site_description: PhyloNetworks is a Julia package for the manipulation, visualization and inference of phylogenetic networks.
site_author:      Claudia Sol&iacute;s-Lemus

theme:
  name: 'material'
  palette:
    primary: 'teal'
    accent:  'teal'
  logo: 'snaq_small.png'

extra_css:
  - assets/Documenter.css

extra_javascript:
  - https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/MathJax.js?config=TeX-AMS_HTML
  - assets/mathjaxhelper.js

markdown_extensions:
  - extra
  - tables
  - fenced_code
  - mdx_math
  - codehilite

docs_dir: 'build'

nav:
  - Home: index.md
  - Manual:
    - Installation: man/installation.md
    - Input Data for SNaQ: man/inputdata.md
    - TICR pipeline: man/ticr_howtogetQuartetCFs.md
    - Network estimation and display: man/snaq_plot.md
    - Network comparison and manipulation: man/dist_reroot.md
    - Candidate Networks: man/fixednetworkoptim.md
    - Extract Expected CFs: man/expectedCFs.md
    - Bootstrap: man/bootstrap.md
    - Multiple Alleles: man/multiplealleles.md
    - Continuous Trait Evolution: man/trait_tree.md
    - Parsimony on networks: man/parsimony.md
  - Library:
    - Public: lib/public.md
    - Internals: lib/internals.md
```

to check/install MkDocs:

```shell
pip install --upgrade pip
pip install --upgrade mkdocs
pip install --upgrade mkdocs-material
pip install --upgrade python-markdown-math
```
and check the installed versions:
(in comments are versions that work okay together):
```shell
python --version # Python 3.5.5 :: Anaconda, Inc.
mkdocs --version              # v0.17.4  v1.0.4
pip show mkdocs-material      # v2.9.2   v3.2.0
pip show Pygments             # v2.2.0   v2.3.1
pip show pymdown-extensions   # v4.11    v4.11
pip show python-markdown-math # v0.6     v0.6
```

then use mkdocs to build the site.
this step creates a `site/` directory with html files.
they can be viewed at http://127.0.0.1:8000 (follow instructions)

```shell
mkdocs build
mkdocs serve
```
