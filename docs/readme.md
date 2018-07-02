# notes to maintain documentation

- built with [Documenter](https://juliadocs.github.io/Documenter.jl).
- deployed [here](http://crsl4.github.io/PhyloNetworks.jl/)
  (go to `latest/` or `stable/`)
  using GitHub and files committed to the `gh-pages` branch.

## how it works: overview

- `.travis.yml` asks to run `docs/make.jl` after a successful test & build.
- the julia script `docs/make.jl` has 2 steps:
  1. run `makedocs()` from `Documenter`: make the documentation.
     also runs all `jldoctest` blocks in the source files, to check that
     the output in the blocks matches the actual output.
     This steps also translate the "Documenter md" documentation files
     into vanilla "GitHub md" (see below).
  2. run `deploydocs(...)` also from Documenter. This step calls `mkdocs`,
     which turns the markdown files in `docs/.../*.md` into html files.

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
and then calls the script `docs/make.sh`, that installs `PhyloPlots` before
starting the build in itself.
Note that, for an unknown reason, `R` must be installed *outside* of `make.sh`,
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

Note that in previous versions of the doc, `Weave` was used to format the
documentation pages (until version 0.7.0), and the saving of the plots on the
Git repository was handled with an extra Travis environment variable DRAW_FIG.
Instructions about this previous previous setup can be found in this very
readme file, in version 0.7.0.

## to make a local version of the website

```shell
cd ~/.julia/v0.6/PhyloNetworks/docs
julia --color=yes make.jl
```

first line: adapt to where the package lives.
second line:
- tests the `jldoctest` blocks of examples in the docstrings
- creates or updates a `build/` directory with markdown files.
- does *not* convert the markdown files into html files.

To do this html conversion, use [MkDocs](http://www.mkdocs.org) directly,
and the mkdocs-material package (for the "material" theme).
First check/install MkDocs:

```shell
pip install --upgrade pip
pip install --upgrade mkdocs
pip install --upgrade mkdocs-material
pip install --upgrade python-markdown-math
```
and check the installed versions:
```shell
python --version
mkdocs --version
pip show mkdocs-material
pip show Pygments
pip show pymdown-extensions
pip show python-markdown-math
```

then use mkdocs to build the site.
this step creates a `site/` directory with html files.
they can be viewed at http://127.0.0.1:8000 (follow instructions)

```shell
mkdocs build
mkdocs serve
```

## references

big difference:

    [blabla](@ref)
    [`blabla`](@ref)

The first version will look for a *section* header "blabla", to link to that section.
The secon version will look for a *function* named "blabla",
to link to the documentation for that function.
