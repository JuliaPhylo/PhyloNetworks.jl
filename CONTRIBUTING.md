# Contributing to PhyloNetworks

The following guidelines are designed for contributors of `PhyloNetworks`.

## Reporting Issues and Questions

For reporting a bug, a failed function or requesting a new feature,
please open an issue in [the issue tracker](https://github.com/juliaphylo/PhyloNetworks.jl/issues).
First, search through existing issues (open or closed) that might have the answer to your question.
If you are reporting a bug, please also include a minimal code example or all relevant information for us to replicate the issue.

For general questions, make sure to check out the [google user group](https://groups.google.com/g/phylonetworks-users). If you cannot find answers to your question, please post a new question. We do our best to reply in a timely fashion, but we are undermanned so we appreciate your patience.

## Contributing Code

To make contributions to `PhyloNetworks`, you need to set up a [GitHub](https://github.com/) account (if you do not have one) and request your change(s) or contribution(s) via a pull request against the `master` branch of `PhyloNetworks` from a non-master branch in your fork.
Using a non-master branch on your end will give developers push access to your
branch to make edits to it (in case we want to work collaboratively on the new code).

Please use the following steps:

1. Fork the `PhyloNetworks` repository to your GitHub account
2. Clone your fork locally with `git clone`
3. Create a new branch with a name that describes your contribution.
   For example, if your contribution is fixing a bug in `readnewick`,
   your new branch can be named `fix-bug-readnewick`.
   You can create it and switch with:
   ```
   git checkout -b fix-bug-readnewick
   ```
4. Make your changes on this new branch.
   Write a [docstring](https://docs.julialang.org/en/v1/manual/documentation/)
   and add one or more tests for each function.
   Make sure that your code passes all the
   [tests](https://pkgdocs.julialang.org/v1/managing-packages/#Testing-packages).
5. Push your changes to your fork
6. [Submit a pull request](https://github.com/juliaphylo/PhyloNetworks.jl/pulls) against the `master` branch in `PhyloNetworks`. Make sure that your code passes all the automatic tests and that it is not in conflict with the current status of `master`

Please make sure to follow the Julia package guidelines and conventions on your code. `PhyloNetworks` was created before these conventions were catalyzed, but we are attempting to follow them going forward.

To learn more about the Julia conventions, check out the following links:

- [Julia package development](http://web.mit.edu/julia_v0.6.0/julia/share/doc/julia/html/en/manual/packages.html#Package-Development-1)
- [Creating Julia packages](https://pkgdocs.julialang.org/v1/creating-packages/)
- [Julia style guide](https://docs.julialang.org/en/v1/manual/style-guide/)
