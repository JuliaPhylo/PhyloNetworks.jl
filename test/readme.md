# Tests functions

- all in `runtests.jl`, which calls other test files
- to see deprecation warnings when running things locally, start julia with
  ```shell
  julia --depwarn=yes
  ```
  and possibly other options (like --project).
- generally, code in file `src/x.jl` is tested by `test/test_x.jl`,
  but see below for what older test files do (related to SNaQ).  
  checkout PhyloNetworks v0.9.1 or older to see those older files.
