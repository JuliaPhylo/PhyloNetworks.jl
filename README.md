CFnetworks
==========

Dissertation work on julia

Claudia August 2014

### Steps to setup git:

1. git clone https://github.com/crsl4/CFnetworks (only once)

** in Mac**

1. git pull origin master

2. make changes

3. git add .

4. git commit -m "message"

5. git push origin master (or other branch)

**in stat**

1. Run JULIA/compare_jl.sh and/or compare date_sync_public.txt with date in git_work

2. Check if there are differences between classes_public.jl and classes.jl 
   (and other .jl files)

3. git pull origin master in git_work/CFnetworks to update classes.jl et al

4. Add differences to classes.jl in git_work from public 
(consider using another branch)

5. git add .
   git commit -m ""
   git push origin master

6. Run JULIA/sync_public.sh to copy updated versions in public

### Create branches

git branch bla

git checkout bla

## Implementation of pseudolikelihood in Julia
## need to log into desk00 for julia

```julia
include("types.jl")
include("functions.jl")
include("case_f_example.jl")
```
