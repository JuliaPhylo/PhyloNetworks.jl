CFnetworks
==========

Dissertation work on julia

Claudia August 2014

### Steps to setup git:

git clone https://github.com/crsl4/CFnetworks (only once)

**in Mac:git_laptop**

1. git pull (origin master (or other branch)):

better to do: git fetch
              git status (to check if remote is ahead of local)
              git merge FETCH_HEAD

2. make changes

3. git add .

4. git commit -m "message"

5. git push (origin master (or other branch))

**in stat:git_work**

need to compare to version in ane/public/quartetNetwork

2. run compare_jl.sh (with argument "classes",...) to check if two versions are equal, or check date file

1. git pull (origin master (or other branch))

better to do: git fetch
              git status (to check if remote is ahead of local)
              git merge FETCH_HEAD


3. adapt changes from public, consider different branches

4. make new changes if needed

5. git add .

6. git commit -m "message"

7. git push (origin master (or other branch))

8. run sync_public.sh to have same version in ane/public, add date in date file

### Create branches

git branch bla

git checkout bla

## Implementation of pseudolikelihood in Julia

Need to log into desk00 for julia

```julia
include("types.jl")
include("functions.jl")
include("case_f_example.jl")
```

## Procedure to create hybrid network:

1. create edges defined as hybrid or not

2. create nodes defined as hybrid, leaves or tree with such edges

3. setNode! to add nodes into edges

4. create hybrid network

5. updateInCycle! updateGammaz! updateGamma2z! updateContainRoot!
