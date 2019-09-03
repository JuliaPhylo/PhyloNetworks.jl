# TICR pipeline

PhyloNetworks' [wiki](https://github.com/crsl4/PhyloNetworks.jl/wiki/TICR:-from-alignments-to-quartet-concordance-factors)
has a step-by-step tutorial,
to go from multiple sequence alignments
to a table of quartet gene frequencies (concordance factors: CFs),
through BUCKy (to integrate out gene tree uncertainty) or through RAxML.
To get the `raxml.pl` perl script to run RAxML on each gene,
download the content of that wiki with

`git clone https://github.com/crsl4/PhyloNetworks.jl.wiki.git`

then go to the `script/` folder.  
Full information and code is [here](https://github.com/nstenz/TICR).
Below is more detail to insert data into the pipeline at various stages,
using one of **two pipelines**, depending on your machine configuration:

- "no-scheduler" pipeline: original set of Perl scripts from the
  TICR pipeline, well suited for a machine or a cluster of machines
  without a job scheduler.
  The scripts automatically parallelize the work across the available cores.

- "slurm" pipeline: well suited for a cluster where users submit jobs
  via a job scheduler like [SLURM](https://slurm.schedmd.com/) or
  [SGE](https://en.wikipedia.org/wiki/Oracle_Grid_Engine). The job scheduler
  does the work of parallelizing the work across available cores.
  The scripts, in this second pipeline, were created to take full advantage
  of job scheduler capabilities. They were developed for a cluster running SLURM.
  Adjustments to the submit scripts will be needed, to adapt to your own
  SLURM configuration or to the syntax that your job scheduler wants.

## To run MrBayes: we already have alignments

### no-scheduler pipeline

We don't need to run `mdl.pl` if we already have aligned gene sequences
from separate loci. To run MrBayes on each locus, we can simply
create a tarball of the Nexus files we wish to use (fasta won't work at this stage).
The command below assumes that we want to use all the files ending with ".nex"
in the current directory, one file per locus:

```bash
tar czf my-genes.tar.gz *.nex
```

Once the tarball has been successfully generated, we can then specify this
file as input for [mb.pl](https://github.com/nstenz/TICR/blob/master/scripts/mb.pl)
assuming we have a valid MrBayes block located in the file "bayes.txt":

```bash
mb.pl my-genes.tar.gz -m bayes.txt -o my-genes-mb
```

If we get an error message like `mb.pl: Command not found`, it might be because
`mb.pl` has no execute permission, or the current directory is not in our "path".
An easy fix is to run this command instead:

```bash
perl mb.pl my-genes.tar.gz -m bayes.txt -o my-genes-mb
```

The resulting output tarball would now be located in `my-genes-mb/my-genes.mb.tar`,
and can be used normally with
[bucky.pl](https://github.com/nstenz/TICR/blob/master/scripts/bucky.pl),
that is, like this:

```bash
bucky.pl my-genes-mb/my-genes.mb.tar -o mygenes-bucky
```

The output, with the table of concordance factors for all sets of 4 taxa,
will be in a file named `my-genes.CFs.csv` inside a directory named `mygenes-bucky`.
That's the file containing the quartet concordance factors to give to SNaQ as input.
There is *no* need to do any of the steps below: they are already done by `bucky.pl`.

### slurm pipeline

SLURM will parallelize the MrBayes runs across genes.

1. Navigate in some "working" directory where you place:
   - a folder containing all nexus files, which we will call "nexusfolder" below
   - a text file named `mb-block.txt` with the MrBayes block to be used
   for all the genes (containing the options for MrBayes: model of sequence
   evolution, number of generations etc.).
   If we want a different MrBayes block for different genes, step 2 should be skipped,
   and we should instead find some other way to put the specific MrBayes block
   at the end of each nexus file.

2. In the "working" directory above, run the julia script
   [`paste-mb-block.jl`](https://github.com/nstenz/TICR/blob/master/scripts-cluster/paste-mb-block.jl)
   with "nexusfolder" as argument, to tell the script where to find all the nexus files:

   ```bash
   julia path/to/paste-mb-block.jl nexusfolder
   ```
   This script will read all the nexus files in the directory `nexusfolder`,
   will create a new directory `nexusfolder-block`,
   and will create new nexus files (containing the MrBayes block found in file `mb-block.txt`)
   as `1.nex, 2.nex, ...` in the new directory. A `translate.txt` file will also be created
   to map the original gene file names to the new (numbered) file names.
   If we named our MrBayes block file differently: we can edit the script and modify it
   to replace `mb-block.txt` by our actual file name for the MrBayes block.

3. Modify the submit script
   [`mb-slurm-submit.sh`](https://github.com/nstenz/TICR/blob/master/scripts-cluster/mb-slurm-submit.sh),
   which will parallelize all the individual-gene MrBayes runs with SLURM:

   - change `--array` to the correct number of genes
   - change `--mail-user` to the user's email (if this is an option for your job scheduler)
   - replace the `/workspace/software/bin` in `PATH="/workspace/software/bin:$PATH"`
     to the path where the `mb` executable is located or put the whole path in the command:
     `/s/mrbayes-3.2.6-1/bin/mb`

   In slurm, we can then submit the MrBayes array job with:

   ```bash
   sbatch mb-slurm-submit.sh
   ```

  With this slurm pipeline, the steps below are needed: keep reading.

## To run mbsum on the output of MrBayes for each gene

If we have the output of MrBayes and want to run BUCKy,
we must first run `mbsum` on the output from MrBayes, separately for each gene.

For a gene with output tree files named `gene1.run1.t`, `gene1.run2.t` and `gene1.run3.t`,
and a desired burnin of 1000 trees per tree file, we do this:

```bash
mbsum -n 1000 -o gene1.in gene1.run1.t gene1.run2.t gene1.run3.t
```

This `mbsum` command will need to be executed for each gene.
Then we can continue to the next section to run bucky.

Alternatively, we can use the julia script
[`mbsum-t-files.jl`](https://github.com/nstenz/TICR/blob/master/scripts-cluster/mbsum-t-files.jl),
and give it as argument the directory that has the output tree files from MrBayes,
to run mbsum for *all* the genes.
`mbsum` is fast, so there is no attempt to parallelize the various mbsum commands.

```bash
julia mbsum-t-files.jl mbfolder outputfolder burnin                # or
julia --color=yes -- mbsum-t-files.jl mbfolder outputfolder burnin # for colorized messages to the screen
```
where `burnin` is replaced by the number of trees to ignore in each tree file
for burnin. This `burnin` argument is optional (default: 2501).
The `outputfolder` will contain the output of `mbsum`.

## To run bucky on all 4-taxon sets: we already have the mbsum output

### no-scheduler pipeline

If we already ran `mbsum` on the output from MrBayes, for each individual gene,
we can simply create a tarball containing all the mbsum output files.
So if we had mbsum output in files named `gene1.in`, `gene2.in`, ... , `gene100.in`,
we would want to run something similar to the following command to create the tarball:

```bash
tar czf my-genes-mbsum.tar.gz gene*.in
```

We can now use this tarball along with the `-s` option in
[bucky.pl](https://github.com/nstenz/TICR/blob/master/scripts/bucky.pl) like this:

```bash
bucky.pl my-genes-mbsum.tar.gz -s -o mygenes-bucky
```

Again, if we get an error like `bucky.pl: Command not found`, we could run instead

```bash
perl bucky.pl my-genes-mbsum.tar.gz -s -o mygenes-bucky
```

The output, with the table of concordance factors for all sets of 4 taxa,
will be in a file named `my-genes.CFs.csv` inside directory `mygenes-bucky`.
That's the file containing the quartet concordance factors to give to SNaQ as input.

### slurm pipeline

We want to run `bucky` on every 4-taxon set.
SLURM will parallelize these jobs with the submit script
[`bucky-slurm-submit.sh`](https://github.com/nstenz/TICR/blob/master/scripts-cluster/bucky-slurm-submit.sh),
which calls the perl script
[`bucky-slurm.pl`](https://github.com/nstenz/TICR/blob/master/scripts-cluster/bucky-slurm.pl).

The perl script
[`bucky-slurm.pl`](https://github.com/nstenz/TICR/blob/master/scripts-cluster/bucky-slurm.pl)
runs `bucky` on a single 4-taxon set.
It takes the following arguments, which must be modified in the submit script
[`bucky-slurm-submit.sh`](https://github.com/nstenz/TICR/blob/master/scripts-cluster/bucky-slurm-submit.sh):
- name of the folder containing the `mbsum` output files (one per locus) from previous step.
  This folder is named `mbsum` in the submit script: adapt if needed.
- output name: `-o` or `--out-dir` name of the directory to store output files in.
  This option is not used in the default submit script
- bucky arguments: `-a` or `--alpha` for the prior alpha value,
  and `-n` or `--ngen` number of generations. These options are not used either,
  in the script: the defaults are used then (α=1, 1 million generations)
- integer for the given quartet, via option `-q`.
  The quartet ID is specified by SLURM with its own array ID: `$SLURM_ARRAY_TASK_ID`.

In the submit script that gives instructions to the job scheduler:
- adapt the name of the `$SLURM_ARRAY_TASK_ID` variable, which captures the task number
  in the array of tasks, to your scheduler syntax
- change `--array` to the correct number of 4-taxon sets.
  For example, if there are 15 taxa in the dataset, there are `1365` 4-taxon sets.
  To get this number, if you are unsure, use `choose(15,4)` in R or `binomial(15,4)` in Julia,
  but replace 15 by your actual number of individuals.
- change `--mail-user` to the user's email (if this is an option for your job scheduler)
- replace the `/workspace/software/bin` in `PATH="/workspace/software/bin:$PATH"`
  by the path where the `bucky` executable is located.
  Also, replace `/workspace/claudia/software/TICR/scripts/` by the full path where the
  `bucky-slurm.pl` script is located.


In slurm, we would submit the BUCKy array job with:

```bash
sbatch bucky-slurm-submit.sh
```

At the end, the array job will produce
- a `.concordance` file for every 4-taxon set
- a `.cf` file with the parsed output for that same 4-taxon set,
  in the format needed for the final CF table.

The `.cf` files can be concatenated to produce the file containing
the quartet concordance factors across all 4-taxon sets, to give to SNaQ as input:

```bash
cat *.cf > CFtable.csv
```

Alternatively, if the list of `.cf` files is not easily captured by `*.cf`
(because the list is too long for a shell command), the following julia script
can do the concatenation. Just copy-paste the commands below within a Julia session,
started from the directory that contains the `.cf` files:

```julia
files = String[] # empty vector of strings: will contain the .cf file names later
for f in filter(x -> endswith(x, ".cf"), readdir())
    push!(files,f)
end
println("found $(length(files)) cf files") # to check how many .cf output files were found
open("CFtable.csv","w") do f_out
  # write the header:
  write(f_out, "taxon1,taxon2,taxon3,taxon4,CF12_34,CF12_34_lo,CF12_34_hi,CF13_24,CF13_24_lo,CF13_24_hi,CF14_23,CF14_23_lo,CF14_23_hi,ngenes\n")
  for file in files
    @show file # to see the .cf file name: comment this out if that's too much screen output
    open(file) do f_in
        line = read(f_in, String)
        write(f_out, string(line,"\n"))
    end # closes "file" safely
  end
end # closes "CFtable.csv" safely
```
When this is done, we will have a file `CFtable.csv` containing the
quartet concordance factors, to give to SNaQ as input :smiley:
