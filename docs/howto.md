# Inserting data into the TICR pipeline at various stages
## To run MrBayes: You already have alignments
If you don't need to run mdl.pl because you already have aligned gene sequences which you would like to run through MrBayes, you can simply create a tarball of the Nexus files (fasta won't work at this stage) you wish to use. This command assumes that you want to use all the files ending with ".nex" in the current directory, one file per locus:

`
tar czf my-genes.tar.gz *.nex
`

Once the tarball has been successfully generated, you can then specify this file as input for [mb.pl](https://github.com/nstenz/TICR/blob/master/scripts/mb.pl) assuming you have a a valid MrBayes block located in the file "bayes.txt":

`
mb.pl my-genes.tar.gz -m bayes.txt -o my-genes-mb
`

If you get an error message like `mb.pl: Command not found`, it might be because
`mb.pl` has no execute permission or your current directory is not in your path.
An easy fix is to run this command instead:

    perl mb.pl my-genes.tar.gz -m bayes.txt -o my-genes-mb

The resulting output tarball would now be located in my-genes-mb/my-genes.mb.tar, and can be used normally with [bucky.pl](https://github.com/nstenz/TICR/blob/master/scripts/bucky.pl), that is, like this:

`
bucky.pl my-genes-mb/my-genes.mb.tar -o mygenes-bucky
`

The output, with the table of concordance factors for all sets of 4 taxa, will be in a file named `my-genes.CFs.csv` inside a directory named `mygenes-bucky`. That's the file containing the quartet concordance factors to give to SNaQ as input.

## To run BUCKy: You already have MrBayes output

### To run mbsum on the output of MrBayes for each gene

You must now run `mbsum` separately on each gene's MrBayes output. For a gene with output tree files named gene1.run1.t gene1.run2.t gene1.run3.t, and a desired burnin of 1000 trees per tree file:

`
mbsum -n 1000 -o gene1.in gene1.run1.t gene1.run2.t gene1.run3.t
`

This `mbsum` command will need to be executed for each gene. Now continue to the next section.

### To run bucky on all 4-taxon sets: you already have the mbsum output
If you have already run mbsum on each individual gene's MrBayes output, you can simply create a tarball containing all the mbsum output files. So if you had mbsum output in files named gene1.in, gene2.in, ... , gene100.in, you would want to run something similar to the following command:

`
tar czf my-genes-mbsum.tar.gz gene*.in
`

You can now use this tarball along with the -s option in [bucky.pl](https://github.com/nstenz/TICR/blob/master/scripts/bucky.pl) like this:

`
bucky.pl my-genes-mbsum.tar.gz -s -o mygenes-bucky
`

Again, if you get an error like `bucky.pl: Command not found`, run instead

    bucky.pl my-genes-mbsum.tar.gz -s -o mygenes-bucky

The output, with the table of concordance factors for all sets of 4 taxa, will be in a file named `my-genes.CFs.csv` inside directory `mygenes-bucky`. That's the file containing the quartet concordance factors to give to SNaQ as input.
