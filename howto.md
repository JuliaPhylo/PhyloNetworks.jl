# Inserting data into the pipeline at various stages
## Before running MrBayes / You already have alignments
If you don't need to run mb.pl, and you already have aligned gene sequences which you would like to run through MrBayes, you can simply create a tarball of the Nexus files (fasta won't work at this stage) you wish to use:

`
tar czf my-genes.tar.gz *.nex
`

Once the tarball has been successfully generated, you can then specify this file as input for mb.pl assuming you have a a valid MrBayes block located in the file "bayes.txt":

`
mb.pl my-genes.tar.gz -m bayes.txt -o my-genes-mb
`

The resulting output tarball would now be located in my-genes-mb/my-genes.mb.tar, and can be used normally with bucky.pl.

## Before running BUCKy / You already have MrBayes output

### Before running mbsum

You must now run mbsum separately on each gene's MrBayes output, i.e. for a gene with output tree files named gene1.run1.t gene1.run2.t gene1.run3.t, and a desired burnin of 1000 trees per tree file:

`
mbsum gene1.run1.t gene1.run2.t gene1.run3.t -n 1000
`

Now continue to the next section.

### After running mbsum
If you have already run mbsum on each individual gene's MrBayes output, you can simply create a tarball containing each mbsum output file. So if you had mbsum output in three files named gene1.in, gene2.in, and gene3.in, you would want to run somethingn similar to the following command:

`
tar czf my-genes-mbsum.tar.gz gene1.in gene2.in gene3.in
`

You can now use this tarball along with the -s option in bucky.pl like so:

`
bucky.pl my-genes-mbsum.tar.gz -s -o mygenes-bucky
`

