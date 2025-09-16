Small example of usage of the tigre CLI.

In this example folder we have 3 subfolder, one for each command of the CLI, inside of each subfolder we have two folder, one for multiple or single modes.

All the commands are suposed to be run from the examples folder, as their paths are relative to it.

Suppose your are inside the root of the tigre repository, then run

```bash
cd examples
```

### `clean` command

#### Multiple mode

```bash
tigre clean multiple -v --tsv clean/multiple/example_dataset.tsv --gdict clean/multiple/plants_mit.gdict --outdir clean/multiple/output
```
#### Single mode

```bash
tigre clean single -v --gff-in clean/single/NC_007982.1.gff3 --gff-out clean/single/NC_007982.1_clean.gff3
```

### `extract` command

#### Multiple mode

```bash
tigre extract multiple -v --tsv extract/multiple/example_dataset.tsv --an-column "random_col_name"
```

#### Single mode

```bash
tigre extract single -v --gff-in extract/single/NC_007982.1_clean.gff3 --gff-out extract/single/NC_007982.1_intergenic.gff3
```

### `getfasta` command

#### Multiple mode

```bash
tigre getfasta multiple -v --tsv getfasta/multiple/example_dataset.tsv
```
#### Single mode

```bash
tigre getfasta single -v --gff-in getfasta/single/NC_007982.1_intergenic.gff3 --fasta-in getfasta/single/NC_007982.1.fasta --fasta-out getfasta/single/NC_007982.1_intergenic.fasta
```

### Comparing results with expected_outcome folder

Each AN folder has an `expected_outcome` folder with the expected results of each command. You can use `diff` to compare your results with the expected ones. For example, to compare the results of the `clean` command in multiple mode, of the AN NC_035963.1, you can run:

```bash
diff -r clean/multiple/NC_035963.1/expected_outcome/NC_035963.1_clean.gff3 clean/multiple/output/NC_035963.1/NC_035963.1_clean.gff3
```