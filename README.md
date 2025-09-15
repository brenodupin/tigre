<div align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://github.com/brenodupin/tigre/releases/download/v0.1.0-alpha/tigre_logo_dark_mode.png">
    <source media="(prefers-color-scheme: light)" srcset="https://github.com/brenodupin/tigre/releases/download/v0.1.0-alpha/tigre_logo_light_mode.png">
    <img src="https://github.com/brenodupin/tigre/releases/download/v0.1.0-alpha/tigre_logo_light_mode.png" width="50%" alt="TIGRE Logo">
  </picture>

![Build Status](https://img.shields.io/badge/tests-in_development-yellow)
[![Checked with mypy](https://www.mypy-lang.org/static/mypy_badge.svg)](https://mypy-lang.org/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Linting: Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/charliermarsh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![License](https://img.shields.io/badge/license-MIT-purple)](https://github.com/brenodupin/gdt/blob/master/LICENSE)
</div>

# TIGRE - Tool for InterGenic Region Extraction

TIGRE is a Python CLI tool that extracts intergenic regions from GFF3 file(s). It is designed to handle overlapping features and provides options for customizing the extraction process. It has 3 main commands: `clean`, `extract`, and `getfasta`, which cleans and prepares GFF3 files, extracts intergenic regions (or just regions between features), and retrieves sequences from a FASTA file, respectively.

In `clean` command, TIGRE optionally supports the use of a [GDT](https://github.com/brenodupin/gdt) .gdict file to standardize feature names.

##  Features

- **Clean & Standardize**: Process and prepares GFF3 files with optional support for [GDT](https://github.com/brenodupin/gdt) to standardize gene names.
- **Extract Regions**: Extract intergenic regions with customizable feature filtering.
- **Generate Sequences**: Retrieve FASTA sequences from extracted regions.
- **Parallel Processing**: On demand multi-thread/multi-process execution for batch operations.

## Table of Contents
 - [Quick Start](#quick-start)
   - [Installation](#installation)
     - [Dependencies](#dependencies)

 - [Commands Overview](#commands-overview)
   - [`clean`](#1-clean---prepare-gff3-files)
   - [`extract`](#2-extract---extract-intergenic-regions)
   - [`getfasta`](#3-getfasta---generate-sequences)
 
 - [Usage Modes](#usage-modes)
   - [Single Mode](#single-mode)
   - [Multiple Mode](#multiple-mode)
     - [TSV format](#tsv-format)
     - [Input organization requirement](#input-organization-requirement)

 - [Detailed Usage](#detailed-usage)
   - [Global Options](#global-options)
   - [`tigre clean`](#tigre-clean)
   - [`tigre extract`](#tigre-extract)
   - [`tigre getfasta`](#tigre-getfasta)
   - [Example Usage](#example-usage)

## Quick Start

### Installation

```bash
# Basic installation
pip install tigre

# With optional dependencies
pip install tigre[bio]    # for biopython (getfasta command)
pip install tigre[gdt]    # for gdt (GDT support in clean command)
pip install tigre[all]    # for all optional dependencies
```

### Dependencies

**Dependencies:**
- [Python](https://www.python.org/) `(>=3.10)`
- [pandas](https://pandas.pydata.org/) `(>=1.5.3,<3.0.0)`

**Optional Dependencies:**
- [biopython](https://biopython.org) `(>=1.80)`
- [gdt](https://github.com/brenodupin/gdt) `(>=1.0.1)`

## Commands Overview

TIGRE has the following three commands:

### 1. `clean` - Prepare GFF3 Files
Processes GFF3 files by sorting, solving overlaps and format attributes for further steps, optionally using a GDT file for gene name standardization.

```bash
# Single file
tigre clean single --gff-in input.gff3 --gff-out clean.gff3

# Multiple files from TSV
tigre clean multiple --tsv samples.tsv --gff-in-suffix "" --gff-out-suffix "_clean"
```

### 2. `extract` - Extract Intergenic Regions  
Identifies and extracts regions between features in processed GFF3s.

```bash
# Single file
tigre extract single --gff-in clean.gff3 --gff-out intergenic.gff3

# Multiple files
tigre extract multiple --tsv samples.tsv --gff-in-suffix "_clean" --gff-out-suffix "_intergenic"
```

### 3. `getfasta` - Generate Sequences
Retrieves FASTA sequences of extracted regions.

```bash
# Single file
tigre getfasta single --gff-in intergenic.gff3 --fasta-in genome.fasta --fasta-out sequences.fasta

# Multiple files  
tigre getfasta multiple --tsv samples.tsv --gff-in-suffix "_intergenic" --fasta-out-suffix "_intergenic"
```

## Usage Modes

### Single Mode
Processes individual files by specifying direct file paths:

```bash
tigre <command> single --gff-in input.gff3 --gff-out output.gff3 [options]
```
<div align="center">

| Flag | Description | `clean` | `extract` | `getfasta` |
| :---: | :---: | :---: | :---: | :---: |
| `--gff-in`    | Input GFF3 file path   | ✅ | ✅ | ✅ |
| `--gff-out`   | Output GFF3 file path  | ✅ | ✅ | ❌ |
| `--fasta-in`  | Input FASTA file path  | ❌ | ❌ | ✅ |
| `--fasta-out` | Output FASTA file path | ❌ | ❌ | ✅ |

</div>

### Multiple Mode

The `multiple` mode allows TIGRE to process batches of files indexed by the input TSV file containing accession numbers.
For each accession, TIGRE builds the expected file paths using suffixes and extensions.
Each command has its own default suffixes and extensions, which are applied automatically unless the user overrides them.

#### TSV format:
The TSV must have a column with accession numbers.
By default, the column name is `AN`, but it can be changed with the `--an-column` flag.
Example:
```
AN          other_columns
NC_123456.1   ...
NC_999999.1   ...
```

#### Input organization requirement:
The TSV file and the accession-number folders must be located in the **same directory**.
Each folder must be named exactly after the corresponding accession number from the TSV.

The file paths are constructed as:
```
<AN>/<AN><suffix><extension>
```

Using the example above, the corresponding folder structure should look like this:
```
your_project/
├──example.tsv
├── NC_123456.1/
│   ├── NC_123456.1.gff3
│   └── NC_123456.1.fasta (if applicable)
└── NC_999999.1/
    ├── NC_999999.1.gff3
    └── NC_999999.1.fasta (if applicable)
```

<div align="center">

File extension default values for each command

| Flag | Default |`clean` | `extract` | `getfasta` |
| :---: | :---: | :---: | :---: | :---: |
| `--gff-in-ext`    | ".gff3"  | ✅ | ✅ | ✅ |
| `--gff-out-ext`   | ".gff3"  |✅ | ✅ | ❌ |
| `--fasta-in-ext`  | ".fasta" |❌ | ❌ | ✅ |
| `--fasta-out-ext` | ".fasta" |❌ | ❌ | ✅ |

File suffix default values for each command

| Flag | `clean` | `extract` | `getfasta` |
| :---: | :---: | :---: | :---: |
| `--gff-in-suffix`    | "" | "_clean" | "_intergenic" |
| `--gff-out-suffix`   | "_clean" | "_intergenic" | ❌ |
| `--fasta-in-suffix`  | ❌ | ❌ | "" |
| `--fasta-out-suffix` | ❌ | ❌ | "_intergenic" |

</div>

## Detailed Usage

### Global Options
Available for all commands:

#### Logging Options

TIGRE has extensive logging capabilities, allowing users to customize the verbosity and destination of log messages, by default logs are saved to `tigre_<command>-<timestamp>.log` in the current working directory.

- `-v, --verbose`: Increase verbosity level. Use multiple times for more verbose output: -v (INFO console, TRACE file), -vv (DEBUG console, TRACE file), -vvv (TRACE console, TRACE file). Default: INFO console + DEBUG file.
- `--log PATH`: Custom log file path.
- `--no-log-file`: Disable file logging.
- `--quiet`: Suppress console output.

#### Output Options
- `--overwrite`: Overwrite existing output files. By default, TIGRE will not even execute if one of the output files already exists, to prevent accidental overwrites.

#### Multiple Mode Options
- `--tsv PATH`: TSV file with accession numbers and other columns.
- `--an-column STR`: Column name for accession numbers in TSV (default: `"AN"`).
- `--workers INT`: Number of workers for parallel processing (default: `0`, which uses all available CPU cores).

### `tigre clean`

The `clean` command processes GFF3 files by retaining features that match the `--query-string` criteria, resolving overlapping annotations, and optionally normalizing gene names using GDT.

**Options:**
- `--gdict PATH`: GDT .gdict file for standardizing gene names.
- `--query-string STR`: pandas query string for feature filtering (default: `"type in ('gene', 'tRNA', 'rRNA', 'region')"`).
- `--keep-orfs`: Keep ORF sequences in output.
- `--extended-filtering`: Enable extended filtering to more aggressively remove ORFs and their related features. This may remove more features than intended in some cases. Default: False.
- `--overwrite`: Overwrite existing output files. By default, TIGRE will not even execute if one of the output files already exists, to prevent accidental overwrites.

**Single Mode Options:**
- `--gff-in PATH`: Input GFF3 file path.
- `--gff-out PATH`: Output GFF3 file path.

**Multiple Mode Options:**
- `--tsv PATH`: TSV file with accession numbers and other columns.
- `--an-column STR`: Column name for accession numbers in TSV (default: `"AN"`).
- `--workers INT`: Number of workers for parallel processing (default: `0`, which uses all available CPU cores).
- `--gff-in-suffix STR`: Suffix for input GFF3 files (default: `""`).
- `--gff-in-ext STR`: Extension for input GFF3 files (default: `".gff3"`).
- `--gff-out-suffix STR`: Suffix for output GFF3 files (default: `"_clean"`).
- `--gff-out-ext STR`: Extension for output GFF3 files (default: `".gff3"`).

When `clean` is called with a GDT file in `multiple` mode, TIGRE automatically uses a server-client architecture to optimize memory usage. A dedicated dictionary server process loads and queries the gene dictionary (gdict), while worker processes communicate with it via inter-process queues and pipes to perform gene name lookups. This design prevents each worker from loading its own copy of the potentially large gdict file, significantly reducing overall memory overhead.

#### Processing Details
TIGRE handles complex genomic scenarios including overlapping features and circular genome boundaries. For detailed information about overlap resolution, circular genome handling, and the underlying algorithms, see [PROCESSING.md](https://github.com/brenodupin/tigre/blob/master/PROCESSING.md#tigre-clean)

### `tigre extract`

The `extract` command identifies and extracts intergenic regions (spaces between features) from cleaned GFF3 files, creating new annotations for these intervening sequences.

**Options:**
- `--add-region`: Add region line to the output GFF3 file.
- `--feature-type STR`: Feature type name for intergenic regions (default: `"intergenic_region"`).
- `--overwrite`: Overwrite existing output files. By default, TIGRE will not even execute if one of the output files already exists, to prevent accidental overwrites.

**Single Mode Options:**
- `--gff-in PATH`: Input GFF3 file path.
- `--gff-out PATH`: Output GFF3 file path.

**Multiple Mode Options:**
- `--tsv PATH`: TSV file with accession numbers and other columns.
- `--an-column STR`: Column name for accession numbers in TSV (default: `"AN"`).
- `--workers INT`: Number of workers for parallel processing (default: `0`, which uses all available CPU cores).
- `--gff-in-suffix STR`: Suffix for input GFF3 files (default: `"_clean"`).
- `--gff-int-ext STR`: Extension for input GFF3 files (default: `".gff3"`).
- `--gff-out-suffix STR`: Suffix for output GFF3 files (default: `"_intergenic"`).
- `--gff-out-ext STR`: Extension for output GFF3 files (default: `".gff3"`).

#### Intergenic Region Creation

TIGRE identifies and extracts intergenic regions (spaces between annotated features) from GFF3 files, creating annotations for these intervening sequences. For detailed information about the algorithm and how circular genomes are handled, see [PROCESSING.md](https://github.com/brenodupin/tigre/blob/master/PROCESSING.md#tigre-extract).


### `tigre getfasta`

`tigre getfasta` retrieves the nucleotide sequences of the intergenic regions extracted by `tigre extract`. The nucleotide sequences are saved in a multi-FASTA file (one per Accession Number).

**Requirements:**
- Requires [biopython](https://biopython.org) (`pip install tigre[bio]`)

**Options:**
- `--bedtools-compatible`: Use 0-based indexing in FASTA headers for bedtools compatibility
- `--overwrite`: Overwrite existing output files. By default, TIGRE will not even execute if one of the output files already exists, to prevent accidental overwrites.

**Single Mode Options:**
- `--gff-in PATH`: Input GFF3 file path.
- `--fasta-in PATH`: Input FASTA file path.
- `--fasta-out PATH`: Output FASTA file path.

**Multiple Mode Options:**
- `--tsv PATH`: TSV file with accession numbers and other columns.
- `--an-column STR`: Column name for accession numbers in TSV (default: `"AN"`).
- `--workers INT`: Number of workers for parallel processing (default: `0`, which uses all available CPU cores).
- `--gff-in-suffix STR`: Suffix for input GFF3 files (default: `"_intergenic"`).
- `--gff-in-ext STR`: Extension for input GFF3 files (default: `".gff3"`).
- `--fasta-in-suffix STR`: Suffix for input FASTA files (default: `""`).
- `--fasta-in-ext STR`: Extension for input FASTA files (default: `".fasta"`).
- `--fasta-out-suffix STR`: Suffix for output FASTA files (default: `"_intergenic"`).
- `--fasta-out-ext STR`: Extension for output FASTA files (default: `".fasta"`).

#### FASTA Header Format
The FASTA headers are formatted as follows:
`>type::seqid:start-end`

Where:
- `type`: Type of the feature (e.g., `intergenic_region`, `intergenic_region_merged`).
- `seqid`: `seqid` column from the GFF3 file (mostly the accession number).
- `start`: Start position of the intergenic region (1-based index).
- `end`: End position of the intergenic region (1-based index).

By default, TIGRE uses 1-based indexing in sequence headers matching GFF3 conventions. When `--bedtools-compatible` is enabled, headers use 0-based indexing for seamless integration with bedtools workflows, while the actual sequences remain unchanged.

### `tigre combine`

`tigre combine` merges two GFF3 files into a single GFF3 file. It follows the same input/output conventions as the other commands, when calcutaling paths.

**Options:**
- `--overwrite`: Overwrite existing output files. By default, TIGRE will not even execute if one of the output files already exists, to prevent accidental overwrites.

**Single Mode Options:**
- `--gff-in-1 PATH`: First input GFF3 file path.
- `--gff-in-2 PATH`: Second input GFF3 file path.
- `--gff-out PATH`: Output GFF3 file path.

**Multiple Mode Options:**
- `--tsv PATH`: TSV file with accession numbers and other columns.
- `--an-column STR`: Column name for accession numbers in TSV (default: `"AN"`).
- `--workers INT`: Number of workers for parallel processing (default: `0`, which uses all available CPU cores).
- `--gff-in-suffix-1 STR`: Suffix for first input GFF3 files (default: `""`).
- `--gff-in-ext-1 STR`: Extension for first input GFF3 files (default: `".gff3"`).
- `--gff-in-suffix-2 STR`: Suffix for second input GFF3 files (default: `"_intergenic"`).
- `--gff-in-ext-2 STR`: Extension for second input GFF3 files (default: `".gff3"`).
- `--gff-out-suffix STR`: Suffix for output GFF3 files (default: `"_combined"`).
- `--gff-out-ext STR`: Extension for output GFF3 files (default: `".gff3"`).

This combination follows the procedure:
1. Read both GFF3 files as DataFrames.
2. Copy region line (first row of GFF3) from the `gff-in-1` file.
3. Remove region lines (`type == 'region'`) from both DataFrames.
4. Check for duplicated IDs between both files and log a warning if any is found.
5. Concatenate both DataFrames, with the region line extracted earlier as the first row.
6. Write the combined DataFrame to a `gff-out` path, preserving the header from the `gff-in-1` file.


### Example Usage

The input files used in the examples below can be found in the `examples` folder of this repository, and should be run **inside the `example` folder**.

Single mode example:
```bash
# Clean a GFF3 file
tigre clean single -vv --gff-in clean/single/NC_007982.1.gff3 --gff-out clean/single/NC_007982.1_clean.gff3 --query-string "type in ('gene', 'tRNA', 'rRNA', region')" --overwrite

# Extract intergenic regions from the cleaned GFF3 file
tigre extract single -vv --gff-in extract/single/NC_007982.1_clean.gff3 --gff-out extract/single/NC_007982.1_intergenic.gff3 --feature-type "intergenic_region"

# Get FASTA sequences of the intergenic regions
tigre getfasta single -vv --gff-in getfasta/single/NC_007982.1_intergenic.gff3 --fasta-in getfasta/single/NC_007982.1.fasta --fasta-out getfasta/single/NC_007982.1_intergenic.fasta --bedtools-compatible
```

Multiple mode example:
```bash
# Clean multiple GFF3 files using a TSV file, and a GDT file for gene name standardization
tigre clean multiple -vv --tsv clean/multiple/example_dataset.tsv --gdict clean/multiple/plants_mit.gdict --server 

# Extract intergenic regions from the cleaned GFF3 files
tigre extract multiple -vv --tsv extract/multiple/example_dataset.tsv --gff-in-suffix "_clean" --feature-type "intergenic_region"

# Get FASTA sequences of the intergenic regions
tigre getfasta multiple -vv --tsv getfasta/multiple/example_dataset.tsv --bedtools-compatible --gff-in-suffix "_intergenic"
```