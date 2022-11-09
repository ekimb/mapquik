
`mapquik`: Efficient low-divergence mapping of long reads in minimizer space
=========

`mapquik` is an ultra-fast read mapper based on $k$-min-mers (matches of $k$ consecutively-sampled minimizers). It aligns long and accurate reads such as [PacBio HiFi](https://www.pacb.com/smrt-science/smrt-sequencing/hifi-reads-for-highly-accurate-long-read-sequencing/) to a reference genome.

## Rationale

## Limitations

## Installation

Pre-requisites: a working Rust environment (https://rustup.rs/), `llvm`, and `clang` (for `libwfa`).

Clone the repository, and run 

```
rustup install nightly
export CARGO_NET_GIT_FETCH_WITH_CLI=true
cargo +nightly build --release
```

The nightly version of `cargo` is required because `mapquik` uses experimental language features (such as SIMD and intrinsics). The `export` command is needed because `libwfa` has WFA2 as a `git` submodule with a SSH url. It fixes the error "`failed to authenticate when downloading repository`".

## Quick start

`target/release/mapquik <reads.fq> --reference <reference.fa>`

## Input

`mapquik` takes a single FASTA/FASTQ input (`gzip`-compressed or not) as input. Multi-line sequences, and sequences with lowercase characters, are not supported. 

## Output

The output of `mapquik` is a regular PAF file.

## Running an example

An example reference genome, and a script to simulate reads using `pbsim` are provided in the `example/` folder. To run `mapquik` on a small set of 100 reads, type:

`cd example && bash run_ecoli.sh`

which will run both `mapquik` and `minimap2` on the simulate reads, and return the output of `paftools.js mapeval` on both PAF files.

To simulate a larger set of reads using pbsim and map, type:

`bash simulate_pbsim.sh && bash run_ecoli_full.sh`


## Parameters

For further information on usage and parameters, run

`target/release/mapquik -h`

for a one-line summary of each flag, or run

`target/release/mapquik --help`

for a lengthy explanation of each flag.

## Data Availability

All scripts used to generate the figures and tables in the paper can be found in the `experiments/` folder. Specifically, the `simulate_chm13.sh` and `simulate_maize.sh` scripts can be used similarly to simulate reads. 

In order to obtain and map DeepConsensus reads, first run

```
wget https://storage.googleapis.com/brain-genomics-public/research/deepconsensus/data/v0.3/assembly_analysis/fastqs/HG002_24kb_2SMRT_cells.dc.v0.3.q20.fastq.gz
gunzip -c HG002_24kb_2SMRT_cells.dc.v0.3.q20.fastq.gz | grep -v TOTAL > dc.hg002.fastq
```

and map to a reference genome `reference.fa` in your directory with `mapquik` using

`target/release/mapquik dc.hg002.fastq --reference reference.fa -p mapquik-dc`

## Performance

## License

`mapquik` is freely available under the [MIT License](https://opensource.org/licenses/MIT).

## Developers

* [Barış Ekim](http://people.csail.mit.edu/ekim/), supervised by [Bonnie Berger](http://people.csail.mit.edu/bab/) at the Computer Science and Artificial Intelligence Laboratory (CSAIL) at Massachusetts Institute of Technology (MIT)
* [Rayan Chikhi](http://rayan.chikhi.name) at the Department of Computational Biology at Institut Pasteur


## Citation

`mapquik` is not yet published. For now, please cite our original mdBG article: [Minimizer-space de Bruijn graphs: Whole-genome assembly of long reads in minutes on a personal computer](https://www.sciencedirect.com/science/article/pii/S240547122100332X) (2021).

```
@article {mdbg,
	author = {Ekim, Bar{\i}{\c s} and Berger, Bonnie and Chikhi, Rayan},
	title = {Minimizer-space de Bruijn graphs: Whole-genome assembly of long reads in minutes on a personal computer},
	year = {2021},
	doi = {10.1016/j.cels.2021.08.009},
	journal = {Cell Systems}
	volume={12},
  	number={10},
  	pages={958--968},
  	year={2021},
  	publisher={Elsevier}
}
```

## Contact

Should you have any inquiries, please contact [Barış Ekim](http://people.csail.mit.edu/ekim/) at baris [at] mit [dot] edu, or [Rayan Chikhi](http://rayan.chikhi.name) at rchikhi [at] pasteur [dot] fr.
