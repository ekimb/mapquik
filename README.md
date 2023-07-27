
`mapquik`: Efficient mapping of accurate long reads in minimizer space
=========

`mapquik` is an ultra-fast read mapper based on $k$-min-mers (matches of $k$ consecutively-sampled minimizers). It aligns long and accurate reads such as [PacBio HiFi](https://www.pacb.com/smrt-science/smrt-sequencing/hifi-reads-for-highly-accurate-long-read-sequencing/) to a reference genome.

## Rationale

The underlying seed constructs ($k$-mers) in state-of-the-art long-read mappers are tailored to noisy reads, and small seed sizes induce longer computation times due to multiple potential mapping locations. Recent advances in short-read alignment methods have demonstrated that [98\% of many organisms' genomes are non-repetitive and can be uniquely aligned to with longer seeds](https://peerj.com/articles/9338/). Therefore, we explore the use of longer, non-exact seeds ($k$-min-mers) in accurate long reads. See [our manuscript](https://genome.cshlp.org/content/early/2023/06/29/gr.277679.123) for details.

## Limitations
The mapping performance of `mapquik` degrades markedly when identity between reads and the reference is lower than $97$\%, and less than $1$\% of the reads are mapped at $Q60$ for identities below $93$\%. Therefore, `mapquik` is not suitable for mapping PacBio CLR reads, and potentially also Oxford Nanopore reads until base-calling consistently reaches identity levels above $98$\%. 


## Installation

Pre-requisites: [A working Rust environment](https://rustup.rs/).

Clone the repository, and run 

```
rustup install nightly
cargo +nightly build --release
```

The nightly version of `cargo` is required because `mapquik` uses experimental language features (such as SIMD and intrinsics).

## Quick start

`target/release/mapquik <reads.fq> --reference <reference.fa>`

## Input

`mapquik` takes a single FASTA/FASTQ input (`gzip`-compressed or not) as input. Multi-line sequences, and sequences with lowercase characters, are not supported. 

## Output

The output of `mapquik` is a regular PAF file.

## Running an example

An example reference genome, and a script to simulate reads using `pbsim` are provided in the `example/` folder. To run `mapquik` on a small set of 100 reads, type:

`cd example && bash run_ecoli.sh`

which will run both `mapquik` and `minimap2` on 100 simulated reads, and return the output of `paftools.js mapeval` on both PAF files.

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

mapquik significantly accelerates the seeding and chaining steps for both the human and maize genomes with $>96$\% sensitivity and near-perfect specificity. On the human genome, for both real and simulated reads, mapquik achieves a $37\times$ speed-up over `minimap2`, and on the maize genome, a $410\times$ speed-up over `minimap2`. 

`mapquik` indexing is $9\times$ faster than `minimap2`, which is of independent interest.

## License

`mapquik` is freely available under the [MIT License](https://opensource.org/licenses/MIT).

## Developers

* [Barış Ekim](http://people.csail.mit.edu/ekim/), supervised by [Bonnie Berger](http://people.csail.mit.edu/bab/) at the Computer Science and Artificial Intelligence Laboratory (CSAIL) at Massachusetts Institute of Technology (MIT)
* [Rayan Chikhi](http://rayan.chikhi.name) at the Department of Computational Biology at Institut Pasteur


## Citation

```
@article{mapquik,
  title={Efficient mapping of accurate long reads in minimizer space with mapquik},
  author={Ekim, Bar{\i}{\c{s}} and Sahlin, Kristoffer and Medvedev, Paul and Berger, Bonnie and Chikhi, Rayan},
  journal={Genome Research},
  pages={gr--277679},
  year={2023},
  publisher={Cold Spring Harbor Laboratory}
}
```

## Contact

Should you have any inquiries, please contact [Barış Ekim](http://people.csail.mit.edu/ekim/) at baris [at] mit [dot] edu, or [Rayan Chikhi](http://rayan.chikhi.name) at rchikhi [at] pasteur [dot] fr.
