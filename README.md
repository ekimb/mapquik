<<<<<<< HEAD
`hifimap`: A fast HiFi read mapper
=========

## Rationale

=======
`hifimap`: fast HiFi read mapper using k-min-mers
=========

`hifimap` is an ultra-fast minimizer-space read mapper based on the de Bruijn graph (mdBG). It aligns long and accurate reads such as [PacBio HiFi](https://www.pacb.com/smrt-science/smrt-sequencing/hifi-reads-for-highly-accurate-long-read-sequencing/) to a reference genome.

## Rationale


>>>>>>> 40e15b2c44dd857d43c1ac45b5186b991fdd1510
## Limitations

## Installation

Clone the repository (make sure you have a working Rust environment), and run 

`cargo build --release`

<<<<<<< HEAD
## Quick start

## Overview

## Input

=======

## Quick start

```
cargo build --release
target/release/hifimap reads-0.00.fa.gz reference.fa
```

## Input

`hifimap` takes a single FASTA/FASTQ input (gzip-compressed or not) as input. Multi-line sequences, and sequences with lowercase characters, are not supported. 

If you have [seqtk](https://github.com/lh3/seqtk) installed, you can use

`seqtk seq -AU reads.unformatted.fq > reads.fa`

to format reads accordingly.

>>>>>>> 40e15b2c44dd857d43c1ac45b5186b991fdd1510
## Output data 

## Running an example

<<<<<<< HEAD
## Parameters

## Performance

=======
A sample set of reads is provided in the `example/` folder. Run

`target/release/hifimap reads-0.00.fa.gz reference.fa`


## Parameters


For further information on usage and parameters, run

`target/release/hifimap -h`

for a one-line summary of each flag, or run

`target/release/hifimap --help`

for a lengthy explanation of each flag.

## Performance


>>>>>>> 40e15b2c44dd857d43c1ac45b5186b991fdd1510
## License

`hifimap` is freely available under the [MIT License](https://opensource.org/licenses/MIT).

## Developers

* [Barış Ekim](http://people.csail.mit.edu/ekim/), supervised by [Bonnie Berger](http://people.csail.mit.edu/bab/) at the Computer Science and Artificial Intelligence Laboratory (CSAIL) at Massachusetts Institute of Technology (MIT)
<<<<<<< HEAD
* [Kristoffer Sahlin](https://sahlingroup.github.io/) at the Department of Mathematics at Stockholm University
=======
* Kristoffer Sahlin
>>>>>>> 40e15b2c44dd857d43c1ac45b5186b991fdd1510
* [Rayan Chikhi](http://rayan.chikhi.name) at the Department of Computational Biology at Institut Pasteur


## Citation

<<<<<<< HEAD
=======
For now, please cite out other mdBG article: [Minimizer-space de Bruijn graphs](https://www.biorxiv.org/content/10.1101/2021.06.09.447586v1) (2021) BiorXiv
```
@article {mdbg,
	author = {Ekim, Bar{\i}{\c s} and Berger, Bonnie and Chikhi, Rayan},
	title = {Minimizer-space de Bruijn graphs},
	year = {2021},
	doi = {10.1101/2021.06.09.447586},
	publisher = {Cold Spring Harbor Laboratory},
	journal = {bioRxiv}
}
```


>>>>>>> 40e15b2c44dd857d43c1ac45b5186b991fdd1510
## Contact

Should you have any inquiries, please contact [Barış Ekim](http://people.csail.mit.edu/ekim/) at baris [at] mit [dot] edu, [Kristoffer Sahlin](https://sahlingroup.github.io/) at ksahlin [at] math [dot] su [dot] se, or [Rayan Chikhi](http://rayan.chikhi.name) at rchikhi [at] pasteur [dot] fr.


