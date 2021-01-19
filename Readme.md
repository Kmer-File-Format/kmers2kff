# kmers2kff 🧬 💻

[![License](https://img.shields.io/badge/license-AGPL--3.0-green)](https://github.com/Kmer-File-Format/kmers2kff/blob/master/LICENSE)
[![CI](https://github.com/Kmer-File-Format/kmers2kff/workflows/CI/badge.svg)](https://github.com/Kmer-File-Format/kmers2kff/actions)
[![CodeCov](https://codecov.io/gh/Kmer-File-Format/kmers2kff/branch/master/graph/badge.svg)](https://codecov.io/gh/Kmer-File-Format/kmers2kff)

A tools to convert a list of kmer with count (in csv format) as kff file with minimizer compression

## Dependencie

A [rust tool chain (> 1.47)](https://rustup.rs/)

## Install

```
cargo install --git https://github.com/Kmer-File-Format/kmers2kff.git
```

## Usage

```
mkdir tmp/
kmers2kff -i <input.csv> -o <output.kff> -k <kmer-size> -m <minimizer-size> -p tmp
rm tmp/
```

Kmer size must be lower than 65, minimizer size must be lower than kmer size.
By default csv delimiter is comma if you want use tabulation call with `-d $'\t'`.

**Warning** kmers2kff writes temporary files with a name like `{prefix}{a number}` if you run several kmers2kff at the same time be careful that their prefix is different to avoid collision between temporary files.

If `RUST_LOG` is set with `info` value kmers2kff write some log in stderr.


Complete cli:
```
kmers2kff 0.1
Pierre Marijon <pierre.marijon@hhu.de>
Convert a raw kmer count in kff format with a minimizer compression

USAGE:
    kmers2kff [OPTIONS] --input <input> --output <output> --kmer-size <k> --minimizer-size <m> --prefix <prefix>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -d, --delimiter <delimiter>    Set delimiter between kmer and count in input [default: ,]
    -i, --input <input>            Path of kmers counts in csv format
    -k, --kmer-size <k>            Kmer size
    -m, --minimizer-size <m>       Minimizer size
    -o, --output <output>          Path of the kff file
    -p, --prefix <prefix>          Prefix add before temporary file
```

## Build and contribute

```
git clone https://github.com/Kmer-File-Format/kmers2kff.git
cd kmers2kff
cargo build [--release]
```

You can create pull request but before please run this command:

```
cargo test
cargo fmt
cargo clippy -- -D warnings
cargo test
```

All this command must finish without error.

## Minimum supported Rust version

Currently the minimum supported Rust version is 1.47.0.
