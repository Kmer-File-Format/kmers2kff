/* std use */
use std::io::Write;
use std::str::FromStr;

/* crate use */
use anyhow::{Context, Result};

/* local use */
use crate::seq2bits;

pub fn create(
    input: &str,
    prefix: &str,
    k: u8,
    m: u8,
    delimiter: u8,
) -> Result<rustc_hash::FxHashSet<u128>> {
    let mut bob = rustc_hash::FxHashSet::<u128>::default();

    let input =
        std::io::BufReader::new(std::fs::File::open(input).with_context(|| "Open input file")?);
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(delimiter)
        .has_headers(false)
        .from_reader(input);
    let mut iter = reader.records();

    while let Some(Ok(record)) = iter.next() {
        let mut kmer = seq2bits::seq2bit(record[0].as_bytes());
        let count = u8::from_str(&record[1])?;

        let (minimizer, _, forward) = seq2bits::get_minimizer(kmer, k, m);
        kmer = if forward {
            kmer
        } else {
            seq2bits::revcomp(kmer, k)
        };

        if seq2bits::multiple_mini(kmer, minimizer, k, m) {
            write_kmer(kmer, count, &format!("{}multiple", prefix))?;
        } else {
            write_kmer(kmer, count, &format!("{}{}", prefix, minimizer))?;
        }

        bob.insert(minimizer);
    }

    Ok(bob)
}

pub fn write_kmer(kmer: u128, count: u8, filename: &str) -> Result<()> {
    let mut bucket = if std::path::Path::new(filename).exists() {
        std::fs::OpenOptions::new()
            .append(true)
            .create(true)
            .open(filename)
            .with_context(|| "Open bucket file")?
    } else {
        std::fs::File::create(filename).with_context(|| "Open bucket file")?
    };

    writeln!(bucket, "{},{}", kmer, count).with_context(|| "Write kmer")?;

    Ok(())
}

pub fn read(path: &str) -> Result<rustc_hash::FxHashMap<u128, u8>> {
    let mut res = rustc_hash::FxHashMap::default();

    let input = std::io::BufReader::new(
        std::fs::File::open(path).with_context(|| format!("Read bucket {}", path))?,
    );
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b',')
        .has_headers(false)
        .from_reader(input);
    let mut iter = reader.records();

    while let Some(Ok(record)) = iter.next() {
        res.insert(u128::from_str(&record[0])?, u8::from_str(&record[1])?);
    }

    Ok(res)
}
