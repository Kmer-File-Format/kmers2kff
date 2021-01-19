/* std use */
use std::io::Read;
use std::process::{Command, Stdio};
use std::str::FromStr;

use kff::seq2bits::Bits2Nuc;

fn seq2bit(subseq: &[u8]) -> u128 {
    let mut kmer: u128 = 0;

    for n in subseq {
        kmer <<= 2;
        kmer |= nuc2bit(*n);
    }

    kmer
}

fn canonical(kmer: u128, k: u8) -> (u128, bool) {
    let rev = revcomp(kmer, k);

    if kmer < rev {
        (kmer, true)
    } else {
        (rev, false)
    }
}

fn nuc2bit(nuc: u8) -> u128 {
    (nuc as u128 >> 1) & 0b11
}

fn revcomp(kmer: u128, k: u8) -> u128 {
    rev(speed_comp(kmer), k)
}

fn speed_comp(kmer: u128) -> u128 {
    kmer ^ 0xAAAA_AAAA_AAAA_AAAA_AAAA_AAAA_AAAA_AAAA
}

fn rev(mut kmer: u128, k: u8) -> u128 {
    // Thank to needtail people ! :)
    kmer = (kmer >> 2 & 0x3333_3333_3333_3333_3333_3333_3333_3333)
        | (kmer & 0x3333_3333_3333_3333_3333_3333_3333_3333) << 2;
    kmer = (kmer >> 4 & 0x0F0F_0F0F_0F0F_0F0F_0F0F_0F0F_0F0F_0F0F)
        | (kmer & 0x0F0F_0F0F_0F0F_0F0F_0F0F_0F0F_0F0F_0F0F) << 4;
    kmer = (kmer >> 8 & 0x00FF_00FF_00FF_00FF_00FF_00FF_00FF_00FF)
        | (kmer & 0x00FF_00FF_00FF_00FF_00FF_00FF_00FF_00FF) << 8;
    kmer = (kmer >> 16 & 0x0000_FFFF_0000_FFFF_0000_FFFF_0000_FFFF)
        | (kmer & 0x0000_FFFF_0000_FFFF_0000_FFFF_0000_FFFF) << 16;
    kmer = (kmer >> 32 & 0x0000_0000_FFFF_FFFF_0000_0000_FFFF_FFFF)
        | (kmer & 0x0000_0000_FFFF_FFFF_0000_0000_FFFF_FFFF) << 32;
    kmer = (kmer >> 64 & 0x0000_0000_0000_0000_FFFF_FFFF_FFFF_FFFF)
        | (kmer & 0x0000_0000_0000_0000_FFFF_FFFF_FFFF_FFFF) << 64;

    kmer >> (128 - k * 2)
}

fn read_kmer_list(path: &str) -> Vec<(u128, u8)> {
    let input = std::io::BufReader::new(std::fs::File::open(path).unwrap());

    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b',')
        .has_headers(false)
        .from_reader(input);
    let mut iter = reader.records();

    let mut res = Vec::new();

    while let Some(Ok(record)) = iter.next() {
        let cano = canonical(seq2bit(record[0].as_bytes()), 11).0;

        res.push((cano, u8::from_str(&record[1]).unwrap()));
    }

    res.sort();

    res
}

fn read_kff(path: &str) -> Vec<(u128, u8)> {
    let mut input = std::io::BufReader::new(std::fs::File::open(path).unwrap());

    let mut reader = kff::Reader::new(&mut input).unwrap();
    let rev_encoding = reader.rev_encoding();

    let mut res = Vec::new();

    while let Ok(section) = reader.next_section() {
        let mut it = section.into_iter();
        while let Some(Ok(kmer)) = it.next() {
            let cano = canonical(seq2bit(&kmer.seq().into_nuc(rev_encoding)), 11).0;

            res.push((cano, kmer.data()[0]));
        }
    }

    res.sort();

    res
}

#[test]
fn all_kmer_is_present() {
    let mut child = Command::new("./target/debug/kmers2kff")
        .args(&[
            "-i",
            "tests/data/test.csv",
            "-o",
            "tests/test.kff",
            "-k",
            "11",
            "-m",
            "6",
            "-p",
            "tests",
        ])
        .stderr(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .expect("Couldn't create kmers2kff subprocess");

    if !child.wait().expect("Error durring kmers2kff run").success() {
        let mut stdout = String::new();
        let mut stderr = String::new();

        child.stdout.unwrap().read_to_string(&mut stdout).unwrap();
        child.stderr.unwrap().read_to_string(&mut stderr).unwrap();

        println!("stdout: {}", stdout);
        println!("stderr: {}", stderr);
        panic!();
    }

    assert_eq!(
        read_kff("tests/test.kff"),
        read_kmer_list("tests/data/test.csv")
    );
}
