/* std use */
use std::io::Write;
use std::str::FromStr;

/* crate use */
use anyhow::{anyhow, Context, Result};
use clap::Clap;

/* mod declaration */
mod cli;
mod error;
mod seq2bits;

fn main() -> Result<()> {
    let params = cli::check_params(cli::Command::parse()).with_context(|| "Check parameter")?;

    // generate bucket
    let bob = create_bucket(&params.input, &params.prefix, params.k, params.m)?;

    // create kff
    let mut writer = kff::Writer::new(std::fs::File::create(params.output)?, 0b00011011, b"")?;

    writer.variables().insert("k".to_string(), params.k as u64);
    writer.variables().insert("m".to_string(), params.m as u64);
    writer.variables().insert("max".to_string(), 255);
    writer.variables().insert("data_size".to_string(), 1);

    writer.write_variables()?;

    // iterate over bucket
    for b_id in bob {
        let bucket = read_bucket(format!("{}{}", params.prefix, b_id))?;
        let mut seens = std::collections::HashSet::new();

        let mut mini_poss = Vec::new();
        let mut sequences = Vec::new();
        let mut datas = Vec::new();

        let mut fusion = Vec::new();
        for kmer in bucket.keys() {
            let mut counts: Vec<u8> = Vec::new();

            if seens.contains(kmer) {
                continue;
            }
            seens.insert(*kmer);

            let mut current: u128 = *kmer;
            while let Some((pred, nuc)) = predecessor(current, params.k, &bucket, &mut seens) {
                current = pred;
                fusion.push(nuc);
                seens.insert(current);
                counts.extend(
                    &bucket
                        .get(&current)
                        .ok_or_else(|| anyhow!("counts conversion"))?
                        .to_le_bytes(),
                );
            }

            fusion.reverse();
            counts.reverse();

            current = *kmer;
            fusion.extend(seq2bits::kmer2seq(current, params.k).bytes());
            counts.extend(
                &bucket
                    .get(&current)
                    .ok_or_else(|| anyhow!("counts conversion"))?
                    .to_le_bytes(),
            );

            while let Some((succ, nuc)) = successor(current, params.k, &bucket, &mut seens) {
                current = succ;
                fusion.push(nuc);
                seens.insert(current);
                counts.push(
                    *bucket
                        .get(&current)
                        .ok_or_else(|| anyhow!("counts conversion"))?,
                );
            }

            let mini_pos = String::from_utf8(fusion.clone())?
                .find(&seq2bits::kmer2seq(b_id, params.m))
                .unwrap() as u64;
            mini_poss.push(mini_pos);

            let mut tmp = Vec::new();
            tmp.extend(&fusion[0..mini_pos as usize]);
            tmp.extend(&fusion[(mini_pos as usize + params.m as usize)..]);

            sequences.push(tmp);
            datas.push(counts);

            fusion.clear();
        }

        writer.write_minimizer_seq_section(
            &seq2bits::kmer2seq(b_id, params.m).into_bytes(),
            &mini_poss[..],
            &sequences[..],
            &datas,
        )?;
    }

    Ok(())
}

fn predecessor(
    kmer: u128,
    k: u8,
    set: &std::collections::HashMap<u128, u8>,
    seens: &mut std::collections::HashSet<u128>,
) -> Option<(u128, u8)> {
    let sub = kmer >> 2;

    for nuc in 0..4 {
        let pred = ((nuc as u128) << ((k - 1) * 2)) ^ sub;

        if !seens.contains(&pred) && set.contains_key(&pred) {
            return Some((pred, seq2bits::bit2nuc(nuc)));
        }
    }

    None
}

fn successor(
    kmer: u128,
    k: u8,
    set: &std::collections::HashMap<u128, u8>,
    seens: &mut std::collections::HashSet<u128>,
) -> Option<(u128, u8)> {
    let mask = (2_u128.pow(k as u32 * 2) - 1) >> 2;
    let sub = (kmer & mask) << 2;

    for nuc in 0..4 {
        let succ = sub ^ nuc as u128;

        if !seens.contains(&succ) && set.contains_key(&succ) {
            return Some((succ, seq2bits::bit2nuc(nuc)));
        }
    }

    None
}

fn create_bucket(
    input: &str,
    prefix: &str,
    k: u8,
    m: u8,
) -> Result<std::collections::HashSet<u128>> {
    let mut bob = std::collections::HashSet::<u128>::new();

    let input =
        std::io::BufReader::new(std::fs::File::open(input).with_context(|| "Open input file")?);
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b',')
        .has_headers(false)
        .from_reader(input);
    let mut iter = reader.records();

    while let Some(Ok(record)) = iter.next() {
        let mut kmer = seq2bits::seq2bit(record[0].as_bytes());

        let (minimizer, _, forward) = seq2bits::get_minimizer(kmer, k, m);

        kmer = if forward {
            kmer
        } else {
            seq2bits::revcomp(kmer, k)
        };

        let mut bucket = if bob.contains(&minimizer) {
            std::fs::OpenOptions::new()
                .append(true)
                .create(true)
                .open(format!("{}{}", prefix, minimizer))
                .with_context(|| "Open bucket file")?
        } else {
            bob.insert(minimizer);

            std::fs::File::create(format!("{}{}", prefix, minimizer))
                .with_context(|| "Open bucket file")?
        };

        writeln!(bucket, "{},{}", kmer, &record[1]).with_context(|| "Write kmer")?;
    }

    Ok(bob)
}

fn read_bucket(path: String) -> Result<std::collections::HashMap<u128, u8>> {
    let mut res = std::collections::HashMap::new();

    let input = std::io::BufReader::new(std::fs::File::open(path)?);
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
