/* crate use */
use anyhow::{anyhow, Context, Result};
use clap::Clap;

/* mod declaration */
mod bucket;
mod cli;
mod error;
mod seq2bits;

fn main() -> Result<()> {
    env_logger::init();

    let params = cli::check_params(cli::Command::parse()).with_context(|| "Check parameter")?;

    log::info!("Start of bucket creation");

    // generate bucket
    let bob = bucket::create(
        &params.input,
        &params.prefix,
        params.k,
        params.m,
        params.delimiter as u8,
    )?;

    log::info!("End of bucket creation");

    // create kff
    let mut writer = kff::Writer::new(std::fs::File::create(params.output)?, 0b00011011, b"")?;

    writer.variables().insert("k".to_string(), params.k as u64);
    writer.variables().insert("m".to_string(), params.m as u64);
    writer.variables().insert("max".to_string(), 255);
    writer.variables().insert("data_size".to_string(), 1);

    writer.write_variables()?;

    // iterate over bucket
    for b_id in bob.iter() {
        log::info!("Compress bucket {}", b_id);

        let bucket = bucket::read(&format!("{}{}", params.prefix, b_id))?;
        let mut seens = rustc_hash::FxHashSet::default();

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
                .find(&seq2bits::kmer2seq(*b_id, params.m))
                .unwrap() as u64;
            mini_poss.push(mini_pos);

            let mut tmp = Vec::new();
            tmp.extend(&fusion[0..mini_pos as usize]);
            tmp.extend(&fusion[(mini_pos as usize + params.m as usize)..]);

            sequences.push(tmp);
            datas.push(counts);

            fusion.clear();
        }

        log::info!("Write bucket {}", b_id);
        writer.write_minimizer_seq_section(
            &seq2bits::kmer2seq(*b_id, params.m).into_bytes(),
            &mini_poss[..],
            &sequences[..],
            &datas,
        )?;
    }

    log::info!("Write bucket of multiple minimizer");
    if std::path::Path::new(&format!("{}multiple", params.prefix)).exists() {
        let bucket = bucket::read(&format!("{}multiple", params.prefix))?;
        let mut sequences = Vec::new();
        let mut datas = Vec::new();

        for kmer in bucket.keys() {
            sequences.push(seq2bits::kmer2seq(*kmer, params.k));
            datas.push([*bucket
                .get(&kmer)
                .ok_or_else(|| anyhow!("counts not present"))?])
        }

        writer.write_raw_seq_section(&sequences[..], &datas[..])?;
    }
    clean_temp_file(bob, &params.prefix).with_context(|| "clean temporary file")?;

    Ok(())
}

fn clean_temp_file(minis: rustc_hash::FxHashSet<u128>, prefix: &str) -> Result<()> {
    for mini in minis {
        std::fs::remove_file(format!("{}{}", prefix, mini))?;
    }

    if std::path::Path::new(&format!("{}multiple", prefix)).exists() {
        std::fs::remove_file(format!("{}multiple", prefix))?;
    }

    Ok(())
}

fn predecessor(
    kmer: u128,
    k: u8,
    set: &rustc_hash::FxHashMap<u128, u8>,
    seens: &mut rustc_hash::FxHashSet<u128>,
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
    set: &rustc_hash::FxHashMap<u128, u8>,
    seens: &mut rustc_hash::FxHashSet<u128>,
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
