//! A set of function to convert small sequence (less than 32 nucleotide) in 2 bit representation.
//! - A or a -> 00
//! - C or c -> 01
//! - T or t -> 10
//! - G or g -> 11
//!
//! We use the second and thrid bit of each value provide, if you provide no ACTG value this function silently convert to A, C, T or G, for exemple N or n is convert in G.
//!
//! With this coding and if kmer size is odd, if the popcount of forward is odd the popcount of reverse is even. In this library if a kmer have even popcount is the canonical kmer.
//!
//! If we work only with canonical kmer, we can remove one bit at any extremity. To reconstruct lost bit, if result have even popcount we add a 0, if it's ood we add 1.
//!
//! This 2bit coding is inspired by https://cs.stackexchange.com/questions/82644/compact-mapping-from-an-involuted-set

use fasthash::FastHash;

/// Convert a sequence in 2 bit representation if suseq is larger than 32 only the last 32 nuc is store
#[inline(always)]
pub fn seq2bit(subseq: &[u8]) -> u128 {
    let mut kmer: u128 = 0;

    for n in subseq {
        kmer <<= 2;
        kmer |= nuc2bit(*n);
    }

    kmer
}

/// Convert a nucleotide in 2bit representation, by use conversion present in [seq2bit](seq2bit)
#[inline(always)]
pub fn nuc2bit(nuc: u8) -> u128 {
    (nuc as u128 >> 1) & 0b11
}

/// Convert a 2 bit repersentation in String.
#[inline(always)]
pub fn kmer2seq(mut kmer: u128, k: u8) -> String {
    let mut buffer: [u8; 64] = [0; 64];

    for i in (0..k).rev() {
        buffer[i as usize] = bit2nuc(kmer & 0b11);

        kmer >>= 2;
    }

    unsafe { String::from_utf8_unchecked((&buffer[..k as usize]).to_vec()) }
}

/// Convert the 2bit representation of a nucleotide in nucleotide
#[inline(always)]
pub fn bit2nuc(bit: u128) -> u8 {
    match bit {
        0 => b'A',
        1 => b'C',
        2 => b'T',
        3 => b'G',
        _ => b'G',
    }
}

/// Take a kmer and return the canonical form
#[inline(always)]
pub fn canonical(kmer: u128, k: u8) -> (u128, bool) {
    let rev = revcomp(kmer, k);

    if kmer < rev {
        (kmer, true)
    } else {
        (rev, false)
    }
}

/// Return true if the kmer parity is even
#[inline(always)]
pub fn parity_even(kmer: u128) -> bool {
    kmer.count_ones() % 2 == 0
}

/// Return the reverse complement of kmer
#[inline(always)]
pub fn revcomp(kmer: u128, k: u8) -> u128 {
    rev(speed_comp(kmer), k)
}

#[inline(always)]
fn speed_comp(kmer: u128) -> u128 {
    kmer ^ 0xAAAA_AAAA_AAAA_AAAA_AAAA_AAAA_AAAA_AAAA
}

/// Return the complement of kmer
#[inline(always)]
pub fn comp(kmer: u128, k: u8) -> u128 {
    speed_comp(kmer) & ((1 << (2 * k)) - 1)
}

/// Return the reverse of kmer
#[inline(always)]
pub fn rev(mut kmer: u128, k: u8) -> u128 {
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

/// return minimizer of kmer and position
pub fn get_minimizer(mut kmer: u128, k: u8, m: u8) -> (u128, usize, bool) {
    let max_len = (k - m + 1) as usize;
    let mask = (1 << (m * 2)) - 1;

    let mut score = u128::max_value();
    let mut minimizer = 0;
    let mut position = 0;
    let mut forward = true;

    for i in 0..max_len {
        let rb_index = (max_len - i - 1) as usize;

        let (mini, local_forward) = canonical(kmer & mask, m);

        let local_score = fasthash::murmur3::Hash128_x64::hash(mini.to_be_bytes());

        if local_score < score {
            score = local_score;
            position = rb_index;
            minimizer = mini;
            forward = local_forward;
        }

        kmer >>= 2;
    }

    (minimizer, position, forward)
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn seq2bit_() {
        // TAGGC -> 1000111101
        assert_eq!(seq2bit(b"TAGGC"), 0b1000111101);

        // GCCTA -> 110101100
        assert_eq!(seq2bit(b"GCCTA"), 0b1101011000);
    }

    #[test]
    fn bit2seq_() {
        // 1000111101 -> TAGGC
        assert_eq!(kmer2seq(0b1000111101, 5), "TAGGC");

        // 110101100 -> GCCTA
        assert_eq!(kmer2seq(0b1101011000, 5), "GCCTA");

        assert_eq!(
            kmer2seq(0b1101011000, 31),
            "AAAAAAAAAAAAAAAAAAAAAAAAAAGCCTA"
        );
    }

    #[test]
    fn canonical_() {
        // TAGGC -> 1000111101 canonical TAGGC -> 1000111101
        assert_eq!(canonical(0b1000111101, 5), (0b1000111101, true));

        // GCCTA -> 1101011000 canonical TAGGC -> 1000111101
        assert_eq!(canonical(0b1101011000, 5), (0b1000111101, false));
    }

    #[test]
    fn parity_even_() {
        assert_eq!(parity_even(0b1111), true);
        assert_eq!(parity_even(0b1110), false);
    }

    #[test]
    fn revcomp_() {
        // TAGGC -> 1000111101 revcomp GCCTA -> 1101011000
        assert_eq!(0b1000111101, revcomp(0b1101011000, 5))
    }

    #[test]
    fn comp_() {
        // TAGGC -> 1000111101 comp 0001001011
        assert_eq!(comp(0b1000111101, 5), 0b0010010111);
    }

    #[test]
    fn rev_() {
        // TAGGC -> 1000111101 rev CGGAT -> 0111110010
        let var = 0b1000111101;

        assert_eq!(498, rev(var, 5));
    }

    #[test]
    fn large_kmer() {
        assert_eq!(
            "CGGAGAGCAGAACAGTCTTACTTTCGTGCACAG",
            kmer2seq(seq2bit(b"CGGAGAGCAGAACAGTCTTACTTTCGTGCACAG"), 33)
        );
        assert_eq!(
            "CTGTGCACGAAAGTAAGACTGTTCTGCTCTCCG",
            kmer2seq(
                revcomp(seq2bit(b"CGGAGAGCAGAACAGTCTTACTTTCGTGCACAG"), 33),
                33
            )
        );

        assert_eq!(
            "AAGCGTTTGCGGGGCGCGGTCCAAAGCATTACACAGGGAGACGTTCGCATGATCTCCCGAATGG",
            kmer2seq(
                seq2bit(b"AAGCGTTTGCGGGGCGCGGTCCAAAGCATTACACAGGGAGACGTTCGCATGATCTCCCGAATGG"),
                64
            )
        );

        assert_eq!(
            "CCATTCGGGAGATCATGCGAACGTCTCCCTGTGTAATGCTTTGGACCGCGCCCCGCAAACGCTT",
            kmer2seq(
                revcomp(
                    seq2bit(b"AAGCGTTTGCGGGGCGCGGTCCAAAGCATTACACAGGGAGACGTTCGCATGATCTCCCGAATGG"),
                    64
                ),
                64
            )
        );
    }

    #[test]
    fn canonical_by_minimizer() {
        let k = 7 as u8;
        let m = 5;

        let seq = "ACCTATGAACTTACG";

        let mut canos = Vec::new();
        for i in 0..(seq.len() + 1 - k as usize) {
            let kmer = seq2bit(seq[i..(i + k as usize)].as_bytes());

            let (mini, pos, forward) = get_minimizer(kmer, k, m);

            if forward {
                canos.push(kmer);
            } else {
                canos.push(revcomp(kmer, k));
            }

            println!(
                "kmer {} rev {} mini {} pos {} mini forward {} min forward {} popcount forward {}",
                kmer2seq(kmer, k),
                kmer2seq(revcomp(kmer, k), k),
                kmer2seq(mini, m),
                pos,
                forward,
                kmer < revcomp(kmer, k),
                parity_even(kmer)
            );
        }

        assert_eq!(
            &[
                seq2bit(b"CATAGGT"),
                seq2bit(b"TCATAGG"),
                seq2bit(b"TTCATAG"),
                seq2bit(b"TATGAAC"),
                seq2bit(b"ATGAACT"),
                seq2bit(b"TGAACTT"),
                seq2bit(b"TAAGTTC"),
                seq2bit(b"AACTTAC"),
                seq2bit(b"CGTAAGT"),
            ],
            &canos[..]
        );
    }
}
