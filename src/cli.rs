/* crate use */
use crate::error::Error;

#[derive(clap::Clap, Debug)]
#[clap(
    version = "0.1",
    author = "Pierre Marijon <pierre.marijon@hhu.de>",
    about = "Convert a raw kmer count in kff format with a minimizer compression"
)]
pub struct Command {
    #[clap(
        short = 'i',
        long = "input",
        about = "Path of kmers counts in csv format"
    )]
    pub input: String,

    #[clap(short = 'o', long = "output", about = "Path of the kff file")]
    pub output: String,

    #[clap(short = 'k', long = "kmer-size", about = "Kmer size")]
    pub k: u8,

    #[clap(short = 'm', long = "minimizer-size", about = "Minimizer size")]
    pub m: u8,

    #[clap(
        short = 'd',
        long = "delimiter",
        about = "Set delimiter between kmer and count in input",
        default_value = ","
    )]
    pub delimiter: char,

    #[clap(
        short = 'p',
        long = "prefix",
        about = "Prefix add before temporary file"
    )]
    pub prefix: String,
}

pub fn check_params(params: Command) -> Result<Command, Error> {
    if params.k > 64 {
        Err(Error::CliKUpperThan64)
    } else if params.m > 32 {
        Err(Error::CliMUpperThan32)
    } else if params.m >= params.k {
        Err(Error::CliMUpperOrEqualThanK)
    } else {
        Ok(params)
    }
}
