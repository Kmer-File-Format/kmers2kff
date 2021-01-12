/* crate use */
use thiserror::Error;

#[derive(Debug, Error)]
pub enum Error {
    #[error("Minimizer size is upper or equal than kmer size")]
    CliMUpperOrEqualThanK,

    #[error("Kmer size is upper than 64")]
    CliKUpperThan64,

    #[error("Minimizer size is upper than 32")]
    CliMUpperThan32,
}
