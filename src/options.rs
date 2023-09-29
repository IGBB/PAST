use std::path::Path;

use clap::{Parser, ValueEnum};

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum Mode {
    /// Find pathways associated with an increase in the trait
    Increasing,
    /// Find pathways associated with a decrease in the trait
    Decreasing,
    /// Run both analyses
    Both,
}

impl Mode {
    pub fn mode(&self) -> Vec<String> {
        match &self {
            Mode::Increasing => vec!["increasing".into()],
            Mode::Decreasing => vec!["decreasing".into()],
            Mode::Both => vec!["increasing".into(), "decreasing".into()],
        }
    }
}

#[derive(Parser)]
#[command(
    author,
    version,
    about,
    long_about = "Pathways Analysis Study Tool\nReport bugs at https://github.com/IGBB/PAST/issues.",
    max_term_width = 80
)]
pub struct Options {
    /// a file containing annotations formatted as GFF3
    #[arg(short, long, value_name = "FILE", value_parser = validate_file)]
    pub annotations: String,

    /// the attribute that contains gene names that match those in the pathways file
    #[arg(short = 'x', long, value_name = "ATTRIBUTE", default_value_t = String::from("ID"))]
    pub attribute: String,

    /// a file containing GWAS data
    #[arg(short, long, value_name = "FILE", value_parser = validate_file)]
    pub gwas: String,

    /// a comma-separated value providing column numbers of required GWAS data; column order is marker, sequence name, position, p-value, effect
    #[arg(short = 'c', long, value_name = "1,2,3,4,5", value_delimiter = ',')]
    pub gwas_columns: Vec<usize>,

    /// a file containing linkage disequilibrium data
    #[arg(short, long, value_name = "FILE", value_parser = validate_file)]
    pub linkage_disequilibrium: String,

    /// a comma-separated value providing column numbers of required linkage disequilibrium data; column order is first sequence name, first position, second sequence name, second position, R^2
    #[arg(short = 'k', long, value_name = "1,2,3,4,5", value_delimiter = ',')]
    pub linkage_columns: Vec<usize>,

    /// the value of R^2 at which two SNPs are considered linked [0.0 - 1.0]
    #[arg(short, long, value_parser = validate_r_squared)]
    pub r_squared_cutoff: f64,

    /// drop linkages between positions with different loci
    #[arg(short, long)]
    pub drop_different_loci: bool,

    /// a file containing pathways data
    #[arg(short, long, value_name = "FILE", value_parser = validate_file)]
    pub pathways: String,

    /// analysis mode
    #[arg(short, long, value_name = "MODE")]
    pub mode: Mode,

    /// the number of permutations to determine pathway significance
    #[arg(short = 'n', long)]
    pub permutations: usize,

    /// keep only pathways with this fraction of the genes or higher linked to SNPs in the GWAS data [0.0 - 1.0]
    #[arg(short = 'f', long, value_parser = validate_membership_filter)]
    pub membership_filter: Option<f64>,

    /// the directory to which results will be written (will be created if it does not exist)
    #[arg(short, long, value_name = "DIRECTORY", default_value_t = String::from("./"))]
    pub output: String,
}

fn validate_r_squared(option: &str) -> Result<f64, String> {
    let r_squared: f64 = option
        .parse()
        .map_err(|_| format!("`{option}` cannot be parsed a number"))?;
    if (0.0..=1.0).contains(&r_squared) {
        Ok(r_squared)
    } else {
        Err("R^2 must be between 0.0 and 1.0 inclusive".to_string())
    }
}

fn validate_membership_filter(option: &str) -> Result<f64, String> {
    let membership_filter: f64 = option
        .parse()
        .map_err(|_| format!("`{option}` cannot be parsed a number"))?;
    if (0.0..=1.0).contains(&membership_filter) {
        Ok(membership_filter)
    } else {
        Err("Membership filter must be between 0.0 and 1.0 inclusive".to_string())
    }
}

fn validate_file(option: &str) -> Result<String, String> {
    if Path::new(option).exists() {
        Ok(option.to_string())
    } else {
        Err(format!("The file {option} does not exist"))
    }
}
