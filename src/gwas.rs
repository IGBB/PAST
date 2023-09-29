use std::{
    collections::HashMap,
    fmt, fs,
    hash::{Hash, Hasher},
    io::{BufRead, BufReader},
};

use anyhow::{bail, Result};

#[derive(Clone, Debug)]
pub struct Snp {
    pub effect: f64,
    pub p_value: f64,
    pub marker: String,
    pub linkages: Vec<(String, f64)>,
}

impl Hash for Snp {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.marker.hash(state);
    }
}

impl Eq for Snp {}

impl PartialEq for Snp {
    fn eq(&self, other: &Self) -> bool {
        self.marker == other.marker
    }
}

impl fmt::Display for Snp {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{}",
            self.effect,
            self.p_value,
            self.linkages.len()
        )
    }
}

// GWAS data is stored as a HashMap of sequence_position keys with SNPs as values.
#[derive(Debug)]
pub struct Data {
    pub snps: HashMap<String, Snp>,
}

impl fmt::Display for Data {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:#?}", self.snps)
    }
}

impl Data {
    pub fn new<S: AsRef<str>>(filename: S, columns: &[usize], header: bool) -> Result<Self> {
        // Initialize variables to store data.
        let mut record_data: Vec<&str>;
        let mut snps: HashMap<String, Snp> = HashMap::new();

        // Open file and get lines.
        let buff_reader = BufReader::new(fs::File::open(filename.as_ref())?);
        let mut lines = buff_reader.lines();

        // Skip the first line if the file has a header.
        if header {
            lines.next();
        }

        // Loop over all lines.
        for record in lines.flatten() {
            // Split the line on TAB characters.
            record_data = record.split('\t').collect::<Vec<&str>>();

            // Get the reference sequence name and position to create a key.
            let reference_sequence_name = record_data
                .get(columns[1] - 1)
                .expect("could not get reference_sequence_name from record")
                .to_string();
            let position = record_data
                .get(columns[2] - 1)
                .expect("could not get position from record")
                .parse::<i32>()
                .expect("could not parse position as integer");
            let key = format!("{}_{}", reference_sequence_name, position);

            // Store the SNP by key.
            snps.insert(
                key,
                Snp {
                    p_value: record_data
                        .get(columns[3] - 1)
                        .expect("could not get p_value from record")
                        .parse::<f64>()
                        .expect("could not parse p-value as float"),
                    effect: record_data
                        .get(columns[4] - 1)
                        .expect("could not get effect from record")
                        .parse::<f64>()
                        .expect("could not parse effect as float"),
                    linkages: Vec::new(),
                    marker: record_data
                        .get(columns[0] - 1)
                        .expect("could not get marker from record")
                        .to_string(),
                },
            );
        }

        Ok(Self { snps })
    }

    pub fn link<S: AsRef<str>>(
        &mut self,
        filename: S,
        r_squared_cutoff: f64,
        drop_different_loci: bool,
        columns: &[usize],
        header: bool,
    ) -> Result<()> {
        // Initialize variables to store data.
        let mut record_data: Vec<&str>;

        // Open file and get lines.
        let buff_reader = BufReader::new(fs::File::open(filename.as_ref())?);
        let mut lines = buff_reader.lines();

        // Skip the first line if the file has a header.
        if header {
            lines.next();
        }

        // Loop over all lines.
        for record in lines.flatten() {
            // Split the line on TAB characters.
            record_data = record.split('\t').collect::<Vec<&str>>();

            // Get the reference sequence names and positions to create keys.
            let locus_1 = record_data
                .get(columns[0] - 1)
                .expect("could not get locus");
            let position_1 = record_data
                .get(columns[1] - 1)
                .expect("could not get position1")
                .parse::<i32>()
                .expect("could not parse position1 as integer");
            let locus_2 = if drop_different_loci {
                record_data
                    .get(columns[2] - 1)
                    .expect("could not get locus")
            } else {
                locus_1
            };
            let position_2 = if drop_different_loci {
                record_data
                    .get(columns[2] - 1)
                    .expect("could not get position2")
                    .parse::<i32>()
                    .expect("could not parse position2 as integer")
            } else {
                record_data
                    .get(columns[3] - 1)
                    .expect("could not get position2")
                    .parse::<i32>()
                    .expect("could not parse position2 as integer")
            };

            // Get the R^2 value.
            let r_squared = if drop_different_loci {
                record_data
                    .get(columns[3] - 1)
                    .expect("could not get R^2")
                    .parse::<f64>()
            } else {
                record_data
                    .get(columns[4] - 1)
                    .expect("could not get R^2")
                    .parse::<f64>()
            };

            // If the R^2 value isn't null, check whether the two SNPs should
            // be linked.
            if let Ok(r_squared) = r_squared {
                // Link the two SNPs if the loci aren't checked or if they are
                // the same and the R^2 value is greater than the user-provided
                // cutoff.
                if (!drop_different_loci || locus_1 == locus_2) && r_squared >= r_squared_cutoff {
                    self.link_snps(locus_1, position_1, locus_2, position_2, r_squared)?;
                    self.link_snps(locus_2, position_2, locus_1, position_1, r_squared)?;
                }
            }
        }
        Ok(())
    }

    // Link two SNPs.
    fn link_snps(
        &mut self,
        reference_sequence_name: &str,
        position: i32,
        reference_sequence_name_2: &str,
        position_2: i32,
        r_squared: f64,
    ) -> Result<()> {
        // Create keys and add them to the SNP's linkages vector.
        let first_key = format!("{}_{}", reference_sequence_name, position);
        let second_key = format!("{}_{}", reference_sequence_name_2, position_2);
        self.snps
            .entry(first_key)
            .and_modify(|snp| snp.linkages.push((second_key, r_squared)));
        Ok(())
    }

    // Get SNP by key.
    pub fn get(&self, snp: &str) -> Result<Snp> {
        if let Some(snp) = self.snps.get(snp) {
            Ok(snp.clone())
        } else {
            bail!("could not get SNP")
        }
    }
}
