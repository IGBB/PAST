use std::{
    collections::{HashMap, HashSet},
    ffi::OsStr,
    fmt, fs,
    fs::File,
    hash::{Hash, Hasher},
    io::{BufRead, BufReader, Lines},
    iter::Peekable,
    path::Path,
    str::FromStr,
};

use anyhow::{bail, Result};
use regex::Regex;

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
    pub fn new<S: AsRef<str>>(
        filenames: Vec<S>,
        columns_strings: &[String],
        tassel: bool,
    ) -> Result<Self> {
        // Initialize variables to store data.
        let mut record_data: Vec<&str>;
        let mut snps: HashMap<String, Snp> = HashMap::new();

        // Open file and get lines.
        let buff_reader = BufReader::new(fs::File::open(filenames[0].as_ref())?);
        let mut lines = buff_reader.lines().peekable();

        let delimiter = Self::get_delimiter(filenames[0].as_ref(), &mut lines)?;
        let regex_split = Regex::new(&delimiter).expect("invalid regex");

        let columns: Vec<usize> = columns_strings[0]
            .split(',')
            .map(|column| column.parse().expect("could not parse column as usize"))
            .collect();

        // Loop over all lines.
        for record in lines.flatten() {
            // Split the line on the delimiter.
            record_data = regex_split.split(&record).collect::<Vec<&str>>();

            // Get the reference sequence name and position to create a key.
            let reference_sequence_name = record_data
                .get(columns[1] - 1)
                .expect("could not get reference_sequence_name from record")
                .to_string();

            let position: i32 = match Self::convert("position", &record_data, &columns, 2) {
                Ok(value) => value,
                Err(err) => bail!(err),
            };

            let key = format!("{}_{}", reference_sequence_name, position);

            // Get the effect or initialize it for later.
            let effect: f64 = if filenames.len() == 1 {
                match Self::convert("effect", &record_data, &columns, 4) {
                    Ok(value) => value,
                    Err(err) => bail!(err),
                }
            } else {
                0.0
            };

            // Store the SNP by key.
            snps.insert(
                key,
                Snp {
                    p_value: match Self::convert("p-value", &record_data, &columns, 3) {
                        Ok(value) => value,
                        Err(err) => bail!(err),
                    },
                    effect,
                    linkages: Vec::new(),
                    marker: record_data
                        .get(columns[0] - 1)
                        .expect("could not get marker from record")
                        .to_string(),
                },
            );
        }

        if filenames.len() == 2 {
            let mut keys: HashSet<String> = HashSet::new();

            // Open file and get lines.
            let buff_reader = BufReader::new(fs::File::open(filenames[1].as_ref())?);
            let mut lines = buff_reader.lines().peekable();

            let delimiter = Self::get_delimiter(filenames[1].as_ref(), &mut lines)?;
            let regex_split = Regex::new(&delimiter).expect("invalid regex");

            let columns: Vec<usize> = columns_strings[1]
                .split(',')
                .map(|column| column.parse().expect("could not parse column as usize"))
                .collect();

            // Loop over all lines.
            for record in lines.flatten() {
                // Split the line on TAB characters.
                record_data = regex_split.split(&record).collect::<Vec<&str>>();

                // Get the reference sequence name and position to create a key.
                let reference_sequence_name = record_data
                    .get(columns[0] - 1)
                    .expect("could not get reference_sequence_name from record")
                    .to_string();
                let position: i32 = match Self::convert("position", &record_data, &columns, 1) {
                    Ok(value) => value,
                    Err(err) => bail!(err),
                };

                let key = format!("{}_{}", reference_sequence_name, position);

                // Get the effect.
                let effect: f64 = match Self::convert("effect", &record_data, &columns, 2) {
                    Ok(value) => value,
                    Err(err) => bail!(err),
                };

                // Update the SNP's effect by key.
                if !tassel || !keys.contains(&key) {
                    snps.entry(key.clone())
                        .and_modify(|snp| snp.effect = effect);
                }

                keys.insert(key);
            }
        }

        Ok(Self { snps })
    }

    pub fn link<S: AsRef<str>>(
        &mut self,
        filename: S,
        r_squared_cutoff: f64,
        drop_different_loci: bool,
        columns: &[usize],
    ) -> Result<()> {
        // Initialize variables to store data.
        let mut record_data: Vec<&str>;

        // Open file and get lines.
        let buff_reader = BufReader::new(fs::File::open(filename.as_ref())?);
        let mut lines = buff_reader.lines().peekable();

        let delimiter = Self::get_delimiter(filename.as_ref(), &mut lines)?;
        let regex_split = Regex::new(&delimiter).expect("invalid regex");

        // Loop over all lines.
        for record in lines.flatten() {
            // Split the line on the delimiter.
            record_data = regex_split.split(&record).collect::<Vec<&str>>();

            // Get the reference sequence names and positions to create keys.
            let locus_1 = record_data
                .get(columns[0] - 1)
                .expect("could not get locus");
            let position_1: i32 = match Self::convert("position_1", &record_data, columns, 1) {
                Ok(value) => value,
                Err(err) => bail!(err),
            };
            let locus_2 = if drop_different_loci {
                record_data
                    .get(columns[2] - 1)
                    .expect("could not get locus")
            } else {
                locus_1
            };
            let position_2: i32 = if drop_different_loci {
                match Self::convert("position_2", &record_data, columns, 3) {
                    Ok(value) => value,
                    Err(err) => bail!(err),
                }
            } else {
                match Self::convert("position_2", &record_data, columns, 2) {
                    Ok(value) => value,
                    Err(err) => bail!(err),
                }
            };

            // Get the R^2 value.
            let r_squared: f64 = if drop_different_loci {
                match Self::convert("r_squared", &record_data, columns, 4) {
                    Ok(value) => value,
                    Err(err) => bail!(err),
                }
            } else {
                match Self::convert("r_squared", &record_data, columns, 3) {
                    Ok(value) => value,
                    Err(err) => bail!(err),
                }
            };

            // If the R^2 value isn't null, check whether the two SNPs should
            // be linked.
            // if let Ok(r_squared) = r_squared {
            // Link the two SNPs if the loci aren't checked or if they are
            // the same and the R^2 value is greater than the user-provided
            // cutoff.
            if (!drop_different_loci || locus_1 == locus_2) && r_squared >= r_squared_cutoff {
                self.link_snps(locus_1, position_1, locus_2, position_2, r_squared)?;
                self.link_snps(locus_2, position_2, locus_1, position_1, r_squared)?;
            }
            // }
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

    fn convert<T: FromStr>(
        field: &str,
        record_data: &[&str],
        columns: &[usize],
        index: usize,
    ) -> Result<T, anyhow::Error> {
        if let Some(returned_value) = record_data.get(columns[index] - 1) {
            if let Ok(parsed_value) = returned_value.parse::<T>() {
                Ok(parsed_value)
            } else {
                bail!("could not parse {field} ({returned_value}) as numeric")
            }
        } else {
            bail!("could not get {field} from record")
        }
    }

    // Get delimiter based on extension or first two lines.
    pub fn get_delimiter(
        filename: &str,
        lines: &mut Peekable<Lines<BufReader<File>>>,
    ) -> Result<String> {
        // Get the header line by iterating.
        let header_line = lines
            .next()
            .expect("could not read line")
            .expect("could not read header");

        // Get the first line by peeking.
        let first_line = lines
            .peek()
            .expect("could not read line")
            .as_ref()
            .expect("could not read first line");

        // Get the delimiter.
        //   .tsv/.csv imply the delimiter.
        //   .txt could be TSV, CSV, or fixed width.
        let delimiter = match Path::new(&filename)
            .extension()
            .and_then(OsStr::to_str)
            .expect("could not get OsStr")
        {
            "tsv" => "\t",
            "csv" => ",",
            _ => {
                // Count tabs and commas in header and first line.
                let tab = (
                    header_line.matches('\t').count(),
                    first_line.matches('\t').count(),
                );
                let comma = (
                    header_line.matches(',').count(),
                    first_line.matches(',').count(),
                );

                // Check if they have the same number in both lines.
                let tab_equal = tab.0 == tab.1;
                let comma_equal = comma.0 == comma.1;

                // Assume <TAB> if unequal number of commas.
                if tab_equal && !comma_equal {
                    "\t"
                // Assume , if unequal number of <TAB>s.
                } else if comma_equal && !tab_equal {
                    ","
                // Assume <TAB> if there are more tabs than ,
                } else if tab_equal && comma_equal && tab.0 > comma.0 {
                    "\t"
                // Assume , if there are more , than <TAB>s
                } else if tab_equal && comma_equal && comma.0 > tab.0 {
                    ","
                // Otherwise, use 1 or more spaces as the delimiter.
                } else {
                    r"\s+"
                }
            }
        };

        Ok(delimiter.to_string())
    }
}
