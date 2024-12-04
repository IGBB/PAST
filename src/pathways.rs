use std::{
    cmp::Ordering,
    collections::{HashMap, HashSet},
    fmt,
    fs::{self, File},
    io::{BufRead, BufReader, Write},
    path::Path,
};

use anyhow::Result;
use plotters::prelude::*;
use rand::Rng;
use statrs::{
    distribution::{ContinuousCDF, Normal},
    statistics::Statistics,
};

use crate::{annotations, gwas::Snp, options::Mode};
#[derive(Clone, Debug)]
pub struct Pathway {
    pub id: String,
    name: String,
    genes: Vec<String>,
    p_value: Option<f64>,
    q_value: Option<f64>,
    percent_genes_in_gwas: f64,
    enrichment_score: Option<f64>,
    enrichment_scores: Option<Vec<f64>>,
}

impl Pathway {
    pub fn plot(
        &self,
        path: &str,
        assignments: &[(String, usize, (String, i32, Snp))],
    ) -> Result<()> {
        let path = Path::new(path).join("plots");
        fs::create_dir_all(&path)?;

        let ranks: Vec<f32> = assignments
            .iter()
            .filter(|(gene, _rank, _snp)| self.genes.contains(gene))
            .map(|(_gene, rank, _snp)| *rank as f32)
            .collect();

        // Generate the image name and caption from the pathway ID and pathway name
        let image = format!(
            "{}/{}.png",
            path.to_str().expect("could not convert path to string"),
            self.id
        );
        let caption = format!("{}: {}", self.id, self.name);

        // Create the new image with a white background
        let root = BitMapBackend::new(&image, (1280, 960)).into_drawing_area();
        root.fill(&WHITE)?;

        // Split the image into two sections: rug marks and line chart
        let areas = root.split_vertically(60);

        // Draw the rug marks.
        // This part of the chart contains the caption.
        // Left and right margins must be set the same in both sections.
        // The rectangles are expanded by 15 on each side to make sure they appear.
        let mut chart = ChartBuilder::on(&areas.0)
            .caption(caption, ("sans-serif", 30, FontStyle::Bold).into_font())
            .margin_left(10)
            .margin_right(20)
            .margin_top(10)
            .y_label_area_size(60)
            .build_cartesian_2d(0f32..25664f32, 0f32..1f32)?;
        chart.draw_series((0..ranks.len()).map(|index| {
            Rectangle::new(
                [
                    (ranks[index] as f32 - 15_f32, 1f32),
                    (ranks[index] as f32 + 15_f32, 0f32),
                ],
                BLACK.filled(),
            )
        }))?;

        // Draw the line chart.
        // Space is added at the bottom and left sides for the axes labels.
        let mut chart = ChartBuilder::on(&areas.1)
            .x_label_area_size(60)
            .y_label_area_size(60)
            .margin_left(10)
            .margin_right(20)
            .margin_bottom(20)
            .build_cartesian_2d(0f32..25664f32, 0f32..1f32)?;

        // Draw the grid lines and set up the formatting for the x-axis.
        chart
            .configure_mesh()
            .x_max_light_lines(5)
            .x_label_formatter(&|v| format!("{:.0}", v))
            .y_max_light_lines(2)
            .label_style(("sans-serif", 20, FontStyle::Normal, &BLACK))
            .x_desc("Gene Rank")
            .y_desc("Running Enrichment Score")
            .draw()?;

        // Draw points for the line chart.
        chart.draw_series((0..ranks.len()).map(|index| {
            Circle::new(
                (
                    ranks[index] as f32,
                    self.enrichment_scores.as_ref().unwrap()[index] as f32,
                ),
                3,
                RED.filled(),
            )
        }))?;

        // Draw the line.
        chart.draw_series(LineSeries::new(
            (0..ranks.len()).map(|index| {
                (
                    ranks[index] as f32,
                    self.enrichment_scores.as_ref().unwrap()[index] as f32,
                )
            }),
            &RED,
        ))?;
        root.present()?;

        Ok(())
    }

    fn score(
        &mut self,
        assignments: &[(String, usize, (String, i32, Snp))],
        permutations: usize,
    ) -> Result<()> {
        let observed_score = {
            // These indices are the indices of the pathway's genes in the
            // vector of genes.
            let indices: Vec<usize> = assignments
                .iter()
                .filter(|(gene, _rank, _snp)| self.genes.contains(gene))
                .map(|(_gene, rank, _snp)| rank - 1)
                .collect();

            // Enrichment scores are calculated by passing the indices and genes
            // to process_permutation.
            let enrichment_scores = Self::process_permutation(&indices, assignments);

            // These enrichment scores are stored for use in plotting.
            self.enrichment_scores = enrichment_scores;

            // The final observed score is the maximum enrichment score.
            // Sometimes, enrichment_scores contains no values. In that case,
            // -999.99 is returned. Since these pathways have no genes (and thus
            // fewer genes than the cutoff), they'll be dropped later.
            self.enrichment_scores.as_ref().map(|enrichment_scores| {
                enrichment_scores
                    .iter()
                    .max_by(|a, b| a.partial_cmp(b).unwrap())
                    .unwrap_or_else(|| {
                        panic!(
                            "could not get maximum value from {:?}",
                            self.enrichment_scores
                        )
                    })
            })
        };

        let total_genes_in_pathway = self.genes.len() as f64;

        // Filter the genes in the pathway to only contain genes found in the
        // vector of genes.
        self.genes = assignments
            .iter()
            .filter_map(|(gene, _, _)| {
                if self.genes.contains(gene) {
                    Some(gene.clone())
                } else {
                    None
                }
            })
            .collect();

        self.percent_genes_in_gwas = self.genes.len() as f64 / total_genes_in_pathway;

        let permuted_scores: Vec<_> = (0..permutations)
            .collect::<Vec<usize>>()
            .iter()
            .filter_map(|_permutation| {
                let mut indices = HashSet::new();
                let mut rng = rand::thread_rng();

                while indices.len() != self.genes.len() {
                    indices.insert(rng.gen_range(0..assignments.len()));
                }
                let mut indices: Vec<usize> = indices.iter().copied().collect();

                if self.genes.len() == 1 {
                    let mut index = indices[0];
                    while assignments[index].2 .2.effect == 0.0 {
                        index = rng.gen_range(0..assignments.len())
                    }
                    indices[0] = index;
                }

                indices.sort_unstable();
                let enrichment_scores = Self::process_permutation(&indices, assignments);

                if let Some(enrichment_scores) = &enrichment_scores {
                    let max_score = *enrichment_scores
                        .iter()
                        .max_by(|a, b| a.partial_cmp(b).unwrap())
                        .unwrap_or_else(|| {
                            panic!(
                                "could not get maximum value from {:?}",
                                self.enrichment_scores
                            )
                        });
                    Some(max_score)
                } else {
                    None
                }
            })
            .collect();

        // Calculate mean and standard deviation of the enrichment scores and
        // store them.
        let permutation_mean = permuted_scores.clone().mean();
        let permutation_sd = permuted_scores.std_dev();

        if let Some(observed_score) = observed_score {
            // Normalize the observed enrichment score and store it.
            let normalized_observed_score = (observed_score - permutation_mean) / permutation_sd;

            // Create a normal distribution and calculate the p-value of the pathway
            // using the normalized observed score. Store the enrichment score.
            let normal = Normal::new(0.0, 1.0).unwrap();
            self.p_value = Some(1.0 - normal.cdf(normalized_observed_score));
            self.enrichment_score = Some(*observed_score);
        }

        Ok(())
    }

    // This function calculates the enrichment scores, provides a set of indices
    // and genes.
    fn process_permutation(
        indices: &[usize],
        assignments: &[(String, usize, (String, i32, Snp))],
    ) -> Option<Vec<f64>> {
        // Get ranks of genes in permutation.
        let ranks: Vec<_> = indices
            .iter()
            .map(|index| assignments[*index].1 as i64)
            .collect();

        // Get effects of genes in permutation.
        let effects: Vec<_> = indices
            .iter()
            .map(|index| assignments[*index].2 .2.effect.abs())
            .collect();

        // Sum the effects.
        let effects_sum: f64 = effects.iter().sum();

        // Calculate the first stage of factors. Each factor is equal to the
        // rank at an index minus the preceding rank - 1.
        let factors: Vec<i64> = ranks
            .iter()
            .enumerate()
            .map(|(index, rank)| {
                if index == 0 {
                    rank - 1
                } else {
                    rank - ranks[index - 1] - 1
                }
            })
            .collect();

        // To calculate the final stage of factors, cumulatively sum the first
        // stage.
        let factors: Vec<i64> = factors
            .iter()
            .enumerate()
            .map(|(index, factor)| {
                if index == 0 {
                    *factor
                } else {
                    let sum: i64 = factors[0..index].iter().sum();
                    *factor + sum
                }
            })
            .collect();

        // Calculate the pmiss, which reflects a penalty to the factors for
        // genes that aren't in the pathway.
        let pmiss: Vec<f64> = factors
            .iter()
            .map(|factor| *factor as f64 / (assignments.len() - ranks.len()) as f64)
            .collect();

        // Calculate the proportion of each effect, then cumulatively sum these
        // values.
        let phits: Vec<f64> = effects.iter().map(|effect| effect / effects_sum).collect();

        let phits: Vec<f64> = phits
            .iter()
            .enumerate()
            .map(|(index, phit)| {
                if index == 0 {
                    *phit
                } else {
                    let sum: f64 = phits[0..index].iter().sum();
                    phit + sum
                }
            })
            .collect();

        // Calculate the scores by subtracting pmiss from phits.
        let scores: Vec<f64> = phits
            .iter()
            .enumerate()
            .map(|(index, phit)| *phit - pmiss[index])
            .collect();

        if !scores.is_empty() {
            Some(scores)
        } else {
            None
        }
    }
}

impl fmt::Display for Pathway {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}\t{}", self.name, self.genes.len())
    }
}

// Pathways data is stored as
// - a HashMap of pathway IDs with Pathway values
// - a HashMap of mode keys with Vec<Pathway> values.
// - a HashMap of mode keys with Vec<gene name, rank, (reference sequence, position, SNP)
//   values
pub struct Data {
    pub pathways: HashMap<String, Pathway>,
    pub results: HashMap<String, Vec<Pathway>>,
    pub assignments: HashMap<String, Vec<(String, usize, (String, i32, Snp))>>,
}

impl fmt::Display for Data {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:#?}", self.pathways)
    }
}

impl Data {
    pub fn new<S: AsRef<str>>(filename: S) -> Result<Self> {
        // Initialize variables.
        let mut record_data: Vec<&str>;
        let mut pathway_id: String;
        let mut pathways: HashMap<String, Pathway> = HashMap::new();

        // Open pathways file and get lines.
        let buff_reader = BufReader::new(fs::File::open(filename.as_ref())?);
        let mut lines = buff_reader.lines();
        let _header = lines.next();

        // Iterate over lines.
        for record in lines.flatten() {
            // Split lines by TAB cahracters.
            record_data = record.split('\t').collect::<Vec<&str>>();

            // Get pathway data.
            pathway_id = record_data
                .first()
                .expect("could not get pathway ID from record")
                .to_string();
            let pathway_name = record_data
                .get(1)
                .expect("could not get pathway name from record")
                .to_string();
            let gene = record_data
                .get(2)
                .expect("could not get gene from record")
                .to_string();

            // Store pathway data as Pathway.
            pathways
                .entry(pathway_id.clone())
                .and_modify(|pathway| {
                    if !pathway.genes.contains(&gene) {
                        pathway.genes.push(gene.clone());
                    }
                })
                .or_insert(Pathway {
                    id: pathway_id,
                    name: pathway_name,
                    genes: vec![gene],
                    percent_genes_in_gwas: 0.0,
                    p_value: None,
                    q_value: None,
                    enrichment_score: None,
                    enrichment_scores: None,
                });
        }

        Ok(Self {
            pathways,
            results: HashMap::new(),
            assignments: HashMap::new(),
        })
    }

    pub fn score(
        &mut self,
        annotations: &annotations::Data,
        mode: &Mode,
        permutations: usize,
        membership_filter: Option<f64>,
    ) -> Result<()> {
        // Select single SNP to represent gene.
        let assignments: Vec<(String, (String, i32, Snp))> = annotations
            .genes
            .clone()
            .into_iter()
            .filter_map(|(gene, data)| {
                if !data.is_empty() {
                    Some((
                        gene,
                        Self::choose_snp(&data).expect("could not choose SNP to represent gene"),
                    ))
                } else {
                    None
                }
            })
            .collect();

        // Iterate over analysis modes.
        for mode in mode.mode() {
            let mut assignments = assignments.clone();
            // Sort genes by effect based on mode, using p-value and number of
            // linkages as tie-breakers.
            assignments.sort_by(|(name, (_reference_sequence_name, _position, snp)), (other_name, (_other_reference_sequence_name, _other_position, other_snp))| {
                if snp.effect == other_snp.effect {
                    if snp.p_value == other_snp.p_value {
                        if snp.linkages.len() == other_snp.linkages.len() {
                            name.cmp(other_name)
                        } else {
                        snp.linkages
                            .len()
                            .partial_cmp(&other_snp.linkages.len())
                            .unwrap()
                        }
                    } else {
                        snp.p_value.partial_cmp(&other_snp.p_value).unwrap()
                    }
                } else if mode == "increasing" {
                    other_snp.effect.partial_cmp(&snp.effect).unwrap()
                } else if mode == "decreasing" {
                    snp.effect.partial_cmp(&other_snp.effect).unwrap()
                } else {
                    panic!("invalid mode")
                }
            }
            );

            // Rank the assignments and stored the ranked assignments.
            let assignments: Vec<(String, usize, (String, i32, Snp))> = assignments
                .into_iter()
                .enumerate()
                .map(|(index, assignment)| (assignment.0, index + 1, assignment.1))
                .collect();
            self.assignments.insert(mode.clone(), assignments.clone());

            // Iterate over the pathways, scoring and filtering them.
            let mut pathways: Vec<Pathway> = self
                .pathways
                .clone()
                .into_iter()
                .filter_map(|(_pathway_id, mut pathway)| {
                    pathway
                        .score(&assignments, permutations)
                        .expect("could not score pathway");
                    if let Some(membership_filter) = membership_filter {
                        if pathway.percent_genes_in_gwas >= membership_filter {
                            Some(pathway)
                        } else {
                            None
                        }
                    } else if pathway.percent_genes_in_gwas > 0.0 {
                        Some(pathway)
                    } else {
                        None
                    }
                })
                .collect();

            // Extract p-values from the pathways.
            let mut assigned_pvalues: Vec<f64> = pathways
                .iter()
                .map(|pathway| pathway.p_value.expect("no p-value"))
                .collect();

            // Calculate q-values from the p-values and add to pathways.
            let qvalues = Self::calculate_qvalue(&mut assigned_pvalues);
            for (index, pathway) in pathways.iter_mut().enumerate() {
                pathway.q_value = Some(qvalues[index]);
            }

            // Sort the pathways by q-value, using p-value as a tie-breaker.
            pathways.sort_by(|a, b| {
                if a.q_value == b.q_value {
                    a.p_value.partial_cmp(&b.p_value).unwrap()
                } else {
                    a.q_value.partial_cmp(&b.q_value).unwrap()
                }
            });

            self.results.insert(mode, pathways);
        }

        Ok(())
    }

    fn choose_snp(snps: &HashSet<(String, i32, Snp)>) -> Result<(String, i32, Snp)> {
        let mut snps: Vec<(String, i32, Snp)> = snps.clone().into_iter().collect();

        // Count the number of negative and positive SNPs.
        let negative: i64 = snps.iter().filter(|record| record.2.effect < 0_f64).count() as i64;
        let mut positive: i64 = snps.iter().filter(|record| record.2.effect > 0_f64).count() as i64;

        if positive == 0 && negative == 0 {
            positive = snps
                .iter()
                .filter(|record| record.2.effect >= 0_f64)
                .count() as i64;
        }

        Ok(match positive.cmp(&negative) {
            // If there are equal numbers of SNPs with positive effects and SNPs
            // with positive effects, and the numbers for each category are not 0,
            // get the SNP with the largest negative effect and the SNP with the
            // larget positive effect. Select the SNP with the largest absolute
            // value of its effect.
            Ordering::Equal => {
                if positive == 0 && negative == 0 {
                    snps.sort_by(|a, b| {
                        if a.2.effect == b.2.effect {
                            if a.2.p_value == b.2.p_value {
                                a.1.partial_cmp(&b.1).unwrap()
                            } else {
                                a.2.p_value.partial_cmp(&b.2.p_value).unwrap()
                            }
                        } else {
                            b.2.effect.partial_cmp(&a.2.effect).unwrap()
                        }
                    });
                    snps[0].clone()
                } else {
                    let mut positive: Vec<_> = snps
                        .iter()
                        .filter(|record| record.2.effect > 0_f64)
                        .collect();
                    positive.sort_by(|a, b| {
                        if a.2.effect == b.2.effect {
                            if a.2.p_value == b.2.p_value {
                                a.1.partial_cmp(&b.1).unwrap()
                            } else {
                                a.2.p_value.partial_cmp(&b.2.p_value).unwrap()
                            }
                        } else {
                            b.2.effect.partial_cmp(&a.2.effect).unwrap()
                        }
                    });
                    let positive = positive[0].clone();
                    let mut negative: Vec<_> = snps
                        .iter()
                        .filter(|record| record.2.effect < 0_f64)
                        .collect();
                    negative.sort_by(|a, b| {
                        if a.2.effect == b.2.effect {
                            if a.2.p_value == b.2.p_value {
                                a.1.partial_cmp(&b.1).unwrap()
                            } else {
                                a.2.p_value.partial_cmp(&b.2.p_value).unwrap()
                            }
                        } else {
                            a.2.effect.partial_cmp(&b.2.effect).unwrap()
                        }
                    });
                    let negative = negative[0].clone();
                    if positive.2.effect > negative.2.effect.abs() {
                        positive
                    } else {
                        negative
                    }
                }
            }
            // If there are more SNPs with positive effects, choose the SNP with
            // the largest positive effect.
            Ordering::Greater => {
                snps.sort_by(|a, b| {
                    if a.2.effect == b.2.effect {
                        if a.2.p_value == b.2.p_value {
                            a.1.partial_cmp(&b.1).unwrap()
                        } else {
                            a.2.p_value.partial_cmp(&b.2.p_value).unwrap()
                        }
                    } else {
                        b.2.effect.partial_cmp(&a.2.effect).unwrap()
                    }
                });
                snps[0].clone()
            }
            // If there are more SNPs with negative effects, choose the SNP with
            // the largest negative effect.
            Ordering::Less => {
                snps.sort_by(|a, b| {
                    if a.2.effect == b.2.effect {
                        if a.2.p_value == b.2.p_value {
                            a.1.partial_cmp(&b.1).unwrap()
                        } else {
                            a.2.p_value.partial_cmp(&b.2.p_value).unwrap()
                        }
                    } else {
                        a.2.effect.partial_cmp(&b.2.effect).unwrap()
                    }
                });
                snps[0].clone()
            }
        })
    }

    fn calculate_qvalue(p: &mut [f64]) -> Vec<f64> {
        let number_of_p_values = p.len() as f64;

        let mut sorted_indices = (0..p.len()).collect::<Vec<_>>();
        sorted_indices.sort_by(|a, b| p[*b].partial_cmp(&p[*a]).unwrap());
        let mut original_indices = sorted_indices.clone();
        original_indices
            .sort_by(|a, b| sorted_indices[*a].partial_cmp(&sorted_indices[*b]).unwrap());

        p.sort_by(|a, b| b.partial_cmp(a).unwrap());

        let p: Vec<_> = p
            .iter()
            .enumerate()
            .map(|(index, p_value)| {
                p_value * number_of_p_values / (number_of_p_values - index as f64)
            })
            .collect();

        let q_values: Vec<_> = p
            .iter()
            .enumerate()
            .map(|(index, _)| {
                *[
                    *p[0..=index]
                        .iter()
                        .min_by(|a, b| a.partial_cmp(b).unwrap())
                        .unwrap_or_else(|| panic!("could not get minimum value")),
                    1.0,
                ]
                .iter()
                .min_by(|a, b| a.partial_cmp(b).unwrap())
                .unwrap_or_else(|| panic!("could not get maximum value"))
            })
            .collect();

        original_indices
            .iter()
            .map(|index| q_values[*index])
            .collect::<Vec<f64>>()
    }

    // Write genes and pathways data as CSV by mode.
    pub fn create_output(&self, output: &str, mode: &Mode) -> Result<()> {
        for mode in mode.mode() {
            let path = Path::new(output).join(mode.clone());
            fs::create_dir_all(&path)?;

            let mut file = File::create(path.join("genes.csv"))?;

            writeln!(
                file,
                "chromosome,position,effect,p_value,linked_snps,marker,name,rank"
            )?;

            for gene in self
                .assignments
                .get(&mode)
                .expect("could not find assignments results for mode")
            {
                let gene_name = &gene.0;
                let rank = gene.1;
                let reference_sequence_name = &gene.2 .0;
                let position = gene.2 .1;
                let snp = &gene.2 .2;
                writeln!(
                    file,
                    "{},{},{},{},{},{},{},{}",
                    reference_sequence_name,
                    position,
                    snp.effect,
                    snp.p_value,
                    snp.linkages.len(),
                    snp.marker,
                    gene_name,
                    rank,
                )?;
            }

            let mut file = File::create(path.join("pathways.csv"))?;

            writeln!(file, "pathway_id,pathway_name,p_value,q_value,percent_membership_hits,enrichment_score,gene,running_enrichment_score")?;
            for pathway in self
                .results
                .get(&mode)
                .expect("could not find pathways results for mode")
            {
                pathway
                    .plot(
                        path.to_str().expect("could not convert path to string"),
                        self.assignments
                            .get(&mode)
                            .expect("could not find assignments for mdoe"),
                    )
                    .expect(&format!("could not plot pathway {}", pathway));
                for (index, gene) in pathway.genes.iter().enumerate() {
                    writeln!(
                        file,
                        "{},\"{}\",{},{},{},{},{},{}",
                        pathway.id,
                        pathway.name,
                        if let Some(p_value) = pathway.p_value {
                            p_value.to_string()
                        } else {
                            "".to_string()
                        },
                        if let Some(q_value) = pathway.q_value {
                            q_value.to_string()
                        } else {
                            "".to_string()
                        },
                        pathway.percent_genes_in_gwas,
                        if let Some(enrichment_score) = pathway.enrichment_score {
                            enrichment_score.to_string()
                        } else {
                            "".to_string()
                        },
                        gene,
                        if let Some(enrichment_scores) = &pathway.enrichment_scores {
                            enrichment_scores[index].to_string()
                        } else {
                            "".to_string()
                        },
                    )?;
                }
            }
        }

        Ok(())
    }
}
