use std::{
    collections::{HashMap, HashSet},
    fmt,
    fs::File,
    io::BufReader,
};

use anyhow::Result;
use coitrees::{COITree, IntervalNode};
use noodles::gff::{self as gff, Directive, Line};

use crate::gwas::{self, Snp};

// Annotations data is stored as a COITree with one tree per sequence ID and
// a HashMap of gene names and SNPs linked to those genes.
pub struct Data {
    tree: HashMap<String, COITree<String, usize>>,
    pub genes: HashMap<String, HashSet<(String, i32, Snp)>>,
}

impl fmt::Display for Data {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:#?}", self.genes)
    }
}

impl Data {
    pub fn new<S: AsRef<str>>(
        filename: S,
        attribute: S,
        window: usize,
        feature: S,
    ) -> Result<Self> {
        // Create a GFF reader.
        let mut reader = File::open(filename.as_ref())
            .map(BufReader::new)
            .map(gff::Reader::new)?;

        // Initialize variables.
        let mut reference_sequence_name: String;
        let mut interval: IntervalNode<String, usize>;
        let mut annotations: HashMap<String, Vec<IntervalNode<String, usize>>> = HashMap::new();
        let mut genes: HashMap<String, HashSet<(String, i32, Snp)>> = HashMap::new();

        // Iterate over the annotations.
        for result in reader.lines() {
            let line = result?;

            // Only process GFF records that have the correct feature.
            match line {
                Line::Directive(Directive::StartOfFasta) => break,
                Line::Record(record) => {
                    if record.ty() == feature.as_ref() {
                        // Get the attribute that contains the name as specified
                        // by the user.
                        let name = record
                            .attributes()
                            .get(attribute.as_ref())
                            .expect("could not get attribute")
                            .to_string();
                        // Create an interval representing the gene and 1kb
                        // window around it.
                        interval = IntervalNode::<String, usize>::new(
                            (record.start().get() - window) as i32,
                            (record.end().get() + window) as i32,
                            name.clone(),
                        );

                        // Get the reference sequence name and add the interval
                        // to its values.
                        reference_sequence_name = record.reference_sequence_name().into();
                        annotations
                            .entry(reference_sequence_name)
                            .and_modify(|intervals| intervals.push(interval.clone()))
                            .or_insert(vec![interval]);
                        genes.insert(name, HashSet::new());
                    }
                }
                _ => {}
            }
        }

        // Return the data, converting the annotations into a HashMap with
        // reference sequence names as keys and COITrees as values.
        Ok(Self {
            tree: HashMap::from_iter(
                annotations
                    .into_iter()
                    .map(|(reference_sequence_name, intervals)| {
                        (reference_sequence_name, COITree::new(intervals.to_vec()))
                    })
                    .collect::<Vec<(String, COITree<String, usize>)>>(),
            ),
            genes,
        })
    }

    pub fn link(&mut self, data: &gwas::Data) -> Result<()> {
        // Iterate over the keys and SNPs from the GWAS data.
        for (key, snp) in &data.snps {
            // Split the key into a reference sequence name and position, then get
            // the COITree for that reference sequence name.
            let (reference_sequence_name, position) =
                key.split_once('_').expect("could not split key");
            let position = position.parse::<i32>()?;
            let annotations_tree = self
                .tree
                .get(reference_sequence_name)
                .expect("could not get tree for representative sequence name");

            // Find all genes that overlap with the SNP's position.
            annotations_tree.query(position, position, |interval| {
                self.genes
                    .entry(interval.metadata.clone())
                    .and_modify(|snps| {
                        snps.insert((reference_sequence_name.to_string(), position, snp.clone()));
                    });
            });

            // Iterate over the SNPs linked to the primary SNP.
            for linked_snp in &snp.linkages {
                // Split the key into a reference sequence name and position, then get
                // the COITree for that reference sequence name.
                let (linked_reference_sequence_name, linked_position) =
                    linked_snp.0.split_once('_').expect("could not split key");
                let linked_position = linked_position
                    .parse::<i32>()
                    .expect("could not parse linked position to integer");
                let annotations_tree = self
                    .tree
                    .get(linked_reference_sequence_name)
                    .expect("could not get tree for representative sequence name");

                // Find all genes that overlap with the linked SNP's position.
                annotations_tree.query(linked_position, linked_position, |interval| {
                    // Choose which SNP should be associated with the gene.
                    let (reference_sequence_name, position, snp_to_link) =
                        // If the linked SNP is not in the GWAS data, use the
                        // primary SNP.
                        if let Ok(linked_snp) = &data.get(&linked_snp.0) {
                            // If the SNPs don't have the same sign and the
                            // primary SNP has a smaller p-value, or they
                            // have the same p-value, use the primary SNP.
                            if (linked_snp.effect * snp.effect) > 0.0
                                || linked_snp.p_value == snp.p_value
                            {
                                // If the primary SNP has an effect with a larger
                                // absolute value, use the primary SNP.
                                if snp.effect.abs() > linked_snp.effect.abs() {
                                    (reference_sequence_name, position, snp.clone())
                                } else {
                                    (
                                        linked_reference_sequence_name,
                                        linked_position,
                                        linked_snp.clone(),
                                    )
                                }
                            } else if snp.p_value > linked_snp.p_value {
                                (
                                    linked_reference_sequence_name,
                                    linked_position,
                                    linked_snp.clone(),
                                )
                            } else {
                                (reference_sequence_name, position, snp.clone())
                            }
                        } else {
                            (reference_sequence_name, position, snp.clone())
                        };

                    // Add the selected SNP to the gene it or its linked SNP overlaps.
                    self.genes
                        .entry(interval.metadata.clone())
                        .and_modify(|snps| {
                            snps.insert((
                                reference_sequence_name.to_string(),
                                position,
                                snp_to_link.clone(),
                            ));
                        });
                });
            }
        }
        Ok(())
    }
}
