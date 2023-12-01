use anyhow::Result;
use clap::Parser;
use past::{annotations, gwas, options, pathways};

fn main() -> Result<()> {
    // Get command-line options.
    let options = options::Options::parse();

    // Read annotations file, expanding genes by 1kb.
    let mut annotations = annotations::Data::new(
        options.annotations,
        options.attribute,
        1000,
        String::from("gene"),
    )?;

    // Read GWAS data.
    let mut gwas = gwas::Data::new(options.gwas, &options.gwas_columns, options.tassel)?;

    // Link SNPs using linkage disequilibrium data.
    gwas.link(
        options.linkage_disequilibrium,
        options.r_squared_cutoff,
        options.drop_different_loci,
        &options.linkage_columns,
    )?;

    // Link SNPs to genes.
    annotations.link(&gwas)?;

    // Read pathways data.
    let mut pathways = pathways::Data::new(options.pathways)?;

    // Score pathways.
    pathways.score(
        &annotations,
        &options.mode,
        options.permutations,
        options.membership_filter,
    )?;

    // Write CSVs for pathways and ranked genes.
    pathways.write_csv(&options.output, &options.mode)?;

    Ok(())
}
