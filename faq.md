---
layout: page
permalink: /faq/
title: Frequently Asked Questions
---

## Why did PAST crash during the SNP to genes step?

You may be using an older version of PAST. Due to mistakes with the release process, PAST was released without a critical bug fix that vastly increases the speed of the package. If you've been using PAST and finding it slower than you expected, (i.e., it takes more than an hour), please make sure you have PAST 1.2.2 or higher and try your analysis again.

## Do I need to download PAST and PAST Shiny both in order to run PAST Shiny?

Yes. PAST Shiny is a user interface to PAST, so you do need both. The instructions for installing PAST can be found [here](/)

## Where can I get the genes and pathways files?

Genes and pathways files can likely be found at the online database for your organism or a closely related model species. Be sure to choose the correct file based on the reference sequence version that you used to align your SNPs in the initial phase before GWAS. You may need to modify your pathways file to the format needed for PAST.

## Can I load data as a .ZIP file?

No. PAST accepts files that have been compressed using GZIP, BZIP, or XZ compression. ZIP-formatted files are archived and compressed, and PAST doesn't take archives as input. If you are running PAST Shiny locally, using compressed versus uncompressed files is unlikely to matter much, but if you are using an online version of PAST Shiny, you will likely want to compress your files first to save time during the upload process. On Windows, [7-Zip](https://www.7-zip.org/) can be used to created these types of compressed files. On Linux, an archive manager can likely create these files, or they can be compressed from the command-line.