# Red Butte Canyon Methane-Cycling Community 16S Analysis

This repository contains R scripts used for the analyses and figure generation associated with:

> Beals DG, Munn JJ, Puri AW. 2024. *Methane-oxidizing bacterial community dynamics in sub-alpine forest soil*. Microbiology Spectrum 12:e00834-24. https://doi.org/10.1128/spectrum.00834-24

## Project Overview

Methane-oxidizing bacteria play an important role in mitigating methane emissions from soil, but the environmental factors shaping these communities and their activity remain poorly understood. This study compared microbial communities from riparian and upland soils in Red Butte Canyon (Salt Lake City, Utah, USA) using 16S rRNA gene amplicon sequencing of both genomic DNA (gDNA) and complementary DNA (cDNA) libraries.

The primary goal was to compare the composition of the total microbial community (gDNA) with the potentially active microbial community (cDNA) and determine how these communities varied across soil types and seasons. The study also examined methane-cycling microbial populations, including methanotrophs, methylotrophs, and associated heterotrophic bacteria, in relation to environmental measurements and methane flux.

## Main Findings

Key findings from the study included:

* Stable microbial communities (gDNA) differed primarily between riparian and upland soils.
* Potentially active microbial communities (cDNA) differed primarily between sampling months.
* Methanotrophs were detected in all soil samples but were generally present at low relative abundance.
* Potentially active methanotrophs and non-methanotrophic methylotrophs exhibited positive correlations in cDNA libraries.
* Methane flux did not significantly correlate with methanotroph abundance or activity across the sampled sites.

## Repository Contents

| Analysis Category | Script | Purpose |
|------------------|--------|----------|
| Data Processing | `June_Oct_16S_MERGE_processing.R` | Processing, filtering, merging, and organization of 16S rRNA gene amplicon sequencing data. |
| Community Composition | `June_Oct_16S_MERGE_analysis.Rmd` | Analysis of microbial community composition across soil types, months, and nucleic acid templates. |
| Community Composition | `June_Oct_16S_MERGE_analysis_relabund.Rmd` | Relative abundance analyses and taxonomic summaries of microbial communities. |
| Alpha Diversity | `alpha_diversity.Rmd` | Calculation and comparison of Hill-number diversity metrics among sample groups. |
| Beta Diversity | `UNIFRAC.Rmd` | UniFrac-based phylogenetic distance analyses and community comparisons. |
| Beta Diversity | `ADONIS-BETADISPER.Rmd` | PERMANOVA and beta-dispersion analyses of microbial community structure. |
| Phylogenetics | `phylogenetics.Rmd` | Phylogenetic analyses and tree generation for microbial taxa. |
| Differential Abundance | `DESeq2_gDNA_cDNA.Rmd` | Differential abundance testing between gDNA and cDNA libraries using DESeq2. |
| Differential Abundance | `upset.Rmd` | Visualization of shared and unique taxa among sample groups using UpSet plots. |
| Methane-Cycling Communities | `MOB-visualization.Rmd` | Visualization of methanotrophs, methylotrophs, and associated heterotrophic community members. |
| Methane-Cycling Communities | `MOB-visualization-family.Rmd` | Family-level analyses and visualizations of methane-cycling community structure. |
| Taxonomic Composition | `family.Rmd` | Family-level taxonomic composition and abundance analyses. |
| Environmental Analyses | `abiotic.Rmd` | Analysis of methane flux measurements and environmental variables. |
| Visualization | `HEATMAP.Rmd` | Generation of taxonomic heatmaps and abundance visualizations. |
| Visualization | `summary_figures.Rmd` | Creation of manuscript figures and summary visualizations. |
| Documentation | `citations.Rmd` | Reference management and manuscript citation tracking. |

## Data Availability

Raw sequencing data are available through the NCBI Sequence Read Archive under BioProject accession **PRJNA1049057**.

This repository contains analysis scripts only. File paths and software environments reflect the original analysis environment and may require modification before reuse. For additional information regarding the underlying biological questions, experimental methods, or original datasets, please contact [Aaron Puri](https://www.chemistry.utah.edu/faculty/aaron-w-puri/) at the University of Utah or see the associated [publication](https://doi.org/10.1128/spectrum.00834-24).
