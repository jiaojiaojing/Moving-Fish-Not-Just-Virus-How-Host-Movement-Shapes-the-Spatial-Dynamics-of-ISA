# R code accompanying the ISA two-patch manuscript

**Manuscript:** *Moving Fish, Not Just Virus: How Host Movement Shapes the Spatial Dynamics of Infectious Salmon Anaemia*  
**Authors:** Jing Jiao and Yuanming Lu

## Purpose

This repository collects the recovered R analysis scripts associated with the manuscript's main and appendix figures. The model represents two salmon-production patches with susceptible fish, infected fish, and environmental pathogen; directed human-mediated fish movement and waterborne pathogen transport are modeled separately. The numerical parameters are **mechanistic scenarios**, not empirically calibrated estimates for a particular farm pair.

**All files are in this repository's top-level directory.** There are no per-figure subfolders and no repeated copies of the same large historical script.

## Figure-to-code guide

| Manuscript figure | Scientific content | Corresponding file |
|---|---|---|
| **Figure 1** | Conceptual schematic of the two-patch S–I–P model, showing host and environmental pathogen movement. | No R script; the diagram was made separately. |
| **Figure 2** | Interaction between one-way host movement and asymmetric waterborne pathogen transport; reproduction number and exporting-patch infection prevalence. | `Figures_02_to_05.R` (Experiment 20) |
| **Figure 3** | Interaction between host movement and external stocking intensity. | `Figures_02_to_05.R` (Experiment 22) |
| **Figure 4** | Interaction between host movement and environmental transmission heterogeneity. | `Figures_02_to_05.R` (Experiment 21) |
| **Figure 5** | Combined effects of stocking, transmission heterogeneity, host movement and waterborne transport on the reproduction number; four-panel analysis. | `Figures_02_to_05.R` (Experiment 23; `p_multi_R0`) |
| **Appendix Figure 6** | Mechanistic explanation of host-movement effects in the exporting patch. | `Figure_06.R` |
| **Appendix Figure 7** | Additional outcomes for the host-movement × waterborne-transport experiment (related to Figure 2). | `Figure_07.R` (Experiment 20) |
| **Appendix Figure 8** | Additional outcomes for the host-movement × stocking experiment (related to Figure 3). | `Figure_08.R` (Experiment 22) |
| **Appendix Figure 9** | Additional outcomes for the host-movement × transmission-heterogeneity experiment (related to Figure 4). | `Figure_09.R` (Experiment 21) |

## Files

- **`Figures_02_to_05.R`** — One combined script for manuscript Figures 2–5. The manuscript Figure 5 display/export block is explicitly labeled after the Figure 4 block; its calculations and plotting object (`p_multi_R0`) are defined earlier in Experiment 23.
- **`Figure_06.R`** — Separate recovered Appendix Figure 6 script.
- **`Figure_07.R`** — Separate extracted Experiment 20 analysis script for Appendix Figure 7.
- **`Figure_08.R`** — Separate extracted Experiment 22 analysis script for Appendix Figure 8.
- **`Figure_09.R`** — Separate extracted Experiment 21 analysis script for Appendix Figure 9.

## Running the code

The scripts use R; the combined analysis uses packages including `deSolve` and `ggplot2`. Open the desired `.R` file in R or RStudio, set the working directory to the repository root, and run the file. Some scripts create an `ISA_two_patch_outputs/` directory and export multiple individual plots and/or numerical results. Check each script's `ggsave()` calls for the exact filenames. **Do not assume that a single script automatically creates the final multi-panel composite shown in the PDF.**

## Verification and limitations

This is a **clearly labeled recovered-code package, not yet a fully validated reproduction release**. The scripts have been mapped by experiment and separated into the five files above, but they have **not been executed end-to-end and their generated images have not been compared pixel-by-pixel or panel-by-panel with the final manuscript PDF**. In particular, the Appendix Figure 7–9 scripts export multiple individual outcomes; they do not currently guarantee automatic assembly into the exact final three-panel manuscript layouts. The historical source includes provisional comments about parameter provenance; manuscript parameter choices and figure captions should be checked against actual runs before archival publication.

For a definitive reproducibility release, run all five scripts, verify packages and parameter values, inspect generated outputs against the manuscript, document R/package versions, and then publish a tagged GitHub release for Zenodo archiving. Add a suitable `LICENSE` and citation metadata after the authors agree on reuse terms.

## Suggested citation

Until a repository DOI exists, cite the manuscript and refer to this code repository by its final GitHub/Zenodo URL. Once Zenodo assigns a DOI, update this README and the manuscript's Code Availability statement.
