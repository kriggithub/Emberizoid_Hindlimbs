# Emberizoid Hindlimbs

Analysis code for **ecomorphology and phylogenetic signal of hindlimbs in North American sparrows and blackbirds (Passerellidae & Icteridae).** The repository currently contains the blackbird (Icteridae) analyses.

## Questions

1. What are the main axes of hindlimb (tarsus, midtoe and hallux) variation across blackbird species, in terms of size and shape?
2. How much phylogenetic signal do individual hindlimb measurements and PC axes carry?
3. Is hindlimb size (PC1) or shape (PC2, PC3) associated with habitat, foraging, locomotion (hop vs. walk), singing or nesting behaviour?
4. Do these associations still hold once shared ancestry is accounted for?
5. How have habitat use and locomotion mode evolved across the Icteridae phylogeny?

## Methods

- **Data reshaping**: per-specimen hindlimb measurements reformatted from a wide, species-by-column sheet into a tidy table and merged with species-level ecological scores (habitat, foraging, locomotion, singing, nesting, brood parasitism)
- **PCA** (`prcomp`): on 9 natural-log-transformed hindlimb measurements; PC1 is a size axis, PC2 and PC3 are length-vs-width shape axes
- **Phylogenetic signal** (`phytools::phylosig`): Pagel's lambda and Blomberg's K, with significance tests, for each measurement and PC1–PC3 on species means, using a pruned Icteridae MCC species tree (`ape`)
- **Univariate trait–ecology tests**: phylogenetic ANOVA (`phytools::phylANOVA`) and PGLS (`caper::pgls`), compared with non-phylogenetic `lm` / `TukeyHSD`
- **Multivariate tests**: phylogenetic MANOVA (`geiger::aov.phylo`, averaged over 10 runs of 1,000 simulations), multivariate regression (`car::Anova`) and phylogenetic multivariate GLS under Brownian motion (`mvMORPH::mvgls`, compared with a null model by GIC)
- **Ancestral state estimation** (`phytools::make.simmap`): stochastic character mapping of habitat and locomotion, comparing ER, SYM and ARD models by AIC
- **Figures**: PC1 violin plots (`vioplot`) and PC2–PC3 convex-hull plots (`ggplot2`) for each ecological variable, plus ancestral-state trees

## Repository contents

| File | Description |
|---|---|
| `EmberizoidHindlimbs.Rproj` | RStudio project |
| `Blackbird_analyses/0_BlackBirdData Reformatting.R` | Reformats the raw hindlimb sheet and merges it with ecological data |
| `Blackbird_analyses/1_BlackbirdPCA and Exploration.R` | Full analysis: PCA, phylogenetic signal, pANOVA/PGLS, MANOVA, multivariate GLS, ancestral states and figures |
| `Blackbird_analyses/Blackbird Hindlimb Dataset.csv` | Raw hindlimb measurements (up to 10 specimens per species, one column per species) |
| `Blackbird_analyses/Blackbird Ecology Data.csv` | Species-level ecological descriptions and scores |
| `Blackbird_analyses/Icteridae_SpeciesTree_MCC.tre` | Icteridae maximum clade credibility species tree (NEXUS) |
| `Blackbird_analyses/Blackbirdhindlimb_tidy.csv`, `AllBlackbirdData.csv`, `AllBlackbirdData_PCA.csv` | Cleaned per-specimen datasets (tidy measurements; plus ecology; plus PC1–PC3 scores) |
| `Blackbird_analyses/PCA_Eigenvectorblackbird.csv` | PCA loadings and proportion of variance per axis |
| `Blackbird_analyses/Phylogenetic_Signal_Results_blackbirds.csv` | Pagel's lambda and Blomberg's K results per trait |
| `Blackbird_analyses/*Size.pdf`, `*Shape.pdf` | PC1 violin plots and PC2–PC3 hull plots by habitat, foraging, locomotion, singing and nesting |
| `Blackbird_analyses/Habitat_Tree_Fig.pdf`, `Locomotion_Tree_Fig.pdf` | Ancestral-state reconstructions on the phylogeny |

## Reproducing the analysis

Open `EmberizoidHindlimbs.Rproj` in RStudio, set the working directory to `Blackbird_analyses/` (the scripts use relative paths), install the packages below and run the scripts in order:

1. `0_BlackBirdData Reformatting.R`
2. `1_BlackbirdPCA and Exploration.R`

Before running script 1, point the `read.nexus()` call at `Icteridae_SpeciesTree_MCC.tre` in this folder; it currently uses an absolute path from the original machine. The `write.csv()` calls in script 1 are commented out, so the provided CSVs are not overwritten.

```r
install.packages(c("ape", "phytools", "dplyr", "phylolm", "geomorph", "caper",
                   "geiger", "car", "mvMORPH", "vioplot", "ggplot2", "plyr"))
```

## Author

**Kurt Riggin**: [GitHub](https://github.com/kriggithub) · [ORCID](https://orcid.org/0009-0004-4700-1251)
