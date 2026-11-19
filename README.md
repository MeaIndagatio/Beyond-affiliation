# Beyond-affiliation: code and output

This repository contains the code and output for a latent profile analysis of religious and spiritual belief and practice among U.S. adults, using the 2023–24 Religious Landscape Study (RLS).

## Data

The data are not included. The analysis uses the public use file of the 2023–24 Religious Landscape Study, available free of charge from Pew Research Center (https://www.pewresearch.org) after registration and acceptance of Pew's terms of use. Save the file in Stata format as `RLS Public Use File_1_original.dta` in a folder of your choice.

## Software

- Stata, for data preparation.
- R, for all analyses, with the packages haven, dplyr, survey, lavaan, semTools, psych, tidySEM, OpenMx and car.

## How to run the analysis

1. In the Stata do-file and in every R script, set the data folder (`data_dir` in R) to the folder that contains the data file. The R scripts write all output to a subfolder of it, `R output`.
2. Run the Stata do-file. It creates `RLS Public Use File_2_Clean.dta`.
3. Run the R scripts in numerical order. Each script reads the files saved by the earlier ones.

The models in Stage 3 take many hours to estimate. Each fitted model is saved to its own file as soon as it is estimated, and a model that has already been saved is skipped if the script is run again.

## Files

### Data preparation (Stata)

`RLS Public Use File_2_cleaning syntax.do`
- Sets refusals and "don't know" answers to missing.
- Reverse-codes the items so that higher values indicate more frequent practice or stronger belief.
- Builds the four indices (religious behavior, religious belief, spiritual behavior, spiritual belief), each computed only for respondents with complete data on all of its items.
- Codes current and raised religious tradition in nine categories, the demographic and political covariates, the childhood religiosity measures, and the spouse or partner measures.
- Flags the analytic sample (complete data on all four indices, n = 31,470) and divides it into two random halves.

### Analysis (R)

| Script | Content | Main output files |
| --- | --- | --- |
| `0. Setup and sample.R` | Analytic sample; item missingness; survey design with the 500 replicate weights; comparison of included and excluded respondents | `stage0_missingness.csv`, `stage0_sample_comparison.csv` |
| `1. Measurement.R` | Item distributions; reliability of the four indices. The script also contains exploratory and confirmatory factor analyses and invariance tests, which are not reported in the article | `stage1_reliability_full_sample.csv` |
| `2. Indicator description.R` | Descriptive statistics, correlations and histograms of the four indices | `stage2_descriptives.csv`, `stage2_weighted_descriptives.csv`, `stage2_correlations_weighted.csv`, `stage2_histograms.png` |
| `3. Class enumeration.R` | Latent profile models with one to eight profiles under specifications A and C; fit statistics | `stage3_fit.csv`, `stage3_bic.png` |
| `4. Profiles.R` | Classification quality and description of the selected five-profile model: average posterior probabilities, odds of correct classification, within-profile correlations, profile means with confidence intervals | `stage4_classification_A5.csv`, `stage4_occ_A5.csv`, `stage4_within_profile_cor_A5.csv`, `stage4_profile_means_ci_A5.csv`, `stage4_floor_ceiling.csv` |
| `5. Prevalence.R` | Weighted share of each profile; profile composition of each current and raised religious tradition; religious composition of each profile | `stage5_weighted_shares.csv`, `stage5_profiles_within_traditions.csv`, `stage5_profiles_within_raised.csv`, `stage5_traditions_within_profiles.csv`, `stage5_raised_within_profiles.csv` |
| `6. Robustness.R` | Split-half replication; estimation on all respondents with any indicator data; survey-weighted estimation; adjacent four- and six-profile models; models with three indices | `stage6_summary.csv`, `stage6_profile_means.csv`, `stage6_adjacent_profiles.csv` |
| `7. Covariates.R` | Multinomial logistic regressions of profile membership on the covariates (Models 1 to 3), by the maximum likelihood three-step method with 500 replicate weights; collinearity diagnostics; average predicted percentages | `stage7_rrr_M1.csv`, `stage7_rrr_M2.csv`, `stage7_rrr_M3.csv`, `stage7_model_fit.csv`, `stage7_average_marginal_effects_M3.csv`, `stage7_gvif.csv`, `stage7_classification_matrix.csv` |
| `8. Outcomes.R` | Reported lifetime change in religiousness and spirituality by profile, corrected for classification error by the BCH method. The script also contains wellbeing analyses, which are not reported in the article | `stage8_lifetime_change.csv` |
| `9. Validation.R` | Bivariate associations of each covariate with profile membership. The script also contains a wellbeing analysis, which is not reported in the article | `stage9_bivariate_categorical.csv`, `stage9_bivariate_numeric.csv`, `stage9_covariates_within_profiles.csv` |

### Output files and the tables in the article

| Article table | Output files |
| --- | --- |
| Table 1 | `stage1_reliability_full_sample.csv` |
| Table 2 | `stage3_fit.csv` |
| Table 3 | `stage4_profile_means_ci_A5.csv`, `stage5_weighted_shares.csv`, `stage4_classification_A5.csv`, `stage4_occ_A5.csv` |
| Table 4 | `stage5_profiles_within_traditions.csv` |
| Table 5 | `stage7_average_marginal_effects_M3.csv` |
| Table 6 | `stage8_lifetime_change.csv` |
| Supplementary Table S1 | `stage3_fit.csv` |
| Supplementary Table S2 | `stage6_summary.csv` |
| Supplementary Table S3 | `stage5_profiles_within_traditions.csv`, `stage5_profiles_within_raised.csv` |
| Supplementary Table S4 | `stage5_traditions_within_profiles.csv`, `stage5_raised_within_profiles.csv` |
| Supplementary Table S5 | `stage7_rrr_M3.csv` |
| Supplementary Table S6 | `stage7_average_marginal_effects_M3.csv` |
| Supplementary Table S7 | `stage8_lifetime_change.csv` |
| Supplementary Table S8 | `stage0_sample_comparison.csv` |

Respondent-level files produced by the scripts (`rls_stage0.rds`, `stage4_posteriors.rds` and the saved models) are not included, because they contain rows of the Pew data.

## Notes on the profile indices

- Belief in God and certainty of that belief are combined into one item: 1 = does not believe; 2 to 5 = believes, from "not at all certain" to "absolutely certain".
- Belief in heaven and belief in hell are combined into one afterlife item: 1 = neither, 2 = one of the two, 3 = both.
- Practicing yoga for spiritual reasons is not included in the spiritual behavior index.
- Each index is the raw sum of its items. Ranges: religious behavior 9 to 46; religious belief 4 to 16; spiritual behavior 5 to 25; spiritual belief 5 to 18.
