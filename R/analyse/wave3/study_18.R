study_18_claim_ids <- function() {"study_18_claim_01"}

load_study_18_data <- function(path = "data/raw/study_18/Psychological_literacy_subject_study_finaldataset.sav") {
  raw <- haven::read_sav(path)
  
  skill_items <- list(
    skill_1 = c("Skill_1_Awareness","Skill_1_Development",
                "Skill_1_Confidence","Skill_1_Importance"),
    skill_2 = c("Skill_2_Awareness","Skill_2_Development",
                "Skill_2_Confidence","Skill_2_Importance"),
    skill_3 = c("Skill_3_Awareness","Skill_3_Development",
                "Skill_3_Confidence","Skill_3_Importance"),
    skill_4 = c("Skill_4_Awareness","Skill_4_Development",
                "Skill_4_Confidence","Skill_4_Importance"),
    skill_5 = c("Skill_5_Awareness","Skill_5_Development",
                "Skill_5_Confidence","Skill_5_Importance"),
    skill_6 = c("Skill_6_Awareness","Skill_6_Development",
                "Skill_6_Confidence","Skill_6_Importance"),
    skill_7 = c("Skill_7_Awareness","Skill_7_Development",
                "Skill_7_Confidence","Skill_7_Importance"),
    skill_8 = c("Skill_8_Awareness","Skill_8_Development",
                "Skill_8_Confidence","Skill_8_Importance")
    )
  
  data <- raw
  
  for (skill in names(skill_items)) {
    data[[skill]] <- rowMeans(as.data.frame(data[, skill_items[[skill]]]),na.rm = FALSE)
  }
  
  data |>
    dplyr::transmute(
      subject_group = factor(
        haven::zap_labels(STEM_NONSTEM_PSYCH),
        levels = c(1, 2, 3),
        labels = c("psychology","stem","humanities")
      ),
      skill_1,
      skill_2,
      skill_3,
      skill_4,
      skill_5,
      skill_6,
      skill_7,
      skill_8
    ) |>
    dplyr::filter(
      stats::complete.cases(
        subject_group,
        skill_1,
        skill_2,
        skill_3,
        skill_4,
        skill_5,
        skill_6,
        skill_7,
        skill_8
      )
    )
}


fit_study_18_model <- function(data) {
  stats::lm(
    cbind(
      skill_1,
      skill_2,
      skill_3,
      skill_4,
      skill_5,
      skill_6,
      skill_7,
      skill_8
    ) ~ subject_group,
    data = data
  )
}


study_18_hypothesis <- function() {
  parameters <- c(paste0("subject_groupstem_on_skill_", 1:8),
                  paste0("subject_grouphumanities_on_skill_", 1:8)
                  )
  paste0(parameters, " = 0",collapse = " & ")
}

compute_study_18_bayes_factors <- function(claim,priors = NULL,iter = 100000) {
  data <- load_study_18_data()
  model <- fit_study_18_model(data)
  wave3_check_assumptions(model, claim$claim_id)
  set.seed(123)
  
  result <- BFpack::BF(model,hypothesis = study_18_hypothesis(),iter = iter)
  bf10 <- as.numeric(result$BFmatrix_confirmatory["H2", "H1"])
  
  wave3_row(
    claim = claim,
    bf10 = bf10,
    model_null = "all subject-group coefficients = 0",
    model_alt = "at least one subject-group coefficient != 0",
    bf_family = "generalized_fractional",
    prior_family = "fractional",
    method = "BFpack_generalized_fractional",
    prior_label = "primary"
  )
}

