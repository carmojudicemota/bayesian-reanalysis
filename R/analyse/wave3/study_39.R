study_39_claim_ids <- function() {
  c("study_39_claim_01","study_39_claim_02")
}


load_study_39_raw <- function(path = "data/raw/study_39/Open_Pedagogy_Student_Perceptions.sav") {
  if (!file.exists(path)) {path <- "data/raw/study_39/Untitled3.sav"}
  haven::read_sav(path)
}


study_39_find_columns <- function(raw) {
  labels <- vapply(raw,
    function(x) {
      label <- attr(x, "label")
      if (is.null(label)) {""} else {as.character(label)}
    },
    character(1)
  )
  
  find_column <- function(name_pattern = NULL,label_pattern = NULL) {
    matches <- character()
    
    if (!is.null(name_pattern)) {
      matches <- c(matches,
                   grep(name_pattern,
                        names(raw),
                        value = TRUE,
                        ignore.case = TRUE
                        )
                   )
    }
    
    if (!is.null(label_pattern)) {
      matches <- c(matches,
                   names(raw)[grepl(label_pattern,labels,ignore.case = TRUE)]
                   )
    }
    unique(matches)[[1]]
  }
  
  list(
    motivation = find_column(
      name_pattern = "^HOWMOTIV",
      label_pattern = "motivating.*final product.*openly available"
    ),
    diversity = find_column(
      name_pattern = "^IBELIEV",
      label_pattern = "photographs add diversity"
    ),
    engagement = find_column(
      name_pattern = "^DIDTHISC",
      label_pattern = "course seem more engaging.*project"
    ),
    flickr = find_column(
      name_pattern = "^(V41|V41_A)$",
      label_pattern = "sharing your photos.*world.*Flickr"
    )
  )
}


study_39_numeric <- function(x) {
  x <- haven::zap_labels(x)
  if (is.factor(x)) {x <- as.character(x)}
  suppressWarnings(as.numeric(x))
}


load_study_39_data <- function(claim_id,path = "data/raw/study_39/Open_Pedagogy_Student_Perceptions.sav") {
  raw <- load_study_39_raw(path)
  columns <- study_39_find_columns(raw)
  variables <- switch(claim_id,
                      study_39_claim_01 = c(columns$motivation,
                                            columns$diversity),
                      study_39_claim_02 = c(columns$engagement,
                                            columns$flickr),
    
    stop("Unknown Study 39 claim: ", claim_id, call. = FALSE)
  )
  
  x <- study_39_numeric(raw[[variables[[1]]]])
  y <- study_39_numeric(raw[[variables[[2]]]])
  keep <- stats::complete.cases(x, y)
  tibble::tibble(x = x[keep],y = y[keep])
}


prepare_study_39_rank_data <- function(data) {
  x_levels <- sort(unique(data$x))
  y_levels <- sort(unique(data$y))
  
  list(
    N = nrow(data),
    K_x = length(x_levels),
    K_y = length(y_levels),
    x_category = match(data$x, x_levels),
    y_category = match(data$y, y_levels),
    x_levels = x_levels,
    y_levels = y_levels
  )
}



