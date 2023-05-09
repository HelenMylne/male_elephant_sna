#### set up ####
# load packages
library(tidyverse)
library(dplyr)
library(cmdstanr)
library(bisonR)

# set seed
set.seed(12345)

#### create model inputs ####
### import data for aggregated model (binomial) -- counts of positive associations and total sightings
counts_df <- read_csv('../data_processed/motnp_binomialpairwiseevents_malesonly.csv')
str(counts_df)

### create data frame for edge weight model -- select only nodes, times together and total sightings
counts_df_model <- counts_df[, c('node_1_males','node_2_males','event_count','count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### set priors
priors <- get_default_priors('binary') # obtain structure for bison model priors
prior_check(priors, 'binary')

### create dummy variable indicating if ever seen together or not -- use to select which prior to draw edge weight from
counts_df_model$seen_together <- ifelse(counts_df_model$event > 0, 1, 0)

### new prior for zero-inflated dyads (never seen together)
plot(density( LaplacesDemon::invlogit(rnorm(1000, -5, 3)) ), col = 'blue')      # new -- dyads never sighted together
lines(density(LaplacesDemon::invlogit(rnorm(1000, -2.5, 1.5)) ), col = 'red')   # original -- dyads ≥ 1 seen together 

#### run function codes needed for models below to run -- directly copied from bisonR GitHub ####
### extract_distribution
extract_distribution <- function(prior_string) {
  parameter_string <- stringr::str_replace_all(prior_string, " ", "")
  parameter_split <- stringr::str_split(parameter_string, "\\(|\\)|,")[[1]]
  distribution_name <- parameter_split[1]
  parameter_values <- parameter_split[2:(length(parameter_split) - 1)]
  parameter_values <- as.numeric(parameter_values)
  return(list(distribution_name=distribution_name, parameter_values=parameter_values))
}

### extract_prior_parameters
extract_prior_parameters <- function(priors) {
  prior_parameters <- list()
  for (parameter_name in names(priors)) {
    prior_distribution <- extract_distribution(priors[parameter_name])
    distribution_name <- prior_distribution$distribution_name
    parameter_values <- prior_distribution$parameter_values
    if (distribution_name == "normal") {
      prior_parameters[[paste("prior", parameter_name, "mu", sep="_")]] <- parameter_values[1]
      prior_parameters[[paste("prior", parameter_name, "sigma", sep="_")]] <- parameter_values[2]
    }
    if (distribution_name == "half-normal") {
      prior_parameters[[paste("prior", parameter_name, "sigma", sep="_")]] <- parameter_values[1]
    }
    if (distribution_name == "beta") {
      prior_parameters[[paste("prior", parameter_name, "alpha", sep="_")]] <- parameter_values[1]
      prior_parameters[[paste("prior", parameter_name, "beta", sep="_")]] <- parameter_values[2]
    }
    if (distribution_name == "gamma") {
      prior_parameters[[paste("prior", parameter_name, "alpha", sep="_")]] <- parameter_values[1]
      prior_parameters[[paste("prior", parameter_name, "beta", sep="_")]] <- parameter_values[2]
    }
  }
  # print(prior_parameters)
  prior_parameters
}

### get_bison_model_data
get_bison_model_data <- function(formula, observations, directed, model_type) {
  # Get model specification from formula
  model_spec <- get_bison_model_spec(formula)
  
  # Automatically detect and apply factor levels
  if (!is.null(model_spec$node_1_name)) {
    node_1_names <- dplyr::pull(observations, model_spec$node_1_name)
    node_2_names <- dplyr::pull(observations, model_spec$node_2_name)
    if (!is.factor(node_1_names) || !all(levels(node_1_names) == levels(node_2_names))) {
      unique_node_names <- sort(unique(c(node_1_names, node_2_names)))
      observations[, model_spec$node_1_name] <- factor(node_1_names, levels=unique_node_names)
      observations[, model_spec$node_2_name] <- factor(node_2_names, levels=unique_node_names)
    }
  }
  
  if (!is.null(model_spec$node_1_name)) {
    # Get node to index mapping
    node_to_idx <- 1:length(levels(dplyr::pull(observations, model_spec$node_1_name)))
    names(node_to_idx) <- levels(dplyr::pull(observations, model_spec$node_1_name))
    
    # Get number of nodes
    num_nodes <- length(node_to_idx)
    
    # Get dyad to index mapping
    
    # Get unique dyads in both directions
    unique_dyads_1 <- unique(cbind(
      node_to_idx[dplyr::pull(observations, model_spec$node_1_name)],
      node_to_idx[dplyr::pull(observations, model_spec$node_2_name)]
    ))
    unique_dyads_2 <- unique(cbind(
      node_to_idx[dplyr::pull(observations, model_spec$node_2_name)],
      node_to_idx[dplyr::pull(observations, model_spec$node_1_name)]
    ))
    if (directed == FALSE) {
      dyad_to_idx <- unique_dyads_1
      rownames(dyad_to_idx) <- 1:nrow(dyad_to_idx) # Assign dyad IDs
      
      # Generate dyad names
      dyad_names <- sapply(
        1:nrow(dyad_to_idx),
        function(x) paste0(names(node_to_idx)[dyad_to_idx[x, ]], collapse=" <-> ")
      )
    } else {
      dyad_to_idx = unique(rbind(unique_dyads_1, unique_dyads_2))
      rownames(dyad_to_idx) <- 1:nrow(dyad_to_idx) # Assign dyad IDs
      
      # Generate dyad names
      dyad_names <- sapply(
        1:nrow(dyad_to_idx),
        function(x) paste0(names(node_to_idx)[dyad_to_idx[x, ]], collapse=" -> ")
      )
    }
    num_dyads <- nrow(dyad_to_idx)
  } else {
    dyad_ids = rep(1, nrow(observations))
    node_1_names = NULL
    node_2_names = NULL
    node_to_idx = NULL
    dyad_to_idx = NULL
    dyad_names = NULL
    num_nodes = 0
    num_dyads = 0
  }
  
  # Get events
  event <- dplyr::pull(observations, model_spec$event_var_name)
  
  # Create empty design matrices for fixed (design_fixed) and random (design_random) effects
  design_fixed <- data.frame(empty_col = rep(0, nrow(observations)))
  design_random <- data.frame(empty_col = rep(0, nrow(observations)))
  
  # Add intercept term if needed
  if (model_spec$intercept) {
    design_fixed[, "intercept"] <- 1
  }
  
  # If using dyad-level weights, populate design matrix
  if (!is.null(model_spec$node_1_name)) {
    # Get dyad IDs in the correct order
    node_1_names <- dplyr::pull(observations, model_spec$node_1_name)
    node_2_names <- dplyr::pull(observations, model_spec$node_2_name)
    
    dyad_ids = as.factor(get_dyad_ids(
      node_to_idx[node_1_names],
      node_to_idx[node_2_names],
      dyad_to_idx,
      directed=directed
    ))
  }
  
  # Variable grouping for random effects
  random_group_index <- c()
  
  # Get additional fixed effects
  if (!is.null(model_spec$fixed)) {
    for (term_name in model_spec$fixed) {
      # If it's a factor, create a column for each level.
      if (is.factor(dplyr::pull(observations, term_name))) {
        var_group <- paste0("fixed_", term_name)
        term_levels <- levels(dplyr::pull(observations, term_name))
        for (term_level in term_levels) {
          new_term_name <-  paste0("fixed_", term_name, term_level)
          design_fixed[, new_term_name] <- 1 * (dplyr::pull(observations, term_name) == term_level)
        }
      } else {
        # Otherwise, create a single column:
        new_term_name <- paste0(c("fixed_", term_name), collapse="")
        design_fixed[, new_term_name] <- dplyr::pull(observations, term_name)
      }
    }
  }
  
  # Get additional random effects
  if (!is.null(model_spec$random)) {
    for (term_name in model_spec$random) {
      var_group <- paste0("random_", term_name)
      term_levels <- levels(as.factor(dplyr::pull(observations, term_name)))
      for (term_level in term_levels) {
        new_term_name <-  paste0("random_", term_name, term_level)
        design_random[, new_term_name] <- 1 * (as.factor(dplyr::pull(observations, term_name)) == term_level)
        random_group_index[length(random_group_index) + 1] <- var_group
      }
    }
  }
  
  # Check if divisor is a real
  if (!is.na(str_match(model_spec$divisor_var_name, "\\d+"))) {
    divisor <- rep(as.numeric(model_spec$divisor_var_name), nrow(observations))
  } else {
    divisor <- dplyr::pull(observations, model_spec$divisor_var_name)
  }
  
  # If divisor is a single value, convert it to a list.
  if (length(divisor) == 1) {
    divisor <- rep(divisor, nrow(observations))
  }
  
  model_data <- list(
    num_rows=nrow(design_fixed),
    event=event,
    divisor=divisor,
    dyad_ids=dyad_ids,
    design_fixed=data.matrix(design_fixed[, -1]),
    design_random=data.matrix(design_random[, -1]),
    num_edges = num_dyads,
    num_fixed=ncol(as.matrix(design_fixed[, -1])),
    num_random=ncol(as.matrix(design_random[, -1])),
    num_random_groups=length(unique(random_group_index)),
    random_group_index=as.integer(as.factor(random_group_index))
  )
  
  obj <- list(
    model_data=model_data,
    node_to_idx=node_to_idx,
    dyad_to_idx=dyad_to_idx,
    dyad_names=dyad_names,
    num_nodes=num_nodes,
    num_dyads=num_dyads,
    row_dyad_ids=dyad_ids
  )
  
  return(obj)
}

### get_bison_model_spec
get_bison_model_spec <- function(formula) {
  model_spec <- list()
  
  x <- str_split(deparse1(formula), "~")[[1]]
  lhs <- x[1]
  rhs <- x[2]
  
  # Process left hand side
  lhs_split <- str_split(lhs, "\\|")[[1]]
  event_var_name <- lhs_split[1]
  event_var_name <- str_replace_all(event_var_name, "\\(", "")
  event_var_name <- str_replace_all(event_var_name, " ", "")
  model_spec$event_var_name <- event_var_name
  
  divisor_var_name <- lhs_split[2]
  divisor_var_name <- str_replace_all(divisor_var_name, "\\)", "")
  divisor_var_name <- str_replace_all(divisor_var_name, " ", "")
  model_spec$divisor_var_name <- divisor_var_name
  
  # Set intercept to false by default
  model_spec$intercept <- FALSE
  
  model_spec$fixed <- c()
  model_spec$random <- c()
  
  rhs_split <- str_split(rhs, "\\+")[[1]]
  for (term in rhs_split) {
    term <- str_replace_all(term, " ", "")
    # Is it an intercept, a dyad, a fixed effect, or a random effect?
    if (!is.na(str_match(term, "^0|1$"))) {
      # Intercept
      if (term == "0") {
        model_spec$intercept <- FALSE
      } else {
        model_spec$intercept <- TRUE
      }
    } else if (!is.na(str_match(term, "^dyad\\(.*,.*\\)$")[[1]])) {
      # dyad(,) term
      node_names <- str_split(term, "\\(|\\)")[[1]][2]
      node_names <- str_replace_all(node_names, " ", "")
      node_names_split <- str_split(node_names, ",")[[1]]
      model_spec$node_1_name <- node_names_split[1]
      model_spec$node_2_name <- node_names_split[2]
    } else if (is.na(str_match(term, "[^a-zA-Z0-9_]"))) {
      # Fixed effect
      model_spec$fixed[length(model_spec$fixed) + 1] <- term
    } else if (!is.na(str_match(term, "^\\(1\\|.*\\)$"))) {
      # Random intercept
      term_name <- str_split(term, "\\(|\\||\\)")[[1]][3]
      model_spec$random[length(model_spec$random) + 1] <- term_name
    } else {
      warning(paste0("Formula term \"", term, "\" not supported by bisonR. Check the formula is correctly specified."))
    }
  }
  
  model_spec
}

### get_dyad_ids
get_dyad_ids <- function(node_id_1, node_id_2, dyad_to_idx, directed) {
  dyad_ids <- rep(NA, length(node_id_1))
  for (i in 1:length(node_id_1)) {
    if (directed == TRUE) {
      dyad_ids[i] <- which(dyad_to_idx[, 1] == node_id_1[i] & dyad_to_idx[, 2] == node_id_2[i])[[1]]
    }
    if (directed == FALSE) {
      dyad_ids[i] <- which(
        (dyad_to_idx[, 1] == node_id_1[i] & dyad_to_idx[, 2] == node_id_2[i]) |
          (dyad_to_idx[, 1] == node_id_2[i] & dyad_to_idx[, 2] == node_id_1[i])
      )[[1]]
    }
  }
  dyad_ids
}

#### obtain Stan code for working model ####
### run OLD edge weight model -- COMMENTED OUT TO AVOID RUNNING UNECESSARILY, BUT HIGHLIGHTS WHERE I GOT STAN CODE FROM TO ADAPT
#motnp_edge_weights_strongpriors <- bison_model(
#  ( event | duration ) ~ dyad(node_1_id, node_2_id)
#  data = counts_df_model,
#  model_type = "binary",
#  priors = priors
#)
#motnp_edge_weights_strongpriors$stan_model$code()

#### obtain and adapt prior parameters ####
### extract_prior_parameters() -- function code copied directly from bisonR GitHub
### extract prior parameters
prior_parameters <- extract_prior_parameters(priors)

### adapt prior parameters so there are two possible distributions:
prior_parameters$prior_edge_mu_0 <- -5
prior_parameters$prior_edge_mu_1 <- -2.5
prior_parameters$prior_edge_sigma_0 <- 3
prior_parameters$prior_edge_sigma_1 <- 1.5
prior_parameters <- prior_parameters[c(10:13,3:9)]
prior_parameters

#### adapt bison_model() function code to use new Stan code ####
bison_model_hkm_2priors <- function(formula, data, model_type=c("binary", "count"),
                                    directed=FALSE, partial_pooling=FALSE, zero_inflated=FALSE,
                                    prior_parameters=NULL,                        # INPUT PRIOR PARAMETERS DIRECTLY RATHER THAN AS FULL PRIOR LIST
                                    refresh=0, mc_cores=4, iter_sampling=1000,
                                    iter_warmup=1000, priors_only=FALSE) {
  
  # If user-specified priors haven't been set, use the defaults
  if (is.null(prior_parameters)) {
    message("No priors set by user, using default priors instead. We recommend setting and checking priors explicitly for reliable inference.")
    priors <- get_default_priors(model_type)
    prior_parameters <- extract_prior_parameters(priors)                  # ALTERED INPUT SO ONLY USES DEFAULT PRIOR PARAMETERS IF NONE ARE PROVIDED
  }
  
  # If the model is a conjugate model, change the data type
  use_conjugate_model <- FALSE
  if (length(str_match_all(model_type, "conjugate")[[1]]) > 0) {
    use_conjugate_model <- TRUE
    model_type <- str_split(model_type, "_")[[1]][1]
    if (!(model_type %in% c("binary", "count"))) {
      stop("Conjugate models only available for binary and count data")
    }
  }
  
  model_spec <- get_bison_model_spec(formula)
  
  model_info <- get_bison_model_data(formula, data, directed, model_type)
  
  # Set the priors in model data
  model_data <- c(model_info$model_data, prior_parameters)
  # Set whether only the priors should be sampled
  model_data$priors_only <- priors_only
  model_data$partial_pooling <- (partial_pooling * 1)
  model_data$zero_inflated <- zero_inflated
  
  # Select sampling method
  if (use_conjugate_model) {
    # Fit conjugate model
    model <- fit_conjugate_model(model_type, model_data, num_samples=iter_sampling*mc_cores, priors_only=priors_only)
    chain <- model$chain
    event_preds <- model$event_preds
    fit <- NULL
    log_lik <- NULL
  } else {
    # Build Stan model
    model <- build_stan_model_hkm(model_type)       # ALTERED THE CODE FOR THIS TO ACCEPT NEW STAN CODE
    
    # Fit model
    fit <- model$sample(
      data=model_data,
      refresh=refresh,
      chains=4,
      parallel_chains=mc_cores,
      step_size=0.1,
      iter_sampling=iter_sampling,
      iter_warmup=iter_warmup
    )
    
    # Extract edge weights from fitted edge model.
    if (model_info$num_dyads > 0) {
      chain <- fit$draws("edge_weight", format="matrix")
    } else {
      chain <- NULL
    }
    
    event_preds <- fit$draws("event_pred", format="matrix")
    
    log_lik <- fit$draws("log_lik", format="matrix")
  }
  
  # Extract edge samples
  edge_samples <- chain
  
  # Prepare output object.
  obj <- list(
    chain = chain,
    event_preds = event_preds,
    edge_samples=edge_samples,
    num_nodes = model_info$num_nodes,
    num_dyads = model_info$num_dyads,
    node_to_idx = model_info$node_to_idx,
    dyad_to_idx = model_info$dyad_to_idx,
    dyad_names = model_info$dyad_names,
    model_info = model_info,
    directed = directed,
    fit = fit,
    formula = formula,
    model_type = model_type,
    model_data = model_data,
    input_data = data,
    stan_model = model,
    conjugate = use_conjugate_model,
    log_lik = log_lik
  )
  class(obj) <- "bison_model"
  return(obj)
}

build_stan_model_hkm <- function(model_name) {
  model_filepath <- "models/two_prior_edge_model.stan"      # USE CUSTOM STAN CODE INSTEAD OF WRITING NEW CODE -- REPLACES: model_filepath <- system.file("stan", paste0(model_name, ".stan"), package="bisonR")
  model <- cmdstanr::cmdstan_model(model_filepath, compile=FALSE, stanc_options = list("O1"))
  model$compile(dir=tempdir())
  return(model)
}

#### run new bison model using 2 priors ####
motnp_edge_weights_strongpriors <- bison_model_hkm_2priors(
  ( event | duration ) ~ dyad(node_1_id, node_2_id),
  data = counts_df_model,
  model_type = "binary",
  prior_parameters = prior_parameters
)

#### alternative attempt where set prior with Stan code ####
### adapt prior for different types of pair
priors$edge
#priors$edge <- c('normal(-2.5, 1.5', 'normal(-5, 3)') # have also tried various combinations of c('normal(c(-2.5, 5), c(1.5, 3))) etc. but same issues of bisonR can't work out where to split it and also even if it could there is nothing in here about when to use and when to use the other
priors$edge <- brms::prior_string("
real conditional_prior(real mu, int seen_together) {
    if (seen_together == 1) {
      return normal(mu | -2.5, 1.5);
    } else {
      return normal(mu | -5, 3);
    }
  }
")

## run new bison model using 2 priors
motnp_fit_edgeweights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id)*seen_together,
  data = counts_df_model,
  model_type = "binary",
  priors = priors
)

