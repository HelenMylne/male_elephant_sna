library(cmdstanr)
library(igraph)
library(brms)
library(purrr)

# METHOD 1: Extract posterior edge weights and do multiple imputation regression

# Compile your Stan model
model <- cmdstan_model("your_model.stan")

# Assume `data_list` contains your data for the Stan model
fit <- model$sample(data = data_list)

# Extract posterior samples
posterior_samples <- fit$draws()

# Convert the array to a matrix
edge_weights_matrix <- posterior_samples$variables()$edge_weight

# Initialize a list to store the centrality measures
centrality_list <- list()

# Initialize a list to store your datasets for the brms model
brms_data_list <- list()

# Loop over each set of edge weights in the posterior
for (i in 1:nrow(edge_weights_matrix)) {
  # Get the edge weights for this set
  edge_weights <- edge_weights_matrix[i, ]

  # Construct a social network based on the dyads and edge weights
  # You'll need to fill in the details of how to do this based on your data
  g <- graph_from_data_frame(dyads, directed=FALSE, vertices=NULL)
  E(g)$weight <- edge_weights

  # Calculate the eigenvector centrality
  centrality <- evcent(g, weights = E(g)$weight)$vector

  # Store the centrality measures in the list
  centrality_list[[i]] <- centrality

  # Create a dataset for the brms model
  # Assume `age` is a vector of the ages of the individuals
  brms_data_list[[i]] <- data.frame(centrality = centrality, age = age)
}

# Define the formula
formula <- bf(centrality ~ age)

# Fit the models and combine the results
brms_fit <- brm_multiple(formula, brms_data_list, chains = 4, cores = 4)

# Print the results
print(brms_fit)

# Method 2: Fit Gaussians to the logit edge weights
# See example here: https://github.com/JHart96/bison_examples/blob/main/examples/dyadic_regression_stan.md
# Advantage: Quicker to fit and might solve the problem you're experiencing on Viking
# Disadvantage: Maybe not as accurate as the multiple imputation approach (but should be fine)







# And to do a dyadic regression
mm_model <- brms_multiple(
  y ~ 1 + ( 1|mm(id1, id2) ),
  data = df,
  family = gaussian()
)
