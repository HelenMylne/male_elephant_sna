lapply(c("foreach", "doParallel"), require, character.only = TRUE)

library(parallel)
cl <- makeCluster(detectCores()) # can replace detectCores() with a number if known
registerDoParallel(cl)

# create empty data frame
nrows <- (length(unique(ate$id))-1) * nrow(ate)
gbi_df <- data.frame(node_1 = rep(NA, nrows),
                     node_2 = rep(NA, nrows),
                     social_event = rep(NA, nrows),
                     obs_id = rep(NA, nrows),
                     block = rep(NA, nrows))


# ============== UNPARALLELISED VERSION ==================

# loop through each block and extract the corresponding data ####
for (block in sort(unique(ate$block))) {
  block_data <- ate[ate$block %in% block,]
  for (obs_id in min(block_data$observation_number):max(block_data$observation_number)) {
    for (i in which(gbi_matrix[obs_id, ] == 1)) {
      for (j in 1:ncol(gbi_matrix)) {
        if (i != j) {
          node_1 <- ifelse(i < j, i, j)
          node_2 <- ifelse(i < j, j, i)
          
          empty_row <- which(is.na(gbi_df$node_1) == TRUE)[1]
          gbi_df[empty_row, ]$node_1 <- node_1
          gbi_df[empty_row, ]$node_2 <- node_2
          gbi_df[empty_row, ]$social_event <- ifelse(gbi_matrix[obs_id, i] == gbi_matrix[obs_id, j], 1, 0)
          gbi_df[empty_row, ]$obs_id <- obs_id
          gbi_df[empty_row, ]$block <- block
        }
      }
    }
    #if(obs_id %% 10 == 0) {print(obs_id)}
    #if(obs_id %% 10 == 0) {print(Sys.time())}
    print(paste0('obs_id ', obs_id, ' in block ', block, ' finished at time ', Sys.time()))
  }
}
gbi_df

# ============PARALLELISED VERSION ====================

# Parallelise the blocks

# Function to be parallelized
process_blocks <- function(block_data)  {
  num_rows <- max(block_data$observation_number) - min(block_data$observation_number)
  # Pre-allocate list (might need to make dataframe)
  rows <- list( node_1 = rep(NA, num_rows),
                node_2 = rep(NA, num_rows),
                social_event = rep(NA, num_rows),
                obs_id = rep(NA, num_rows),
                block = rep(NA, num_rows)
                )

  for (obs_id in min(block_data$observation_number):max(block_data$observation_number)) {
    row_index <- 1
    for (i in which(gbi_matrix[obs_id, ] == 1)) {
      for (j in 1:ncol(gbi_matrix)) {
        if (i != j) {
          node_1 <- ifelse(i < j, i, j)
          node_2 <- ifelse(i < j, j, i)
          
          empty_row <- which(is.na(gbi_df$node_1) == TRUE)[1]
          rows[row_index] <- list(node_1 = node_1,
                                  node_2 = node_2,
                                  social_event = ifelse(gbi_matrix[obs_id, i] == gbi_matrix[obs_id, j], 1, 0),
                                  obs_id = obs_id,
                                  block = block)
        }
      }
    }
    row_index <- row_index + 1
    #if(obs_id %% 10 == 0) {print(obs_id)}
    #if(obs_id %% 10 == 0) {print(Sys.time())}
    #print(paste0('obs_id ', obs_id, ' in block ', block, ' finished at time ', Sys.time()))
  }
  return(rows)
}


sorted_unique_blocks <- sort(unique(ate_test$block))
num_unique_blocks <- length(sorted_unique_blocks)
foreach(block_id = 1:num_unique_blocks) %dopar% {
  block <- sorted_unique_blocks[block_id]
  block_data <- ate_test[ate_test$block %in% block,]
  process_blocks(block_data)
}
gbi_df
