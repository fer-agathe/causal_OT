#' Define the function to split the dataset into train/test sets
#'
#' @param data Dataset.
#' @param seed Seed number to use for replication, if provided.
#' @param train_ratio Ratio of data to use in the train set (default to 0.7).
#'
#' @returns A list with two elements:
#' * `data_train`: The train set.
#' * `data_test`: The test set.
#' @export
#' @md
#' @seealso [simul_dataset()]
#' @examples
#' df <- simul_dataset()
#' sets <- split_dataset(df)
#' data_train <- sets$data_train
#' data_test <- sets$data_test
split_dataset <- function(data, seed = NULL, train_ratio = 0.7) {
  if (!is.null(seed)) set.seed(seed)
  train_size <- floor(train_ratio * nrow(data))
  train_indices <- sample(seq_len(nrow(data)), replace = TRUE, size = train_size)
  # Split the data into train and test sets
  train_set <- data[train_indices, ]
  test_set <- data[-train_indices, ]
  # Return the train and test sets as a list
  list(
    data_train = train_set,
    data_test = test_set
  )
}
