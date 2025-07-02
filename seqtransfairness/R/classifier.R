#' Fit a Logistic Regression Model to a Binary Outcome With or Without a
#' Sensitive Attribute.
#'
#' @description
#' The aim is to predict a binary outcome (Y) using a set of predictors (X)
#' that may include (aware model) or may not include (unaware model) a sensitive
#' binary attribute (S). A logistic regression is used to to so. It is fitted on
#' a train set. Predicted scores are then returned on both the train and the
#' test set.
#'
#' WARNING
#' This function is used to illustrate the method. We do not recommend
#' using it for other purposes.
#'
#' @param train_data Train set.
#' @param test_data Test set.
#' @param s Name of the sensitive attribute.
#' @param y Name of the target variable.
#' @param type If `"type=aware"`, the model includes the sensitive attributes,
#'        otherwise, if `type=unaware`, it does not.
#'
#' @returns A list with three elements:
#' * `model`: The estimated logistic regression model.
#' * `pred_train`: Estimated scores on the train set.
#' * `pred_test`: Estimated scores on the test set.
#'
#' @export
#'
#' @importFrom dplyr select
#' @importFrom rlang !!
#' @importFrom stats glm predict as.formula binomial
#' @seealso [simul_dataset()], [split_dataset()]
#' @examples
#' # WARNING: this function is used to illustrate the method. We do not recommend
#' # using it for other purposes.
#'
#' # Simulate data
#' df <- simul_dataset()
#' # Split: train (70%) / test (30%)
#' sets <- split_dataset(df, train_ratio = 0.7)
#' data_train <- sets$data_train
#' data_test <- sets$data_test
#' # Logistic regression to predict Y, making the model blind to S
#' pred_unaware <- log_reg_train(
#'   data_train, data_test, s = "S", y = "Y", type = "unaware"
#' )
#' # Fitted model:
#' pred_unaware$model
#' # Predicted values on train set
#' head(pred_unaware$pred_train)
#' # Predicted values on test set
#' head(pred_unaware$pred_test)
log_reg_train <- function(train_data,
                          test_data,
                          s,
                          y,
                          type = c("aware", "unaware")) {
  if (type == "unaware") {
    train_data_ <- train_data |> select(-!!s)
    test_data_ <- test_data |> select(-!!s)
  } else {
    train_data_ <- train_data
    test_data_ <- test_data
  }
  # Train the logistic regression model
  form <- paste0(y, "~.")
  model <- glm(as.formula(form), data = train_data_, family = binomial)
  # Predictions on train and test sets
  pred_train <- predict(model, newdata = train_data_, type = "response")
  pred_test <- predict(model, newdata = test_data_, type = "response")
  list(
    model = model,
    pred_train = pred_train,
    pred_test = pred_test
  )
}
