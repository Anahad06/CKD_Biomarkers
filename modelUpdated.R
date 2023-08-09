# Removes 'e'
options(scipen = 999)

library(GEOquery)
library(preprocessCore)
library(sva)

# Get DataSet from GEO Database
gse <- getGEO("GSE66494", GSEMatrix = TRUE)

# Get the Labels (possible future use)
gseData <- pData(gse[[1]])

# creating labels for the train dataset 
statustrain <- c(
  rep("CKD", 47),
  rep("N", 5)
)

# creating labels for the test dataset
statustest <- c(
  rep("CKD", 5),
  rep("N", 3)
)

# Import Dataset and Save
Expr_GSM1623299 <- exprs(gse[[1]])

# Set column names for Train, Control, and Test datasets
sample_names <- colnames(Expr_GSM1623299)
Train_cols <- sample_names[2:53]
Test_cols <- sample_names[54:ncol(Expr_GSM1623299)]

# Splice Data
Train <- Expr_GSM1623299[, Train_cols, drop = FALSE]
Test <- Expr_GSM1623299[, Test_cols, drop = FALSE]

# Normalize Data
Train_normalized <- normalize.quantiles(Train)
Test_normalized <- normalize.quantiles(Test)

# Set the row and column names for the normalized data
rownames(Train_normalized) <- rownames(Train)
colnames(Train_normalized) <- colnames(Train)

rownames(Test_normalized) <- rownames(Test)
colnames(Test_normalized) <- colnames(Test)

# Log2 transform the normalized data
Train_normalized_log2 <- log2(Train_normalized)
Test_normalized_log2 <- log2(Test_normalized)

# Add disease status as row names for the Train_normalized_log2 dataset
colnames(Train_normalized_log2) <- paste(colnames(Train_normalized_log2), statustrain, sep = "")

# Add disease status as row names for the Test_normalized_log2 dataset
colnames(Test_normalized_log2) <- paste(colnames(Test_normalized_log2), statustest, sep = "")

# shuffles the dataset 
set.seed(123)  # For reproducibility
random_order_train <- sample(1:ncol(Train_normalized_log2))
Train_normalized_log2_randomized <- Train_normalized_log2[, random_order_train]

# randomizing the train labels with the dataset (future use)
status_Trainshuffled = statustrain[random_order_train]

# shuffles the dataset (column)
set.seed(456)  # For reproducibility
random_order_test <- sample(1:ncol(Test_normalized_log2))
Test_normalized_log2_randomized <- Test_normalized_log2[, random_order_test]

# randomizing the test labels with the dataset (future use)
status_Testshuffled = statustest[random_order_test]

# comparing
View(Train_normalized_log2)
View(Train_normalized_log2_randomized)
print(status_Trainshuffled)


