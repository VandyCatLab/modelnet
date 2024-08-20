# Libraries ----
library(tidyverse)
library(psych)

# Helper functions ----
long2Matrix <- function(data, idCol, corrCol, trialNCol, standardize = FALSE) {
    # Get subjects and trials
    sbjs <- data[idCol] |>
        unique() |>
        pull() |>
        sort()
    trials <- unique(data[trialNCol]) |>
        pull() |>
        as.character()

    # Create a tibble with subject id row
    taskMatrix <- tibble(!!idCol := sbjs)

    for (trial in trials) {
        trialData <- data[data[trialNCol] == trial, c(idCol, corrCol)]
        # Rename corr to the trial
        names(trialData)[2] <- trial

        taskMatrix <- suppressMessages(taskMatrix |>
            left_join(
                trialData
            ))
    }

    # Set the SbjID column as row names
    taskMatrix <- column_to_rownames(taskMatrix, var = idCol)

    # Save raw matrix for later
    rawInfo <- taskMatrix

    # Add mean row
    taskMatrix <- rbind(taskMatrix, mean = apply(taskMatrix, 2, mean))

    # Add participant mean row
    taskMatrix <- cbind(taskMatrix, Sum = apply(taskMatrix, 1, sum))

    # Add correlation row
    totals <- taskMatrix[seq_len(nrow(taskMatrix) - 1), "Sum"]
    tmp <- c()
    for (col in seq_len(ncol(rawInfo))) {
        tmp <- c(tmp, cor(rawInfo[col], totals - rawInfo[col]))
    }
    taskMatrix["cor", ] <- c(tmp, 1)

    # Add variance row
    sbjs <- as.character(t(unique(data[idCol])))
    taskMatrix <- rbind(taskMatrix, var = apply(taskMatrix[sbjs, ], 2, var))

    return(taskMatrix)
}

matrix2Long <- function(matrix) {
    # Convert humanLE Matrices to data format
    longData <- tibble(SbjID = character(), Trial = numeric(), Corr = logical())
    # Loop through matrix
    for (sbj in rownames(matrix)) {
        for (trial in colnames(matrix)) {
            longData <- longData |>
                add_row(
                    SbjID = sbj,
                    Trial = as.numeric(trial),
                    Corr = matrix[sbj, trial]
                )
        }
    }
    return(longData)
}

matrix2Rel <- function(mat, check.keys = FALSE, covar = FALSE) {
    # Check if matrix has a cor row
    if ("cor" %in% rownames(mat)) {
        tmp <- mat[, !is.na(mat["cor", ])]
        tmp <- tmp[seq_len(nrow(tmp) - 3), seq_len(ncol(tmp) - 1)]
    } else { # Just process the entire table
        # Remove columns with no variance
        tmp <- mat[, !(colMeans(mat) == 1 | colMeans(mat) == 0)]
    }
    return(splitHalf(tmp, check.keys = check.keys, covar = covar))
}

# Load data ----
# Load human data
humanLEMatrix <- read_csv("./data_storage/humanData/le1.csv") |>
    column_to_rownames("SbjID")
humanMatchMatrix <- read_csv("./data_storage/humanData/match1.csv") |>
    column_to_rownames("SbjID")
humanMOOMatrix <- read_csv("./data_storage/humanData/moo1.csv") |>
    column_to_rownames("SbjID")
# Change the column names to just be the trial number (removing the X)
colnames(humanLEMatrix) <- gsub("X", "", colnames(humanLEMatrix))
colnames(humanMatchMatrix) <- gsub("X", "", colnames(humanMatchMatrix))
colnames(humanMOOMatrix) <- gsub("X", "", colnames(humanMOOMatrix))

# Load model data from grid search ----
modelDir <- "./data_storage/results/v3ModelingDataBackup"

# Handle LE data
modelLEFiles <- list.files(modelDir, pattern = "results_learn_exemp", full.names = TRUE)

humanLEDiff <- colMeans(humanLEMatrix)
LEParamData <- tibble(noise = numeric(), learnAdv = numeric(), MSE = numeric(), Cor = numeric(), Rel = numeric(), NoVar = numeric())
# Loop through files
for (file in modelLEFiles) {
    # Find model parameters in the file
    tmpParams <- str_split(file, "/") |>
        unlist() |>
        tail(1) |>
        str_replace("results_learn_exemp_", "") |>
        str_replace(".csv", "") |>
        strsplit("_") |>
        unlist()

    # LE only noise and learnAdv
    noiseParam <- tmpParams[str_detect(tmpParams, "noise")] |>
        str_replace("noise-", "") |>
        as.numeric()
    learnAdvParam <- tmpParams[str_detect(tmpParams, "learnAdv")] |>
        str_replace("learnAdv-", "") |>
        as.numeric()

    if (length(learnAdvParam) == 0) {
        learnAdvParam <- 0
    }

    # Load data
    tmpData <- read_csv(file, show_col_types = FALSE)

    # Convert to matrix
    tmpMatrix <- long2Matrix(tmpData, "ModelName", "Corr", "Trial")

    # Calculate reliability
    rel <- matrix2Rel(tmpMatrix)$lambda2

    # Remove extra columns and rows
    tmpMatrix <- tmpMatrix |>
        select(-Sum) |>
        filter(!row.names(tmpMatrix) %in% c("mean", "cor", "var"))

    # Calculate difference in difficulty between this model and human
    tmpDiff <- mean((colMeans(tmpMatrix) - humanLEDiff)^2)
    tmpCor <- cor(colMeans(tmpMatrix), humanLEDiff)

    # Calculate ceiling/floor trials
    tmpNoVar <- sum(colMeans(tmpMatrix) == 0 | colMeans(tmpMatrix) == 1)

    # Add to tibble
    LEParamData <- LEParamData |>
        add_row(
            noise = noiseParam,
            learnAdv = learnAdvParam,
            MSE = tmpDiff,
            Cor = tmpCor,
            Rel = rel,
            NoVar = tmpNoVar
        )
}

# Handle Match data
modelMatchFiles <- list.files(modelDir, pattern = "results_threeAFC", full.names = TRUE)

humanMatchDiff <- colMeans(humanMatchMatrix)
matchParamData <- tibble(noise = numeric(), learningAdv = numeric(), encNoise = numeric(), MSE = numeric(), Cor = numeric(), Rel = numeric(), NoVar = numeric())
# Loop through files
for (file in modelMatchFiles) {
    # Find model parameters in the file
    tmpParams <- str_split(file, "/") |>
        unlist() |>
        tail(1) |>
        str_replace("results_threeAFC_", "") |>
        str_replace(".csv", "") |>
        strsplit("_") |>
        unlist()

    # Match has noise, encNoise, learningAdv
    noiseParam <- tmpParams[str_detect(tmpParams, "noise")] |>
        str_replace("noise-", "") |>
        as.numeric()
    encNoiseParam <- tmpParams[str_detect(tmpParams, "encNoise")] |>
        str_replace("encNoise-", "") |>
        as.numeric()
    learnAdvParam <- tmpParams[str_detect(tmpParams, "learnAdv")] |>
        str_replace("learnAdv-", "") |>
        as.numeric()

    if (length(learnAdvParam) == 0) {
        learnAdvParam <- 0
    }
    if (length(encNoiseParam) == 0) {
        encNoiseParam <- 0
    }
    # Load data
    tmpData <- read_csv(file, show_col_types = FALSE)

    # Convert to matrix
    tmpMatrix <- long2Matrix(tmpData, "ModelName", "Corr", "Trial")

    # Calculate reliability
    rel <- matrix2Rel(tmpMatrix)$lambda2

    # Remove extra columns and rows
    tmpMatrix <- tmpMatrix |>
        select(-Sum) |>
        filter(!row.names(tmpMatrix) %in% c("mean", "cor", "var"))


    # Calculate difference in difficulty between this model and human
    tmpDiff <- mean((colMeans(tmpMatrix) - humanMatchDiff)^2)
    tmpCor <- cor(colMeans(tmpMatrix), humanMatchDiff)

    # Calculate ceiling/floor trials
    tmpNoVar <- sum(colMeans(tmpMatrix) == 0 | colMeans(tmpMatrix) == 1)

    # Add to tibble
    matchParamData <- matchParamData |>
        add_row(
            noise = noiseParam,
            learningAdv = learnAdvParam,
            encNoise = encNoiseParam,
            MSE = tmpDiff,
            Cor = tmpCor,
            Rel = rel,
            NoVar = tmpNoVar
        )
}
matchParamData

# Handle MOO data
modelMOOFiles <- list.files(modelDir, pattern = "results_many_odd", full.names = TRUE)

humanMOODiff <- colMeans(humanMOOMatrix)
MOOParamData <- tibble(noise = numeric(), encNoise = numeric(), MSE = numeric(), Cor = numeric(), Rel = numeric(), NoVar = numeric())
# Loop through files
for (file in modelMOOFiles) {
    # Find model parameters in the file
    tmpParams <- str_split(file, "/") |>
        unlist() |>
        tail(1) |>
        str_replace("results_many_odd_", "") |>
        str_replace(".csv", "") |>
        strsplit("_") |>
        unlist()

    # MOO has noise and encNoise
    noiseParam <- str_split(tmpParams[1], "-") |>
        unlist() |>
        tail(1) |>
        as.numeric()
    encNoiseParam <- str_split(tmpParams[2], "-") |>
        unlist() |>
        tail(1) |>
        as.numeric()

    if (length(encNoiseParam) == 0) {
        encNoiseParam <- 0
    }

    # Load data
    tmpData <- read_csv(file, show_col_types = FALSE)

    # Convert to matrix
    tmpMatrix <- long2Matrix(tmpData, "ModelName", "Corr", "Trial")

    # Calculate reliability
    rel <- matrix2Rel(tmpMatrix)$lambda2

    # Remove extra columns and rows
    tmpMatrix <- tmpMatrix |>
        select(-Sum) |>
        filter(!row.names(tmpMatrix) %in% c("mean", "cor", "var"))

    # Calculate difference in difficulty between this model and human
    tmpDiff <- mean((colMeans(tmpMatrix) - humanMOODiff)^2)
    tmpCor <- cor(colMeans(tmpMatrix), humanMOODiff)

    # Calculate ceiling/floor trials
    tmpNoVar <- sum(colMeans(tmpMatrix) == 0 | colMeans(tmpMatrix) == 1)

    # Add to tibble
    MOOParamData <- MOOParamData |>
        add_row(
            noise = noiseParam,
            encNoise = encNoiseParam,
            MSE = tmpDiff,
            Cor = tmpCor,
            Rel = rel,
            NoVar = tmpNoVar
        )
}
MOOParamData

modelLEMatrix <- long2Matrix(modelLEData, "ModelName", "Corr", "Trial")

modelMatchData <- read_csv("./data_storage/results/results_threeAFC_noise-1.0_encNoise-2.0.csv")
modelMatchMatrix <- long2Matrix(modelMatchData, "ModelName", "Corr", "Trial")

modelMOOData <- read_csv("./data_storage/results/results_many_odd_noise-1.0_encNoise-1.0.csv")
modelMOOMatrix <- long2Matrix(modelMOOData, "ModelName", "Corr", "Trial")

# Convert human Matrices to data format
humanLEData <- matrix2Long(humanLEMatrix) |>
    left_join(
        modelLEData |>
            filter(ModelName == modelLEData$ModelName[1]) |>
            select(-ModelName, -Response, -Corr),
        by = "Trial"
    )
humanMatchData <- matrix2Long(humanMatchMatrix) |>
    left_join(
        modelMatchData |>
            filter(ModelName == modelMatchData$ModelName[1]) |>
            select(-ModelName, -Response, -Corr),
        by = "Trial"
    )
humanMOOData <- matrix2Long(humanMOOMatrix) |>
    left_join(
        modelMOOData |>
            filter(ModelName == modelMOOData$ModelName[1]) |>
            select(-ModelName, -Response, -Corr),
        by = "Trial"
    )

# Trial attributes analysis ----
# LE analysis, View and Noise as factors
# Stick human and model data together
allLEData <- humanLEData |>
    mutate(Group = "Human") |>
    bind_rows(
        modelLEData |>
            rename(SbjID = ModelName) |>
            mutate(Group = "Model")
    )

LETrialSummary <- allLEData |>
    group_by(View, Noise, Group) |>
    summarize(
        Corr = mean(Corr),
        .groups = "drop"
    ) |>
    arrange(desc(View), Noise)
LETrialSummary
# Match analysis, Duration as factors
# Stick human and model data together
allMatchData <- humanMatchData |>
    mutate(Group = "Human") |>
    bind_rows(
        modelMatchData |>
            rename(SbjID = ModelName) |>
            mutate(Group = "Model")
    )

matchTrialSummary <- allMatchData |>
    group_by(View, Noise, Duration, Group) |>
    summarize(
        Corr = mean(Corr),
        .groups = "drop"
    ) |>
    arrange(desc(View), Noise, Duration)

matchTrialSummary
# MOO analysis, Duration as factors
# Stick human and model data together
allMOOData <- humanMOOData |>
    mutate(Group = "Human") |>
    bind_rows(
        modelMOOData |>
            rename(SbjID = ModelName) |>
            mutate(Group = "Model")
    )

MOOTrialSummary <- allMOOData |>
    group_by(Duration, Group) |>
    summarize(
        Corr = mean(Corr),
        .groups = "drop"
    ) |>
    arrange(Duration)
MOOTrialSummary
