#' This is simulated data for r4proteomics course
library(tidyverse)
set.seed(42)

simulate_proteomics_data <- function() {
    n_proteins <- 1500
    n_samples <- 24

    sample_metadata <- data.frame(
        sample_id = paste0("S", sprintf("%02d", 1:n_samples)),
        condition = rep(c("Control", "Treatment"), each = 12),
        batch = rep(c("Batch1", "Batch2", "Batch3"), each = 8),
        timepoint = rep(rep(c("Week0", "Week4", "Week8"), each = 4), times = 2),
        patient_id = rep(1:12, times = 2),
        sex = rep(sample(c("M", "F"), 12, replace = TRUE), times = 2),
        age = rep(sample(40:70, 12, replace = TRUE), times = 2),
        stringsAsFactors = FALSE
    )

    protein_annotations <- data.frame(
        protein_id = paste0("P", sprintf("%05d", 1:n_proteins)),
        gene_symbol = paste0("GENE", sprintf("%04d", sample(1:20000, n_proteins))),
        protein_name = paste0("Protein_", sample(1:5000, n_proteins, replace = TRUE)),
        molecular_weight = round(runif(n_proteins, min = 10000, max = 300000), 0),
        peptide_count = sample(1:50, n_proteins, replace = TRUE, prob = 1 / (1:50)),
        confidence_score = round(runif(n_proteins, min = 0, max = 100), 1),
        stringsAsFactors = FALSE
    )

    base_intensities <- rnorm(n_proteins, mean = 1e6, sd = 3e5)

    protein_matrix_raw <- matrix(
        rep(base_intensities, n_samples),
        nrow = n_proteins,
        ncol = n_samples
    )

    # (protein-specific)
    biological_variation <- abs(matrix(
        rnorm(n_proteins * n_samples, mean = 1, sd = 0.25),
        nrow = n_proteins,
        ncol = n_samples
    ))

    protein_matrix_raw <- protein_matrix_raw * biological_variation

    # add DE effect (15% of proteins)
    n_de_proteins <- round(n_proteins * 0.15)
    de_protein_indices <- sample(1:n_proteins, n_de_proteins)

    # varying fold changes
    fc_values <- c(
        runif(round(n_de_proteins * 0.7), min = 1.5, max = 3), # moderate changes
        runif(round(n_de_proteins * 0.2), min = 3, max = 6), # huge changes
        runif(round(n_de_proteins * 0.1), min = 0.1, max = 0.5) # down-regulated
    )

    for (i in 1:n_de_proteins) {
        protein_idx <- de_protein_indices[i]
        treatment_samples <- sample_metadata$condition == "Treatment"

        if (fc_values[i] < 1) {
            # down-regulation
            protein_matrix_raw[protein_idx, treatment_samples] <-
                protein_matrix_raw[protein_idx, treatment_samples] * fc_values[i]
        } else {
            # up-regulation
            protein_matrix_raw[protein_idx, treatment_samples] <-
                protein_matrix_raw[protein_idx, treatment_samples] * fc_values[i]
        }
    }

    # add strong batch effects
    batch_factors <- c(Batch1 = 0.8, Batch2 = 1.0, Batch3 = 1.3)
    for (batch in names(batch_factors)) {
        batch_samples <- sample_metadata$batch == batch
        protein_matrix_raw[, batch_samples] <- protein_matrix_raw[, batch_samples] * batch_factors[[batch]]
    }

    # add timepoint effects (longitudinal component)
    timepoint_factors <- c(Week0 = 1.0, Week4 = 1.1, Week8 = 1.2)
    for (timepoint in names(timepoint_factors)) {
        time_samples <- sample_metadata$timepoint == timepoint
        protein_matrix_raw[, time_samples] <- protein_matrix_raw[, time_samples] * timepoint_factors[[timepoint]]
    }

    # apply log2 transformation (standard in proteomics)
    protein_matrix <- log2(protein_matrix_raw + 1)

    # add missing values (MNAR - Missing Not At Random)
    # Lower abundance proteins are more likely to be missing
    protein_means <- rowMeans(protein_matrix)
    missing_probs <- 1 / (1 + exp(-scale(protein_means))) # Inverse relationship with abundance

    for (i in 1:n_proteins) {
        prob_i <- missing_probs[i] * 0.25
        if (is.na(prob_i) || prob_i <= 0) next

        n_missing <- rbinom(1, size = n_samples, prob = prob_i)

        if (n_missing > 0) {
            missing_samples <- sample(1:n_samples, n_missing)
            protein_matrix[i, missing_samples] <- NA
        }
    }

    # add some completely missing proteins (5%)
    completely_missing <- sample(1:n_proteins, round(n_proteins * 0.05))
    protein_matrix[completely_missing, ] <- NA

    rownames(protein_matrix) <- protein_annotations$protein_id
    colnames(protein_matrix) <- sample_metadata$sample_id

    dataset <- list(
        protein_matrix = protein_matrix,
        sample_metadata = sample_metadata,
        protein_annotations = protein_annotations,
        de_protein_indices = de_protein_indices,
        fc_values = fc_values
    )

    return(dataset)
}

cat("Generating proteomics dataset...\n")
course_dataset <- simulate_proteomics_data()

saveRDS(course_dataset, "data/whole_proteomics_dataset.rds")

# saveRDS(dataset$protein_matrix, "data/protein_matrix.rds")
# saveRDS(dataset$sample_metadata, "data/sample_metadata.rds")
# saveRDS(dataset$protein_annotations, "data/protein_annotations.rds")

write.csv(as.data.frame(course_dataset$protein_matrix),
    "data/protein_matrix.csv",
    row.names = TRUE
)

write.csv(course_dataset$sample_metadata,
    "data/sample_metadata.csv",
    row.names = FALSE
)

write.csv(course_dataset$protein_annotations,
    "data/protein_annotations.csv",
    row.names = FALSE
)

# ground-truth differential expression information
truth_df <- data.frame(
    protein_id = rownames(course_dataset$protein_matrix)[course_dataset$de_protein_indices],
    de_index = course_dataset$de_protein_indices,
    fold_change = course_dataset$fc_values
)

write.csv(truth_df, "data/DE_truth.csv", row.names = FALSE)

cat("=== DATASET SUMMARY ===\n")
cat("Proteins:", nrow(dataset$protein_matrix), "\n")
cat("Samples:", ncol(dataset$protein_matrix), "\n")
cat(
    "Missing values:", sum(is.na(dataset$protein_matrix)),
    "(", round(mean(is.na(dataset$protein_matrix)) * 100, 1), "%)\n"
)
cat("Differential proteins:", length(dataset$de_protein_indices), "\n")
cat("Conditions:", paste(unique(dataset$sample_metadata$condition), collapse = ", "), "\n")
cat("Batches:", paste(unique(dataset$sample_metadata$batch), collapse = ", "), "\n")
cat("Timepoints:", paste(unique(dataset$sample_metadata$timepoint), collapse = ", "), "\n")

cat("\nDataset saved successfully!\n")


# --- Introduce a complete outlier sample ---
# Pick a sample to be an outlier (e.g., randomly or fixed)
set.seed(123)
outlier_sample <- sample(1:ncol(course_dataset$protein_matrix), 1)
cat("Outlier sample:", colnames(course_dataset$protein_matrix)[outlier_sample], "\n")

# Make all protein intensities extremely high for this sample
course_dataset$protein_matrix[, outlier_sample] <- course_dataset$protein_matrix[, outlier_sample] + 3.2
# (Adding 10 in log2 scale creates a strong shift visible in PCA)
write.csv(as.data.frame(course_dataset$protein_matrix),
    "data/protein_matrix_withOutlier.csv",
    row.names = TRUE
)
