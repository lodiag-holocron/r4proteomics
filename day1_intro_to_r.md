# Introduction to R and RStudio {#day1}



## Learning Objectives

By the end of Day 1, you will be able to:

- Install and navigate RStudio effectively
- Understand basic R data structures (vectors, data frames, lists)
- Import and explore simple datasets
- Write basic control flow structures and functions

## Module 1: Setting Up and Getting Started with R {#day1-mod1}

### Introduction

R is a powerful programming language and environment specifically designed for statistical computing and graphics. RStudio is an integrated development environment (IDE) that makes working with R much easier.

### Installing R and RStudio

1. **Install R** (version 4.3.0 or higher)
   - Visit [https://cran.r-project.org/](https://cran.r-project.org/)
   - Download the version appropriate for your operating system
   - Run the installer

2. **Install RStudio Desktop**
   - Visit [https://posit.co/download/rstudio-desktop/](https://posit.co/download/rstudio-desktop/)
   - Download the free Desktop version
   - Run the installer

### RStudio Interface Tour

RStudio has four main panes:

1. **Source Editor** (top-left): Where you write and edit your scripts
2. **Console** (bottom-left): Where code is executed and results appear
3. **Environment/History** (top-right): Shows objects in memory and command history
4. **Files/Plots/Packages/Help** (bottom-right): File browser, plot viewer, package manager, and help documentation

### Scripts vs Console

The **Console** is for:
- Quick calculations
- Testing commands
- Interactive exploration

**Scripts** (.R or .Rmd files) are for:
- Saving your work
- Creating reproducible analyses
- Organizing complex workflows

### Basic Operators


``` r
# Arithmetic operators
5 + 3        # Addition
#> [1] 8
10 - 4       # Subtraction
#> [1] 6
6 * 7        # Multiplication
#> [1] 42
20 / 4       # Division
#> [1] 5
2 ^ 3        # Exponentiation
#> [1] 8
17 %% 5      # Modulo (remainder)
#> [1] 2

# Assignment operator
x <- 10      # Assign 10 to x
y = 5        # Alternative (but <- is preferred)

# Comparison operators
5 == 5       # Equal to
#> [1] TRUE
5 != 3       # Not equal to
#> [1] TRUE
7 > 3        # Greater than
#> [1] TRUE
4 < 8        # Less than
#> [1] TRUE
5 >= 5       # Greater than or equal
#> [1] TRUE
3 <= 10      # Less than or equal
#> [1] TRUE

# Logical operators
TRUE & FALSE  # AND
#> [1] FALSE
TRUE | FALSE  # OR
#> [1] TRUE
!TRUE         # NOT
#> [1] FALSE
```

### Creating Your First Script


``` r
# Create a new R script: File > New File > R Script
# Or use Ctrl+Shift+N (Windows/Linux) or Cmd+Shift+N (Mac)

# Write your code
message("Hello, Proteomics World!")

# Save your script: File > Save
# Run code: Ctrl+Enter (Windows/Linux) or Cmd+Return (Mac)
```

### Exercise 1.1: First Steps

Create a new R script and:

1. Calculate the sum of 123 and 456
2. Assign the result to a variable called `total`
3. Print the value of `total`
4. Calculate what percentage 123 is of the total


```{.r .fold-hide}
# Solution
result1 <- 123 + 456
total <- result1
print(total)
#> [1] 579

percentage <- (123 / total) * 100
print(paste0("123 is ", round(percentage, 2), "% of the total"))
#> [1] "123 is 21.24% of the total"
```

## Module 2: Data Types and Structures {#day1-mod2}

### Vectors

Vectors are the most basic data structure in R. They contain elements of the same type.


``` r
# Numeric vectors
ages <- c(25, 30, 35, 40, 45)
print(ages)
#> [1] 25 30 35 40 45

# Character vectors
names <- c("Alice", "Bob", "Charlie", "Diana", "Eve")
print(names)
#> [1] "Alice"   "Bob"     "Charlie" "Diana"   "Eve"

# Logical vectors
passed_qc <- c(TRUE, TRUE, FALSE, TRUE, TRUE)
print(passed_qc)
#> [1]  TRUE  TRUE FALSE  TRUE  TRUE

# Sequences
seq_1_10 <- 1:10
seq_custom <- seq(from = 0, to = 100, by = 10)
```

### Indexing and Subsetting


``` r
# Access elements by position (1-indexed!)
ages[1]           # First element
#> [1] 25
ages[c(1, 3, 5)]  # Multiple elements
#> [1] 25 35 45
ages[-2]          # All except second element
#> [1] 25 35 40 45

# Logical indexing
ages[ages > 35]   # Elements greater than 35
#> [1] 40 45

# Named vectors
protein_abundance <- c(ACTB = 1500, GAPDH = 2000, MYC = 800)
protein_abundance["ACTB"]
#> ACTB 
#> 1500
```

### Data Frames

Data frames are the most common structure for storing tabular data.


``` r
# Create a data frame
patient_data <- data.frame(
  patient_id = 1:5,
  name = c("Alice", "Bob", "Charlie", "Diana", "Eve"),
  age = c(25, 30, 35, 40, 45),
  treatment = c("A", "B", "A", "B", "A"),
  response = c(TRUE, TRUE, FALSE, TRUE, FALSE),
  stringsAsFactors = FALSE
)

print(patient_data)
#>   patient_id    name age treatment response
#> 1          1   Alice  25         A     TRUE
#> 2          2     Bob  30         B     TRUE
#> 3          3 Charlie  35         A    FALSE
#> 4          4   Diana  40         B     TRUE
#> 5          5     Eve  45         A    FALSE

# View structure
str(patient_data)
#> 'data.frame':	5 obs. of  5 variables:
#>  $ patient_id: int  1 2 3 4 5
#>  $ name      : chr  "Alice" "Bob" "Charlie" "Diana" ...
#>  $ age       : num  25 30 35 40 45
#>  $ treatment : chr  "A" "B" "A" "B" ...
#>  $ response  : logi  TRUE TRUE FALSE TRUE FALSE

# Summary statistics
summary(patient_data)
#>    patient_id     name                age    
#>  Min.   :1    Length:5           Min.   :25  
#>  1st Qu.:2    Class :character   1st Qu.:30  
#>  Median :3    Mode  :character   Median :35  
#>  Mean   :3                       Mean   :35  
#>  3rd Qu.:4                       3rd Qu.:40  
#>  Max.   :5                       Max.   :45  
#>   treatment          response      
#>  Length:5           Mode :logical  
#>  Class :character   FALSE:2        
#>  Mode  :character   TRUE :3        
#>                                    
#>                                    
#> 

# Access columns
patient_data$age
#> [1] 25 30 35 40 45
patient_data[, "name"]
#> [1] "Alice"   "Bob"     "Charlie" "Diana"   "Eve"
patient_data[, 2]
#> [1] "Alice"   "Bob"     "Charlie" "Diana"   "Eve"

# Access rows
patient_data[1, ]           # First row
#>   patient_id  name age treatment response
#> 1          1 Alice  25         A     TRUE
patient_data[1:3, ]         # First three rows
#>   patient_id    name age treatment response
#> 1          1   Alice  25         A     TRUE
#> 2          2     Bob  30         B     TRUE
#> 3          3 Charlie  35         A    FALSE

# Access specific cells
patient_data[2, 3]          # Row 2, Column 3
#> [1] 30
patient_data[2, "age"]      # Same, using column name
#> [1] 30

# Subset by condition
patient_data[patient_data$age > 30, ]
#>   patient_id    name age treatment response
#> 3          3 Charlie  35         A    FALSE
#> 4          4   Diana  40         B     TRUE
#> 5          5     Eve  45         A    FALSE
patient_data[patient_data$treatment == "A", ]
#>   patient_id    name age treatment response
#> 1          1   Alice  25         A     TRUE
#> 3          3 Charlie  35         A    FALSE
#> 5          5     Eve  45         A    FALSE
```

### Lists

Lists can contain elements of different types and structures.


``` r
# Create a list
experiment <- list(
  experiment_id = "EXP001",
  date = "2025-01-15",
  samples = c("S1", "S2", "S3"),
  data = patient_data,
  validated = TRUE
)

# Access list elements
experiment$experiment_id
#> [1] "EXP001"
experiment[[1]]
#> [1] "EXP001"
experiment[["samples"]]
#> [1] "S1" "S2" "S3"
```

### Factors

Factors are used for categorical data.


``` r
# Create factor
treatment_factor <- factor(c("Control", "Drug A", "Drug B", "Control", "Drug A"))
print(treatment_factor)
#> [1] Control Drug A  Drug B  Control Drug A 
#> Levels: Control Drug A Drug B

# Check levels
levels(treatment_factor)
#> [1] "Control" "Drug A"  "Drug B"

# Ordered factors
severity <- factor(
  c("Mild", "Severe", "Moderate", "Mild", "Severe"),
  levels = c("Mild", "Moderate", "Severe"),
  ordered = TRUE
)
print(severity)
#> [1] Mild     Severe   Moderate Mild     Severe  
#> Levels: Mild < Moderate < Severe
```

### Type Coercion


``` r
# Implicit coercion
mixed <- c(1, 2, "three", 4)  # All converted to character
print(mixed)
#> [1] "1"     "2"     "three" "4"

# Explicit coercion
numbers_char <- c("1", "2", "3", "4")
numbers_num <- as.numeric(numbers_char)
print(numbers_num)
#> [1] 1 2 3 4

# Check types
class(mixed)
#> [1] "character"
is.numeric(mixed)
#> [1] FALSE
is.character(mixed)
#> [1] TRUE
```

### Exercise 1.2: Data Structures

Create a data frame for a proteomic experiment with:

- 10 protein IDs (P001 to P010)
- Random abundance values between 100 and 5000
- Random p-values between 0 and 1
- Significance status (TRUE if p-value < 0.05)


```{.r .fold-hide}
# Solution
set.seed(42)  # For reproducibility

proteins <- data.frame(
  protein_id = paste0("P", sprintf("%03d", 1:10)),
  abundance = round(runif(10, min = 100, max = 5000), 2),
  p_value = runif(10, min = 0, max = 1),
  stringsAsFactors = FALSE
)

proteins$significant <- proteins$p_value < 0.05

print(proteins)
#>    protein_id abundance   p_value significant
#> 1        P001   4582.55 0.4577418       FALSE
#> 2        P002   4691.67 0.7191123       FALSE
#> 3        P003   1502.08 0.9346722       FALSE
#> 4        P004   4169.19 0.2554288       FALSE
#> 5        P005   3244.55 0.4622928       FALSE
#> 6        P006   2643.57 0.9400145       FALSE
#> 7        P007   3709.28 0.9782264       FALSE
#> 8        P008    759.87 0.1174874       FALSE
#> 9        P009   3319.26 0.4749971       FALSE
#> 10       P010   3554.82 0.5603327       FALSE

# Summary
cat("\nNumber of significant proteins:", sum(proteins$significant), "\n")
#> 
#> Number of significant proteins: 0
```

## Module 3: Control Flow and Functions {#day1-mod3}

### Conditional Statements


``` r
# if statement
x <- 10

if (x > 5) {
  print("x is greater than 5")
}
#> [1] "x is greater than 5"

# if-else
if (x > 15) {
  print("x is greater than 15")
} else {
  print("x is 15 or less")
}
#> [1] "x is 15 or less"

# if-else if-else
score <- 75

if (score >= 90) {
  grade <- "A"
} else if (score >= 80) {
  grade <- "B"
} else if (score >= 70) {
  grade <- "C"
} else {
  grade <- "F"
}

print(paste("Grade:", grade))
#> [1] "Grade: C"

# Vectorized ifelse
values <- c(1, 5, 10, 15, 20)
categories <- ifelse(values > 10, "High", "Low")
print(categories)
#> [1] "Low"  "Low"  "Low"  "High" "High"
```

### Loops


``` r
# for loop
for (i in 1:5) {
  print(paste("Iteration:", i))
}
#> [1] "Iteration: 1"
#> [1] "Iteration: 2"
#> [1] "Iteration: 3"
#> [1] "Iteration: 4"
#> [1] "Iteration: 5"

# Loop through vector
proteins <- c("ACTB", "GAPDH", "MYC")
for (protein in proteins) {
  print(paste("Processing:", protein))
}
#> [1] "Processing: ACTB"
#> [1] "Processing: GAPDH"
#> [1] "Processing: MYC"

# while loop
counter <- 1
while (counter <= 5) {
  print(paste("Counter:", counter))
  counter <- counter + 1
}
#> [1] "Counter: 1"
#> [1] "Counter: 2"
#> [1] "Counter: 3"
#> [1] "Counter: 4"
#> [1] "Counter: 5"

# Loop with condition
numbers <- 1:10
for (num in numbers) {
  if (num %% 2 == 0) {
    print(paste(num, "is even"))
  }
}
#> [1] "2 is even"
#> [1] "4 is even"
#> [1] "6 is even"
#> [1] "8 is even"
#> [1] "10 is even"
```

### Functions


``` r
# Basic function
greet <- function(name) {
  message <- paste("Hello,", name, "!")
  return(message)
}

greet("Alice")
#> [1] "Hello, Alice !"

# Function with multiple parameters
calculate_fold_change <- function(treatment, control) {
  fc <- treatment / control
  log2_fc <- log2(fc)
  return(log2_fc)
}

calculate_fold_change(treatment = 200, control = 100)
#> [1] 1

# Function with default parameters
normalize_abundance <- function(abundance, method = "median") {
  if (method == "median") {
    normalized <- abundance / median(abundance, na.rm = TRUE)
  } else if (method == "mean") {
    normalized <- abundance / mean(abundance, na.rm = TRUE)
  } else {
    stop("Method must be 'median' or 'mean'")
  }
  return(normalized)
}

values <- c(100, 200, 300, 400, 500)
normalize_abundance(values)
#> [1] 0.3333333 0.6666667 1.0000000 1.3333333 1.6666667
normalize_abundance(values, method = "mean")
#> [1] 0.3333333 0.6666667 1.0000000 1.3333333 1.6666667
```

### Apply Family Functions


``` r
# Create sample data
protein_matrix <- matrix(
  c(100, 150, 200, 250, 
    110, 160, 210, 260,
    120, 170, 220, 270),
  nrow = 3, byrow = TRUE
)
colnames(protein_matrix) <- c("Sample1", "Sample2", "Sample3", "Sample4")
rownames(protein_matrix) <- c("Protein1", "Protein2", "Protein3")

print(protein_matrix)
#>          Sample1 Sample2 Sample3 Sample4
#> Protein1     100     150     200     250
#> Protein2     110     160     210     260
#> Protein3     120     170     220     270

# apply: apply function to rows or columns
row_means <- apply(protein_matrix, 1, mean)  # 1 = rows
col_means <- apply(protein_matrix, 2, mean)  # 2 = columns

print(row_means)
#> Protein1 Protein2 Protein3 
#>      175      185      195
print(col_means)
#> Sample1 Sample2 Sample3 Sample4 
#>     110     160     210     260

# lapply: apply function to list, returns list
my_list <- list(a = 1:5, b = 6:10, c = 11:15)
list_means <- lapply(my_list, mean)
print(list_means)
#> $a
#> [1] 3
#> 
#> $b
#> [1] 8
#> 
#> $c
#> [1] 13

# sapply: simplified version of lapply
vector_means <- sapply(my_list, mean)
print(vector_means)
#>  a  b  c 
#>  3  8 13
```

### Exercise 1.3: Functions and Loops

Write a function that:

1. Takes a vector of protein abundances
2. Calculates the coefficient of variation (CV = sd/mean * 100)
3. Returns "Pass" if CV < 20%, "Fail" otherwise

Apply this function to multiple samples using a loop.


```{.r .fold-hide}
# Solution
calculate_cv_status <- function(abundances) {
  cv <- (sd(abundances, na.rm = TRUE) / mean(abundances, na.rm = TRUE)) * 100
  
  if (cv < 20) {
    status <- "Pass"
  } else {
    status <- "Fail"
  }
  
  return(list(cv = round(cv, 2), status = status))
}

# Create sample data
sample_data <- list(
  sample1 = c(100, 105, 98, 102, 99),
  sample2 = c(100, 150, 90, 200, 80),
  sample3 = c(500, 505, 498, 502, 496)
)

# Apply function
for (sample_name in names(sample_data)) {
  result <- calculate_cv_status(sample_data[[sample_name]])
  cat(sample_name, "- CV:", result$cv, "% - Status:", result$status, "\n")
}
#> sample1 - CV: 2.75 % - Status: Pass 
#> sample2 - CV: 40.56 % - Status: Fail 
#> sample3 - CV: 0.7 % - Status: Pass
```

## Importing and Exploring Data

### Reading CSV Files


``` r
# Read CSV
data <- read.csv("data/proteins.csv")

# Read with tidyverse
library(readr)
data <- read_csv("data/proteins.csv")

# Read tab-delimited
data <- read.delim("data/proteins.txt", sep = "\t")
```

### Basic Data Exploration


``` r
# Create example data
set.seed(123)
protein_data <- data.frame(
  protein_id = paste0("P", 1:100),
  abundance = rnorm(100, mean = 1000, sd = 200),
  condition = rep(c("Control", "Treatment"), each = 50)
)

# Dimensions
dim(protein_data)
#> [1] 100   3
nrow(protein_data)
#> [1] 100
ncol(protein_data)
#> [1] 3

# First and last rows
head(protein_data)
#>   protein_id abundance condition
#> 1         P1  887.9049   Control
#> 2         P2  953.9645   Control
#> 3         P3 1311.7417   Control
#> 4         P4 1014.1017   Control
#> 5         P5 1025.8575   Control
#> 6         P6 1343.0130   Control
tail(protein_data)
#>     protein_id abundance condition
#> 95         P95 1272.1305 Treatment
#> 96         P96  879.9481 Treatment
#> 97         P97 1437.4666 Treatment
#> 98         P98 1306.5221 Treatment
#> 99         P99  952.8599 Treatment
#> 100       P100  794.7158 Treatment

# Summary statistics
summary(protein_data)
#>   protein_id          abundance       condition        
#>  Length:100         Min.   : 538.2   Length:100        
#>  Class :character   1st Qu.: 901.2   Class :character  
#>  Mode  :character   Median :1012.4   Mode  :character  
#>                     Mean   :1018.1                     
#>                     3rd Qu.:1138.4                     
#>                     Max.   :1437.5

# Table for categorical data
table(protein_data$condition)
#> 
#>   Control Treatment 
#>        50        50
```

## Day 1 Summary

Today you learned:

- ✓ How to set up R and RStudio
- ✓ Basic R operators and syntax
- ✓ Data structures: vectors, data frames, lists, factors
- ✓ Indexing and subsetting data
- ✓ Control flow: if/else, loops
- ✓ Writing custom functions
- ✓ Importing and exploring data

### Homework

1. Install all required packages for Day 2
2. Practice writing functions for data manipulation
3. Explore the built-in datasets in R (use `data()` to see available datasets)


``` r
# Install packages for Day 2
install.packages(c("ggplot2", "dplyr", "tidyr", "pheatmap"))

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("limma", "vsn"))
```

## Additional Resources

- [R for Data Science](https://r4ds.had.co.nz/) by Hadley Wickham
- [RStudio Cheat Sheets](https://posit.co/resources/cheatsheets/)
- [Stack Overflow](https://stackoverflow.com/questions/tagged/r) for questions
