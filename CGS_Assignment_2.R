install.packages("BiocManager")
install.packages("DESeq2")
install.packages("ggplot2")
library(DESeq2)
library(ggplot2)

metaData <- readr::read_tsv("C:\\metadata_SRP055573.tsv")
data <- readr::read_tsv("C:\\SRP055573.tsv")

log<-apply(data[,2:69],1,log) 

max<-apply(log[,2:69],2,max)

min<-apply(log[,2:69],2,min)

rng <- maxVals - minVals
rng <- na.omit(rng)
plot(rng)
dens <- density(rng)
plot(dens)

```{r}
# Create the data folder if it doesn't exist
if (!dir.exists("data")) {
  dir.create("data")
}

# Define the file path to the plots directory
plots_dir <- "plots"

# Create the plots folder if it doesn't exist
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Define the file path to the results directory
results_dir <- "results"

# Create the results folder if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}
```

```{r}
# Define the file path to the data directory
# Replace with the path of the folder the files will be in
data_dir <- file.path("data", "SRP055573")

# Declare the file path to the gene expression matrix file
# inside directory saved as `data_dir`
# Replace with the path to your dataset file
data_file <- file.path(data_dir, "SRP055573.tsv")

# Declare the file path to the metadata file
# inside the directory saved as `data_dir`
# Replace with the path to your metadata file
metadata_file <- file.path(data_dir, "metadata_SRP055573.tsv")
```

```{r}
# Check if the gene expression matrix file is at the path stored in `data_file`
file.exists(data_file)

# Check if the metadata file is at the file path stored in `metadata_file`
file.exists(metadata_file)
```

```{r}
if (!("DESeq2" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("DESeq2", update = FALSE)
}
if (!("EnhancedVolcano" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("EnhancedVolcano", update = FALSE)
}
if (!("apeglm" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("apeglm", update = FALSE)
}

# Attach the DESeq2 library
library(DESeq2)

# Attach the ggplot2 library for plotting
library(ggplot2)

# We will need this so we can use the pipe: %>%
library(magrittr)
```

```{r}
set.seed(12345)
```

```{r}
# Read in metadata TSV file
metadata <- readr::read_tsv(metadata_file)

# Read in data TSV file
expression_df <- readr::read_tsv(data_file) %>%
  tibble::column_to_rownames("Gene")
```

```{r}
# Make the data in the order of the metadata
expression_df <- expression_df %>%
  dplyr::select(metadata$refinebio_accession_code)

# Check if this is in the same order
all.equal(colnames(expression_df), metadata$refinebio_accession_code)
```

```{r}
head(metadata$refinebio_title)
```

```{r}
metadata <- metadata %>%
  dplyr::mutate(rot_status = dplyr::case_when(
    stringr::str_detect(refinebio_title, "rotenone") ~ "rotenone",
    stringr::str_detect(refinebio_title, "control") ~ "control"
  ))
```

```{r}
# Let's take a look at the original metadata column's info
# and our new `mutation_status` column
dplyr::select(metadata, refinebio_title, rot_status)
```

```{r}
# Print out a preview of `mutation_status`
str(metadata$rot_status)
```

```{r}
# Make mutation_status a factor and set the levels appropriately
metadata <- metadata %>%
  dplyr::mutate(
    # Here we define the values our factor variable can have and their order.
    rot_status = factor(rot_status, levels = c("control", "rotenone"))
  )
```

```{r}
levels(metadata$rot_status)


# Define a minimum counts cutoff and filter the data to include
# only rows (genes) that have total counts above the cutoff
filtered_expression_df <- expression_df %>%
  dplyr::filter(rowSums(.) >= 10)
```

```{r}
# round all expression counts
gene_matrix <- round(filtered_expression_df)
```

```{r}
ddset <- DESeqDataSetFromMatrix(
  # Here we supply non-normalized count data
  countData = gene_matrix,
  # Supply the `colData` with our metadata data frame
  colData = metadata,
  # Supply our experimental variable to `design`
  design = ~rot_status
)
```

```{r}
deseq_object <- DESeq(ddset)
```

```{r}
deseq_results <- results(deseq_object)
```


```{r}
deseq_results <- lfcShrink(
  deseq_object, # The original DESeq2 object after running DESeq()
  coef = 2, # The log fold change coefficient used in DESeq(); the default is 2.
  res = deseq_results # The original DESeq2 results table
)
```


```{r}
head(deseq_results)
```

```{r}
# this is of class DESeqResults -- we want a data frame
deseq_df <- deseq_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # add a column for significance threshold results
  dplyr::mutate(threshold = padj < 0.05) %>%
  # sort by statistic -- the highest values will be genes with
  # higher expression in RPL10 mutated samples
  dplyr::arrange(dplyr::desc(log2FoldChange))

head(deseq_df)


plotCounts(ddset, gene = "caly", intgroup = "rot_status")


readr::write_tsv(
  deseq_df,
  file.path(
    results_dir,
    "SRP055573_diff_expr_results.tsv" # Replace with a relevant output file name
  )
)

# We'll assign this as `volcano_plot`
volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  deseq_df,
  lab = deseq_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)

# Print out plot here
volcano_plot

ggsave(
  plot = volcano_plot,
  file.path(plots_dir, "SRP055573_volcano_plot.png")
) # Replace with a plot name relevant to your data

install.packages("devtools")
library(devtools)
install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)

deseq_df

filterd_data <- deseq_df[!duplicated(deseq_df$Gene), ]

filterd_data

rownames(filterd_data) <- filterd_data[,1]

filterd_data <- filterd_data[,-1]
mat <- counts(ddset)[rownames(filterd_data), ]
mat <- t(apply(mat, 1, scale))
mat
print(metadata)
coldata <- metadata %>% tibble::column_to_rownames("refinebio_accession_code")
colnames(mat) <- rownames(coldata)

map <- Heatmap(mat, cluster_rows = T, cluster_columns = F,
               column_labels = colnames(mat), name = "Heat Bar")

map

