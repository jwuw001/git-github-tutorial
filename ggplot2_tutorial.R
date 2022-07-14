# View the data we're using in this tutorial in a sheet

View(iris)


## ====== built-in function ===================================================

# Use the built-in plot function to generate a scatter plot

plot(x=iris$Sepal.Length, y=iris$Sepal.Width,
     xlab="Sepal Length", ylab="Sepal Width", main="Sepal Length-Width")


# Use the built-in box plot function to generate a figure

boxplot(Sepal.Length ~ Species, data=iris, xlab="Species", ylab="Sepal Length")


## ====== Basic ggplot step by step ===========================================

# Check whether the package ggplot2 exist and install it if not.

if (!require(ggplot2)) {
  install.packages("ggplot2")
}


# Load ggplot2 library

library(ggplot2)


# Generate plotting space and assign variables to be plotted

fig <- ggplot(iris)
fig


# Assign variables to be plotted

fig <- fig + aes(x=Species, y=Sepal.Length)
fig


# Assign which kind of plot is used to represent the data
# Add point as geometric object to generate a scatter plot

fig + geom_point() + aes(color=Species)


# Generate box plot instead

fig + geom_boxplot() + aes(fill=Species)


# Can even overlap two kind of geometric object together

fig + geom_boxplot() + geom_point() + aes(color=Species)


# Connect all the codes with “+”

ggplot(iris) +
  aes(x=Species, y=Sepal.Length, color=Species) +
  geom_boxplot() +
  geom_point()


## ====== Stat function =======================================================

# Let's try to scatter plot with sepal length in x and sepal width in y

fig <- ggplot(iris)
fig +
  aes(x=Sepal.Length, y=Sepal.Width, color=Species) +
  geom_point()


# However, there are actually lot of overlapped points in the figure.
# Let's try to use the stat function to better describe it.

fig +
  aes(x=Sepal.Length, y=Sepal.Width, color=Species) +
  geom_point() +
  stat_sum(geom = "point")

# You can swap the order of the layers

fig +
  aes(x=Sepal.Length, y=Sepal.Width, color=Species) +
  stat_sum(geom = "point") +
  geom_point()


## ====== Order of the layers =================================================

# First box plot, then scatter plot.

ggplot(iris) +
  aes(x=Species, y=Sepal.Length, color=Species) +
  geom_boxplot() +
  geom_point()


# First scatter plot, then box plot.

ggplot(iris) +
  aes(x=Species, y=Sepal.Length, color=Species) +
  geom_point() +
  geom_boxplot()

## ====== Tips for Histogram ==================================================

# Plot a histogram

fig <- ggplot(iris) +
  aes(x=Sepal.Width, fill=Species) +
  geom_histogram(binwidth=0.2, color="gray") +
  labs(title = "Histogram of Sepal Width", x = "Sepal Width", y = "Frequency")
fig


# If we want to show more x ticks, we can set up this parameter:

fig + scale_x_continuous(breaks = seq(2, 5, 0.2))


# Note that each bin have been centered at the ticks, e.g. the first
# bin have been centered at 2.0, which is not what we expected from a
# histogram. We can specify the argument "boundary = 0" to fix it.

fig <- ggplot(iris) +
  aes(x=Sepal.Width, fill=Species) +
  geom_histogram(binwidth=0.2, boundary = 0, color="gray") +
  scale_x_continuous(breaks = seq(2, 5, 0.2)) +
  labs(title = "Histogram of Sepal Width", x = "Sepal Width", y = "Frequency")
fig


# You can also add a vertical line on the figure.
# e.g. highlight a threshold for the sepal width larger than 3.0.

fig + geom_vline(xintercept=3, color="red", linetype="dotdash")


## ====== Tips for bar plot ===================================================

# Plot a bar plot

fig <- ggplot(iris)
fig +
  aes(x=Species, y=Sepal.Length) +
  geom_bar(stat="identity", position="dodge")


# Calculate the mean and sd by each species. Here we can use the function
# `aggregate` to get the mean and sd for each species.
# For more details about the function `aggregate`, please refer to ?aggregate
# or https://www.datasciencemadesimple.com/aggregate-function-in-r/

iris_mean <- aggregate(Sepal.Length ~ Species, data=iris, FUN=mean)
iris_sd <- aggregate(Sepal.Length ~ Species, data=iris, FUN=sd)
iris_sd <- setNames(iris_sd, c("Species", "sd"))
iris_df <- merge(iris_mean, iris_sd, by="Species")
iris_df


# Now we can use the iris_df instead of the original iris dataset
# to generate the plot. In addition, we can use the sd value we got
# to set up the upper and lower values for the error bar layer.

fig <- ggplot(iris_df)
fig +
  aes(x=Species, y=Sepal.Length) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(ymin=Sepal.Length - sd, ymax=Sepal.Length + sd),
                width=0.2, position=position_dodge(0.9))


## ====== Visualize genomic data by ggbio =====================================

# Check whether the required packages exist and install it if not.

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require(ggbio)) {
  BiocManager::install("ggbio")
}

if (!require(GenomicFeatures)) {
  BiocManager::install("GenomicFeatures")
}


# Load ggbio and GenomicFeatures libraries

library(ggbio)
library(GenomicFeatures)


# Download S.cerevisiae annotation from ensenbl and generate a TxDb object

txdb <- makeTxDbFromBiomart(biomart = "ensembl", dataset = "scerevisiae_gene_ensembl")


# Generate a chromosome level GRanges object, which will later used as one of the layers.

df <- data.frame(chr=seqlevels(txdb), start=1, end=seqlengths(txdb),strand='*')
chroms <- makeGRangesFromDataFrame(df, seqinfo=seqinfo(txdb))


# Get the gene subset from the TxDb object and split the genes by strands.

genes <- genes(txdb)
genes_plus <- genes[strand(genes) == "+"]
genes_minus <- genes[strand(genes) == "-"]


# Plot the circular plot layer by layer. Note that the layers have been placed from inside
# to outside. Here are the definition of each layers. (1) Genes on minus strand; (2) Genes
# on plus strand; (3) chromosome ideogram; (4) chromosome labels

ggbio() +
  circle(genes_minus, geom = "rect", color = "steelblue") +
  circle(genes_plus, geom = "rect", color = "red") +
  circle(chroms, geom = "ideo", fill = "gray70") +
  circle(chroms, geom = "text", aes(label = seqnames), vjust = 0, size = 3) +
  theme(legend.key = element_rect(fill = "white", colour = "black"))


# Important!! If you want to try the bar plot again, please unload the ggbio first.
# It seems like a bug on ggbio will cause a failure when calling `geom_bar` function.

detach("package:ggbio", unload = TRUE)
