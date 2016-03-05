# diversit-tag


R code for invesigating changes in diversity for a massively multi-allelic system. 


## The Model


## Data Format

Raw data should be in files with the first column a set of tags, labelled tag,
other columns then give the counts in a set of organs.  Example 
data sets are given in the files `/data/example1.txt` and `/data/example2.txt`.

The first few lines of `example1.txt` look like


            tag            |  Pop1_A | Pop1_B | Pop1_C 
---------------------------|---------|--------|--------
AAATCACGATGGAAATTGGTTAAACCC|   765   |   61   |     75   
AAATCATGATGCAAAACGGTTCAACAT|   387   |  154   |    156   
AAATCCCGATGCAAACTGGTACAACTC|   1849  |17374   |     20   
AAATCCCGATGCAAATCGGTTTAACTC|   2338  |  888   |    403   
AAATCGCGATGCAAAGCGGTTTAACCC|   1213  |  365   |    400   


These can be read into R in the usual way

````
ex1 <- read.table("data/example1.txt", header=TRUE)
ex2 <- read.table("data/example2.txt", header=TRUE)
```

### Combining data sets

The R function `merge_by_tag` allows you to merge data sets 
together.  Tags that are not present in either of the data sets 
will be given a zero count.

```
source("R/functions.R")
ex <- merge_by_tag(ex1, ex2)
```



## Analysing these data

### Plotting

A simple barplot of these data is providied by the function `plot_tagfreq`. This takes the `n_colours` most frequent
tags across both data sets, sorts by the frequency in the first column, and plots the frequencies.  Tags that are not in the top `n_colours` are labelled grey.

```
plot_tagfreq(ex[, -1],  n_colours=20)
```

![tag plot](fig/plot1.png)


### Estimating Nei's G

We estimate Nei's G using the function `estimate_G`.  This takes a 
pair of parameters, the first is the count, the second, if present gives 
the frequency of tags in a baseline population.  If this is not present then an equal frequency in the baseline is assumed.


```{r}
estimate_G(ex$Pop1_A)
estimate_G(ex$Pop1_A, ex$Pop2_A, min_p=0.0001)
estimate_G(ex$Pop2_C, ex$Pop1_A, min_p=0.0001)
```





