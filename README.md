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


### Shared Tags

We can tabulate the sharing of tags by considering when tags 
have been shared between two populations using the basic R `table` command. 
```
table(A=ex$Pop1_A>0, C=ex$Pop1_C>0)
```

This does not distinguish between rare an common alleles.  We can improve this by working with the tag relative frequencies.  Get the frequencies 
`freq` by sweeping out the sum of tags.  * Note that `ex[, -1]` removes 
returns a `data.frame` without the first column of tag names.

```
freq <- sweep(ex[,-1], 2, colSums(ex[, -1]), "/")
rownames(freq) <- ex[,1]

table(A=cut(freq$Pop1_A, c(0,0.0001,0.005,1)),
      C=cut(freq$Pop1_C, c(0,0.0001,0.005,1)))
```

Simple R functions have been added to simplify the analysis of tag sharing.  



### Distances Between Tags

If tags are coded as DNA then the distances between tags can be simply measured by using the `stringdist` library.

```{r}
stringdist(ex$tag[1], ex$tag[2])
```

The function `tagdistplot` allows you to get plot the frequencies of 
all the tags vs. the distance from the most frequent tag (in addition those tags that are 1 step away from the 2nd, 3rd, ... , `howmany`) tags are colour coded.  This allows you to see how 1-step errors are a major 
reason for low frequency variants in some runs.

```
tagdistplot(ex$tag, freq$Pop1_A, howmany=4)
```


### Estimating Nei's G


$$ 
a = \Sum_{i=1} c_i
$$

We estimate Nei's G using the function `estimate_G`.  This takes a 
pair of parameters, the first is the count, the second, if present gives 
the frequency of tags in a baseline population.  If this is not present then an equal frequency in the baseline is assumed.


```{r}
estimate_G(ex$Pop1_A)
estimate_G(ex$Pop1_A, ex$Pop2_A, min_p=0.0001)
estimate_G(ex$Pop2_C, ex$Pop1_A, min_p=0.0001)
```
![plot of closest tags](fig/closest.png)

