## Anomaly detection for local ancestry painting

This code was used in two papers:

1. [Nelson et al 2017 "Genomewide analysis of admixture and adaptation in the Africanized honeybee"](https://onlinelibrary.wiley.com/doi/10.1111/mec.14122)
2. [Howard-McCombe et al 2023 "Genetic swamping of the critically endangered Scottish wildcat was recent and accelerated by disease"](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4528485)

If you use the MOSAIC code then the appropriate reference is the second. The first is for HAPMIX or ChromoPainter.

The idea is to take any measure of genome-wide painting (such as ancestry estimation by hapmix, or chromosome painting) and "add up" all the individuals probability of being from a particular population. Then perform statistical testing to identify regions of the genome that are "surprisingly high" (or low) for each ancestry. This is based on a Poisson-Binomial model, as the paintings are themselves probabilities.

## How to use these files

The code comes in two parts: painting and estimation. Painting is done via vanilla "fs", or on the direct output of hapmix. The scripts run through this.

### Top level files:
* dopainting_cp.R: uses the provided "readPaintings" functions to read raw fs painting data.
* dopainting_hapmix.R: uses read.table and manually constructs the data in the desired format.
* dopainting_mosaic.R: uses the MOSAIC .RData and converts the data to the desired format.

### Input data generation
To use with **ChromoPainter** you need to run finestructure:

```{sh}
cd subsetdata
./dofs.sh ## Take a look at this file for how to do the painting
cd ..
```

This creates copying probabilities in files called 
`subsetdata/test/stage7/test_stage7_tmp_mainrun.linked_file<X>_ind<Y>.copyprobsperlocus.out.gz`
for each chromosome <X> and individual <Y>.

To use with **HAPMIX** we have provided tha output of HAPMIX

### File format:
we use a list containing:
- for each chromosome:
  - for each target population
	A dataframe containing the painting probability for snps (rows) against haplotypes (columns)

The is assumed to be called "alllist" which can then be processed by sourcing "getpainting.R". Along with some other key parameters:
 
### Getpainting:
This is a script that makes output files of most of what you want.

For example with the root subsetdataout/testcp we get created:

* subsetdataout/testcpPaintingsAFR/AFRPainting_Chr1.png
* subsetdataout/testcpPaintingsAFR/AFRPainting_Chr2.png
* subsetdataout/testcpPaintingsEUR/EURPainting_Chr1.png
* subsetdataout/testcpPaintingsEUR/EURPainting_Chr2.png

which are the visualisation of the whole genome, as well as a mean painting.

There are lots of other things output too...

## paintingfns.R

This is where all the code lives. We might care most about getPvalMatrix which is the implementation of the poisson binomial. We might reimplement it to work in two passes through the data, the first to get the mean painting probabilities for each individual, and the second to compute the p-values.


