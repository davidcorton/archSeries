#archSeries
##Archaeological time series tools

##Nb. THIS DOCUMENT IS NOT 100% UP-TO-DATE. UPDATES COMING SOON.

###Contributors so far
* David Orton, BioArCh, Department of Archaeology, University of York
* James Morris, School of Forensic and Investigative Sciences, University of Central Lancashire

###Overview
This project aims to develop tools for constructing and comparing frequency time series from archaeological data. In particular, we're looking at ways of synthesising ecofactual data from multiple sites and contexts with varying dates and dating precision, using both aoristic and simulation-based approaches (inspired to a great extent by Crema 2012). More experimentally, we're also developing functions for calibrating ecofactual data against time series of research intensity, e.g. based on volumes of processed environmental samples or numbers of excavated contexts.

The files in this repo are designed for use on a pilot dataset of contexts, environmental samples, and zooarchaeological finds supplied by MOLA, one of London's largest archaeological contractors. The data themselves are not included here for obvious reasons).

This is an ongoing project, so this README is intended primarily as a place to update progress and discuss problems. Ultimately, we intend to turn the functions in **date_functions.R** into an R package.

###Main files
1. **date_functions.R** - source file containing the functions under development.
2. **London_analysis.R** - script demonstrating the use of some of these functions to assess variations in research intensity over the course of London's 2000-year existence.
3. **fresh_vs_marine.R** - script expoloratory analysis revisiting (for London) the 'Fish Event Horizon' (FEH) phenomenon, which saw a sudden shift towards marine fish consumption in medieval England at around AD1000 (Barrett et al. 2004).
4. **CPUE_code.R** - code for a forthcoming paper exploring the FEH in London **WORK IN PROGRESS**.
5. **London_prep.R** - script for cleaning and formatting the datasets as supplied by the archaeological contractor. Obviously this is specific to our dataset and unlikely to be of general use.

##Required packages

* data.table
* reshape2

###Analysis functions
The core functions so far are `aorist` and `date.simulate`, which each estimate a chronological distribution from a data table of entities with date ranges, but using two very different approaches. `dummy.simulate` is used to simulate idealised chronological distributions against which empirical results can be compared. `freq.simulate` is a wrapper function that calls both `date.simulate` and `dummy.simulate` on the same data. `comp.simulate` uses a factor variable to compare the results of any of the above simulation functions between subsets of a dataset.

###1. aorist
Calculates the aoristic sum (*sensu* Crema 2012) from given date ranges for a set of entities, for example representing discrete archaeological contexts or samples. Optionally, these can be weighted, e.g. representing the number of items within a context or the volume of a soil sample. The function uses the data.table package for speed, but is nonetheless quite computationally intensive as we weren't able to vectorise all of the calculations.

**Arguments:**

* `data` is a data table with (at least) two numeric columns called Start and End.
* `start.date` is a single numeric value for the start of the time period to be analysed (defaults to 0).
* `end.date` is a single numeric value for the end end of the time period to be analysed (defaults to 2000).
* `bin.width` is a single numeric value setting the resolution of the analysis, in years (defaults to 100).
* `weight` is a numeric vector giving a weight for each context/entity (defaults to a constant of 1).
* `progress` is a logical value specifying whether a progress report should be made (defaults to TRUE)

**Returns:**
A two-column data table with the aoristic sum itself (numeric) and bin labels (character). The number of rows in this data table will be `(end.date-start.date)/bin.width`.

**Also outputs (->>):**
`breaks`: a numeric vector of breaks points.
`params`: a character value summarising the arguments, for use in naming output files.

###2. comp.simulate
Applies `date.simulate`, `dummy.simulate`, or `freq.simulate` to multiple subsets of a data table, combining the results into a single results table (and, optionally, a single summary results table).

**Arguments**
Nb. most of these are simply passed to the specified simulation function, so may or may not be relevant depending on which function that is.

* `data` is a data table with (at least) two numeric columns called Start and End. If additional factor columns are provided, these can be used to filter the data using the `filter.field` and `filter.values` arguments, below.
* `probs` is a vector of relative probabilities from which to sample, e.g. the output of an `aorist` call. These don't need to sum to one -- they will effectively be scaled. The length of this vector should match the desired number of bins in the output. Alternatively, omit this argument to use uniform probabilities (defaults to 1, which will be recycled to the correct number of bins).
* `weight` is a numeric vector giving a weight for each context/entity (defaults to a constant of 1).
* `comp.values` is a vector of values of `comp.field` to be compared. If not provided, all unique values will be compared.
* `comp.field` is a character vector of length one, denoting the name of the coumn by which the data will be subsetted.
* `comp.fun` is the simulation function to be used (defaults to `date.simulate`).
* `filter.field` is a character vector of length one, denoting the name of a column that will be used to filter the data. Defaults to "group"" but ignored unless 'filter.values' is set.
* `filter.values` is a character vector indicating which rows should be included in analysis (defaults to NULL, in which case all values are included).
* `quant.list` is a numeric vector of quantiles to be included in the summary output (defaults to median, quartiles, 2.5th and 97.5th percentiles).
* `start.date` is a single numeric value for the start of the time period to be analysed (defaults to 0).
* `end.date` is a single numeric value for the end end of the time period to be analysed (defaults to 2000).
* `bin.width` is a single numeric value setting the resolution of the analysis, in years (defaults to 100).
* `reps` is the number of times the simulation will be run (defaults to 100).
* `RoC` is a logical value indicating whether rates-of-change should be calculated and appended to the output (defaults to FALSE). Nb. setting to TRUE will make the function *much* slower.
* `summ` is a logical value indicating whether a summary results table should be returned along with the full results (defaults to FALSE).

**Returns:**
A data table giving the sum of weight for each group, for each bin in each repeat, for both the 'real' and 'dummy' simulations (plus columns for rates of change, if RoC==TRUE).
If `summ`=TRUE, a list containing the above and a second data table giving the specified quantiles of the full simulation results for each bin, again for both the 'real' and 'dummy' simulations (and for rates of change, if RoC=TRUE).

###3. date.simulate
Applies a simulation approach to estimate chronological distribution from given date ranges for a set of entities, for example representing discrete archaeological contexts or samples (this is again based on the procedure described by Crema (2011) and has advantages over `aorist` as set out in that paper). Optionally, the entities can be weighted, e.g. representing the number of items within a context or the volume of a soil sample.

The function simulates a date (year) for each entity, assuming a uniform distribution within the limits of its date range, and places the results into chronological bins of a specified resultion. This process is repeated a specified number of times, allowing the number of points falling into each chronological bin to be estimated, with confidence intervals. Weights are applied **after** sampling, such that heavily weighted entities will tend to increase confidence intervals within their date ranges (this will typically be an appropriately conservative approach in an archaeological scenario).

In future this function may be expanded to permit distributions other than uniform, although it is more likely these will end up as separate functions.

**Arguments:**


* `data` is a data table with (at least) two numeric columns called Start and End. If additional factor columns are provided, these can be used to filter the data using the `filter.field` and `filter.values` arguments, below.
* `weight` is a numeric vector giving a weight for each context/entity (defaults to a constant of 1).
* `filter` is a character vector indicating which rows should be included in analysis (defaults to NULL). It will be ignored if no "group" column is provided in 'data'.
* `start.date` is a single numeric value for the start of the time period to be analysed (defaults to 0).
* `end.date` is a single numeric value for the end end of the time period to be analysed (defaults to 2000).
* `bin.width` is a single numeric value setting the resolution of the analysis, in years (defaults to 100).
* `reps` is the number of times the simulation will be run (defaults to 100).

* `RoC` is a logical value indicating whether rates-of-change should be calculated and appended to the output (defaults to FALSE). Nb. setting to TRUE will make the function much slower.

**Returns:**
A long-format data table giving the sum of weight for each bin in each repeat.

**Also outputs (->>):**
`breaks`: a numeric vector of breaks points.
`params`: a character value summarising the arguments, for use in naming output files.

###4. dummy.simulate
Simulates a 'dummy' chronological distribution within specified date limits by sampling from within a distribution defined by an input vector. Designed for use with date.simulate, particularly within wrapper functions like freq.simulate. The idea here is to simulate chronological distributions based on the same number of entities (with the same weights) as a date.simulate call, but unconstrained by the known date ranges.

**Arguments:**

* `weight` is a numeric vector (or data frame/data.table) representing (weighted) instances to be simulated. If a data frame/data table, must have a row named "weight". If given an additional character column called "group", this can be used to select rows for analysis using the `filter` argument. Alternatively, the number of (unweighted) entities to simulate can be set by passing a numeric vector of length one to `weight`.
* `probs` is a vector of relative probabilities from which to sample, e.g. the output of an `aorist` call. These don't need to sum to one -- they will effectively be scaled. The length of this vector should match the desired number of bins in the output. Alternatively, omit this argument to use uniform probabilities (defaults to 1, which will be recycled to the correct number of bins).
* `breaks` is a numeric vector of length one greater than the desired number of bins, indicating the break points between bins. Typically, this might be the `breaks` vector output by the `aorist` or `date_simulate` functions. If omitted, `bin.widths` will be used to set the breaks instead.
* `filter` is a character vector indicating which rows should be included in analysis (defaults to NULL). It will be ignored if no "group" column is provided in `weight`.
* `start.date` is a single numeric value for the start of the time period to be analysed (defaults to 0). Overruled if `breaks` set.
* `end.date` is a single numeric value for the end end of the time period to be analysed (defaults to 2000). Overruled if `breaks` set.
* `bin.width` is a single numeric value setting the resolution of the analysis, in years (defaults to 100). Overruled if `breaks` set.
* `reps` is the number of times the simulation will be run (defaults to 100).
* `RoC` is a logical value indicating whether rates-of-change should be calculated and appended to the output (defaults to FALSE). Nb. setting to TRUE will make the function much slower.

**Returns**
A long-format data table giving the sum of weight for each bin in each repeat.

###5. freq.simulate
Performs both 'real' and dummy simulation (using `date.simulate` and `dummy simulate` respectively) on a set of entities with date ranges and optionally weights so that the two can be compared to detect deviation from a null hypothesis. The dummy set is generated using the same number of entities and the same weights as for the 'real' set. Both the full simulation results and a summary dataset are returned, and optionally also saved as .csv files.

Optionally, the function can also calculate the rate of change (ROC) between each bin and the next in each simulation (for both real and dummy sets). This allows one to test hypotheses concerning increases or decreases at certain points in the time series.

**Arguments:**

* `data` is a data.frame or data.table with columns including Start, End, and Frag. If additional factor columns are provided, these can be used to filter the data using the `filter.field` and `filter.values` arguments, below.
* `probs` is an optional vector of relative probabilities to be passed to dummy.simulate. Defaults to 1, resulting in a uniform dummy set. 
* `filter.field` is a character vector of length one, denoting the name of a column that will be used to filter the data. Defaults to "group"" but ignored unless 'filter.values' is set.
* `filter.values` is a character vector indicating which rows should be included in analysis (defaults to NULL, in which case all values are included).
* `quant.list` is a numeric vector of quantiles to be included in the summary output (defaults to median, quartiles, 2.5th and 97.5th percentiles).
* `start.date` is a single numeric value for the start of the time period to be analysed (defaults to 0).
* `end.date` is a single numeric value for the end end of the time period to be analysed (defaults to 2000).
* `bin.width` is a single numeric value setting the resolution of the analysis, in years (defaults to 100). Overruled if `probs` has length>1.
* `reps` is the number of times the simulation will be run (defaults to 100).
* `RoC` is a logical value indicating whether rates-of-change should be calculated and appended to the output (defaults to FALSE). Nb. setting to TRUE will make the function much slower.
*`save.full` is a logical value indicating whether the full results table should be saved as a csv file (defauls to FALSE).
*`save.summ` is a logical value indicating whether the summary results table should be saved as a csv file (defauls to FALSE).

**Returns:**
A list of length two:

1. A long-format data table giving the sum of weight for each bin in each repeat, for both the 'real' and 'dummy' simulations (plus columns for rates of change, if RoC==TRUE).
2. A data table giving the specified quantiles of the full simulation results for each bin, again for both the 'real' and 'dummy' simulations (and for rates of change, if RoC=TRUE).

**Also outputs:**
Depending on arguments, .csv files for full and/or summary results

###6. sim.summ
This is a utility function that takes the long-format output from `date.simulate` or `dummy.simulate` and creates a wide-format summary dataset based on specified quantiles.

**Arguments**

* `results` is a data table with at least four columns, the first three of which must be called "rep.no", "bin", and "bin.no", while the rest can take any name. Typically this is the output from `date.simulate` or `dummy.simulate`.
* `summ.col` is a numeric vector indicating the columns to be summarised. By default, all columns apart from the first three are summarised.
* `quant.list` is a numeric vector of quantiles to be included in the summary output (defaults to `c(0.025, 0.25, 0.5, 0.75, 0.975)`, i.e. the median, quartiles, 2.5th and 97.5th percentiles).

**Returns**
A wide-format data table with a factor column for bin names and a numeric column for each specified quantile.

##Plotting functions
These functions are all designed to work with the output of the simulate functions above. `box.chron` and `lines.chron` use the full results while `poly.chron` uses summary results, but all will work if passed a list with both full and summary results (i.e. the default output from the simulation functions).

###1. axis.setup
This is a utility function that scans the data and plots an appropriate set of axes. This will typically be called from within one of the other plotting functions, but can also be used manually.

###2. box.chron
This plots 

###3. lines.chron
Plots every single simulation run for each specified variable as a separate semi-transparent line.

###4. poly.chron
Plots a confidence band for each specified variable.

##References

* Barrett, J.H., A.M. Locker & C.M. Roberts (2004) The origins of intensive marine fishing in medieval Europe: the English evidence. *Proceedings of the Royal Society of London* B, **271**, 2417-2421.
* Crema, E. (2012) Modelling temporal uncertainty in archaeological analysis. *Journal of Archaeological Method and Theory*, **19**, 440-461. 

##Current issues to work on
1. Replace comp.sim with extra arguments for date.sim and dummy.sim, that can be passed through freq.sim. In so doing, switch to a single-simulation approach, as for cpue.
2. Add progress reporters for ROC routines.
4. Update freq.simulate and comp.simulate to simulate context date ranges only once (see 1)
5. make plotting function compatible with aorist
6. New function(s) to generate model distributions to feed into dummy.simulate?
