# archSeries
## Tools for Chronological Uncertainty in Archaeology

## Overview

This package is designed to tackle a common problem in archaeology: the integration of information from multiple entities with varying chronological resolution and overlapping date ranges. It includes functions for constructing overall chronological distributions from tables of start dates, end dates, and weights; and for plotting those distributions in various ways. Distributions may be calculated either by basic aoristic analysis or by simulation. In the latter case, 'dummy sets' can also be calculated automatically based on uniformity or on a custom null model, and rates of change between bins can also be simulated. Importantly, archSeries makes use of the data.table package to run simulations extremely quickly.

Borrowing heavily from ideas set out elsewhere (Crema 2012), the code released here was initially developed for meta-
analysis of environmental archaeological data from urban sites (Orton et al. 2017). We hope that its uses within 
archaeology will go beyond this, however -- for example via application to zooarchaeological age-at-death studies (see e.g. Bréhard et al. 2014; the *surv.convert* function is already included to this end) and perhaps eventually to human demographics, though the latter would require a suite of additional features.

At present the package supports aoristic analysis based assumption of uniform probability within date limits, and simulation based on either uniform distribution or beta distribution with specified parameters (see Baxter & Cool 2016). It is set up mainly for frequency data, but the central date.simulate function can also be used with e.g. metrical or isotopic data. In future we hope to implement:

* ability to work with radiocarbon dates (e.g. by interfacing with the Bchron or Rcarbon packages).
* ability to specify null models other than uniform or custom, and to fit these to the data before simulating dummy.
sets (as for example with exponential growth curves in 14C-based demographic studies, e.g. Timpson et al. 2014).
* a global hypothesis tester (see e.g Timpson et al. 2014).
* ability to specify priors for simulations of empirical data, not just for dummy sets (necessary for demographics).
* additional tools for zooarchaeological age-at-death (perhaps interfacing with the zooaRch package).
* faster aorist calculations using data.table::foverlaps() (suggestion from Matt Harris).

The functions fall into two groups:

* **Analysis functions:** take data and produce one or more chronological frequency distributions, potentially including confidence intervals, dummy sets under null models, etc.
* **Plotting functions:** take the output of the analysis functions and plot them in various ways, using R base graphics.

## Analysis functions

There are three main analysis functions: *aorist*, *date.simulate*, and *cpue*. The former is deterministic and produces a single simple table of chronological bins and corresponding probability densities; the latter are both based on simulation and produce output in a standard format. This is a list of two data tables, the first ("full") containing full simulation results and the second ("summary") containing a summary table of quantiles based on the first. In the case of *cpue*, certain arguments can result in the addition of a third component ("small.n"): a list of lists defining areas of low confidence. This standard output format could be thought of as a special class of object, though it is not formally defined as such. In any case, the main plotting functions are designed to work directly with it.

#### **aorist:** calculate the aoristic sum for a group of entities with defined date limits.
Calculates aoristic sum from a table of entities with defined date ranges, based on assumption of uniform probability between limits. Returns a two-column data table listing the calculated total probability density for each chronological bin.

#### **date.simulate:** simulate chronological distribution from a group of entities with defined date limits.
Simulates chronological distributions from a table of entities with defined date ranges, based on specified beta distribution parameters (defaulting to uniformity), then (optionally) simulates a dummy set of the same size drawing from a specified distribution. By default, sums the frequency of entities within specified chronological bins, but user can specify other summary functions (e.g. mean or median) when dealing with non-frequency data. Returns a list with two named elements: "full" is a long-format data table with complete simulation results; "summary" is a data table with summary results by bin.

#### **cpue:** "Catch per Unit Effort".
Estimates the relationship between a numerator ('catch') and a denominator ('effort') variable over time. The idea here is that the former represents  and the latter represents some measure of sampling effort, though there may be other applications. Returns output in the same format as *date.simulate*: a list with elements "full", a long-format data table with complete simulation results, and "summary", a data table giving summary results by bin. Optionally a third element may be added, defining areas where the measure of effort is below a user-defined threshold.

#### **surv.convert:** convert simulated mortality data to survivorship format.
This is the only function so far developed with a view specifically to mortality profiles rather than chronological data (i.e. ages in months or years rather than dates in decades or centuries). It take the output from any of the simulate functions and converts the simulated frequencies into simulated survivorship curves, which are output in the standard format.

#### **RoC.fun:** calculate Rates of Change in simulated variables.
Takes the standard "full" output from one of the archSeries simulation functions and adds columns giving rates of change between bins, for all value columns or a specified subset. Called within simulation functions by setting their 'RoC' arguments to TRUE, but can also be called directly after the fact. Returns the input data table with the addition of one new RoC column for each original value column.

## Plotting functions
There are four main plotting functions in the package: *poly.chron*, *lines.chron*, *box.chron*, and *aorist.plot*. The first three are designed to work with the output of the various simulation functions, and the fourth with output from *aorist*. In addition, two auxiliaries functions called within the above -- *axis.setup* and *grey.zones* -- can also be called on their own if required.

#### **poly.chron:** plot medians and confidence zones for simulation results.
Plots defined confidence intervals (as polygons) and medians (as lines) for output from *date.simulate* or an associated function.

#### **lines.chron:** plot all simulation runs as lines.
Plots every single simulation run for each specified variable as a separate semi-transparent line.

#### **box.chron:** boxplots for simulation output.
Plots a series of boxplots (one per bin) summarising the output from *date.simulate* or an associated function.

#### **aorist.plot:** plot aorist output as a barplot.
Just a wrapper for *barplot* with some tweaks added, e.g. to make bars line up with data plotted by other archSeries functions. NEEDS WORK - DOESN'T YET ALLOW FOR SECONDARY Y-AXIS.

#### **axis.setup:** set up axes for plotting functions
A utility function designed to be used within the various plot functions in this package, but which can also be used alone to set up axes based on simulation data prior to plotting anything.

## References

* Baxter, M.J. & H.E.M. Cool (2016) [Reinventing the wheel? Modelling temporal uncertainty with applications to brooch distributions in Roman Britain](https://doi.org/10.1016/j.jas.2015.12.007). *Journal of Archaeological Science*, **66**, 12-127.
* Bréhard, S., V. Radu, A. Martin, P. Hanot, D. Popovici & A. Bălăşescu (2014) [Food supply strategies in the Romanian Eneolithic: sheep/goat husbandry and fishing activities from Hârşova Tell and Borduşani-Popină (5th Millennium bc)](http://dx.doi.org/10.1179/1461957113Y.0000000051). *European Journal of Archaeology*, **17**, 407–433.
* Crema, E. (2012) [Modelling temporal uncertainty in archaeological analysis](https://doi.org/10.1007/s10816-011-9122-3). *Journal of Archaeological Method and Theory*, **19**, 440-461.
* Orton, D.C., J. Morris & A. Pipe (2017) [Catch Per Unit Research Effort: Sampling Intensity, Chronological Uncertainty, and the Onset of Marine Fish Consumption in Historic London](http://doi.org/10.5334/oq.29). *Open Quaternary*, **2**, 1–26..
* Timpson, A., S. Colledge, E. Crema, K. Edinborough, T. Kerig, K. Manning, M.G. Thomas & S. Shennan (2014) [Reconstructing regional population fluctuations in the European Neolithic using radiocarbon dates: a new case-study using an improved method](https://doi.org/10.1016/j.jas.2014.08.011). *Journal of Archaeological Science*, **52**, 549–557.


