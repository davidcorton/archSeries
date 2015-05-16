#archSeries
##Archaeological time series tools

###Contributors so far
* David Orton, BioArCh, Department of Archaeology, University of York
* James Morris, School of Forensic and Investigative Sciences, University of Central Lancashire

###Overview
This project aims to develop tools for constructing and comparing frequency time series from archaeological data. In particular, we're looking at ways of synthesising ecofactual data from multiple sites and contexts with varying dates and dating precision, using both aoristic and simulation-based approaches (inspired to a great extent by Crema 2011). More experimentally, we're also developing functions for calibrating ecofactual data against time series of research intensity, e.g. based on volumes of processed environmental samples or numbers of excavated contexts.

The files in this repo are designed for use on a pilot dataset of contexts, environmental samples, and zooarchaeological finds supplied by MOLA, one of London's largest archaeological contractors. The data themselves are not included here for obvious reasons).

This is an ongoing project, so this README is intended primarily as a place to update progress and discuss problems. Ultimately, we intend to turn the functions in **date_functions.R** into an R package.

###Main files
1. **date_functions.R** - source file containing the functions under development.
2. **London_analysis.R** - script demonstrating the use of some of these functions to assess variations in research intensity over the course of London's 2000-year existence.
3. **fresh_vs_marine.R** - script expoloratory analysis revisiting the 'Fish Event Horizon' (FEH) phenomenon which saw a sudden shift towards marine fish consumption in medieval England at around AD1000 (Barrett et al. 2004).
4. **CPUE_code.R** - code for a forthcoming paper exploring the FEH in London **WORK IN PROGRESS**.
5. **London_prep.R** - script for cleaning and formatting the datasets as supplied by the archaeological contractor. Obviously this is specific to our dataset and unlikely to be of general use.

###Current issues to work on

