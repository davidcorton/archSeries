#London date distributions
###Contributors so far:
* David Orton, UCL Institute of Archaeology

###Overview
This repo contains several R files related to a project to quantify how archaeological research/sampling effort in London has been distributed across the c.2000 years since the city's foundation. This will allow various archaeolocal datasets - particularly those relating to environmental variables - to be calibrated for possible biases due to varying research intensity.
These files are designed for use on a pilot dataset supplied by MOLA, one of London's largest archaeological contractors. Tthe data themselves are not included here for obvious reasons).
This is an ongoing project, so this README is intended primarily as a place to update progress and discuss problems.
 

###Files:
1. **London_prep.R** - script for cleaning and formatting the context, sample, periods, and finds-dating data files as supplied by the archaeological contractor. This will allow us to repeat analysis quickly when the raw data are updated.
2. **date_functions.R** - source file defining two custom functions (_aorist_ and _date.simulate_) that estimate overall date distributions for multiple archaeological entities with given date ranges, using aoristic analysis and simulation respectively. 
3. **London_analysis.R** - script that applies these two functions to the London context and sample datasets, saves the results, and presents some simple summary plots **WORK IN PROGRESS**

###Current issues to work on

1. We now have estimated date distributions for contexts and samples using *period*-based dating, but it would be good to compare this with context-level *finds*-based dating. This is going to a bit tricky since each context may have multiple finds records with different date ranges. We COULD just take either the overall range (very conservative) or the overlap (more optimistic, and potentially non-existent), or sample within the combined range with some kind of weighting for any overlap periods. In the latter case we might also think about weighting for the assemblage size (an ordinal variable in the finds dating dataset).
Clearly residual/intrusive records have already been excluded at the data cleaning stage using a flag variable from the original dataset, but it would still be good to find some way of identifying and dealing with apparently contradictory finds records from the same context. Suggestions? 

2. Once we have finds dating up and running, we need to decide how to integrate it with period dating to get the best coverage and the best reliability.

3. We need to establish which samples were actually processed/studied. In particular, we need to distinguish between samples with (a) no fish bone, (b) bone noted but not recorded, and (c) never processed in the first place.

4. Once 3 is resolved, we can think about weighting the sample results by processed volume. 