# Constructing multiple, independent analyses in the regression discontinuity design with multiple cutoffs

## Overview 

The regression discontinuity (RD) design is a commonly used non-experimental approach for evaluating policy or program effects. However, this approach heavily relies on the untestable assumption that observations or average potential outcomes near or at the cutoff are comparable. When there are multiple cutoffs that create several discontinuities in the treatment assignments, factors that can lead this assumption to the failure at one cutoff may overlap with those at other cutoffs, invalidating the causal effects from each cutoff.

In this study, we propose a novel method for testing the causal hypothesis of RD treatment effects that remain valid even when the assumption commonly considered in the RD design does not hold at each cutoff. We propose to leverage a set of instrumental variables (IVs) constructed from multiple cutoffs and combine the evidence from multiple IVs with a direct comparison under the local randomization framework. This reinforced design that combines multiple factors from a single data can produce several, nearly independent inferential results that depend on very different assumptions with each other. Our proposed approach can be extended to a fuzzy RD design. We apply our method to evaluate the effect of having access to higher achievement schools on students' academic performances in Romania.

This repository provides the code and sample data that can reproduce the results in the manuscript.

## Data


Our data application study of school admission is based on this [publicly available data](https://www.openicpsr.org/openicpsr/project/112645/version/V1/view?path=/openicpsr/112645/fcr:versions/V1/data/data-AER-1.dta&type=file). The pre-processing procedure is demonstrated in the `Code/preprocess.R`, with the data dictionary available. 

`Data/data-AER-1.dta` is the original data file and `Data/analysis.dat.RData` is the data used for analysis in `Code/realdata.R`. 

- `Y`: Baccalaureate exam score
- `D`: An indicator of entering the best school in town
- `C`: Cutoff
- `Students`: Student ID
- `Towns`: Towns ID
- `W`: Transition score 
- `X1`-`X12`: List of potential cutoff-based IVs 
- `X13`: Indicator of `W` >= `C` 
- `res.Y`: Transformed outcomes (residuals from the third-degree polynomial regression model on `W`)





## Code

- `Code/simple_sim.R`: generates the main simulation studies with five, individual-specific cutoffs and ran simulations while varying the value of a parameter(s) that could create the bias for each analysis. Both of a sharp and a fuzzy regression discontinuity design can be considered. 

- `Code/cluster_sim.R`: generates the simulation studies with cluster-level cutoffs 

- `Code/preprocess.R`: illustrates the pre-processing procedure from the original data

- `Code/realdata.R`: analyzes the real data application studies using the administrative data from Romania. There are 21 towns available each of which has own cutoff for entering the best school. 

