# test
repository so called because I was just testing out how to get my laptop contents up here, committing each change one at a time so I would still have the full version history, but now I can't get it to correct...

First chapter of PhD research -- for two populations of male elephants (MOTNP = Mosi-Oa-Tunya National Park, Zambia; ANP = Amboseli National Park, Kenya; MPNP = Makgadikgadi Pans National Park, Botswana, but this one is no longer in use for this study), estimates dyad edge weights, and performs regression analyses of edge weight on age and age difference.

Process in order:


STEP 1: DATA PROCESSING

1a -- motnp_dataprocessing) Convert all sightings from MOTNP into matrices of dyadic sightings to make data frames of group-by-individual data

1b -- anp_dataprocessing1) Convert all sightings from ANP into matrices of dyadic sightings to make data frames of group-by-individual data

1c -- anp_dataprocessing2) Convert gbi data to counts of frequency with which a dyad were seen together and apart in a 2 year period

DO NOT RUN -- motnp_dataprocessing_bisonR.R) creates the same data but using a different method to check that the method used in data processing scripts is working as it should.

DO NOT RUN -- mpnp_dataprocessing1.R) No longer using MPNP data, takes an exceptionally long time to convert gbi-matrix to dataframe

DO NOT RUN -- mpnp_dataprocessing2.R) No longer using MPNP data


SKIP STEP 2: AGE ESTIMATION -- No longer necessary, as the predictive model was not working properly so we switched back to using ordered categorical age categories instead of attempting a continuous age probability distribution for MOTNP

2a -- anp_survivalcurve.R) Use ANP data to develop a survival curve for elephants

2b -- motnp_bayesianestimation_age.R) Use ANP survival curve to convert categorical age estimations of MOTNP elephants into probability distributions of true age


STEP3: EDGE WEIGHT ESTIMATION

3a -- motnp_bayesianestimation_edgeweight.R) Use counts of dyad sightings produced in step 1a to estimate dyadic edge weights for MOTNP

3b -- anp_bayesianestimation_edgeweight.R) Use counts of dyad sightings produced in steps 1b&c to estimate dyadic edge weights for ANP

3c -- anp_plot.R) Produce nice network plots for ALL network types (ANP short, ANP long and MOTNP)

3d -- summarise_edges.R) Produce reportable values and density plots for each population


STEP 4: NODAL REGRESSION

4a -- motnp_nodalregression.R) calculate eigenvector centrality distribution from edge weight distribution, then perform eigenvector centrality ~ node age category

4b -- anp_nodalregression.R) calculate eigenvector centrality distribution from edge weight distribution, then perform eigenvector centrality ~ node age category

4c -- anp_nodalregression_randomeffects.R) plot random effects output from ANP nodal regression models


STEP 5: DYADIC REGRESSION

5a -- motnp_dyadicregression.R) edge weight ~ age category difference + multimembership
