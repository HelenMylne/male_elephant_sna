# test
repository so called because I was just testing out how to get my laptop contents up here, committing each change one at a time so I would still have the full version history, but now I can't get it to correct...

First chapter of PhD research -- for three populations of male elephants, estimates ages, dyad edge weights, and performs regression analyses of edge weight on age and age difference.

For all: MOTNP = Mosi-Oa-Tunya National Park (Zambia); ANP = Amboseli National Park (Kenya); MPNP = Makgadikgadi National Park (Botswana).

Process:

Step 1a -- data_processing) Convert all sightings into matrices of dyadic sightings to make data frames of group-by-individual data

Step 1b -- data_processing) Convert gbi data to counts of frequency with which a dyad were seen together and apart in a 2 year period 

Step 2a -- age estimation) Use ANP data to develop a survival curve for elephants

Step 2b -- age estimation) Use ANP survival curve to convert categorical age estimations in MOTNP and MPNP into probability distributions of true age

Step 3 -- edge weight estimation) Use counts of dyad sightings produced in step 1 to estimate dyadic edge weights. This folder contains script from when doing this by a frequentist method, estimating everything directly as SRI -- this is no longer a useful script and purely highlights the difference from the initial method.

Step 4 -- nodal regression) Perform nodal regression of eigenvector and betweenness centrality from node age, using age distributions calculated in step 2 as a predictor and edge weight distributions from 3.

Step 5 -- dyadic regression) Perform dyadic regression of edge weight from node age and age difference, using age distributions calculated in step 2 as a predictor and edge weight distributions from 3.

