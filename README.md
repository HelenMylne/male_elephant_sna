# test
repository so called because I was just testing out how to get my laptop contents up here, committing each change one at a time so I would still have the full version history, but now I can't get it to correct...

First chapter of PhD research -- for three populations of male elephants, estimates ages, dyad edge weights, and performs regression analyses of edge weight on age and age difference.

For all: MOTNP = Mosi-Oa-Tunya National Park (Zambia); ANP = Amboseli National Park (Kenya); MPNP = Makgadikgadi National Park (Botswana).

Process:
Step 1a -- data_processing) Convert all sightings into matrices of dyadic sightings to make data frames of group-by-individual data
Step 1b -- data_processing) Convert gbi data to counts of frequency with which a dyad were seen together and apart in a 2 year period 
Step 2a -- age estimation) Use ANP data to develop a survival curve for elephants
Step 2b -- age estimation) Use ANP survival curve to convert categorical age estimations in MOTNP and MPNP into probability dsitributions of true age
Step 3 -- edge weight estimation) Use counts of dyad sightings produced in step 1 to estimate dyadic edge weights.
Step 4 -- data_analysis) Perform nodal and dyadic regressions of edge weight on age/age difference, using age distributions calculated in step 2 as a predictor and edge weight distributions from 3 as a y variable.

