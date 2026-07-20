# SparseHybrids_TCstage1: simulation 1.0
Simulations for stage 1 testcrosses: GCA and SCA

Simulation1.0_testcross:
* 4Loop_testerchoice_dom.R: Script over different dominance levels and number of testers. Runs the main simulation in AlphasimR + analysis, looping over dominance levels and number of testers, then summarises accuracy (e.g., GCA correlations/rankings) and writes outputs to the results folder. This is the updated version with graphs and results presented in the poster session of 20th April 2026 at Roslin :)

* Loop_dom_test_it.R: Script over different dominance levels and number of testers. Runs the main simulation + analysis, looping over dominance levels and number of testers, then summarises accuracy (e.g., GCA correlations/rankings) and writes outputs to the results folder.

* Step1_testcross_dominance.Rmd: Additive + dominance, but still simplified (no full looping; uses a simple dominance setting). Kept as an intermediate/reference.

* Step1_testcross.Rmd: First, simplest version: additive-only trait. Kept as a baseline/reference.

Simulation2.0_theoretical_framework:
*TesterIndex_v2.R: Script to calculate representativeness, redundancy and discrimination metrics within pool. Simulation of QTL and LD for both pools. 

*Tester4Sven_DT_v2_dummyPE.R: Script of the proof-of-concept of aligned SCA for a single tester with completely complementary heterotic pools.

*peskyGCA.R: Script to demonstrate how allele frequencies of the pools are important for short- and long-term selection of parents by their GCA and importance of developing complementarity in the heterotic pools from the basics of QG.
