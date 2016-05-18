This is the Spark version of Lotka Volterra (LV) simulation using particle marginal Metropolis-Hastings (PMMH) algorithm.

The program is based on the Dr. Darren Wilkinson's bayeskit: https://github.com/darrenjw/djwhacks/tree/master/scala/bayeskit.

In this project, a one prey-three predators LV model is develeoped based on the one prey-one predator LV model in Dr. Wilkinson's original bayeskit, and book chapter "Predicting extinction of biological systems with competition" by Branko Ristic and Alex Skvortsov from "Emerging Trends in Computational Biology, Bioinformatics, and Systems Biology: Algorithms and Software Tools".

The project is adapted to run on Apache Spark deployed on the HPC cluster at Scania IT. It will be tested and compared with the MPI version of the same LV simulation, to reveal if Spark will be a good alternative for large scale compute-intensive applications.
