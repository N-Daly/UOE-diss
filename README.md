# UOE-diss
Repo for dissertation related code. The diss looks at modelling two observers in distance sampling scenarios and the effect this has on estimates of spatial intensity.

This has done with R version 4.5.3 and ``renv`` should perserve all the package dependencies.

There are a number of files where and they all depend on each other e.g. by calling functions defined in others or by executing code in others and using the resulting R objects.

# simulation study.R
The file where stuff happens.
This runs the Simulation study where a fictional animal population is generated, fictional observers perform (line transect) distance sampling to estimate the animal density. Five competing statistical models are fitted and compared to the known ground truth. 
This process is repeated a number of times, the actual code for this is in ``one iteration of realisation and fitting.R``.
The true animal density is log gaussian cox process. The repeated iterations may take up to ten minutes each, and model comparisons and fitting times are printed to the console to give a sense of progress.

## one iteration of realisation and fitting.R
performs one iteration of the simulation study just described.
The script cannot be run in isolation as it depends on objects described in the simulation study script, set the number of simulations to one instead.

## sim_results
Where the simulation results get saved to, as .rda files.

## figs
Many scripts save figures to this location.

## observer models.R
This is where the statistical models were defined.
The models differed only in the assumed form of the "detection function", the probability of an observer(s) detecting an animal that is a certain distance away.
There are functions to extract predictions (of say probabilities of detection or animal densities) for each model as well. It can be seen that these prediction extraction functions are repetitive but I decided it was not worth the time to rewrite working code.

## a good ground truth detection function.R
How I eyeballed what choice of parameters to use in the true detection function based of what gave an interesting behaviour to talk about.

## plot figures for report.R
code for presenting simulation results in a reproducible manner.

# function definitions.R
Utility functions and so on.

# simulate 2D ground truth.R
Contains functions for simulating the ground truth of a population whose density is described by a log gaussian cox process, and for the following distance sampling where not all animals are detected by observers.

# mesh construction methods.R
The software used here for spatial statistics, ``inlabru`` uses meshes (triangulations of 2d space) as a computational tool. 
The choice of mesh is critical to the analysis and this script looks at three different methods of constructing meshes around the pedagogical dataset of an ocean survey in the gulf of mexico.
Also explored is the tradeoff between having very "fine"/dense meshes which are good for modelling and the time taken to fit a statistical model as a result.
Note that this code is intense for a regular laptop, taking well over 5 minutes per mesh and using 12Gb of RAM with all other applications closed.

## precomputed stuff
As mentioned, the meshes take a while to make so they are saved in .rda files.

## Prior sensitivity stuff
messing with spline priors.R
various dual hn fits.R
various dual hr fits.R
