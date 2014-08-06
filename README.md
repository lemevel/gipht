TITLE: General Inversion of Phase Technique (GIPhT)

Copyright (c) Kurt Feigl, Cliff Thurber  All rights reserved.
U.S. Patent No. 7,446,705.


OVERVIEW

Synthetic aperture radar (SAR) is an active remote sensing technique used for measuring geophysical activity on the Earth’s surface. It records microwaves transmitted by a sensor (usually aboard a satellite) and reflected by features on the Earth’s surface (usually on land). The reflected signal contains information in the form of amplitude and phase data, and requires sophisticated post-processing.

A technique known as interferometric SAR (InSAR) measures the difference in phase between two images of the same area, which can be used to measure motion and deformation of the ground. In most applications, the interferogram must be “unwrapped” before it can be interpreted. The unwrapped interferogram may be used to monitor geophysical changes on the Earth’s surface associated with earthquakes, volcanoes, landslides or glaciers, or with the withdrawal of oil, gas, water or minerals by extractive industries. Unwrapping requires considerable computational power and time, and may lead to significant mistakes in the unwrapped interferogram and thus in its intepretation.
THE INVENTIONUW-Madison researchers have developed an algorithm for interpreting an interferogram without the need for unwrapping. To do so, the invention interprets the interferogram by estimating parameters in a quantitative model directly from the wrapped phase data. Alternative unwrapping algorithms have been developed, but these can provide inadequate results in areas where the phase data are imperfect, leading to errors in the unwrapped phase values. Likewise, these algorithms rarely, if ever, provide uncertainty estimates, limiting attempts to weight the data in statistical analysis. Implementation of the invention would reduce the time and resources necessary for advanced interpretation of InSAR data products, and would provide a more accurate result that includes an assessment of the uncertainties of the parameter estimates.
 
APPLICATIONS

    InSAR for monitoring hazardous natural phenomena, e.g., landslides
    InSAR for monitoring subsidence due to extraction, e.g., oil, gas, water

KEY BENEFITS

    Validated on real (noisy) data in a peer-reviewed publication
    Provides a more direct path to a quantitative interpretation than existing methods
    Provides a realistic assessment of uncertainty, unlike existing methods
    Potentially applicable to other phase data, e.g., pattern recognition in color images
     

Publications

Feigl, K. L., and C. H. Thurber (2009)
A method for modelling radar interferograms without phase unwrapping: application to the M 5 Fawnskin, California earthquake of 1992 December 4 
Geophys. J. Int., 176, 491-504. http://dx.doi.org/10.1111/j.1365-246X.2008.03881.x 

Interferometric analysis of synthetic aperture radar images (InSAR) measures the phase shifts between two images acquired at two distinct times. These ambiguous 'wrapped' phase values range from -1/2 to +1/2 cycles. The standard approach interprets the phase values in terms of the change in distance between the ground and the radar instrument by resolving the integer ambiguities in a process known as 'unwrapping'. To avoid unwrapping, we have developed, validated and applied a new method for modelling the wrapped phase data directly. The method defines a cost function in terms of wrapped phase to measure the misfit between the observed and modelled values of phase. By minimizing the cost function with a simulated annealing algorithm, the method estimates parameters in a non-linear model. Since the wrapped phase residuals are compatible with a von Mises distribution, several parametric statistical tests can be used to evaluate the fit of the model to the data. The method, named General Inversion for Phase Technique (GIPhT), can handle noisy, wrapped phase data. Applying GIPhT to two interferograms in the area of Fawnskin, California, we estimate a set of model parameters describing a magnitude 5 aftershock of the 1992 Landers earthquake. The resulting simulation fits the data well. The phase final residuals have a circular mean deviation less than 0.15 cycles per datum. Sampling the final residuals, we find the circular standard deviation of a phase measurement to be approximately 0.2 cycle, corresponding to 6 mm in range.

Contact Information For current licensing status, please contact Mark Staudt at mstaudt@warf.org or (608) 265-3084.


Use of this computer program [GIPhT] is governed by a license agreement
between You and the Wisconsin Alumni Research Foundation (the
"Agreement"). Unauthorized use, reproduction, or distribution of this
program or any portion of it may result in severe civil and criminal
penalties, and will be prosecuted to the maximum extent possible under
the law. By using this program, you agree to be bound by the terms of the
Agreement.

Source code:
Copyright (c) Kurt Feigl, Cliff Thurber, Lee Powell, Peter Sobol, Aaron Masters
University of Wisconsin-Madison


