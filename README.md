General Inversion of Phase Technique (GIPhT)

Copyright (c) Kurt Feigl, Cliff Thurber  All rights reserved.
U.S. Patent No. 7,446,705.

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


Modification History below this line
2007-2008: Kurt
    prototyping
2009-MAR-28: Lee
    clean up for demo on Iceland subsidence: 6 pairs
2009-MAR-29: Kurt
    expand to other examples
2009-APR-4: Lee
    add license check and turn off the diary
2009-MAY-10: Kurt
    add option to use coherence
2009-JUN-18:
    Use signed char variables for all phases for speed
    Introduce pselect = 5 for quadtree
2010-JUN: Get gradients to work with pselect == 7
2011-JUN:
    Fix bug in quadtree routines with midpoint of patch
    Paramterize orbits in terms of incidence angle - correctly
    Handle missing data
    Handle correlated parameters
2011-JUL
    Speed up step 5
    Fix bug with gradients
2011-OCT
    Clean up plots
2011-OCT-11
    for pselect == 7, use test_generalized_paretos to estimate critical
    value of cost
2011-NOV-11
    update for MATLAB R2011b
2012-JAN
    gipht_step3: take max-min for parameter uncertainties
    add bootstrap to anneal4
    measure gradients in dimensionless strain
    Taylor approximation
2012-SEP v. 2.5
    print out derived parameters, too
    identify bugs when step size in latitude (DL) is positive
    Add quadtree coordinates to DST
    Add phaseprefix to gipht.in
    Fix vector components
2013-MAY v. 2.5 (for short course)
    write_dst.m: write quad dimensions to dst_sample
    gipht_step1.m: same as above
    gipht_path.m: handle file separators under Windows and DOS
    quad_tree: STILL TO DO same as above
2013-MAY v. 2.7
