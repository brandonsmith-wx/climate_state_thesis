#!/usr/bin/env python3

import numpy as np
from fastkde import fastKDE

def contour_levels_2d(mypdf,v1,v2,prob=0.8,nlevels=4):
    """
    This function determines the level values that are used 
    to plot two-dimensional contour estimates. The function
    determines the levels such that low probability noise is 
    filtered through the value in 'prob' which determines how 
    much data is encompassed by the contours (0.95 seems to work)
    relatively generically but can be set through the keys. 

    NOTE: The function should generally work for any equally spaced
          probability density estimate. However, it has been only been 
          tested extensively with the fastKDE python package.

    :mypdf:     2-d field of probability density function
    :v1:        vector of values in x direction
    :v2:        vector of values in y direction
    :prob:      fraction of data that is encompassed by the outermost
                contour level
    :nlevels:   number of contour levels
    """

    # Define the space of each grid-point for cumulative integration
    area=(v1[1]-v1[0])*(v2[1]-v2[0])
    # Retrieve the unique likelihood values and there respective counts
    unique,counts=np.unique(mypdf,return_counts=True) 
    # Flip the arrays in order to integrate from largest likelihood to lowest
    unique = np.flip(unique)
    counts = np.flip(counts)
    #Cumulative sum
    cum_sum = 0 
    # Integer to loop through unique likelihood values
    k=0
    # Loop until 'prob' percent is encompassed by the contours
    while cum_sum<prob:
        # Integrate probability density function from largest to lowest probablity
        cum_sum+=unique[k]*counts[k]*area
        k+=1

    # Determine the levels for the contours  
    clevels = np.linspace(unique[k],unique[0],nlevels+1)[0:-1]
    
    return clevels
