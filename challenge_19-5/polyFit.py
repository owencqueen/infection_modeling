## module polyFit
''' c = polyFit(xData,yData,m).
    Returns coefficients of the polynomial
    p(x) = c[0] + c[1]x + c[2]x^2 +...+ c[m]x^m
    that fits the specified data in the least
    squares sense.

    sigma = stdDev(c,xData,yData).
    Computes the std. deviation between p(x)
    and the data.
'''    

# This is a modified version of the module 
#  provided by the textbook.

'''
Citation for code:

'''  

import numpy as np
from math import sqrt
from gaussPivot import *

def polyFit(xData,yData,m):
    a = np.zeros((m+1,m+1))
    b = np.zeros(m+1)
    s = np.zeros(2*m+1)
    for i in range(len(xData)):
        temp = yData[i]
        for j in range(m+1):
            b[j] = b[j] + temp
            temp = temp*xData[i]
        temp = 1.0
        for j in range(2*m+1):
            s[j] = s[j] + temp
            temp = temp*xData[i]
    for i in range(m+1):
        for j in range(m+1):
            a[i,j] = s[i+j]
    return gaussPivot(a,b)

def stdDev(c,xData,yData):
    
    def evalPoly(c,x):
        m = len(c) - 1
        p = c[m]
        for j in range(m):
            p = p*x + c[m-j-1]
        return p    
    
    n = len(xData) - 1
    m = len(c) - 1
    sigma = 0.0
    for i in range(n+1):
        p = evalPoly(c,xData[i])
        sigma = sigma + (yData[i] - p)**2
    sigma = sqrt(sigma/(n - m))
    return sigma

# Code taken from Dr. Berry's HW writeup:

def evalPoly(c,x): # c stores coefficients for polynomial
    m = len(c) - 1 # (copied from polyFit module)
    p = c[m]
    for j in range(m):
        p = p*x + c[m-j-1]
    return p



# Original code from OQ:
# ----------------------

def RMSE(c, xData, yData):
    '''
    Gives the root mean squared error of a fitted regression
        model vs. the original data points

    Arguments:
    ----------
    c: coefficients for polynomial
    xData: x data points 
    yData: y data points corresponding to x points
    
    Returns:
    --------
    RMSE: root mean squared error of regression fit
    '''
    variances = 0
    n = len(xData)

    for i in range(0, n):

        x = xData[i]
        y = yData[i]
        y_hat = evalPoly(c, x)
        
        variances += ( y_hat - y ) ** 2

    return np.sqrt( variances / n )

def r2ADJ(c, xData, yData):
    '''
    Returns 1 - r2adj for minimization purposes

    Arguments:
    ---------
    c: coefficients for polynomial
    xData: x data points 
    yData: y data points corresponding to x points
    
    Returns:
    --------
    (1 - r2adj): this is 1 - (r^2 adjusted statistic) 
    '''

    n = len(xData)
    SSR = 0 # Sum squared of regression error
    SST = 0 # Sum squared of total error

    y_bar = np.average(yData) # Get average y values

    # Calculate SSR and SST needed for r2 calculation
    for i in range(0, n):

        x = xData[i]
        y = yData[i]

        y_hat = evalPoly(c, x)

        SSR += (y - y_hat) ** 2
        SST += (y - y_bar) ** 2

    r2 = 1 - (SSR / SST) # Calculate r2

    deg = len(c) - 1 # Get degree of polynomial

    # Calculate r2adj
    r2adj = 1 - ( ((1 - r2) * (n - 1)) / (n - deg - 1) )

    return (1 - r2adj)







        
             

