import numpy as np

#Please refer to Kreidberg (2015) for explanation and formulas of these stellar profiles.

def uniform(points):
    '''
    Function - calculates the flux weightage of the lightray point that lies that
               a certain radial distance from the center of the star 
               which has a "UNIFORM" limb darkening profile.
              
    Input - coordinates the points/lightrays from the star in R_star units.       
              
    
    Output - Flux weights of the lightrays from the star.
    
    '''
    
    return np.ones(len(points))

def linear(coeffs, points):
    '''
    Function - calculates the flux weightage of the lightray point that lies that
               a certain radial distance from the center of the star 
               which has a "LINEAR" limb darkening profile.
              
    Input - coordinates the points/lightrays from the star in R_star units.       
              
    
    Output - Flux weights of the lightrays from the star.
    
    '''
    r_list = np.ones(len(points))
    for index,item in enumerate(points):
        r_list[index] = np.sqrt((item[0]**2 + item[1]**2))

    mu = np.sqrt(1-r_list**2)
    a1 = coeffs[0] * (1 - mu)
    a2 = np.ones(len(points)) - a1
    
    return a2


def nonlinear(coeffs, points):
    
    '''
    Function - calculates the flux weightage of the lightray point that lies that
               a certain radial distance from the center of the star 
               which has a "NONLINEAR" limb darkening profile.
              
    Input - coordinates the points/lightrays from the star in R_star units.       
              
    
    Output - Flux weights of the lightrays from the star.
    
    '''
    r_list = np.ones(len(points))
    for index,item in enumerate(points):
        r_list[index] = np.sqrt((item[0]**2 + item[1]**2))
 
    mu = np.sqrt(1-r_list**2)
    
    a1 = coeffs[0] * (1-np.sqrt(mu))
    
    a2 = coeffs[1] * (1-mu)
    
    a3 = coeffs[2] * (1-mu**1.5)
    
    a4 = coeffs[3] * (1-mu**2)
    
    foo = np.ones(len(points))-a1-a2-a3-a4
    
    return foo


def quadratic(coeffs, points):

    '''
    Function - calculates the flux weightage of the lightray point that lies that
               a certain radial distance from the center of the star 
               which has a "QUADRATIC" limb darkening profile.
              
    Input - coordinates the points/lightrays from the star in R_star units.       
              
    
    Output - Flux weights of the lightrays from the star.
    
    '''
    r_list = np.ones(len(points))
    
    for index,item in enumerate(points):
        r_list[index] = np.sqrt((item[0]**2 + item[1]**2))
 
    mu = np.sqrt(1-r_list**2)
    
    a1 = coeffs[0] * (1- mu)
    
    a2 = coeffs[1] * (1-mu)**2
    
    a3 = np.ones(len(points)) - a1 - a2
    
    return a3

def sqroot(coeffs, points):
    
    '''
    Function - calculates the flux weightage of the lightray point that lies that
               a certain radial distance from the center of the star 
               which has a "SQUARE ROOT" limb darkening profile.
              
    Input - coordinates the points/lightrays from the star in R_star units.       
              
    
    Output - Flux weights of the lightrays from the star.
    
    '''
    r_list = np.ones(len(points))
    
    for index,item in enumerate(points):
        r_list[index] = np.sqrt((item[0]**2 + item[1]**2))
 
    mu = np.sqrt(1-r_list**2)
    
    a1 = coeffs[0] * (1- mu)
    
    a2 = coeffs[1] * (1-np.sqrt(mu))
    
    a3 = np.ones(len(points)) - a1 - a2
    
    return a3  


def logarithmic(coeffs, points):
    
    '''
    Function - calculates the flux weightage of the lightray point that lies that
               a certain radial distance from the center of the star 
               which has a "LOGARITHMIC" limb darkening profile.
              
    Input - coordinates the points/lightrays from the star in R_star units.       
              
    
    Output - Flux weights of the lightrays from the star.
    
    '''
    r_list = np.ones(len(points))
    
    for index,item in enumerate(points):
        r_list[index] = np.sqrt((item[0]**2 + item[1]**2))/radius_of_star
 
    mu = np.sqrt(1-r_list**2)
    
    a1 = coeffs[0] * (1- mu)
    
    a2 = coeffs[1] * (1-mu*np.log(mu))
    
    a3 = np.ones(len(points)) - a1 - a2
    
    return a3

def exponential(coeffs, points):
    
    '''
    Function - calculates the flux weightage of the lightray point that lies that
               a certain radial distance from the center of the star 
               which has a "EXPONENTIAL" limb darkening profile.
              
    Input - coordinates the points/lightrays from the star in R_star units.       
              
    
    Output - Flux weights of the lightrays from the star.
    
    '''
    r_list = np.ones(len(points))
    
    for index,item in enumerate(points):
        r_list[index] = np.sqrt((item[0]**2 + item[1]**2))
    mu = np.sqrt(1-r_list**2)
    
    a1 = coeffs[0] * (1- mu)
    
    a2 = coeffs[1] / (1-np.exp(mu))
    
    a3 = np.ones(len(points)) - a1 - a2
    
    return a3 
