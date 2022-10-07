import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as const
import math
import time
from limb_darkening import *

from IPython.display import clear_output
from scipy.special import wofz

plt.rcParams['figure.figsize'] = (10, 7)
plt.rcParams['legend.frameon'] = True
plt.rcParams['legend.fontsize'] = 15
plt.rcParams['legend.borderpad'] = 0.2
plt.rcParams['legend.labelspacing'] = 0.2
plt.rcParams['legend.handletextpad'] = 0.2
plt.rcParams['font.family'] = 'stixgeneral'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.size'] = 16



def gen_rays(rs_val, c_x, c_y, n_points): #no change
    '''
    Function: generates points inside a circle of radius "rs_val". 
              These points are considered as the lightrays emerging from a star.
    
    Inputs:
        rs_val - Radius of the star.
        c_x - x-coordinate of the center of the star.
        c-y - y-coordinate of the center of the star.
        n_points - number of gridpoints that will be used to make the star.
        
    Outputs: Coordinates of the lightrays/points on the star.
    '''
    x = np.linspace(-rs_val,rs_val,n_points)
    y = np.linspace(-rs_val,rs_val,n_points)

    xx, yy = np.meshgrid(x, y)


    xx = xx.reshape((np.prod(xx.shape),))
    yy = yy.reshape((np.prod(xx.shape),))

    test = np.stack((xx, yy), axis = 1)
    
    print("Generating star ...")

    indices_star = np.where((xx-c_x)**2+(yy-c_y)**2 <= rs_val**2)[0]
    coord = test[indices_star]
    
    return coord


def gen_atm(cx, cy, R_atm):
    '''
    Function: stores the coordinates of stellar lightrays falling inside the
              atmopshere as the planet transits.
              
    Inputs:
        R_atm - Radius of the atmosphere.
        cx - x-coordinate of the center of the planet.
        cy - y-coordinate of the center of the planet.

    Outputs: Coordinates of the lightrays/points coming through the atmosphere.
    '''
    
    rays_on_atm  = star_coord[np.where(((starx-cx)**2+(stary-cy)**2>R_pl**2)&((starx-cx)**2+(stary-cy)**2<=R_atm**2))[0]] 

    weights_atm = globals()[stellar_model](coeffs,rays_on_atm/R_star)
    
    print("Generating atmopshere ...")
    
    return rays_on_atm, weights_atm

def voigt(del_nu,a,g):

    '''
    Function: returns the Voigt line shape at x with Lorentzian component HWHM 'g'
        and Gaussian component HWHM 'a'.
        
    Inputs:
        del_nu - Frequency offset of the voigt function from the peak value.
        a - or alpha, FWHM of the gaussian part of the voigt function.
        g - or gamms, FWHM of the lorentzian part of the voigt function.

    Outputs: Value of the voigt function at a given x-value (here, frequency) offset.
    '''

    sigma = a / np.sqrt(2 * np.log(2))

    return np.real(wofz((del_nu + 1j*g)/sigma/np.sqrt(2.0))) / sigma/np.sqrt(2.0*np.pi)

###############################################################################
#Parameters of the planetary system

r_wasp69 = 0.813 * 6.96e+10 #centimeters
r_wasp69b = 0.1028*r_wasp69 #centimeters
M_wasp69b = 0.29 * 1.898 * 10**30

R_pl = r_wasp69b
R_star = r_wasp69
M_pl = M_wasp69b

omega = 1.88002696e-5 #rotation frequency of the planet

transit_imp_parm = 0.686 * r_wasp69 #impact parameter of the transit chord

orb_period = 3.8681382 #days

au = const.au.cgs.value

semi_maj_axis =  0.04525 * au #semi-major axis

orb_vel = 2*np.pi*semi_maj_axis/(orb_period*86400) #cmsec^-1

coeff_quad = [0.27168999, 0.26061000] #limb darkening coeffs from EXOFAST website

stellar_model = 'quadratic' #limb darkening model of the star

coeffs = coeff_quad

cx_star = 0

cy_star = 0

star_coord = gen_rays(R_star, cx_star,cy_star, 100)

weights_allrays = globals()[stellar_model](coeffs,star_coord/R_star)

starx = star_coord[:,0]
stary = star_coord[:,1]

G = const.G.cgs.value


###############################################################################
#Parameters of the atmosphere + model

R_atm = 5.41*R_pl 

rest_frame = "planetary"

points_along_ray = 100

###############################################################################
#1D isothermal Parker model file

test_model = "parker_test_model.txt" # a test parker model, T_atm = 7700K and M_dot = 10^11.8 g/sec

T_atm = 7700

he_dens_file = np.genfromtxt(test_model)

altitude = he_dens_file[:,0]
velocity = he_dens_file[:,2] #outflowing velocity.
n_he = he_dens_file[:,6] # helium number density.
he_trip_frac = he_dens_file[:,8] # fraction of metastable helium triplet.

###############################################################################


e = 4.8032e-10
m_e = const.m_e.cgs.value
c = const.c.cgs.value
k_b = const.k_B.cgs.value
m_he = 6.6464731e-24 #grams

he_trip = np.array([10829.09114e-8,  10830.25010e-8,  10830.33977e-8])

f = np.array([5.9e-02, 1.7974e-01, 2.9958e-01]) #oscillator strengths of the He triplet lines, values taken from NIST database

sigma_o = math.pi*e**2*f/(m_e*c)

gamma = 1.0216e7/(4*math.pi) 

lam  = np.linspace(10828,10832,70)

nu = c*10**8/lam




def gen_spec(cx, cy, t_o, sr_val, Ext_dnw, v_dnw, rest_frame):
    
    '''
    Function: generates the transmission spectrum in the band of the helium 
              triplet with given atmospheric parameters.
        
    Inputs:
        
        cx - x-coordinate of the center of the planet.
        cy - y-coordinate of the center of the planet.
        t_o - Temperature of the atmosphere/parker model.
        sr_val - value of the super rotation constant.
        Ext_dnw - Extent of the day-to-night side winds, given in radius of the
                  planet units.
        v_dnw - velocity of the day-to-night side winds in cm/sec.
        rest_frame - rest frame of the spectrum (planetary/stellar).

    Outputs: 
        
        f_norm - 1D array of normalised spectrum calculated as f_sum/f_tot. 
                 Here, f_sum is the total flux output when the planet is 
                 in transit and the atmosphere is absorbing radiation. f_tot is
                 when the planet is still in transit but the atmosphere is not 
                 absorbing any radiation.
                 
        f_inout - 1D array of normalised spectrum calculated as f_in/f_out. 
                  Here, f_in is the total in-transit flux. f_out is the total 
                  out-of-transit flux i.e, just the star's flux.
    '''
    
    t_start = time.time()
    
    if rest_frame == 'planetary':
        v_corr = 0
    if rest_frame == 'stellar':
        v_corr = orb_vel * np.sin(-cx/semi_maj_axis)

    omega_val = sr_val * omega
    
    rays_pl = star_coord[np.where(((starx-cx)**2+(stary-cy)**2<=R_pl**2))[0]] #lightrays blocked by the planet
    weights_pl = globals()[stellar_model](coeffs, rays_pl/R_star) #weights array of those rays blocked behind the planet
    
    rays_atm, weights_atm = gen_atm(cx, cy, R_atm)
    
    print("Number of lightrays in the atmosphere - ", len(weights_atm))

    weights_left_on_star = np.sum(weights_allrays)- (np.sum(weights_atm) + np.sum(weights_pl)) #amount of flux that is remaining on the star.
                                                 
    rays_atm = rays_atm - [cx, cy] #shifting coordinate frame to center of planet

    b_atm = np.sqrt((rays_atm[:,0]**2+rays_atm[:,1]**2)) #impact params of all the rays in the atmosphere

    v_rot = -omega_val*rays_atm[:,0]  #rotation velocity of the longitude of the atmosphere
    
    segments = np.linspace(b_atm,R_atm,points_along_ray) #segments in the lightrays
    
    #optical depth, and there by output flux, is calculated at each of these segments

    midpoints = (segments[1:,:]+segments[:-1,:])/2 #midpoints of these segments, this step is used to 
    #calculate the integral using the midpoint rule
    
    dr = np.diff(midpoints, axis = 0)[0]

    theta = np.arccos(np.divide(b_atm,midpoints)) # angle subtended by the tangent at these segments
    v_out = np.interp(midpoints/R_pl, altitude, velocity) #outlfowing velocity of these segments
    n3_r_val = np.interp(midpoints/R_pl, altitude, n_he*he_trip_frac)  #helium triplet fractions at these segments
    
    dnw_indices = np.where(b_atm<=Ext_dnw)[0] # indices of the b_atm array elements which 
    #fall inside the extend of the day-to-night winds' region of the atmosphere
 
    v_los1 = np.multiply(v_rot, np.cos(theta)) - np.multiply(v_out, np.sin(theta)) # line of sight velocity of the -z hemisphere of the atmosphere
    v_los2 = np.multiply(v_rot, np.cos(theta)) + np.multiply(v_out, np.sin(theta)) # line of sight velocity of the +z hemisphere of the atmsophere
    
    v_los1[:,dnw_indices] = v_los1[:,dnw_indices] + v_dnw #same arrays as before but day-to-night winds magnitude added
    v_los2[:,dnw_indices] = v_los2[:,dnw_indices] + v_dnw

    alpha = np.sqrt(2*np.log(2)*k_b*t_o/m_he)*(1/he_trip) #FWHM of the gaussian part of voigt profile
    
    flux_output = []
    
    print("Generating spectrum ...")

    for i,nu_i in enumerate(nu):

        tau_b = 0
        tau_b1 = 0
        tau_b2 = 0

        for triplet_index in range(3):

                linecenter1 = (c/he_trip[triplet_index])*(1+(v_los1+v_corr)/c) #linecenter of spectral line from left half of the ray(array of 1000 elements)
                linecenter2 = (c/he_trip[triplet_index])*(1+(v_los2+v_corr)/c) #(array of 1000 elements)

                del_nu1 = nu_i*np.ones(np.shape(midpoints)) - linecenter1         
                del_nu2 = nu_i*np.ones(np.shape(midpoints)) - linecenter2

                tau_b1 = np.sum(n3_r_val*sigma_o[triplet_index]*voigt(del_nu1,alpha[triplet_index],gamma)* dr *midpoints/(np.sqrt(midpoints**2-(b_atm)**2)), axis = 0)
                tau_b2 = np.sum(n3_r_val*sigma_o[triplet_index]*voigt(del_nu2,alpha[triplet_index],gamma)* dr * midpoints/(np.sqrt(midpoints**2-(b_atm)**2)), axis = 0)
                tau_b += tau_b1 + tau_b2

        flux_output.append(weights_atm * np.exp(-1*tau_b))  
        
        # print("Progress = {}/{}".format(i,len(nu)))
        # clear_output(wait=True)

    f_atm = np.sum(flux_output,axis = 1) #output flux from the atmosphere

    f_rem = np.ones(len(nu)) * np.sum(weights_left_on_star)

    f_sum = f_atm + f_rem #in transit flux

    f_tot = np.sum(weights_atm) + np.sum(weights_left_on_star) # in transit flux when atmosphere is transparent

    f_norm  = f_sum/f_tot
     
    f_out = np.sum(weights_allrays) #out of transit flux
    
    f_inout = f_sum/f_out
    
    t_end = time.time()
    
    print("Time taken - {} sec".format(np.round(t_end - t_start,2)))
    
    return f_norm, f_inout

cy = transit_imp_parm

f_norm, f_inout = gen_spec(0, cy, T_atm, sr_val = 1, Ext_dnw = 1 * R_atm, v_dnw = 0e5, rest_frame='planetary')

plt.plot(lam, f_norm, color = "k")
# plt.plot(lam, f_inout, color = "k")

for item in he_trip:
    plt.axvline(item*1e8, linestyle = "--", color = "g")
    
plt.xlabel("Wavelength ($\AA$)")
plt.ylabel("Normalised flux")
plt.show()
