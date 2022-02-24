# -*- coding: utf-8 -*-
import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
import os
import sys
import shutil
import configparser
import datetime



# Begin function definitions.
def northwardDist(planet_re, planet_rp, lat1, lat2):
    '''
    Calculate dy between point of interest and previous latitude
    
    Returns the "northward distance" from one latitude on a oblate spheroid
    to another. Does not consider longitude.
     
    Parameters
    ----------
    planet_re : float
        Equatorial radius of the spheroid.
    planet_rp : float
        Polar radius of the spheroid.
    lat1 : float
        Current latitude point.
    lat2 : float
        Previous latitude point.
    
    Returns
    -------
    float
        The northward distance between `lat1` and `lat2`.
        
    Notes
    -----
    The northward distance is given by:
    
    .. math::
        dy = (\phi_{1}-\phi_{2})(\\frac{r_{e}\\cos^{-1}\phi_{mid}}{\sqrt{ \\ 
        1+[\\frac{r_{p}}{r_{e}}]^{2}\\tan^{2}\phi_{mid}}}) \\ 
        (\\sin^{2}\phi_{mid}+[\\frac{r_{e}}{r_{p}}]^{2} \\ 
        \\cos^{2}\phi_{mid})^{-1}
            
    where :math:`\phi_{1}`, :math:`\phi_{2}`, and :math:`\phi_{mid}` are the 
    first, second, and midpoint latitudes, respectively; :math:`r_{e}` is the 
    planet's equatorial radius; and :math:`r_{p}` is the planet's polar radius.
    
    History
    -------
    2/18/2022, planetary@garland.run : v1.0
        Initial published version. Based on IDL code by Dr. Kunio Sayanagi.
    '''      
    lat1_rad = np.pi * (lat1 / 180.0)
    lat2_rad = np.pi * (lat2 / 180.0)
    dlat_rad = lat1_rad - lat2_rad

    lat_mid = 0.5 * (lat1_rad + lat2_rad)
    map_lon = planet_re / np.sqrt((1 + ((planet_rp / planet_re)**2) 
                                   * (np.tan(lat_mid)**2)))
    map_lat = ((map_lon / np.cos(lat_mid)) / ((np.sin(lat_mid)**2) + 
                ((planet_re / planet_rp)**2) * (np.cos(lat_mid)**2)))

    return(dlat_rad * map_lat)


 
def geoTW(p, pPrev, t, tPrev, lat, R, omega, dy, uPrev, eqTuning):
    '''
    Extrapolate cloud-top wind with the geostrophic thermal wind equation.
    
    Uses the geostophic balance version of the thermal wind equation (TWE) to
    calculate winds at depth based on measured cloud-top wind speeds and 
    horizontal (latitudinal) temperature gradients. Includes a term to 
    "turn off" contribution from the coriolis parameter, f, to avoid division
    by zero near the equator.
     
    Parameters
    ----------
    p : float
        Current pressure.
    pPrev : float
        Previous pressure, in the same units as `p`.
    t : float
        Current temperature.
    tPrev : float
        Previous latitude's temperature at the same pressure `p` and in the 
        same units as `t`. 
    lat : float
        Current latitude.
    R : float
        The planet's specific gas constant R, in J/kg*K.
    omega : float
        The planet's rotation rate, in rad/s.
    dy : float
        The "northward distance" between `lat` and the 
        latitude of `tPrev` in meters.
    uPrev : float
        The wind speed at `pPrev`.
    eqTuning : float
        A tuning parameter to "turn off" the coriolis parameter near the 
        equator in the TWE. Lower values result in lower horizontal temperature
        gradients near the equator.
        
    Returns
    -------
    u : float
        The calculated wind speed at `p`.
    
    See Also
    --------
    northwardDist : Calculation of dy.
    
    Notes
    -----
    The discretized form of the TWE implimented is:
    
    .. math::
        u_{n} = u_{n-1}+(\\frac{-R[p_{n-1}-p_{n}]}{fp_{n}}) \\ 
        (\\frac{t_{y-1}-t_{y}}{dy})(\\frac{\\arctan[x\phi ]}{\\frac{\pi }{2}})
        
    where :math:`u_{n}` and :math:`p_{n}` are the current wind and pressure, 
    :math:`u_{n-1}` and :math:`p_{n-1}` are the previous, vertically-displaced 
    wind and pressure, :math:`t_{y}` is the current latitude's temperature, 
    :math:`t_{y-1}` is the previous latitude's temperature, :math:`R` is the 
    specific gas constant of the atmosphere.and :math:`x` is a tuning 
    parameter that determines how fast the third "equatorial" term grows. 
    The third term serves to "turn off" the contribution from :math:`f`, the
    Coriolis parameter, near the equator to prevent :math:`u` blowing up to 
    infinity as :math:`f` approaches zero. :math:`dy` is the 
    "northwards distance", or the north/south distance between the latitudes 
    that :math:`t_{y}` and :math:`t_{y-1}` are measured at.
    
    History
    -------
    2/18/2022, planetary@garland.run : v1.0
        Initial published version. 
    '''  
    # "Textbook" TWE with a tuning term to make the contribution 
    # from f vanish at the equator.
    f          = (2 * omega * np.sin(lat))
    term1      = ((-1 * R * (pPrev - p)) / (f * p))
    term2      = ((tPrev - t) / dy)
    termTuning = (np.arctan(eqTuning * lat) / (0.5 * np.pi))
 
    # Handle the sign flip across the equator
    if lat <= 0.0:
        u = uPrev + (term1 * term2 * termTuning)
    else:
        u = uPrev - (term1 * term2 * termTuning)
    
    return u
  
  

def rEllipsoid(re, f, phi):
    '''
    Calculate radius of an ellipsoid at a particular point on the surface.
    
    Returns the length of the radius to the surface of an ellipsoid with 
    given flattening at a particular latitude.
     
    Parameters
    ----------
    re : float
        Equatorial radius of the ellipsoid.
    f : float
        Flattening of the ellipsoid, given by (`re` - polar radius) / `re`
    phi : float
        Current latitude.
        
    Returns
    -------
    r : float
        The calculated radius at `r`.
    
    History
    -------
    2/18/2022, planetary@garland.run : v1.0
        Initial published version. 
    '''  
    # Convert phi to radians, then calculate r at phi.
    phiRad = np.pi * (phi / 180.0)
    r      = np.sqrt(((re**2) / (1 + (((1 / ((1 - f)**2))) - 1) 
                                 * (np.sin(phiRad)**2))))
    return r



def gAtr(M, re, r, phi, J2, omega):
    '''
    "g at r", calculates gravitational acceleration on at a particular r
    
    Used to calculate local g at each latitude on the planet.
     
    Parameters
    ----------
    M : float
        Mass of the planet, in kg.
    re : float
        Equatorial radius of the planet.
    r :  float
        Radius of interest.
    phi : float
        Current latitude.
    J2 : float
        The 2nd zonal harmonic coeffecient of the planet.
    omega : float
        The planet's rotation rate, in rad/s.
        
    Returns
    -------
    float
        The local gravitational acceleration at `r` and `phi`.
        
    See Also
    --------
    rEllipsoid : Find the `r` value at a particular `phi`.
    
    History
    -------
    2/18/2022, planetary@garland.run : v1.0
        Initial published version. 
    '''  
    # Define universal gravitational constant. 
    G = 6.674e-11 # m^3 kg^-1 s^-2
    
    # A fairly long equation here. Break up into terms for ease of calculation.
    term1 = ((G * M) / (r**2))  
    term2 = (1 - ((9/2) * ((re**2) / (r**2)) * J2 * (np.sin(phi)**2)) + 
             ((3/2) * ((re**2) / (r**2)) * J2))
    term3 = ((omega**2) * r * (np.cos(phi)**2))

    return ((term1 * term2) - term3)



def tFromConstBV(t_0, p, R, cp, g, bvFreq):
    '''
    Extrapolate cloud-top observed temperatures to depth.
    
    Assuming a constant BV frequency, solve the BV equation for temperature
    to extrapolate measured cloud-top temperatures downwards.
     
    Parameters
    ----------
    t_0 : float
        Measured, cloud-top temperature.
    p : array_like
        Pressure levels of interest to calculate t at.
    R : float
        The planet's specific gas constant R, in J/kg*K.
    cp : float
        The planet's specific heat at constant pressure cp, in J/kg*K.
    g : float
        The local gravitational acceleration at latitude of `p`.
    bvFreq : float
        The assumed Brunt-Vaisala frequency.
        
    Returns
    -------
    t : array_like
        The calculated temperature at each input value of `p`.
        
    See Also
    --------
    gAtR: Find the `g` value at a particular latitude.
    
    Notes
    -----
    The BV equation implimented to solve for t is:
    
    .. math:: 
        t_{n} = \\frac{t_{n-1}(\\frac{\kappa g^{2}}{RN^{2}}) \\ 
        (\\frac{p_{n}}{p_{n-1}})^{\kappa}}{\\frac{\kappa g^{2}}{RN^{2}}+ \\ 
        (t_{n-1}[\\frac{p_{n}}{p_{n-1}}]^{\kappa})-t_{n-1}}
        
    where :math:`t_{n}` is the current temperature at pressure 
    :math:`p_{n}`, :math:`t_{n-1}` is the previous temperature at pressure 
    :math:`p_{n-1}`, :math:`g` is the local gravity, :math:`R` is the specific 
    gas constant of the atmosphere, :math:`N` is the BV frequency, 
    and :math:`\kappa` is the Poisson constant = :math:`\\frac{R}{c_{p}}` 
    where :math:`c_{p}` is the specific heat at constant pressure 
    of the atmosphere. 
    
    History
    -------
    2/18/2022, planetary@garland.run : v1.0
        Initial published version. Based on Excel spreadsheet equations by 
        Dr. Kunio Sayanagi.
    '''
    # Initialize t array the same length as the p array. t_0 should come from 
    # the cloudtop temperature map.
    t    = np.empty(p.size)
    t[0] = t_0
    
    # Calculate TN term.
    TN = (((R/cp) * ((g)**2)) / (R * (bvFreq)**2))
    
    # Loop over t array, filling with correct values.
    for i in range(1, t.size):
        numer = (t[i-1] * TN * ((p[i] / p[i-1])**(R/cp)))
        denom = (TN + (t[i-1] * ((p[i] / p[i-1])**(R/cp))) - t[i-1])
        t[i]  = (numer / denom)
    
    return t


 
def tempDigitize(temp, phiStart, phiStop, phiSteps):
    '''
    Interpolates observed latitudinal, cloud-top temperature profiles to a 
    given latitude range and step size.
    
    Parameters
    ----------
    temp : array_like
        Measured, cloud-top temperatures.
    phiStart : float
        Starting latitude of model.
    phiStop : float
        Final latitude of model
    phiSteps : int
        Amount of points between `phiStart` and `phiStop` in model.
        
    Returns
    -------
    tempInterp : array_like
        The interpolated, cloud-top temperatures in the same shape as our
        model's latitude array.
    
    History
    -------
    2/18/2022, planetary@garland.run : v1.0
        Initial published version. 
    '''
    # Temp is an ndarray of arbitrary length with lats and temps. Interpolate  
    # to the size of the phi array.
    xVals       = np.linspace(phiStart, phiStop, phiSteps)
    tempInterp  = np.interp(xVals, temp[:,0], temp[:,1])

    return(tempInterp)



def tempDigitizeMean(low, high, phiStart, phiStop, phiSteps):
    '''
    Interpolates observed latitudinal, cloud-top temperature profiles to a 
    given latitude range and step size.
    
    Also performs a mean between two input temperature observations.
    
    Parameters
    ----------
    low : array_like
        Measured, cloud-top temperatures for "low" observation.
    high : array_like
        Measured, cloud-top temperatures for "high" observation. 
    phiStart : float
        Starting latitude of model.
    phiStop : float
        Final latitude of model
    phiSteps : int
        Amount of points between `phiStart` and `phiStop` in model.
        
    Returns
    -------
    array_like
        The averaged, interpolated, cloud-top temperatures in the same shape 
        as our model's latitude array.
    
    See Also
    --------
    tempDigitize : The base, single observation function.
    
    Notes
    -----
    `low` and `high`, while named as such, are not required to be lower or
    higher than each other. The mean will be performed successfully regardless.
    The parameters are named as they are in the main config.ini.
    
    History
    -------
    2/18/2022, planetary@garland.run : v1.0
        Initial published version. 
    '''
    # Low and high are ndarrays of arbitrary length with lats and temps. 
    # Interpolate them to  size of phi array and 
    # average to get mean temperature.
    xVals      = np.linspace(phiStart, phiStop, phiSteps)
    lowInterp  = np.interp(xVals, low[:,0], low[:,1])
    highInterp = np.interp(xVals, high[:,0], high[:,1])
    
    return(np.mean(np.array([lowInterp, highInterp]), axis=0))
 

 
def windDigitize(wind, phiStart, phiStop, phiSteps):
    '''
    Interpolates observed zonal wind profiles to a given latitude range 
    and step size.
    
    Parameters
    ----------
    wind : array_like
        Measured, cloud-top zonal winds.
    phiStart : float
        Starting latitude of model.
    phiStop : float
        Final latitude of model
    phiSteps : int
        Amount of points between `phiStart` and `phiStop` in model.
        
    Returns
    -------
    windInterp : array_like
        The interpolated, cloud-top zonal winds in the same shape as our
        model's latitude array.
    
    History
    -------
    2/18/2022, planetary@garland.run : v1.0
        Initial published version. 
    '''
    # Wind is a ndarray of arbitrary length with lats and winds.  
    # Interpolate them to size of phi array.
    xVals       = np.linspace(phiStart, phiStop, phiSteps)
    windInterp  = np.interp(xVals, wind[:,0], wind[:,1])
    
    return(windInterp)



# Function to 'smooth' the temperatures with depth. Towards our lower bound,
# temperatures become homogenous with a tanh function.
def tanhSmoothing(x, intersect, tuning):
    '''
    Tanh-based function to "smooth" temperatures with depth.
    
    Towards the lower pressure bound of the model, the horizontal temperature
    gradient approaches zero. Results in homogeneous temperatures when 
    in the deep troposphere. 
    
    Parameters
    ----------
    x : int
        Current index in the model pressure array.
    intersect : int
        Index of the point to begin making temperatures homogeneous.
    tuning : float
        Controls how fast temperatures become homogeneous. Higher values 
        result in a slower transition.
        
    Returns
    -------
    float
        Number from 0.0 -> 1.0 determining the contribution from the actual
        temperature vs. the average temperature at a pressure level.
    
    History
    -------
    2/18/2022, planetary@garland.run : v1.0
        Initial published version. 
    '''
    # Return the input function, smoothed using a tanh.
    return(0.5 + 0.5*np.tanh(((x - intersect)/tuning)))



def main(example=False, exampleName=None):
    # Begin main program.
    # Set user-defined parameters from config. Defined in the config file.
    currentTime = str(datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
    config      = configparser.ConfigParser()
    
    if not(example):
        config.read(os.path.join(os.path.dirname(os.path.abspath(__file__)), 
                                 'config.ini'))
    else:
        config.read(os.path.join(os.path.dirname(os.path.abspath(__file__)), 
                             'examples', exampleName, 'config.ini'))
    
    # INPUT FILE PARAMETERS
    try:
        lowTemps  = config['USER']['lowTemps']
        highTemps = config['USER']['highTemps']
        wind      = config['USER']['wind']
        
        # PLANET PARAMETERS
        planet = config['USER']['planet'].lower()
        omega  = float(config['USER']['omega'])
        R      = float(config['USER']['R'])
        cp     = float(config['USER']['cp'])
        re     = float(config['USER']['re']) 
        rp     = float(config['USER']['rp']) 
        M      = float(config['USER']['M']) 
        J2     = float(config['USER']['J2'])
        
        # MODEL PARAMETERS
        p_0            = float(config['USER']['p_0']) 
        pFinal         = float(config['USER']['pFinal']) 
        pSteps         = int(config['USER']['pSteps']) 
        phiStart       = float(config['USER']['phiStart']) 
        phiStop        = float(config['USER']['phiStop']) 
        phiSteps       = int(config['USER']['phiSteps']) 
        bvFreq         = float(config['USER']['bvFreq']) 
        eqTuning       = float(config['USER']['eqTuning']) 
        tempType       = config['USER']['tempType'].lower()
        intersectDepth = float(config['USER']['intersectDepth']) 
        tuning         = float(config['USER']['tuning'])
        svFilter       = config['USER'].getboolean('svFilter')
        windowLength   = int(config['USER']['windowLength'])
        polyOrder      = int(config['USER']['polyOrder'])
        
        # PLOTTING PARAMETERS
        plotInputTemps   = config['USER'].getboolean('plotInputTemps')
        plotInputWind    = config['USER'].getboolean('plotInputWind')
        plot2DTemps      = config['USER'].getboolean('plot2DTemps')
        plot2DDensity    = config['USER'].getboolean('plot2DDensity')
        plot2DWind       = config['USER'].getboolean('plot2DWind')
        plotWindProfiles = config['USER'].getboolean('plotWindProfiles')
        plot1DStepSize   = int(config['USER']['plot1DStepSize']) 
        dpi              = int(config['USER']['dpi'])
        savePlots        = config['USER'].getboolean('savePlots')
        
        # RESULTS PARAMETERS
        saveData = config['USER'].getboolean('saveData')
    except KeyError:
        print('config.ini not found! Use and rename the configTemplate.ini ' + 
              'in /docs and move to the same directory as tweModel.py')
    
    # Ensure windowLength and polyOrder are set properly.
    if svFilter:
        svFilterCorrect = True
        if (windowLength % 2) == 0:
            print('ERROR: windowLength must be odd. Currently: ' + 
                  str(windowLength))
            svFilterCorrect = False
        if windowLength <= polyOrder:
            print('ERROR: polyOrder (' + str(polyOrder) + ') must be less ' +
                  'than windowLength (' + str(windowLength) + ')')
            svFilterCorrect = False
        if not(svFilterCorrect):
            sys.exit(1)
    
    
    
    # Set data directory.
    if not(example):
        dataPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), 
                                'data', str(planet))
    else:
        dataPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), 
                                'examples', exampleName, 'data', str(planet))
        
    print('Loading data from directory: ' + dataPath + '\n')
    
    
    
    # Set results directory.
    if not(example):
        resultsPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), 
                                   'results')
    else:
        resultsPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), 
                                   'examples', exampleName, 'results')
        
    if not os.path.exists(resultsPath):
        os.mkdir(resultsPath)
    resultsPath = os.path.join(resultsPath, str(planet))
    if not os.path.exists(resultsPath):
        os.mkdir(resultsPath)
    
    # If we're saving anything, make the run directory. 
    if saveData or savePlots:
        if not(example):
            runPath = os.path.join(resultsPath, currentTime)
        else:
            runPath = os.path.join(resultsPath, 'runResults')
        if not os.path.exists(runPath):
            os.mkdir(runPath)
        
        
     
    # Create an array of latitudes to invesitgate.
    print('Building cloud-top temperature / wind structure. \n')
    phi = np.linspace(phiStart, phiStop, phiSteps) # degrees
      
    # Create an array of pressures to investigate, logspace.
    p = np.geomspace(p_0, pFinal, pSteps, endpoint=True)
    
    
    
    # Read temperture array(s) and interpolate them to our phi values. 
    if tempType == 'single':
        low  = np.loadtxt(os.path.join(dataPath, lowTemps), delimiter=',')   
        t    = tempDigitize(low, phiStart, phiStop, phiSteps)
    elif tempType == 'mean':    
        low  = np.loadtxt(os.path.join(dataPath, lowTemps), delimiter=',')
        high = np.loadtxt(os.path.join(dataPath, highTemps), delimiter=',')
        t    = tempDigitizeMean(low, high, phiStart, phiStop, phiSteps)
    else:
        print('tempType parameter is unknown value: ' + tempType)
    
    # Plot temperatures to ensure data has been read 
    # and interpolated correctly.
    if plotInputTemps:
        plt.plot(phi, t)
        plt.title('Input Cloud-Top Temperatures ' + lowTemps)
        plt.ylabel('Temperature (K)')
        plt.xlabel('Latitude (deg.)')
        if savePlots:
            plt.savefig(os.path.join(runPath, planet) + '_inputTemps.png', 
                        dpi=dpi)
        plt.show(block=False)
        plt.clf()
    
    
    
    # Read wind array and interpolate them to our phi values.
    windInput   = np.loadtxt(os.path.join(dataPath, wind), delimiter=',')
    windProfile = windDigitize(windInput, phiStart, phiStop, phiSteps)
    
    # Plot wind profile to ensure data has been read 
    # and interpolated correctly.
    if plotInputWind:
        plt.plot(phi, windProfile)
        plt.title('Input Wind Profile ' + wind)
        plt.ylabel('Wind (m/s)')
        plt.xlabel('Latitude (deg.)')
        if savePlots:
            plt.savefig(os.path.join(runPath, planet) + 
                        '_inputWindProfile.png', dpi=dpi)
        plt.show(block=False)
        plt.clf()
         
    
    
    # Create storage arrays for t, u, and rho values. 
    # These will hold our outputs. UnfilteredU stores the wind before
    # a Savitzky-Golay filter is applied.
    finalT           = np.ones((p.size, phi.size))
    finalU           = np.ones((p.size, phi.size))
    finalUnfilteredU = np.ones((p.size, phi.size))
    finalRho         = np.ones((p.size, phi.size))
    
    
    
    # Enter radius loop, over all phi. This calculates the radius on an
    # ellipsoidal planet based on flattening and the eq. radius.
    r = np.ones(phi.size)
    f = (re - rp) / re
    
    for i in np.arange(0, phi.size):
        rCalc = rEllipsoid(re, f, phi[i])
        r[i]  = rCalc
    
    
    
    # Enter gravity loop, over all phi. This calculates g at each r in our
    # calcuated r array.
    g = np.ones(phi.size)
    
    for i in np.arange(0, phi.size):
        gCalc = gAtr(M, re, r[i], phi[i], J2, omega)
        g[i]  = gCalc
    
    # Perform a polynomial fit of g to eliminate sinusodal signal.
    g = np.polyval(np.polyfit(phi, g, 6), phi)
    
    
    
    # Enter temperature loop, over all phi. This performs the basic propagation
    # of the observed, cloud-top temps down to our highest pressure using an
    # assumed constant BV frequency.
    print('Propagating cloud-top temperatures down to ' + str(pFinal) +
          ' bar with constant BV freq. ' + str(bvFreq) + '\n')
    for i in np.arange(0, phi.size):
        # Calculate t array at this phi and put in finalT output array.
        tCalc       = tFromConstBV(t[i], p, R, cp, g[i], bvFreq)
        finalT[:,i] = tCalc
        
        
        
    # Adjust the temperature array, according to user parameters. Described in
    # detail in the config file.
    intersect = int(p.size / 2) - intersectDepth
    
    for i in np.arange(0, p.size):
        smoothing   = tanhSmoothing(i, intersect, tuning)
        finalT[i,:] = (smoothing * finalT[i,:]) + ((1 - smoothing) * 
                                                   np.average(finalT[i,:]))
    
    print('Temperature smoothing is ' + 
          str(np.round(tanhSmoothing((p.size-1), intersect, tuning), 3)) + 
          ' of real value.\n')
    
    
    
    # Calculate the final output density array. 
    print('Calculating densities. \n')
    finalRho = ((np.tile(p.reshape(-1,1), finalT[0].size) * 100000) / 
                (R * finalT))
    
    
    
    # Enter wind loop, over all phi and r.
    print('Calculating wind at depth via the thermal wind equation. \n')
    for i in np.arange(0, phi.size):    
        if i == 0:
            uCalc                 = np.ones(p.size) 
            finalUnfilteredU[:,i] = uCalc * windProfile[i]
            continue
        
        # Calculate northward distance
        dy = northwardDist(re, rp, phi[i], phi[i-1])
        
        # Create array to temporarily store this phi's u values
        uCalc    = np.zeros(p.size) 
        uCalc[0] = windProfile[i]
    
        # For each pressure, calculate the vertical wind sheer and store the 
        # current pressure's thermal-wind extrapolated wind speed.
        for j in np.arange(1, p.size):
            uCalc[j] = geoTW(p[j], p[j-1], finalT[j,i], finalT[j,i-1], 
                             np.radians(phi[i]), R, omega, dy, uCalc[j-1], 
                             eqTuning) 
        
        # Store the final,  unfiltered result.
        finalUnfilteredU[:,i] = uCalc
    
    # Second wind loop. Use a Savitzky-Golay filter to smooth noise in 
    # final wind profiles at each pressure.
    if svFilter:
        for i in np.arange(0, p.size):
            finalU[i,:] = signal.savgol_filter(finalUnfilteredU[i,:], 
                                               windowLength, polyOrder)
    else:
        finalU = finalUnfilteredU
    
    
    
    # Save pressure, temperature, wind, and density results to .csv files.
    if saveData:
        print('Saving results to ' + resultsPath + '\n')
        shutil.copy('config.ini', os.path.join(runPath, 'config.ini'))
        np.savetxt(os.path.join(runPath, planet) + '_pressures.csv', 
                   p, delimiter=',')
        np.savetxt(os.path.join(runPath, planet) + '_latitudes.csv', 
                   phi, delimiter=',')
        np.savetxt(os.path.join(runPath, planet) + '_temperatures.csv', 
                   finalT, delimiter=',')
        np.savetxt(os.path.join(runPath, planet) + '_wind.csv', 
                   finalU, delimiter=',')
        np.savetxt(os.path.join(runPath, planet) + '_densities.csv', 
                   finalRho, delimiter=',')
        np.savetxt(os.path.join(runPath, planet) + '_gravity.csv', 
                   g, delimiter=',')
    
    
    
    # Plotting snippits!
    if plot2DTemps or plot2DDensity or plot2DWind or plotWindProfiles:
        print('Plotting... \n')
    
    # 2D temp plot
    if plot2DTemps:
        plt.imshow(finalT, cmap='inferno', 
                   extent=[phiStart, phiStop, pFinal, p_0])
        cb1 = plt.colorbar(fraction=0.035, pad=0.035)
        cb1.ax.tick_params(labelsize=9)
        cb1.ax.invert_yaxis()
        cb1.ax.set_title('K', fontsize=10)
        plt.contour(finalT, levels=15, colors='white', alpha=0.4, 
                    extent=[phiStart, phiStop, pFinal, p_0])
        plt.title(str(planet).capitalize() + 
                  ' Temperature from BV Freq: ' + str(bvFreq))
        plt.ylabel('Pressure (bar)')
        plt.xlabel('Latitude (deg.)')
        if savePlots:
            plt.savefig(os.path.join(runPath, planet) + '_temperature2D.png', 
                        dpi=dpi)
        plt.show(block=False)
        plt.clf()
    
    # 2D density plot
    if plot2DDensity:
        plt.imshow(finalRho, cmap='cividis', 
                   extent=[phiStart, phiStop, pFinal, p_0])
        cb2 = plt.colorbar(fraction=0.035, pad=0.035)
        cb2.ax.tick_params(labelsize=9)
        cb2.ax.invert_yaxis()
        cb2.ax.set_title('$kg/m^3$', fontsize=10)
        plt.contour(finalRho, levels=15, colors='white', alpha=0.4, 
                    extent=[phiStart, phiStop, pFinal, p_0])
        plt.title(str(planet).capitalize() + ' Density')
        plt.ylabel('Pressure (bar)')
        plt.xlabel('Latitude (deg.)')
        if savePlots:
            plt.savefig(os.path.join(runPath, planet) + '_density2D.png', 
                        dpi=dpi)
        plt.show(block=False)
        plt.clf()
    
    # 2D wind plot, smooth
    if plot2DWind:
        plt.imshow(finalU, extent=[phiStart, phiStop, pFinal, p_0])
        plt.clim(-150, 400)
        plt.title(str(planet).capitalize() + ' Wind from TWE')
        plt.ylabel('Pressure (bar)')
        plt.xlabel('Latitude (deg.)')
        cb3 = plt.colorbar(fraction=0.035, pad=0.035)
        cb3.ax.tick_params(labelsize=9)
        cb3.ax.set_title('$m/s$', fontsize=10)
        if savePlots:
            plt.savefig(os.path.join(runPath, planet) + '_wind2D.png', dpi=dpi)
        plt.show(block=False)
        plt.clf()
    
    # 2D wind plotting, indv. profiles.
    if plotWindProfiles:
        ax = plt.subplot(111)
        numSteps = int(finalU.shape[0] / plot1DStepSize) + 1
        ax.set_prop_cycle('color', plt.cm.jet(np.linspace(0, 1, numSteps)))
        for bar in np.arange(0, finalU.shape[0], plot1DStepSize):
            plt.plot(phi, finalU[bar], lw=0.75, 
                     label = str(np.round(p[bar], 2)))
        plt.plot(phi, finalU[-1], lw=0.75, label = str(np.round(p[-1], 2)))    
        plt.title(str(planet).capitalize() + 
                  ' Zonal Wind Profiles with Depth')
        plt.ylabel('Wind (m/s)')
        plt.xlabel('Latitude (deg.)')
        plt.legend(ncol=1, fontsize=10, loc='center left', 
                   title='Pressure (bar)', bbox_to_anchor=(1, 0.5))
        if savePlots:
            plt.savefig(os.path.join(runPath, planet) + '_windProfiles.png', 
                        dpi=dpi, bbox_inches='tight')
        plt.show(block=False)
        plt.clf()
        
    # 2D wind plotting, indv. profiles. No Savitzky-Golay filter.
    if plotWindProfiles and svFilter:
        ax = plt.subplot(111)
        numSteps = int(finalUnfilteredU.shape[0] / plot1DStepSize) + 1
        ax.set_prop_cycle('color', plt.cm.jet(np.linspace(0, 1, numSteps)))
        for bar in np.arange(0, finalUnfilteredU.shape[0], plot1DStepSize):
            plt.plot(phi, finalUnfilteredU[bar], lw=0.75, 
                     label = str(np.round(p[bar], 2)))
        plt.plot(phi, finalUnfilteredU[-1], lw=0.75, 
                 label = str(np.round(p[-1], 2)))
        plt.title(str(planet).capitalize() + 
                  ' Zonal Wind Profiles with Depth, no S-V Filter')
        plt.ylabel('Wind (m/s)')
        plt.xlabel('Latitude (deg.)')
        plt.legend(ncol=1, fontsize=10, loc='center left', 
                   title='Pressure (bar)', bbox_to_anchor=(1, 0.5))
        if savePlots:
            plt.savefig(os.path.join(runPath, planet) + 
                        '_windProfilesNoSVFil.png', 
                        dpi=dpi, bbox_inches='tight')
        plt.show(block=False)
        plt.clf()
    
        
    
    # End program.
    endTime = str(datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
    print('Model complete. Runtime:')
    print(datetime.datetime.strptime(endTime, '%Y-%m-%d_%H-%M-%S') - 
          datetime.datetime.strptime(currentTime, '%Y-%m-%d_%H-%M-%S'))



if __name__ == '__main__':
    if len(sys.argv) == 1:
        main()
    elif len(sys.argv) > 1:
        main(example=sys.argv[1], exampleName=sys.argv[2])