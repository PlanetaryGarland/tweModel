# -*- coding: utf-8 -*-
import numpy as np
import os
import sys
import datetime
import matplotlib.pyplot as plt
import configparser



def main(example=False, exampleName=None):
    # Begin main program.
    # Read planet and run date from the validation config.
    currentTime = str(datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
    valConfig   = configparser.ConfigParser()
    if not(example):  
        valConfig.read(os.path.join(os.path.dirname(os.path.abspath(__file__)), 
                                 'validationConfig.ini'))
    else:
        valConfig.read(os.path.join(os.path.dirname(os.path.abspath(__file__)), 
                             'examples', exampleName, 'validationConfig.ini'))
    
    try:
        planet = valConfig['USER']['planet'] 
        date   = valConfig['USER']['date']
        
        # PLOTTING PARAMETERS
        plot2DBV         = valConfig['USER'].getboolean('plot2DBV')
        plot2DPotTemp    = valConfig['USER'].getboolean('plot2DPotTemp')
        plot2DPotTempDel = valConfig['USER'].getboolean('plot2DPotTempDel')
        dpi              = int(valConfig['USER']['dpi'])
        savePlots        = valConfig['USER'].getboolean('savePlots')
    
        # RESULTS PARAMETERS
        saveData = valConfig['USER'].getboolean('saveData')
    except KeyError:
        print('validationConfig.ini not found! ' +
              'Use and rename the valConfigTemplate.ini in ' + 
              '/docs and move to the same directory as validation.py')
    
    
    
    # Set saved results directory.
    if not(example):
        resultsPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), 
                                   'results', str(planet), str(date))
    else:
        resultsPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), 
                                   'examples', exampleName, 'results', 
                                   str(planet), 'runResults')
                                   
    print('Loading data from directory: ' + resultsPath + '\n')
    
    # Set user-defined parameters from saved config. 
    # Parameters defined in the config file.
    config = configparser.ConfigParser()
    config.read(os.path.join(resultsPath, 'config.ini'))
    
    try:
        # PLANET PARAMETERS
        R         = float(config['USER']['R'])
        cp        = float(config['USER']['cp'])
        
        # MODEL PARAMETERS
        p_0       = float(config['USER']['p_0']) 
        pFinal    = float(config['USER']['pFinal']) 
        phiStart  = float(config['USER']['phiStart']) 
        phiStop   = float(config['USER']['phiStop'])
    except KeyError:
        print('Saved config.ini not found! Are you sure the ' +
              'run date string is correct?')
    
    
    
    # Load saved temperatures, pressures, and gravities.
    temperatures = np.genfromtxt(os.path.join(resultsPath, planet) +
                                              '_temperatures.csv', 
                                              delimiter=',')
    pressures    = np.genfromtxt(os.path.join(resultsPath, planet) + 
                                              '_pressures.csv', 
                                              delimiter=',')
    gravities    = np.genfromtxt(os.path.join(resultsPath, planet) +
                                              '_gravity.csv', 
                                              delimiter=',')
    
    
    
    # Begin BV validation.
    print('Starting BV validation.')
    
    # Calculate kappa, which is the same everywhere.
    kappa = (R / cp)
    
    # Initialize N2 array for storing BV freq. results
    BVSquared = np.ones(temperatures.shape)
    
    # Enter main Brunt Vaisala validation loops, over latitude and pressure.
    count = 0
    # Latitude loop
    for i in np.arange(0, temperatures.shape[1]):
        # Set current gravity at r.
        currentGrav = gravities[i]
        
        # Pressure loop
        for j in np.arange(0, temperatures.shape[0]-1):
            # Set current and next pressure and temperature.
            currentPress = pressures[j]
            nextPress    = pressures[j + 1]
            currentTemp  = temperatures[j, i]
            nextTemp     = temperatures[j + 1, i]
            
            # Calculate scale height, H, and layer height, dz.
            H  = ((R * currentTemp) / currentGrav)
            dz = (H * (np.log((nextPress / currentPress))))
        
            # Calculate the Brunt-Vaisala frequency squared, a check on the
            # stability of the atmosphere. Should be positive everywhere, 
            # ensuring the BV frequency is real.
            term1           = (R / H)
            term2           = ((currentTemp - nextTemp) / dz)
            term3           = ((kappa * currentTemp) / H)
            N2              = term1 * (term2 + term3)
            BVSquared[j, i] = N2
    
            # If N2 is negative, warn of an unstable value after validation.
            if N2 < 0:
                count +=1
    
    if count > 0:           
        print('VALIDATE FAILED! Brunt Vaisala freq. is imaginary at: ' +
              str(count) + ' points in the lat. / pressure array! \n')
    else:
        print('VALIDATE SUCCESS! Brunt Vaisala freq. is real at all points ' +
              'in the lat. / pressure array. \n')
    
    
    
    # Begin potential temperature validation.
    print('Starting potential temperature validation.')
    
    # Initialize pot. temp array for storing results.
    potTemp = temperatures.copy()
    
    # Pressure loop. Do this pythonicly with numpy to avoid two loops.
    for i in np.arange(0, temperatures.shape[0]):
        # Calculate the potential temperature at all points in our temperature
        # array. Potential temperature should increase with height or 
        # our model atmosphere is unstable.
        potTemp[i,:] = (potTemp[i,:] * ((pressures[0] / pressures[i])**(R/cp)))
        
    # Initialize pot. temp delta array for storing results.
    potTempDel = potTemp.copy()
    
    # Calculate d potTemp / d P to see if it's positive or negative.
    # It should be negative at all points, indicating that potential temp.  
    # increases with height.
    for i in np.arange(0, temperatures.shape[0]-1):
        potTempDel[i,:] = (np.subtract(potTempDel[i+1,:], potTempDel[i,:]) / 
                           np.subtract(pressures[i+1], pressures[i]))
    
    if np.where(potTempDel[0:-1,:] > 0)[0].size > 0:
        print('VALIDATE FAILED! Pot. temperature increases with height at: ' +
              str(np.where(potTempDel[0:-1,:] > 0)[0].size) + 
              ' points in the lat. / pressure array! \n')
    else:
        print('VALIDATE SUCCESS! Pot. temperature decreases with height ' + 
              'everywhere in the lat. / pressure array! \n')
    
    
    
    # Save validation results.
    if saveData:
        print('Saving results to ' + resultsPath + '\n')
        np.savetxt(os.path.join(resultsPath, planet) + 
                   '_BVSquared_validation.csv', BVSquared, delimiter=',')
        np.savetxt(os.path.join(resultsPath, planet) + 
                   '_potTemp_validation.csv', potTemp, delimiter=',')
        np.savetxt(os.path.join(resultsPath, planet) + 
                   '_delPotTemp_validation.csv', potTempDel, delimiter=',')
    
    
    
    # Plotting snippets!
    if plot2DBV or plot2DPotTemp or plot2DPotTempDel:
        print('Plotting... \n')
    
    # BV plot, Just to display the variation and that we're not negative
    if plot2DBV:
        BV    = np.sqrt(BVSquared[0:-1,:])
        bvMin = np.abs(np.amin(BV))
        bvMax = np.abs(np.amax(BV))
        
        if bvMax >= bvMin:
            bvClip = bvMax * 1.1
        else:
            bvClip = bvMin * 1.1
    
        plt.imshow(BV, cmap='RdBu_r', clim=((-1 * bvClip), bvClip), 
                   extent=[phiStart, phiStop, pFinal, p_0])
        cb1 = plt.colorbar(fraction=0.035, pad=0.035)
        cb1.ax.tick_params(labelsize=9)
        cb1.ax.set_title('$s^{-1}$', fontsize=10)
        plt.title(str(planet).capitalize() + ' Brunt-Vaisala Frequency')
        plt.ylabel('Pressure (bar)')
        plt.xlabel('Latitude (deg.)')
        if savePlots:
            plt.savefig(os.path.join(resultsPath, planet) + 
                        '_BVSquared_validation.png', dpi=dpi)
        plt.show()
        plt.clf()
    
    # 2D pot. temp. plot
    if plot2DPotTemp:
        plt.imshow(potTemp, cmap='inferno', 
                   extent=[phiStart, phiStop, pFinal, p_0])
        cb2 = plt.colorbar(fraction=0.035, pad=0.035)
        cb2.ax.tick_params(labelsize=9)
        cb2.ax.set_title('K', fontsize=10)
        plt.contour(potTemp, levels=15, colors='white', alpha=0.4, 
                    extent=[phiStart, phiStop, pFinal, p_0])
        plt.title(str(planet).capitalize() + ' 2D Potential Temperature')
        plt.ylabel('Pressure (bar)')
        plt.xlabel('Latitude (deg.)')
        if savePlots:
            plt.savefig(os.path.join(resultsPath, planet) + 
                        '_potTemp_validation.png', dpi=dpi)
        plt.show()
        plt.clf()
    
    # Pot. temp. delta.
    if plot2DPotTempDel:
        ptdMin = np.abs(np.amin(potTempDel[0:-1,:]))
        ptdMax = np.abs(np.amax(potTempDel[0:-1,:]))
        if ptdMax >= ptdMin:
            ptdClip = ptdMax * 1.1
        else:
            ptdClip = ptdMin * 1.1
            
        plt.imshow(potTempDel[0:-1,:], clim=((-1 * ptdClip), ptdClip), 
                   cmap='RdBu_r', extent=[phiStart, phiStop, pFinal, p_0])
        cb3 = plt.colorbar(fraction=0.035, pad=0.035)
        cb3.ax.tick_params(labelsize=9)
        cb3.ax.set_title('$K/bar$', fontsize=10)
        plt.title(str(planet).capitalize() + ' $d\\theta / dp$')
        plt.ylabel('Pressure (bar)')
        plt.xlabel('Latitude (deg.)')
        if savePlots:
            plt.savefig(os.path.join(resultsPath, planet) + 
                        '_delPotTemp_validation.png', dpi=dpi)
        plt.show()
        plt.clf()
    
    
    
    # End program.
    endTime = str(datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
    print('Validation complete. Runtime:')
    print(datetime.datetime.strptime(endTime, '%Y-%m-%d_%H-%M-%S') - 
          datetime.datetime.strptime(currentTime, '%Y-%m-%d_%H-%M-%S'))



if __name__ == '__main__':
    if len(sys.argv) == 1:
        main()
    elif len(sys.argv) > 1:
        main(example=sys.argv[1], exampleName=sys.argv[2])