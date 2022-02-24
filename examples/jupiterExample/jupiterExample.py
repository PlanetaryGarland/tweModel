import os
import sys
sys.path.insert(1, os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
import numpy as np
import datetime
import tweModel
import validation



def checkRun(examplePath, runPath, planet, parameter):
    '''
    Checks an output .csv file against a stored, expected output.
    
    Parameters
    ----------
    examplePath : string
        Path to the directory where the expected output is stored.
    runPath : string
        Path to the directory where the run output is stored.
    planet : string
        Name of the planet of interest.
    parameter : string
        Name of parameter to be checked.
        
    Returns
    -------
    int
        Value indicating success (0), values withing rounding error (1) or
        failure (2).
    
    Notes
    -----
    The example and run files will be searched for with the following directory
    and filename structure:
        
    <examplePath>/<planet>_<parameter>.csv
    
    <runPath>/<planet>_<parameter>.csv
    
    History
    -------
    2/22/2022, planetary@garland.run : v1.0
        Initial published version. 
    '''
    example = np.genfromtxt(os.path.join(examplePath, planet) + 
                            ('_' + str(parameter) + '.csv'), delimiter=',')
    run     = np.genfromtxt(os.path.join(runPath, planet) + 
                            ('_' + str(parameter) + '.csv'), delimiter=',')
    
    check   = np.array_equal(example, run)
    close   = np.isclose(example, run)
    
    if check:
        print('INSTALL SUCCESSFUL: ' + str(parameter).capitalize() + 
              ' results match exactly.\n')
        return 0
    elif close.all():
        print('INSTALL SUCCESSFUL: ' + str(parameter).capitalize() + 
              ' results are extremely close, but do ' +
              'not match exactly.\n')
        return 1
    else:
        print('INSTALL FAILURE! ' + str(parameter).capitalize() + 
              ' results do NOT match.\n')
        return 2
        
        
        
# Begin example.
currentTime = str(datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
planet      = 'jupiter'
exampleName = 'jupiterExample'
print('Beginning example: ' + exampleName)
print('Checking installed version of tweModel.py and validation.py against ' + 
      'expected values.\n')



# Running tweModel.py in example mode.
tweModel.main(example=True, exampleName=exampleName)

# Running validation.py in example mode.
validation.main(example=True, exampleName=exampleName)
print()



# Set paths to the expected results and those from our run.
examplePath = os.path.join(os.path.dirname(os.path.abspath(__file__)), 
                           'results', 'exampleResults', exampleName)
runPath     = os.path.join(os.path.dirname(os.path.abspath(__file__)), 
                           'results', planet, 'runResults')



# Check run and validation output vs. stored expected results.
exitCode = 0
for param in ['pressures', 'temperatures', 'gravity', 'densities', 'wind', 
              'BVSquared_validation', 'potTemp_validation', 
              'delPotTemp_validation']:
    code = checkRun(examplePath, runPath, planet, param)

    if code == 1:
        exitCode = 1
    if code == 2:
        exitCode = 2
        
    

# End example and print a final pass/fail.
endTime = str(datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))

if exitCode == 0:
    print('INSTALL FULLY SUCCESSFUL. Runtime:')
if exitCode == 1:
    print('INSTALL SUCCESSFUL. At least 1 parameter was ' +
          'different but within machine error. Runtime:')
if exitCode == 2:
    print('INSTALL FAILURE! At least 1 parameter did not match ' +
          'the expected result. Runtime:')

print(datetime.datetime.strptime(endTime, '%Y-%m-%d_%H-%M-%S') - 
      datetime.datetime.strptime(currentTime, '%Y-%m-%d_%H-%M-%S'))