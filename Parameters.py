import numpy as np

class parameters:
    def __init__(self):
        self.kAtm = {'u': 0.50,
                     'g': 0.21,
                     'r': 0.13,
                     'i': 0.10,
                     'z': 0.07,
                     'y': 0.18} 
        
        
        self.msky = {'u': 22.95,
                     'g': 22.24,
                     'r': 21.20,
                     'i': 20.47,
                     'z': 19.60,
                     'y': 18.63} 
        
        self.Cm = {'u':22.94,
                   'g':24.46,
                   'r':24.48,
                   'i':24.34,
                   'z':24.18,
                   'y':23.73}
        
        self.dCm_infinity = {'u':0.56,
                             'g':0.12,
                             'r':0.06,
                             'i':0.05,
                             'z':0.03,
                             'y':0.02}
        
        
       #FWHM_500 = seeing at 500 nm
       # FWHM_Sys_Zenith = sqrt(telSeeing**2 + opticalDesSeeing**2 + cameraSeeing**2)
       # Filter_Wavelength_Correction = (500 nm / Filter_Effective_Wavelength)**0.3
       # Airmass_Correction = airmass**0.6
       # FWHM_Sys = FWHM_Sys_Zenith * Airmass_Correction
       # FWHM_Atm = FWHM_500 * Filter_Wavelength_Correction * Airmass_Correction
       # FWHM_Eff = scaleToNeff * sqrt(FWHM_Sys**2 + atmNeffFactor * FWHM_Atm**2)
       # FWHM_Eff is the value in ObsHistory.finSeeing for the observations filter
       #
       # Units = unitless, Format = float, no default
       #
        
        self.telSeeing = 0.250 # design goal
        self.opticalDesSeeing = 0.08
        self.cameraSeeing = 0.300
       # Scaling factors for above seeing calculation
        self.scaleToNeff = 1.16
        self.atmNeffFactor = 1.04
        self.FWHM_Sys_Zenith = np.sqrt(self.telSeeing**2 + self.opticalDesSeeing**2 + self.cameraSeeing**2)
        
        self.filterWave = {'u': 367.0, 'g': 482.5, 'r': 622.2, 'i': 754.5, 'z': 869.1, 'y': 971.0}
        
