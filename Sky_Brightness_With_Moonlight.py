import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
import astropy.coordinates as coord
import astropy.units as u
import math
import healpy as hp
from astropy.table import Table
from optparse import OptionParser
import palpy as pal
from mpl_toolkits.mplot3d import Axes3D
from lsst.sims.photUtils import Bandpass
from lsst.sims.photUtils import Sed
#from Throughputs import Throughputs
#from Parameters import parameters
from lsst.sims.photUtils import SignalToNoise
from lsst.sims.photUtils import PhotometricParameters

DEG2RAD = math.pi / 180.    # radians = degrees * DEG2RAD
RAD2DEG = 180. / math.pi    # degrees = radians * RAD2DEG
TWOPI = 2 * math.pi
DAY = 86400.

class SkyBrightness:
    def __init__(self,ra,dec,mjd,date,filtre,airmass_opsim,extinction=0.172, skyBrightness=21.587):
        
        self.ra_RAD=ra
        self.dec_RAD=dec
        self.mjd=mjd
        self.k=extinction
        self.skyBrightness=skyBrightness
        self.date=date
        self.simEpoch=59580
        #self.simEpoch=53371
        #if date<0.:
        ddate=self.date
        mjd_new=float(ddate)/float(DAY)+float(self.simEpoch)
        self.date= int((self.mjd-self.simEpoch)*float(DAY))
        self.filter=filtre
        self.airmass_opsim=airmass_opsim
        self.lon_RAD = -70.7494 *DEG2RAD #obs site long (from sims_operations/conf/system/SiteCP.conf) 
        self.lat_RAD = -30.2444* DEG2RAD #obs site lat (from sims_operations/conf/system/SiteCP.conf)
        self.lst_RAD  = self.Get_Sidereal_Time_at_Site() #local sidereal time at site (radians)
        # TWILIGHT LIMITS
        # Altitude of the Sun in degrees that define the twilight.
        # When the sun is above this limit and below the night limit, a special
        # twilight factor is included in the sky brightness model
        self.SunAltitudeTwilightLimit = -18.0
        self.SunAltitudeNightLimit = -12.0 
        self.n1dp92104 = (1. / 0.92104)
        self.log34p08 = math.log(34.08)
        
        self.twilightBrightness = 17.3
        self.fluxTwilight = 10. ** ((-self.twilightBrightness) / 2.5)

        self.MoonProfile()
        self.Moon_Distance()
        alpha=self.Get_alpha()
        self.distance2moon_RAD = pal.dsep(self.moonRA_RAD, self.moonDec_RAD, self.ra_RAD, self.dec_RAD)
        self.distance2moon_DEG = self.distance2moon_RAD * RAD2DEG
        lha_RAD = self.lst_RAD - self.ra_RAD
        (az_RAD, d1, d2, alt_RAD, d4, d5, pa_RAD, d7, d8) = pal.altaz(lha_RAD, self.dec_RAD,self.lat_RAD)
        # Altitude -> Zenith distance (degrees)
        targetZD_DEG = 90. - (alt_RAD* RAD2DEG)
        # Altitude -> Zenith distance (radian)
        zd_RAD = 1.5707963 - alt_RAD
        # Airmass
        #am = slalib.sla_airmas (zd_RAD)
        am = pal.airmas(zd_RAD)
        
        
        # Compute Raileigh scattering
        rs = (10. ** 5.36) * (1.06 + (math.cos(self.distance2moon_RAD)) ** 2.)
        #rs = 2.27E5 * (1.06 + (math.cos(self.distance2moon_RAD)) ** 2.)

        # Compute Mie scattering
        #print 'distance',self.distance2moon_DEG
        if self.distance2moon_DEG > 10.:
            ms = 10. ** (6.15 - (self.distance2moon_DEG / 40.))
        else:
            ms = 6.2E7 / (self.distance2moon_DEG ** 2.)

        # Compute illumninace of the Moon
        i = 10. ** (-0.4 * (3.84 + 0.026 * abs(alpha) + 4.E-9 * (alpha) ** 4.))

        # Compute optical pathlength (in units of airmass) as a
        # function of the Zenith Distance (ZD)
        x = lambda z: (1. - 0.96 * (math.sin(z * DEG2RAD)) ** 2.) ** -0.5
        #xa = lambda z: pal.airmas(z* DEG2RAD)

        #print 'airmass here',x(self.moonZD_DEG),xa(self.moonZD_DEG),x(targetZD_DEG),xa(targetZD_DEG),self.distance2moon_RAD,self.moonRA_RAD, self.moonDec_RAD,lha_RAD,zd_RAD,alt_RAD,self.lst_RAD,targetZD_DEG,self.mjd

         # Put all together (nanoLamberts)!

        moonBr = (rs + ms) * i * 10. ** (-0.4 * self.k * x(self.moonZD_DEG))
        moonBr *= (1. - 10. ** (-0.4 * self.k * x(targetZD_DEG)))

        self.moon_airmass=x(self.moonZD_DEG)
        self.targetZD_DEG=targetZD_DEG
        self.target_airmass=x(targetZD_DEG)
        #print 'airmasses',self.target_airmass,self.airmass_opsim,pal.airmas(zd_RAD)

        #  Taper brightness to 0, at moon altitude = self.sunAltitudeNightLimit
        if self.moonAlt_DEG < self.SunAltitudeNightLimit:
            moonBr = 0.
        elif self.moonAlt_DEG < 0. and self.moonAlt_DEG >= self.SunAltitudeNightLimit:
            moonBr *= (1 - self.moonAlt_DEG / self.SunAltitudeNightLimit)
        

        # Now, compute the dark sky brightness (in nanoLamberts) at
        # the Zenith from the same value in mag/arcsec^s
        skyBrightness = 34.08 * math.exp(20.7233 - 0.92104 * skyBrightness)

        # ... and at our position in the sky (targetZD_DEG)
        skyBr = skyBrightness * x(targetZD_DEG) * 10. ** (-0.4 * self.k * (x(targetZD_DEG) - 1.))

        # Add the brightness of the moon to that of the dark sky (nanoLamberts)
        totBr = moonBr + skyBr

        # Transform it into mag/arcsec^2
        #totBr = (1. / 0.92104) * math.log (34.08 * math.exp (20.7233) / totBr)
        totBr = self.n1dp92104 * (self.log34p08 + 20.7233 - math.log(totBr))

        #print 'totbr',totBr

        #Get the Twilight profile

        self.TwilightProfile()
      
        #if self.date < self.sunSetTwil or self.date > self.sunRiseTwil:
        (yy, mm, dd,hh,min,sec) = self.mjd2gre(mjd)[:6]
        if (hh < 20 and self.date > self.sunRiseTwil) or (hh>20  and self.date < self.sunSetTwil):
            totBr = -2.5 * math.log10(10. ** ((-totBr) / 2.5) + self.fluxTwilight)
            
        #print 'after twilight',totBr,self.date,self.sunSetTwil,self.sunRiseTwil

        #print 'filterskybrightness',self.FilterSkyBrightness(totBr)
        
        self.skybrightness_moon=self.FilterSkyBrightness(totBr)
        #self.skybrightness_moon=totBr

    def new_skybrightness(self):

        return self.skybrightness_moon

    def Get_Sidereal_Time_at_Site(self):

        lst_RAD = pal.gmst(self.mjd) + self.lon_RAD
        if lst_RAD < 0:
            lst_RAD += TWOPI

        return lst_RAD

    def MoonProfile(self):
        (self.moonRA_RAD, self.moonDec_RAD, self.moonDiam) = pal.rdplan(self.mjd, 3, self.lon_RAD, self.lat_RAD)
        #print 'hello moon',self.moonRA_RAD, self.moonDec_RAD, self.moonDiam
        
        #getting the moon phase

        v6moon = pal.dmoon(self.mjd)
        
        Rmatrix = pal.prenut(2000.0, self.mjd)
        
        xyzMoon2000 = pal.dmxv(Rmatrix, v6moon)
       
        (moonra, moondec) = pal.dcc2s(xyzMoon2000)
        moonra = pal.dranrm(2.0 * math.pi + moonra)
        sun12 = pal.evp(self.mjd, 2000.0)
        sun3heliocentric = sun12[3]
        for i in range(3):
            sun3heliocentric[i] *= -1
        (sunra, sundec) = pal.dcc2s(sun3heliocentric)
        sunra = pal.dranrm(2.0 * math.pi + sunra)
        moonsunsep = pal.dsep(sunra, sundec, moonra, moondec)
        self.moonPhase = ((1.0 - math.cos(moonsunsep)) / 2.0) * 100

        #print 'and the phase',self.moonPhase

    def Moon_Distance(self):
        # Compute moon altitude in radians
        self.moonha_RAD = self.lst_RAD - self.moonRA_RAD

        (moonAz_RAD, d1, d2, self.moonAlt_RAD, d4, d5, d6, d7, d8) = \
            pal.altaz(self.moonha_RAD, self.moonDec_RAD, self.lat_RAD)
        self.moonAlt_DEG = self.moonAlt_RAD * RAD2DEG
        
        # compute moon's zenith distance
        if self.moonAlt_RAD < 0.:
            moonZD_RAD = 1.5707963
        else:
            moonZD_RAD = 1.5707963 - self.moonAlt_RAD
            
        self.moonZD_DEG = moonZD_RAD * RAD2DEG

    def Get_alpha(self):
        # Compute the phase angle given the illuminated fraction
        # formula from http://astro.ft.uam.es/TJM/tjm/webpaginas/practicas/ephemeris/moon.js
        alpha = math.acos(2. * self.moonPhase  / 100. - 1.) * RAD2DEG
        alpha = self.normalize(alpha, min=0., max=180., degrees=True)
        
        return alpha

    def normalize(self,angle, min=0., max=None, degrees=True):
        """
        Given an angle, make sure that its value is within the range
        [min, max].
        
        If degrees=True the ammount to add/subtract is 360.
        If degrees=False the ammount to add/subtract is 2*pi
        
        If min/max is not specified, it reverts to its default value
        (which depends on the value of degrees):
        min = 0
        max = 360./2*pi
        
        Return the normalized angle in its original units
        """
        angle = float(angle)
        
        if min is None:
            min = 0.
        if degrees:
            if max is None:
                max = 360.
            addit = 360.
        else:
            if max is None:
                max = TWOPI
            addit = TWOPI

        while angle < min:
            angle += addit
        while angle > max:
            angle -= addit
        return angle


    def TwilightProfile(self):
         (sunRise, sunSet, sunRiseMJD, sunSetMJD, self.sunRiseTwil,
            self.sunSetTwil) = self.getTwilightSunriseSunset(self.date)

         #print 'TwilightProfile','sunRiseMJD=',sunRiseMJD,'sunSetMJD=',sunSetMJD,'sunRiseTwil=',self.sunRiseTwil,'sunSetTwil=',self.sunSetTwil

    def getTwilightSunriseSunset(self, date):
        """
        Compute the time of the astronomical twilight for Sun raise
        and set times at the location on Earth specified by
        (self.latitude_RAD, self.longitude_RAD).

        Input
        date:   date in seconds from Jan 1 of the simulated year.

        Return
        A four-element array:
        (sunriseTime, sunsetTime, sunriseMJD, sunsetMJD)
        sunriseTime & sunsetTime are in seconds from Jan 1 of simulated year;
        sunriseMJD and sunsetMJD are MJD equivalents.
        """
        # Convert date in MJD
        mjd = (float(self.date) / float(DAY)) + self.simEpoch

        #mjd=49353.
        # MJD -> calendar date
        #mjd=self.mjd
        (yy, mm, dd,hh,min,sec) = self.mjd2gre(mjd)[:6]

        #print 'hello',yy,mm,dd,hh,min,sec
        # Compute sunset and sunrise at twilight. Return values are in
        # decimal hours.
        import Sun
        s = Sun.Sun()
        #(sunRise, sunSet) = s.sunRiseSet(yy, mm, dd,
        #                        self.longitude_RAD * RAD2DEG, self.latitude_RAD * RAD2DEG)
        # Following set to nauticalTwilight
        (sunRise, sunSet) = s.__sunriset__(yy, mm, dd, self.lon_RAD * RAD2DEG,
                                           self.lat_RAD * RAD2DEG, self.SunAltitudeNightLimit, 0)
        # -12.0, 0)
        (sunRiseTwil, sunSetTwil) = s.__sunriset__(yy, mm, dd, self.lon_RAD * RAD2DEG,
                                                   self.lat_RAD * RAD2DEG, self.SunAltitudeTwilightLimit,
                                                   0)
        
        # -18.0, 0)

        # Compute MJD values for sunrise and sunset
        sunSetMJD = int(mjd) + (sunSet / 24.)
        sunRiseMJD = int(mjd) + (sunRise / 24.)
        self.sunSetTwilMJD = int(mjd) + (sunSetTwil / 24.)
        self.sunRiseTwilMJD = int(mjd) + (sunRiseTwil / 24.)
        #print 'ici', sunSetTwilMJD,sunRiseTwilMJD
        # MJD -> simulated seconds
        sunsetDate = (sunSetMJD - self.simEpoch) * float(DAY)
        sunriseDate = (sunRiseMJD - self.simEpoch) * float(DAY)
        sunsetTwilDate = (self.sunSetTwilMJD - self.simEpoch) * float(DAY)
        sunriseTwilDate = (self.sunRiseTwilMJD - self.simEpoch) * float(DAY)
        
        #print 'twilight',self.mjd2gre(sunSetTwilMJD)[:6],self.mjd2gre(sunRiseTwilMJD)[:6]
        return (sunriseDate, sunsetDate, sunRiseMJD, sunSetMJD, sunriseTwilDate, sunsetTwilDate)

    def mjd2gre(self,mjd):
        # Use SLALIB to convert from MJD to (year, month, day)
    #(year, month, day, fraction, error) = slalib.sla_djcl (mjd)
        (year, month, day, fraction) = pal.djcl(mjd)
        
    #if (error):
    #    msg = 'Fatal error: MJD must correspond to a date later than 4701BC March 1'
    #    raise (SyntaxError, msg)
        
    # Now take care of the fraction of day
        hh = math.floor(24. * fraction)
        mm = math.floor(1440. * (fraction - hh / 24.))
        ss = 86400. * (fraction - hh / 24. - mm / 1440.)
        return (year, month, day, int(hh), int(mm), ss)

    def Corrections_Filters_for_moonPhase(self):
        
        self.skyBrightKeys = [0, 18, 50, 80, 100]
        self.filterOffset = {}
        
        # Corrections for moonPhase = 0 percent (new moon)
        self.filterOffset['u', 0.] = 0.66
        self.filterOffset['g', 0.] = 0.41
        self.filterOffset['r', 0.] = -0.28
        self.filterOffset['i', 0.] = -1.36
        self.filterOffset['z', 0.] = -2.15
        
        # Corrections for moonPhase = 18 percent
        self.filterOffset['u', 18.] = 0.28
        self.filterOffset['g', 18.] = 0.30
        self.filterOffset['r', 18.] = -0.19
        self.filterOffset['i', 18.] = -1.17
        self.filterOffset['z', 18.] = -1.99
        
        # Corrections for moonPhase = 50 percent
        self.filterOffset['u', 50.] = -1.05
        self.filterOffset['g', 50.] = 0.03
        self.filterOffset['r', 50.] = 0.02
        self.filterOffset['i', 50.] = -0.96
        self.filterOffset['z', 50.] = -1.78
        
        # Corrections for moonPhase = 80 percent
        self.filterOffset['u', 80.] = -1.83
        self.filterOffset['g', 80.] = -0.08
        self.filterOffset['r', 80.] = 0.10
        self.filterOffset['i', 80.] = -0.78
        self.filterOffset['z', 80.] = -1.54
        
        # Corrections for moonPhase = 100 percent (full moon)
        self.filterOffset['u', 100.] = -2.50
        self.filterOffset['g', 100.] = -0.35
        self.filterOffset['r', 100.] = 0.31
        self.filterOffset['i', 100.] = -0.47
        self.filterOffset['z', 100.] = -1.16


    def FilterSkyBrightness(self,skyBrightness):
        
        self.Corrections_Filters_for_moonPhase()
        moonPhase_PERCENT=self.moonPhase 

        if self.filter == 'y':
            filterSkyBright = 17.3
        else:      # g,r,i,z,u
            # If moon below horizon, use new moon offset for filter
            # brightness - MM
            if math.degrees(self.moonAlt_RAD) <= -6.0:
                adjustBright = self.filterOffset[self.filter, 0.]

            # Interpolate if needed. Note: moonPhase is a float not int
            elif moonPhase_PERCENT not in self.skyBrightKeys:
                i = 0
                while self.skyBrightKeys[i] < moonPhase_PERCENT:
                    i += 1

                # find upper and lower bound
                upperMoonPhase = self.skyBrightKeys[i]
                lowerMoonPhase = self.skyBrightKeys[i - 1]
                lowerAdjustBright = self.filterOffset[self.filter, lowerMoonPhase]
                upperAdjustBright = self.filterOffset[self.filter, upperMoonPhase]
                # linear interpolation
                diffMoonPhase = upperMoonPhase - lowerMoonPhase
                diffAdjustBright = upperAdjustBright - lowerAdjustBright
                moonPhaseAdj = moonPhase_PERCENT - lowerMoonPhase
                adjustBright = lowerAdjustBright + (((moonPhaseAdj) * (diffAdjustBright)) / (diffMoonPhase))

            else:          # moon not set and moon phase is key
                adjustBright = self.filterOffset[self.filter, moonPhase_PERCENT]
            #print 'adjust',adjustBright
            filterSkyBright = skyBrightness + adjustBright

            # z sky brightness should never be under 17.0
            if self.filter == 'z' and filterSkyBright < 17.0:
                filterSkyBright = 17.0

        # If twilight, set brightness for z and y
        (yy, mm, dd,hh,min,sec) = self.mjd2gre(self.mjd)[:6]
        if (hh < 20 and self.date > self.sunRiseTwil) or (hh>20  and self.date < self.sunSetTwil):
            if self.filter == 'z' or self.filter == 'y':
                filterSkyBright = 17.0

        return filterSkyBright

