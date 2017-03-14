import numpy as np
import matplotlib.pyplot as plt
from lsst.sims.maf.metrics import BaseMetric
import pickle as pkl
from Sky_Brightness_With_Moonlight import SkyBrightness
import palpy as pal
import math

DEG2RAD = math.pi / 180.    # radians = degrees * DEG2RAD
RAD2DEG = 180. / math.pi    # degrees = radians * RAD2DEG
TWOPI = 2 * math.pi
DAY = 86400.
simEpoch=59580

class Fake_Rolling(BaseMetric):
    """
    Measure how many time series meet a given time and filter distribution requirement.Fake 
    """
    
    def __init__(self, metricName='Fake_Rolling', m5Col='fiveSigmaDepth', units='', badval=-666,uniqueBlocks=False,fieldname='DD',fieldID=[290],merge_factor=0.8,**kwargs):
        self.mjdCol ='expMJD'
        self.filterCol='filter'
        self.m5Col = m5Col
        self.fieldID='fieldID'
        self.dateCol = 'expDate'
        self.fieldRA='fieldRA'
        self.fieldDec='fieldDec'
        self.fieldID='fieldID'
        self.ditheredRA='ditheredRA'
        self.ditheredDec='ditheredDec'
        self.visitTime='visitExpTime'
        self.finSeeing='finSeeing'
        self.rawSeeing='rawSeeing'
        self.moonPhase='moonPhase'
        self.airmass='airmass'
        self.filtSkyBrightness='filtSkyBrightness'
        
        self.Fields={}
        self.Fields_Seasons={}

        self.fieldName=fieldname
        self.fieldID_ref=fieldID
        self.merge_factor=merge_factor

        self.bands=['u','g','r','i','z','y']

        super(Fake_Rolling, self).__init__(col=[self.mjdCol,self.m5Col, self.filterCol, self.dateCol,self.fieldRA,self.fieldDec, self.ditheredRA,self.ditheredDec,self.visitTime,self.rawSeeing,self.moonPhase,self.airmass,self.filtSkyBrightness,self.fieldID],metricName=metricName, units=units, badval=badval, **kwargs)

        self.uniqueBlocks = uniqueBlocks
        

    def run(self, dataSlice, slicePoint=None):

        if len(self.fieldID_ref)==1:
            self.Fake_Single(dataSlice)
        else:
            self.Fake_Multiple(dataSlice)

    def Fake_Single(self,dataSlice):

        dataSlice = dataSlice[np.where(dataSlice[self.fieldID]==self.fieldID_ref)]
        if len(dataSlice) > 0:
            dataSlice.sort(order=self.mjdCol)
            
            seasons=self.Get_Seasons(dataSlice)

            rolling_cadence={}
            
            for i in range(len(seasons)):
                rolling_cadence[i]=np.zeros((0,1),dtype=seasons[0].dtype)
                
            #Season 0 : keep the same
            season_keep,season_remain=self.Split_Season(seasons[0])
            print 'allo',len(season_keep),len(season_remain),len(rolling_cadence[0]),season_keep.shape,rolling_cadence[0].shape
            rolling_cadence[0]=np.concatenate((rolling_cadence[0],season_keep))
             
            # Season i in [1,4,7]:
            for iseason in [1,4,7]:

                shift_one=self.Shift(seasons[iseason],seasons[iseason+1])
                shift_two=self.Shift(seasons[iseason],seasons[iseason+2])
                
                for band in self.bands:
                    season_ref=seasons[iseason][np.where(seasons[iseason][self.filterCol]==band)]
                    season_ref_plus_one=seasons[iseason+1][np.where(seasons[iseason+1][self.filterCol]==band)]
                    season_ref_plus_two=seasons[iseason+2][np.where(seasons[iseason+2][self.filterCol]==band)]         
                    season_keep,season_remain=self.Split_Season(season_ref)
                    # keep season i
                    rolling_cadence[iseason]=np.concatenate((rolling_cadence[iseason],season_keep))
                    # from season i+1: keep 20% in i+1, put 80% in i (MJD[i] = MJD[i+1]-365.)
                    season_keep,season_remain=self.Split_Season(season_ref_plus_one,self.merge_factor)
                   
                    rolling_cadence[iseason]=np.concatenate((rolling_cadence[iseason],self.Modify(season_keep,shift=shift_one)))
                   
                    rolling_cadence[iseason+1]=np.concatenate((rolling_cadence[iseason+1],season_remain))
                    #from season i+2: keep 20% in i+2, put 80% in i (MJD[i] = MJD[i+2]-730.)
                    season_keep,season_remain=self.Split_Season(season_ref_plus_two,self.merge_factor)
                    
                    rolling_cadence[iseason]=np.concatenate((rolling_cadence[iseason],self.Modify(season_keep,shift=shift_two)))
                    rolling_cadence[iseason+2]=np.concatenate((rolling_cadence[iseason+2],season_remain))

            
            #Save the final cadence in a pkl file
            final_cadence=np.zeros((0,1),dtype=seasons[0].dtype)
            for key,val in rolling_cadence.items():
                final_cadence=np.concatenate((final_cadence,val))

            pkl_file = open('Rolling_Cadence_'+self.fieldName+'_'+str(self.fieldID_ref[0])+'.pkl','w')
            pkl.dump(final_cadence, pkl_file)
            pkl_file.close()


            #Illustrate results with a plot...
            self.Plot(self.fieldID_ref,seasons,rolling_cadence)
            plt.show()

    def Fake_Multiple(self,dataSlice):
        
        itag=False
        if dataSlice[self.fieldID][0] in self.fieldID_ref:
            itag=True
            dataSlice.sort(order=self.mjdCol)
            self.Fields[dataSlice[self.fieldID][0]]=dataSlice
            self.Fields_Seasons[dataSlice[self.fieldID][0]]=self.Get_Seasons(dataSlice)

        if itag and len(self.Fields)==len(self.fieldID_ref):

            self.Fields_Rolling={}
            
            for key,val in self.Fields_Seasons.items():
                rolling_cadence={}
            
                for i in range(len(val)):
                    rolling_cadence[i]=np.zeros((0,1),dtype=val[0].dtype)
                self.Fields_Rolling[key]=rolling_cadence
        
            #Now let us do the merging
            #First season : no change

            for key, seasons in self.Fields_Seasons.items():
                season_keep,season_remain=self.Split_Season(seasons[0])
                #print 'allo',len(season_keep),len(season_remain),len(rolling_cadence[0]),season_keep.shape,rolling_cadence[0].shape
                self.Fields_Rolling[key][0]=np.concatenate((self.Fields_Rolling[key][0],season_keep))

            #season 1 : field_A = field_A + merge_factor * field_B + merge_factor*field_C
            #           field_B = (1.-merge_factor) * field_B
            #           field_C = (1.-merge_factor) * field_C
            #           .......

            #season 2 : field_ B= field_B + merge_factor * field_A + merge_factor*field_B
            #           field_A = (1.-merge_factor) * field_A
            #           field_C = (1.-merge_factor) * field_C
            # 

            seasons_to_tag=[]
            nseason=1
            while 10.-nseason >= len(self.fieldID_ref):
                seasons_to_tag.append(nseason)
                nseason+=len(self.fieldID_ref)

            for season_tag in seasons_to_tag:

                combi=list(self.fieldID_ref)
                jcount=0
                for i in range(season_tag,season_tag+len(self.fieldID_ref)):
                    jcount+=1
                
                    iseason=i
            
                    self.Concat(combi,iseason)

                    if jcount < len(combi):
                        combi[0],combi[jcount]=combi[jcount],combi[0]
                        

            #and for the remaining periods...
            for idx in self.fieldID_ref:
                for io in range(len(self.Fields_Rolling[idx])):
                    if len(self.Fields_Rolling[idx][io]) == 0.:
                        season_keep,season_remain=self.Split_Season(self.Fields_Seasons[idx][io])
                        self.Fields_Rolling[idx][io]=np.concatenate((self.Fields_Rolling[idx][io],season_keep))

            #Dump Results in pkl files
            for idx in self.fieldID_ref:
                final_cadence=np.zeros((0,1),dtype=self.Fields_Rolling[idx][0].dtype)
                for key,val in self.Fields_Rolling[idx].items():
                    final_cadence=np.concatenate((final_cadence,val.reshape((len(val),1))))

                pkl_file = open('Rolling_Cadence_'+self.fieldName+'_'+str(idx)+'_'+str(len(self.fieldID_ref))+'_'+str(int(100.*self.merge_factor))+'.pkl','w')
                #print 'size',idx,len(final_cadence)
                pkl.dump(final_cadence, pkl_file)
                pkl_file.close()


            for idx in self.fieldID_ref:
                self.Plot(idx,self.Fields_Seasons[idx],self.Fields_Rolling[idx])
            plt.show()

    def Concat(self,combi,iseason):

        percent=[self.merge_factor]*len(self.fieldID_ref)
        percent[0]=1

        shift=[]
        
        ref_season=self.Fields_Seasons[combi[0]][iseason]

        for j,val in enumerate(combi):
            season=self.Fields_Seasons[val][iseason]
            shift.append(self.Shift(ref_season,season))

        for band in self.bands:
            season_keep=[]
            season_remain=[]
            for j,val in enumerate(combi):
                season=self.Fields_Seasons[val][iseason]
                season_band=season[np.where(season[self.filterCol]==band)]
                season_k,season_rem=self.Split_Season(season_band,percent[j])
                season_keep.append(season_k)
                season_remain.append(season_rem)

            for j,val in enumerate(combi):
                if j ==0:
                    self.Fields_Rolling[val][iseason]=np.concatenate((self.Fields_Rolling[val][iseason],season_keep[j]))
                    fieldRa=self.Fields_Rolling[val][iseason][0][self.fieldRA]
                    fieldDec=self.Fields_Rolling[val][iseason][0][self.fieldDec]
                    for kk in range(1,len(self.fieldID_ref)):
                    #shift=self.Shift(self.Fields_Rolling[val][iseason],season_keep[j+kk])
                        self.Fields_Rolling[val][iseason]=np.concatenate((self.Fields_Rolling[val][iseason],self.Modify(season_keep[j+kk],fieldRa[0],fieldDec[0],shift[j+kk])))               
                else:
                    self.Fields_Rolling[val][iseason]=np.concatenate((self.Fields_Rolling[val][iseason],season_remain[j]))

                

    def Modify(self, array,fieldRA=None, fieldDec=None,shift=None):
        
        array_copy=array.copy()
        
        #print 'before',array[self.fieldRA][0],array[self.fieldDec][0],array[self.mjdCol][0],array[self.airmass][0],array[self.filtSkyBrightness][0]
        
        if fieldRA is not None :
            array_copy[:][self.fieldRA]=fieldRA
        if fieldDec is  not  None:   
            array_copy[:][self.fieldDec]=fieldDec
        if shift is not None:
            array_copy[:][self.mjdCol]+=shift

            
        #An adjustment may need to be applied here
        #Since we are making a raw translation of mjd, it is not guaranteed that
        # shifted observations will keep a reasonable airmass (lower than 1.5; Opsim cut)
        #Thus need to tune mjd so as to get a reasonable airmass and retain observation

        for i in range(len(array_copy)):
            data = array_copy[i]
            ra=data[self.fieldRA][0]
            dec=data[self.fieldDec][0]
            mjd_ref=data[self.mjdCol][0]
            airmass=self.Get_airmass(mjd_ref,ra,dec)
            array_copy[i][self.airmass]=airmass
            array_copy[i][self.dateCol]=int((mjd_ref-simEpoch)*float(DAY))

            if airmass > 1.5:
                for tirage in range(100):
                    mjd_rand=np.random.uniform(mjd_ref-0.5,mjd_ref+0.5)
                    airmass_rand= self.Get_airmass(mjd_rand,ra,dec)
                    if airmass_rand < 1.5:
                        array_copy[i][self.mjdCol]=mjd_rand
                        array_copy[i][self.airmass]=airmass_rand
                        array_copy[i][self.dateCol]=int((mjd_rand-simEpoch)*float(DAY))
                        break
            mysky=SkyBrightness(data[self.fieldRA][0],data[self.fieldDec][0],data[self.mjdCol][0],data[self.dateCol][0],data[self.filterCol][0],data[self.airmass][0])
            array_copy[i][self.filtSkyBrightness]=mysky.new_skybrightness()

        return array_copy
            
    def Shift(self,array_ref,array_add):

        """
        array_ref.sort(order=self.mjdCol)
        array_add.sort(order=self.mjdCol)

        min_ref=np.min(array_ref[self.mjdCol])
        min_add=np.min(array_add[self.mjdCol])
        t=array_ref[self.mjdCol]
        mean_diff_ref=np.mean([t[i+1]-t[i] for i in range(len(t)-1)])
        
        if min_ref >=min_add:
            resu=min_ref-min_add+mean_diff_ref/2.
        else:
            resu=min_ref-min_add-mean_diff_ref/2.
        
        """
        mean_ref=np.mean(array_ref[self.mjdCol])
        mean_add=np.mean(array_add[self.mjdCol])

        resu = int(mean_ref - mean_add)
            
        return resu

    def Get_Seasons(self,filtc):
 
        #Get the seasons of a given set of obs
        #A season is defined by a set of observations 
        # for which the consecutive difference in time is lower than 100 days
        #input : fitc = set of observations
        #output: dict of seasons; key = seson number (starting at 0), val = corresponding set of observations

        dict_for_seasons={}
        if len(filtc) > 0:
            inum=0
            dict_for_seasons[inum]=np.zeros((60,1),dtype=filtc.dtype)
           
                                
            iloop=0
            iinside=0
            dict_for_seasons[inum][iinside]=filtc[iloop]
            
            if len(filtc) > 1:
                while iloop < len(filtc)-1: 
                    iinside+=1
                    diff_time_days=filtc[self.mjdCol][iloop+1]-filtc[self.mjdCol][iloop]
                    if diff_time_days > 100.:
                        dict_for_seasons[inum]=np.resize(dict_for_seasons[inum],iinside)
                        inum+=1
                        dict_for_seasons[inum]=np.zeros((60,1),dtype=filtc.dtype)
                        iinside=0
                    if len(dict_for_seasons[inum]) <= iinside:
                        dict_for_seasons[inum]=np.resize(dict_for_seasons[inum],(len(dict_for_seasons[inum])+50,1))

                    dict_for_seasons[inum][iinside]=filtc[iloop+1]
   
                    iloop+=1
                
        #removing zeros (if any):
        outdict={}
        for key, val in dict_for_seasons.items():
            outdict[key]=val[np.where(val[self.fieldID]>0.)]
    
        return outdict

    def Split_Season(self,season,frac=1.):
    
        #Split a given season in two sets, depending on the fraction requested (frac)
        #input : a given season
        #Output:
        #season_copy : copy of (random) frac of season, with MJD->MJD-shift
        #season_remain : copy of (random) (1-frac) of season 


        if frac > 0.99:
            return np.reshape(season,(len(season),1)),np.zeros((0,1),dtype=season.dtype)
        

        season.sort(order=self.mjdCol)
        season_copy=np.zeros((0,1),dtype=season.dtype)
        season_remain=np.zeros((0,1),dtype=season.dtype)

        n_tobe_excluded=(1.-frac)*len(season)

        if n_tobe_excluded < 1 and  n_tobe_excluded > 0.:
            n_tobe_excluded=1
        
        
        excluded=[]
        while len(excluded)< np.rint(n_tobe_excluded):
            aleat=np.random.randint(0,len(season))
            if aleat not in excluded:
                excluded.append(aleat)
    
        inum=-1
        inum_remain=-1

        for i in range(len(season)):
            if i not in excluded:
                season_copy=np.vstack([season_copy,season[i]])
            else:
                season_remain=np.vstack([season_remain,season[i]])

        return season_copy,season_remain

    def Print_Results(self,idx,what,vals):

        print 'Field ID',idx
        print what
        print 'Season Ntot_Obs #u #g #r #i #z #y DT DT_u DT_g DT_r DT_i DT_z DT_y'
        for i in range(len(vals)):
            val=[np.asscalar(theval) for theval in vals[i][self.mjdCol]]
            """
            for theval in vals[i][self.mjdCol]:
                val.append(np.asscalar(theval))
            """
            vala=vals[i]
            tot=''
            val.sort()
            #print 'hello',val
            #print 'helli',[jj-llj for llj, jj in zip(val[:-1],val[1:])]
            mean_visit_av=np.mean([jj-llj for llj, jj in zip(val[:-1],val[1:])])

            for band in self.bands:
                select=vala[np.where(vala[self.filterCol]==band)]
                tot+=str(len(select))+' '

            tot+=str(round(mean_visit_av,2))+' '
            for band in self.bands:
                select=vala[np.where(vala[self.filterCol]==band)]
                if len(select)> 1:
                    bal=[np.asscalar(theval) for theval in select[self.mjdCol]]
                    mean_visit=np.mean([j-ll for ll, j in zip(bal[:-1],bal[1:])])
                else:
                   mean_visit=0 
               
                
                tot+=str(round(mean_visit,2))+' '
    
                    
                #tot+=str(len(select))+' '
            print i,len(val),tot
            
    def Plot(self,fieldid,seasons,rolling_cadence):
        
        figa, axa = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
        fontsize=15.
        
        #axa.scatter(dataSlice['expMJD'],dataSlice['airmass'],facecolors='none', edgecolors='r',marker='*')
       
        what=self.airmass
        for i in range(len(seasons)):
            axa.scatter(seasons[i][self.mjdCol],seasons[i][what],facecolors='none', edgecolors='r',marker='*')


        self.Print_Results(fieldid,'Initial dataSlice',seasons)
        self.Print_Results(fieldid,'Rolling Cadence',rolling_cadence)
       
        for i in range(len(rolling_cadence)):
            axa.scatter(rolling_cadence[i][self.mjdCol],rolling_cadence[i][what],facecolors='none', edgecolors='b')
            axa.set_ylabel(r'airmass',{'fontsize': fontsize})
            axa.set_xlabel(r'MJD',{'fontsize': fontsize})

        figa.suptitle('Field type:'+' Field ID: '+str(fieldid))

        #plt.show()

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
 
    def Get_Sidereal_Time_at_Site(self,mjd):
        
        lon_RAD=-70.7494 *DEG2RAD #obs site long (from sims_operations/conf/system/SiteCP.conf)
        
        lst_RAD = pal.gmst(mjd) + lon_RAD
        if lst_RAD < 0:
            lst_RAD += TWOPI

        return lst_RAD

    def Get_airmass(self,mjd,ra_RAD,dec_RAD):
        
        lha_RAD = self.Get_Sidereal_Time_at_Site(mjd)-ra_RAD

        lat_RAD= -30.2444* DEG2RAD #obs site lat (from sims_operations/conf/system/SiteCP.conf)

        (az_RAD, d1, d2, alt_RAD, d4, d5, pa_RAD, d7, d8) = pal.altaz(lha_RAD,dec_RAD,lat_RAD)

        # Altitude -> Zenith distance (degrees)
        targetZD_DEG = 90. - (alt_RAD* RAD2DEG)
        # Altitude -> Zenith distance (radian)
        zd_RAD = 1.5707963 - alt_RAD
        # Airmass
        #am = slalib.sla_airmas (zd_RAD)
        am = pal.airmas(zd_RAD)
        return am
