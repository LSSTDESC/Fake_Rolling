import numpy as np
import matplotlib.pyplot as plt
from lsst.sims.maf.metrics import BaseMetric
import pickle as pkl
from Sky_Brightness_With_Moonlight import SkyBrightness
import palpy as pal
import math
from Parameters import parameters
from lsst.sims.photUtils import SignalToNoise
#from Write_DB import Write_DB
import sqlite3
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import Table


DEG2RAD = math.pi / 180.    # radians = degrees * DEG2RAD
RAD2DEG = 180. / math.pi    # degrees = radians * RAD2DEG
TWOPI = 2 * math.pi
DAY = 86400.
simEpoch=59580

class Fake_Rolling(BaseMetric):
    """
    Measure how many time series meet a given time and filter distribution requirement.Fake 
    """
    
    def __init__(self, metricName='Fake_Rolling', units='', badval=-666,uniqueBlocks=False,fieldname='DD',fieldID=[[309,310,311]],merge_factor=0.8,n_merger=3,ra_grid=2.,db_name='Rolling.db',**kwargs):

        self.mjdCol ='expMJD'
        self.filterCol='filter'
        self.m5Col = 'fiveSigmaDepth'
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
        self.n_merger=n_merger #taken only when the full sky is reshuffled
        self.ra_grid=ra_grid

        self.bands=['u','g','r','i','z','y']

        self.cols=['obsHistID', 'sessionID', 'propID', 'fieldID', 'fieldRA', 'fieldDec', 'filter', 'expDate', 'expMJD', 'night', 'visitTime', 'visitExpTime', 'finRank', 'FWHMeff', 'FWHMgeom', 'transparency', 'airmass', 'vSkyBright', 'filtSkyBrightness', 'rotSkyPos', 'rotTelPos', 'lst', 'altitude', 'azimuth', 'dist2Moon', 'solarElong', 'moonRA', 'moonDec', 'moonAlt', 'moonAZ', 'moonPhase', 'sunAlt', 'sunAz', 'phaseAngle', 'rScatter', 'mieScatter', 'moonIllum', 'moonBright', 'darkBright', 'rawSeeing', 'wind', 'humidity', 'slewDist', 'slewTime', 'fiveSigmaDepth', 'ditheredRA', 'ditheredDec']

        super(Fake_Rolling, self).__init__(col=self.cols,metricName=metricName, units=units, badval=badval, **kwargs)

        self.uniqueBlocks = uniqueBlocks

        self.db_name=db_name
        self.data={}
        itot=0
        for combi in self.fieldID_ref:
            itot+=len(combi)
        self.tot_combi=itot

        
    def run(self, dataSlice, slicePoint=None):

        if self.tot_combi==1:
            self.Fake_Single(dataSlice)
        else:
            self.Fake_Multiple_new(dataSlice)

    def Fake_Single(self,dataSlice):

        dataSlice = dataSlice[np.where(dataSlice[self.fieldID]==self.fieldID_ref)]
        if len(dataSlice) > 0:

            conn = sqlite3.connect(self.db_name)
            cursor=conn.cursor()
            print "Opened database successfully - removing fieldid ",dataSlice[self.fieldID][0]
            cursor.execute("DELETE FROM Summary WHERE "+self.fieldID+" =="+str(dataSlice[self.fieldID][0])+";")
            conn.commit()
            conn.close()

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

            
            
            final_cadence=np.zeros((0,1),dtype=seasons[0].dtype)
            for key,val in rolling_cadence.items():
                final_cadence=np.concatenate((final_cadence,val))
            
            conn = sqlite3.connect(self.db_name)
            for ipo in range(len(final_cadence)):
                self.insert_in_DB(conn,final_cadence[ipo][0])
            conn.close()

            """

            #Save the final cadence in a pkl file 
            final_cadence=np.zeros((0,1),dtype=seasons[0].dtype)
            for key,val in rolling_cadence.items():
                final_cadence=np.concatenate((final_cadence,val))

            dictout={}
            dictout['dataSlice']=final_cadence

            pkl_file = open('Rolling_Cadence_'+self.fieldName+'_'+str(self.fieldID_ref[0])+'.pkl','w')
            pkl.dump(dictout, pkl_file)
            pkl_file.close()
            """

            #Illustrate results with a plot...
            #self.Plot(self.fieldID_ref,seasons,rolling_cadence)
            #plt.show()


    def Fake_Multiple_new(self,dataSlice):

        self.data[dataSlice[self.fieldID][0]]=dataSlice
        
        combi_tot=self.fieldID_ref
        
        if dataSlice[self.fieldID][0]==2776:
            #print 'ok pal'
            if len(combi_tot)==0:
                combi_tot=self.Get_Combis(self.data)

                print combi_tot
            for combi in combi_tot:
                self.Fake_Multiple(combi)

        
    def Fake_Multiple(self,combi_orig=[]):
    
        Fields={}
        Fields_Seasons={}
        Fields_Rolling={}

        for fieldID in combi_orig:
            self.data[fieldID].sort(order=self.mjdCol)
            Fields[fieldID]=self.data[fieldID]
            Fields_Seasons[fieldID]=self.Get_Seasons(self.data[fieldID])
            #open the database and remove this field
            conn = sqlite3.connect(self.db_name, timeout=30)
            cursor=conn.cursor()
            print "Opened database successfully - removing fieldID ",fieldID,'nseasons =',len(Fields_Seasons[fieldID])
            cursor.execute("DELETE FROM Summary WHERE "+self.fieldID+" =="+str(fieldID)+";")
            conn.commit()
            conn.close()
       
        for key,val in Fields_Seasons.items():
            rolling_cadence={}
            
            for i in range(len(val)):
                rolling_cadence[i]=np.zeros((0,1),dtype=val[0].dtype)
            Fields_Rolling[key]=rolling_cadence     

        for key, seasons in Fields_Seasons.items():
            season_keep,season_remain=self.Split_Season(seasons[0])
            Fields_Rolling[key][0]=np.concatenate((Fields_Rolling[key][0],season_keep))

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
        while 10.-nseason >= len(combi_orig):
            seasons_to_tag.append(nseason)
            nseason+=len(combi_orig)

        #print seasons_to_tag
        for season_tag in seasons_to_tag:
            
            combi=list(combi_orig)
            jcount=0
            for i in range(season_tag,season_tag+len(combi)):
                jcount+=1
                
                iseason=i
            
                self.Concat(combi,iseason,Fields_Seasons,Fields_Rolling)

                if jcount < len(combi):
                        combi[0],combi[jcount]=combi[jcount],combi[0]

        #and for the remaining periods...
        for idx in combi_orig:
            for io in range(len(Fields_Rolling[idx])):
                if len(Fields_Rolling[idx][io]) == 0.:
                    season_keep,season_remain=self.Split_Season(Fields_Seasons[idx][io])
                    Fields_Rolling[idx][io]=np.concatenate((Fields_Rolling[idx][io],season_keep))
        #Dump Results in pkl files
        final_cadence={}
        for idx in combi_orig:
            final_cadence[idx]=np.zeros((0,1),dtype=Fields_Rolling[idx][0].dtype)
            for key,val in Fields_Rolling[idx].items():
                final_cadence[idx]=np.concatenate((final_cadence[idx],val.reshape((len(val),1))))

            print 'inserting field',idx 
                
            conn = sqlite3.connect(self.db_name, timeout=50)
            for val in final_cadence[idx]:
                self.insert_in_DB(conn,val[0])
    
            conn.close()

    def Fake_Multiple_old(self,dataSlice):
        
        itag=False
        if dataSlice[self.fieldID][0] in self.fieldID_ref:
            itag=True
            dataSlice.sort(order=self.mjdCol)
            self.Fields[dataSlice[self.fieldID][0]]=dataSlice
            self.Fields_Seasons[dataSlice[self.fieldID][0]]=self.Get_Seasons(dataSlice)
            #open the database and remove this field
            conn = sqlite3.connect(self.db_name)
            cursor=conn.cursor()
            print "Opened database successfully - removing fieldid ",dataSlice[self.fieldID][0]
            cursor.execute("DELETE FROM Summary WHERE "+self.fieldID+" =="+str(dataSlice[self.fieldID][0])+";")
            conn.commit()
            conn.close()

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
            final_cadence={}
            for idx in self.fieldID_ref:
                final_cadence[idx]=np.zeros((0,1),dtype=self.Fields_Rolling[idx][0].dtype)
                for key,val in self.Fields_Rolling[idx].items():
                    final_cadence[idx]=np.concatenate((final_cadence[idx],val.reshape((len(val),1))))

                
                print 'inserting',idx 
                #open the database and remove this field
                
                conn = sqlite3.connect(self.db_name)
                for val in final_cadence[idx]:
                    self.insert_in_DB(conn,val[0])
                """
                if idx == 309:
                    #print set(final_cadence[idx]['fieldID'])
                    select=final_cadence[idx][np.where(np.logical_and(final_cadence[idx]['expMJD'] > 60000, final_cadence[idx]['expMJD']< 60500))]
                    select.sort(order='expMJD')
                    print select['expMJD']
                """
                conn.close()
                
                """
                dictout={}
                dictout['dataSlice']=final_cadence
                pkl_file = open('Rolling_Cadence_'+self.fieldName+'_'+str(idx)+'_'+str(len(self.fieldID_ref))+'_'+str(int(100.*self.merge_factor))+'.pkl','w')
                #print 'size',idx,len(final_cadence)
                pkl.dump(dictout, pkl_file)
                pkl_file.close()
                """
            """
            for idx in self.fieldID_ref:
                #self.Plot(idx,self.Fields_Seasons[idx],self.Fields_Rolling[idx])
                self.Plot(idx,self.Fields_Seasons[idx],final_cadence[idx])
            plt.show()
            """
    def Concat(self,combi,iseason,Fields_Seasons,Fields_Rolling):

        percent=[self.merge_factor]*len(combi)
        percent[0]=1

        shift=[]

        #check whether seasons exist for the considered fields
        #if not -> do not do enything

        for j,val in enumerate(combi):
            #print 'alors?',len(Fields_Seasons[val]),iseason
            if iseason>=len(Fields_Seasons[val]):
                #print 'should be removed ?'
                return
        
        ref_season=Fields_Seasons[combi[0]][iseason]

        for j,val in enumerate(combi):
            season=Fields_Seasons[val][iseason]
            shift.append(self.Shift(ref_season,season))

        for band in self.bands:
            season_keep=[]
            season_remain=[]
            for j,val in enumerate(combi):
                season=Fields_Seasons[val][iseason]
                season_band=season[np.where(season[self.filterCol]==band)]
                season_k,season_rem=self.Split_Season(season_band,percent[j])
                season_keep.append(season_k)
                season_remain.append(season_rem)

            for j,val in enumerate(combi):
                if j ==0:
                    Fields_Rolling[val][iseason]=np.concatenate((Fields_Rolling[val][iseason],season_keep[j]))
                    fieldRa=Fields_Rolling[val][iseason][0][self.fieldRA]
                    fieldDec=Fields_Rolling[val][iseason][0][self.fieldDec]
                    ditheredRa=Fields_Rolling[val][iseason][0][self.ditheredRA]
                    ditheredDec=Fields_Rolling[val][iseason][0][self.ditheredDec]
                    fieldID=Fields_Rolling[val][iseason][0][self.fieldID]


                    for kk in range(1,len(combi)):
                    #shift=self.Shift(self.Fields_Rolling[val][iseason],season_keep[j+kk])
                        Fields_Rolling[val][iseason]=np.concatenate((Fields_Rolling[val][iseason],self.Modify(season_keep[j+kk],fieldRa[0],fieldDec[0],ditheredRa[0],ditheredDec[0],shift[j+kk],fieldID)))               
                else:
                    Fields_Rolling[val][iseason]=np.concatenate((Fields_Rolling[val][iseason],season_remain[j]))

                
    """

    def Concat_old(self,combi,iseason):

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
                    fieldID=self.Fields_Rolling[val][iseason][0][self.fieldID]
                    for kk in range(1,len(self.fieldID_ref)):
                    #shift=self.Shift(self.Fields_Rolling[val][iseason],season_keep[j+kk])
                        self.Fields_Rolling[val][iseason]=np.concatenate((self.Fields_Rolling[val][iseason],self.Modify(season_keep[j+kk],fieldRa[0],fieldDec[0],shift[j+kk],fieldID)))               
                else:
                    self.Fields_Rolling[val][iseason]=np.concatenate((self.Fields_Rolling[val][iseason],season_remain[j]))

    """          

    def Modify(self, array,fieldRA=None, fieldDec=None,ditheredRA=None, ditheredDec=None,shift=None, fieldID=None):
        
        array_copy=array.copy()
        
        #print 'before',array[self.fieldRA][0],array[self.fieldDec][0],array[self.mjdCol][0],array[self.airmass][0],array[self.filtSkyBrightness][0]
        
        if fieldRA is not None :
            array_copy[:][self.fieldRA]=fieldRA
            array_copy[:][self.ditheredRA]=ditheredRA
        if fieldDec is  not  None:   
            array_copy[:][self.fieldDec]=fieldDec
            array_copy[:][self.ditheredDec]=ditheredDec
        if fieldID is  not  None:   
            array_copy[:][self.fieldID]=fieldID

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
            array_copy[i][self.m5Col]=self.Recompute_fiveSigmaDepth(array_copy[i][self.airmass],array_copy[i][self.visitTime],array_copy[i][self.filtSkyBrightness],array_copy[i][self.filterCol][0],array_copy[i][self.rawSeeing])

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

    def Recompute_fiveSigmaDepth(self,airmass,visitexptime,filtskybrightness,filtre,seeing):
        
        param=parameters()
        Filter_Wavelength_Correction = np.power(500.0 / param.filterWave[filtre], 0.3)
        Airmass_Correction = math.pow(airmass,0.6)
        FWHM_Sys = param.FWHM_Sys_Zenith * Airmass_Correction
        FWHM_Atm = seeing * Filter_Wavelength_Correction * Airmass_Correction
        finSeeing = param.scaleToNeff * math.sqrt(np.power(FWHM_Sys,2) + param.atmNeffFactor * np.power(FWHM_Atm,2))
        FWHMeff = SignalToNoise.FWHMgeom2FWHMeff(finSeeing)
         
        Tscale = visitexptime/ 30.0 * np.power(10.0, -0.4*(filtskybrightness - param.msky[filtre]))
        dCm = param.dCm_infinity[filtre] - 1.25*np.log10(1 + np.power(10.,0.8*param.dCm_infinity[filtre]- 1.)/Tscale)

        m5_recalc=dCm+param.Cm[filtre]+0.5*(filtskybrightness-21.)+2.5*np.log10(0.7/finSeeing)-param.kAtm[filtre]*(airmass-1.)+1.25*np.log10(visitexptime/30.)
        
        return m5_recalc

    def insert_in_DB(self,conn,dataSlice):
        cursor = conn.cursor()
        table_name='Summary'                     
        sql = 'insert into %s (obsHistID, sessionID, propID, fieldID, fieldRA, fieldDec, filter, ' %(table_name)
        sql += 'expDate, expMJD, night, visitTime, visitExpTime, finRank, transparency, airmass, vSkyBright, '
        sql += 'filtSkyBrightness, rotSkyPos, rotTelPos, lst, altitude, azimuth, dist2Moon, solarElong, moonRA, moonDec, '
        sql += 'moonAlt, moonAZ, moonPhase, sunAlt, sunAz, phaseAngle, rScatter, mieScatter, moonIllum, '
        sql += 'moonBright, darkBright, rawSeeing, wind, humidity, slewDist, slewTime, fiveSigmaDepth, '
        sql += 'FWHMeff, ditheredRA, FWHMgeom, ditheredDec) values '
        sql += '(%d, %d, %d, %d, %f, %f, "%s", %d, %f, %d, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f)' % (dataSlice['obsHistID'],dataSlice['sessionID'],dataSlice['propID'],dataSlice['fieldID'],dataSlice['fieldRA'],dataSlice['fieldDec'],dataSlice['filter'],dataSlice['expDate'],dataSlice['expMJD'],dataSlice['night'],dataSlice['visitTime'],dataSlice['visitExpTime'],dataSlice['finRank'],dataSlice['transparency'],dataSlice['airmass'],dataSlice['vSkyBright'],dataSlice['filtSkyBrightness'],dataSlice['rotSkyPos'],dataSlice['rotTelPos'],dataSlice['lst'],dataSlice['altitude'],dataSlice['azimuth'],dataSlice['dist2Moon'],dataSlice['solarElong'],dataSlice['moonRA'],dataSlice['moonDec'],dataSlice['moonAlt'],dataSlice['moonAZ'],dataSlice['moonPhase'],dataSlice['sunAlt'],dataSlice['sunAz'],dataSlice['phaseAngle'],dataSlice['rScatter'],dataSlice['mieScatter'],dataSlice['moonIllum'],dataSlice['moonBright'],dataSlice['darkBright'],dataSlice['rawSeeing'],dataSlice['wind'],dataSlice['humidity'],dataSlice['slewDist'],dataSlice['slewTime'],dataSlice['fiveSigmaDepth'],dataSlice['FWHMeff'],dataSlice['ditheredRA'],dataSlice['FWHMgeom'],dataSlice['ditheredDec'])
        cursor.execute(sql)
        conn.commit()


    def Get_Combis(self, data):

         pos_tab= Table(names=('fieldID','fieldRA','fieldDec'), dtype=('i4', 'f8','f8'))

         for key,val in data.items():
             pos_tab.add_row((key,val['fieldRA'][0],val['fieldDec'][0]))

         ra_step=self.ra_grid # degrees
         n_dec_zone=self.n_merger # number of region in dec = number of fields to be merged

         ra_min=np.min(pos_tab['fieldRA'])
    
         ra_dec_strip={}
         istrip=-1
    
         while ra_min < 360.-ra_step:
             istrip+=1
             ra_dec_strip[istrip]={}
             ra_max=ra_min+ra_step
             sel=pos_tab[np.where(np.logical_and(np.rad2deg(pos_tab['fieldRA'])>=ra_min,np.rad2deg(pos_tab['fieldRA'])<ra_max))]
             sel.sort('fieldDec')
             num_per_part=len(sel)/n_dec_zone
             ntag=0
             for count in range(n_dec_zone):
                 if count == n_dec_zone-1:
                     ra_dec_strip[istrip][count]=sel[ntag:]
                 else:
                     ra_dec_strip[istrip][count]=sel[ntag:ntag+num_per_part] 
            #print count, len(sel),num_per_part,len(ra_dec_strip[istrip][count])
                 ntag+=num_per_part

             ra_min+=ra_step
             #break


         ra_dec_final_combi={}
    
         for iv in range(len(ra_dec_strip)):
        #print 'o yes',iv
             ra_dec_final_combi[iv]={}
             strip_copy={}
             icombi=-1
             for i in range(1,n_dec_zone):
                 strip_copy[i]=ra_dec_strip[iv][i].copy()
             for val in ra_dec_strip[iv][0]:
                 icombi+=1
                 restable=Table(names=('fieldID','fieldRA','fieldDec'), dtype=('i4', 'f8','f8'))
                 restable.add_row((val))
                 for  iu in range(1,n_dec_zone):
                     strip_copy[iu],resu=self.Get_Nearest(strip_copy[iu],val)
                     restable.add_row((resu))
                 ra_dec_final_combi[iv][icombi]=restable

         all_combo=[]
         for key,vval in ra_dec_final_combi.items():
             for key,val in vval.items():
                 local=[]
                 for ik in range(len(val)):
                     local.append(val['fieldID'][ik])
                 all_combo.append(local)
         
         return all_combo

    def Get_Nearest(self,orig_table, val):
   
        table=Table(orig_table)
        c = SkyCoord(ra=val['fieldRA']*u.radian, dec=val['fieldDec']*u.radian)  
        catalog = SkyCoord(ra=table['fieldRA']*u.radian, dec=table['fieldDec']*u.radian)  
        idx, d2d, d3d = c.match_to_catalog_sky(catalog)

    #print 'astropy matching',idx,d2d,d3d,len(table),type(idx),table
        theres=[table['fieldID'][int(idx)],table['fieldRA'][int(idx)],table['fieldDec'][int(idx)]]
        table.remove_row(int(idx))
        return table,theres
