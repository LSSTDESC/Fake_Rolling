import numpy as np
import matplotlib.pyplot as plt
from lsst.sims.maf.metrics import BaseMetric
import pickle as pkl

class Fake_Rolling(BaseMetric):
    """
    Measure how many time series meet a given time and filter distribution requirement.Fake 
    """
    
    def __init__(self, metricName='Fake_Rolling', m5Col='fiveSigmaDepth', units='', badval=-666,uniqueBlocks=False,fieldname='DD',fieldID=290,**kwargs):
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


        self.fieldName=fieldname
        self.fieldID_ref=int(fieldID)

        self.filters=['u','g','r','i','z','y']

        super(Fake_Rolling, self).__init__(col=[self.mjdCol,self.m5Col, self.filterCol, self.dateCol,self.fieldRA,self.fieldDec, self.ditheredRA,self.ditheredDec,self.visitTime,self.rawSeeing,self.moonPhase,self.airmass,self.filtSkyBrightness,self.fieldID],metricName=metricName, units=units, badval=badval, **kwargs)

        self.uniqueBlocks = uniqueBlocks
        

    def run(self, dataSlice, slicePoint=None):

        dataSlice = dataSlice[np.where(dataSlice[self.fieldID]==self.fieldID_ref)]
        if len(dataSlice) > 0:
            dataSlice.sort(order=self.mjdCol)
            
            seasons=self.Get_Seasons(dataSlice)

            rolling_cadence={}
            
            for i in range(len(seasons)):
                rolling_cadence[i]=np.zeros((0,1),dtype=seasons[0].dtype)
                
            #Season 0 : keep the same
            season_keep,season_remain=self.Split_Season(seasons[0])
            rolling_cadence[0]=np.concatenate((rolling_cadence[0],season_keep))
             
            # Season i in [1,4,7]:
            for iseason in [1,4,7]:
                # keep season i
                for band in self.filters:
                    season_ref=seasons[iseason][np.where(seasons[iseason][self.filterCol]==band)]
                    season_ref_plus_one=seasons[iseason+1][np.where(seasons[iseason+1][self.filterCol]==band)]
                    season_ref_plus_two=seasons[iseason+2][np.where(seasons[iseason+2][self.filterCol]==band)]         
                    season_keep,season_remain=self.Split_Season(season_ref)
                    rolling_cadence[iseason]=np.concatenate((rolling_cadence[iseason],season_keep))
                    # from season i+1: keep 20% in i+1, put 80% in i (MJD[i] = MJD[i+1]-365.)
                    season_keep,season_remain=self.Split_Season(season_ref_plus_one,365.,0.8)
                    rolling_cadence[iseason]=np.concatenate((rolling_cadence[iseason],season_keep))
                    #print 'ohoh',len(rolling_cadence[iseason]),band,len(season_ref),len(season_ref_plus_one),len(season_ref_plus_two),len(season_remain),len(season_keep)
                    rolling_cadence[iseason+1]=np.concatenate((rolling_cadence[iseason+1],season_remain))
                    #from season i+2: keep 20% in i+2, put 80% in i (MJD[i] = MJD[i+2]-730.)
                    season_keep,season_remain=self.Split_Season(season_ref_plus_two,730.,0.8)
                    rolling_cadence[iseason]=np.concatenate((rolling_cadence[iseason],season_keep))
                    rolling_cadence[iseason+2]=np.concatenate((rolling_cadence[iseason+2],season_remain))

            
            #Save the final cadence in a pkl file
            final_cadence=np.zeros((0,1),dtype=seasons[0].dtype)
            for key,val in rolling_cadence.items():
                final_cadence=np.concatenate((final_cadence,val))

            pkl_file = open('Rolling_Cadence_'+self.fieldName+'_'+str(self.fieldID_ref)+'.pkl','w')
            pkl.dump(final_cadence, pkl_file)
            pkl_file.close()


            #Illustrate results with a plot...
            self.Plot(dataSlice,seasons,rolling_cadence)

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

    def Split_Season(self,season,shift=0.,frac=1.):
    
        #Split a given season in two sets, depending on the fraction requested (frac)
        #input : a given season
        #Output:
        #season_copy : copy of (random) frac of season, with MJD->MJD-shift
        #season_remain : copy of (random) (1-frac) of season 

        season.sort(order=self.mjdCol)
        season_copy=np.zeros((60,1),dtype=season.dtype)
        season_remain=np.zeros((60,1),dtype=season.dtype)

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
                inum+=1
                if len(season_copy) <= inum:
                    season_copy=np.resize(season_copy,(len(season_copy)+50,1))
                for data_name in season.dtype.names:
                    season_copy[data_name][inum]=season[data_name][i]
                    season_copy[self.mjdCol][inum]=season[self.mjdCol][i]-shift
                
            else:
                inum_remain+=1
                if len(season_remain) <= inum_remain:
                    season_remain=np.resize(season_remain,(len(season_remain)+50,1))

                for data_name in season.dtype.names:
                    season_remain[data_name][inum_remain]=season[data_name][i]
                    season_remain[self.mjdCol][inum_remain]=season[self.mjdCol][i] 


        season_copy=np.resize(season_copy,(inum+1,1))
        season_remain=np.resize(season_remain,(inum_remain+1,1))

        return season_copy,season_remain

    def Plot(self,dataSlice,seasons,rolling_cadence):
        
        figa, axa = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
        fontsize=15.
        
        axa.scatter(dataSlice['expMJD'],dataSlice['airmass'],facecolors='none', edgecolors='r',marker='*')

        print 'Initial dataSlice'
        print 'Season Ntot_Obs #u #g #r #i #z #y'
        for key,val in seasons.items():
            tot=''
            for band in ['u','g','r','i','z','y']:
                select=val[np.where(val['filter']==band)]
                tot+=str(len(select))+' '

            print key,len(val),tot

        print 'Rolling cadence'
        print 'Season Ntot_Obs #u #g #r #i #z #y'
        for i,val in seasons.items():
            tot=''
            for band in ['u','g','r','i','z','y']:
                select=rolling_cadence[i][np.where(rolling_cadence[i]['filter']==band)]
                tot+=str(len(select))+' '

            if i == 0:
                rolling_cat=rolling_cadence[i]
            else:
                rolling_cat=np.concatenate((rolling_cat,rolling_cadence[i]))


            print i,len(rolling_cadence[i]),tot

            axa.scatter(rolling_cadence[i]['expMJD'],rolling_cadence[i]['airmass'],facecolors='none', edgecolors='b')
            axa.set_ylabel(r'airmass',{'fontsize': fontsize})
            axa.set_xlabel(r'MJD',{'fontsize': fontsize})

        plt.show()
