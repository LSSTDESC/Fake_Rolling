import sqlite3
import matplotlib.pyplot as plt
import pickle as pkl
from astropy.table import Table
import numpy as np

def append(thelist, what, name,test):
    if name == test:
        thelist.append(what)


conn_origin=sqlite3.connect('/sps/lsst/data/dev/pgris/sims_operations/DB_Files/minion_1016_sqlite.db')

#conn_rolling=sqlite3.connect('Rolling_minion_1016_309_310_311_80.db')
conn_rolling=sqlite3.connect('Rolling_test.db')

selec="SELECT * from Summary WHERE fieldID == 1064 and filter=='y'"
#selec="SELECT * from Summary WHERE propID == 54"

cursor_origin = conn_origin.cursor()
cursor_origin.execute(selec)
names = [description[0] for description in cursor_origin.description]
index_a=names.index('expMJD')
index_b=names.index('airmass')
index_c=names.index('obsHistID')
index_d=names.index('filter')
index_e=names.index('fieldID')
index_f=names.index('fieldRA')
index_g=names.index('fieldDec')
index_h=names.index('filtSkyBrightness')

r=[]
for row in cursor_origin:
    #print 'ici',row[index_a],row[index_c],row[index_h],row[index_d],row[index_b]
    r.append((row[index_a],row[index_c],row[index_h],row[index_d],row[index_b]))

resua=np.rec.fromrecords(r, names = ['expMJD','obsHistID','filtSkyBrightness','filter','airmass'])


cursor_rolling = conn_rolling.cursor()
cursor_rolling.execute(selec)

r=[]
for row in cursor_rolling:
    #print 'la',row[index_a],row[index_c],row[index_h],row[index_d],row[index_b]
    r.append((row[index_a],row[index_c],row[index_h],row[index_d],row[index_b]))

resub=np.rec.fromrecords(r, names = ['expMJD','obsHistID','filtSkyBrightness','filter','airmass'])


resua.sort(order='expMJD')
resub.sort(order='expMJD')

print resua[:20]
print resub[:20]




"""
for row in cursor_origin:
    mjd_orig.append(row[index_a])
    airmass_orig.append(row[index_b])
    print 'ici',row[index_c]
    for rowb in cursor_rolling:
        print 'test',rowb[index_c],row[index_c]
        if rowb[index_c] == row[index_c]:
            print rowb[index_c],rowb[index_h],row[index_h]
"""

"""
        mjd_rolling.append(rowb[index_a])
        airmass_rolling.append(rowb[index_b])
        if rowb[index_e] == 1163:
            fieldids_rolling.append(rowb[index_e])
"""
#print 'Nevts',len(mjd_orig),len(mjd_rolling),len(fieldids_rolling)

"""
cursor_rolling = conn_rolling.cursor()
cursor_rolling.execute(selec)
"""
"""
thelist=[]
for rowb in cursor_rolling:
    print rowb[index_c]
    thelist.append(rowb[index_c])
"""
"""

thedir='/sps/lsst/data/dev/pgris/Make_Cadence/Obs_minion_1016'

prefix='Observations'
fieldname='WFD'
fieldid=309

fichname=prefix+'_'+fieldname+'_'+str(fieldid)+'.pkl'
    
thepkl=pkl.load(open(thedir+'/'+fichname,'rb'))

thedata=thepkl['dataSlice']

print 'hello',len(thedata)

thelist_data=[]

for data in thedata:
    thelist_data.append(data['obsHistID'])

missing=[]
for val in thelist:
    if not val in thelist_data:
        print 'not there',val
        missing.append(val)

print len(missing)
"""

"""
cursor_rolling = conn_rolling.cursor()
cursor_rolling.execute(selec)

thelist=[]

pos_tab= Table(names=('expMJD','airmass','fieldRA','fieldDec'), dtype=('f8','f8','f8','f8'))

for rowb in cursor_rolling:
    if rowb[index_c] in missing:
        pos_tab.add_row((rowb[index_a],rowb[index_b],rowb[index_f],rowb[index_g]))

pos_tab.sort('expMJD')
for val in pos_tab:
    print val['expMJD'],val['fieldRA'],val['fieldDec']
#print set(pos_tab['expMJD'])

figa, axa = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
axa.scatter(mjd_orig,airmass_orig,facecolors='none', edgecolors='r',marker='*')
axa.scatter(mjd_rolling,airmass_rolling,facecolors='none', edgecolors='b')

"""

plt.show()
