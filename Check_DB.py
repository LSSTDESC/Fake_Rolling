import sqlite3
import matplotlib.pyplot as plt
import pickle as pkl
from astropy.table import Table

def append(thelist, what, name,test):
    if name == test:
        thelist.append(what)


conn_origin=sqlite3.connect('/sps/lsst/data/dev/pgris/sims_operations/DB_Files/minion_1016_sqlite.db')

#conn_rolling=sqlite3.connect('Rolling_minion_1016_309_310_311_80.db')
conn_rolling=sqlite3.connect('Rolling.db')

selec="SELECT * from Summary WHERE fieldID == 309"

cursor_origin = conn_origin.cursor()
cursor_origin.execute(selec)

cursor_rolling = conn_rolling.cursor()
cursor_rolling.execute(selec)

names = [description[0] for description in cursor_origin.description]

mjd_orig=[]
airmass_orig=[]

mjd_rolling=[]
airmass_rolling=[]

index_a=names.index('expMJD')
index_b=names.index('airmass')
index_c=names.index('obsHistID')
index_d=names.index('filter')

for row in cursor_origin:
    mjd_orig.append(row[index_a])
    airmass_orig.append(row[index_b])

for rowb in cursor_rolling:
    mjd_rolling.append(rowb[index_a])
    airmass_rolling.append(rowb[index_b])
    
print 'Nevts',len(mjd_orig),len(mjd_rolling)

cursor_rolling = conn_rolling.cursor()
cursor_rolling.execute(selec)

thelist=[]
for rowb in cursor_rolling:
    print rowb[index_c]
    thelist.append(rowb[index_c])

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

cursor_rolling = conn_rolling.cursor()
cursor_rolling.execute(selec)

thelist=[]

pos_tab= Table(names=('expMJD','airmass'), dtype=('f8','f8'))
for rowb in cursor_rolling:
    if rowb[index_c] in missing:
        #print rowb[index_a],rowb[index_b],rowb[index_c],rowb[index_d]
        pos_tab.add_row((rowb[index_a],rowb[index_b]))

pos_tab.sort('expMJD')
print pos_tab['expMJD']
print set(pos_tab['expMJD'])

figa, axa = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
axa.scatter(mjd_orig,airmass_orig,facecolors='none', edgecolors='r',marker='*')
axa.scatter(mjd_rolling,airmass_rolling,facecolors='none', edgecolors='b')




plt.show()
