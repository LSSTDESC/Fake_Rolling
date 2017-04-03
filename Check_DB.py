import sqlite3
import matplotlib.pyplot as plt

def append(thelist, what, name,test):
    if name == test:
        thelist.append(what)


conn_origin=sqlite3.connect('/sps/lsst/data/dev/pgris/sims_operations/DB_Files/minion_1016_sqlite.db')

#conn_rolling=sqlite3.connect('Rolling_minion_1016_309_310_311_80.db')
conn_rolling=sqlite3.connect('Rolling.db')

selec="SELECT * from Summary WHERE fieldID == 310"

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

for row in cursor_origin:
    mjd_orig.append(row[index_a])
    airmass_orig.append(row[index_b])

for rowb in cursor_rolling:
    mjd_rolling.append(rowb[index_a])
    airmass_rolling.append(rowb[index_b])
    
figa, axa = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
axa.scatter(mjd_orig,airmass_orig,facecolors='none', edgecolors='r',marker='*')
axa.scatter(mjd_rolling,airmass_rolling,facecolors='none', edgecolors='b')

plt.show()
