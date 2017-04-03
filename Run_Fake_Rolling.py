import lsst.sims.maf.metricBundles as metricBundles
#import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.db as db
import lsst.sims.maf.utils as utils
from Fake_Rolling import Fake_Rolling
import argparse
from shutil import copyfile


def do_list_of_fields(s):
    out=[]
    for val in s.split(','):
        print val
        outb=[]
        for vval in val.split(' '):
            if vval != '':
                outb.append(int(vval))
        out.append(outb)
    
    #x, y, z = map(int, s.split(','))
    return out


parser = argparse.ArgumentParser()

parser.add_argument("-f", "--fieldname", default='WFD', help="filter [%default]")
parser.add_argument("-i", "--fieldid",default="",help="list of ids [%default]")
parser.add_argument("-d", "--dbFile", default='None', help="dbFile to process [%default]")
parser.add_argument("-m", "--merge_factor", type=float,default=0.8, help="Merge factor between periods [%default]")
parser.add_argument("-s", "--save_DB", type=int,default=1, help="Save in a DB file[%default]")
parser.add_argument("-n", "--new_DB_name", default='Rolling.db', help="output dbFile [%default]")

opts= parser.parse_args()

outDir ='Test'

#dbFile = '/data/pgris/sims_operation/Run_OpSim/enigma_1189_sqlite.db'
#dbFile = '/data/pgris/sims_operation/Run_OpSim/'+opts.opsimrun+'_sqlite.db'
#dbFile = '/sps/lsst/data/dev/pgris/sims_operations/DB_Files/'+opts.opsimrun+'_sqlite.db'
dbFile = opts.dbFile
#dbFile = '/data/pgris/sims_operation/Run_OpSim/clrlsstsrv_1068_sqlite.db'
opsimdb = utils.connectOpsimDb(dbFile)
resultsDb = db.ResultsDb(outDir=outDir)

propinfo, proptags = opsimdb.fetchPropInfo()
print 'hello',proptags,propinfo

#field='DD'
#numID=744
#field='WFD'

copyfile(opts.dbFile,opts.new_DB_name)

combi=[[309,310,311],[312,313,314]]
combi=[[885,1782,2700]]
combi=[]
if opts.fieldid != "":
    combi=do_list_of_fields(opts.fieldid)

metric=Fake_Rolling(fieldname=opts.fieldname,fieldID=combi,merge_factor=opts.merge_factor,db_name=opts.new_DB_name)
#slicer = slicers.HealpixSlicer(nside=256)


slicer=slicers.OpsimFieldSlicer()
sqlconstraint = utils.createSQLWhere(opts.fieldname, proptags)
mb = metricBundles.MetricBundle(metric, slicer, sqlconstraint)

mbD = {0:mb}


mbg =  metricBundles.MetricBundleGroup(mbD, opsimdb,
                                      outDir=outDir, resultsDb=resultsDb)

mbg.runAll()
