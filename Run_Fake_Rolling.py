import lsst.sims.maf.metricBundles as metricBundles
#import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.db as db
import lsst.sims.maf.utils as utils
from Fake_Rolling import Fake_Rolling
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-f", "--fieldname", default='WFD', help="filter [%default]")
parser.add_argument("-i", "--fieldid",  nargs='+',type=int,help="list of ids [%default]")
parser.add_argument("-d", "--dbFile", default='None', help="dbFile to process [%default]")
parser.add_argument("-m", "--merge_factor", type=float,default=0.8, help="Merge factor between periods [%default]")


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


metric=Fake_Rolling(m5Col='fiveSigmaDepth',fieldname=opts.fieldname,fieldID=opts.fieldid,merge_factor=opts.merge_factor)
#slicer = slicers.HealpixSlicer(nside=256)


slicer=slicers.OpsimFieldSlicer()
sqlconstraint = utils.createSQLWhere(opts.fieldname, proptags)
mb = metricBundles.MetricBundle(metric, slicer, sqlconstraint)

mbD = {0:mb}


mbg =  metricBundles.MetricBundleGroup(mbD, opsimdb,
                                      outDir=outDir, resultsDb=resultsDb)

mbg.runAll()
