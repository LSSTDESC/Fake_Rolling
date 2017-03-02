import lsst.sims.maf.metricBundles as metricBundles
#import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.db as db
import lsst.sims.maf.utils as utils
from Fake_Rolling import Fake_Rolling
from optparse import OptionParser

parser = OptionParser()

parser.add_option("-f", "--fieldname", type="string", default='DD', help="filter [%default]")
parser.add_option("-i", "--fieldid", type="float", default=290, help="filter [%default]")
parser.add_option("-o", "--opsimrun", type="string", default='minion_1016', help="filter [%default]")

opts, args = parser.parse_args()

outDir ='Test'

#dbFile = '/data/pgris/sims_operation/Run_OpSim/enigma_1189_sqlite.db'
dbFile = '/data/pgris/sims_operation/Run_OpSim/'+opts.opsimrun+'_sqlite.db'
#dbFile = '/data/pgris/sims_operation/Run_OpSim/clrlsstsrv_1068_sqlite.db'
opsimdb = utils.connectOpsimDb(dbFile)
resultsDb = db.ResultsDb(outDir=outDir)

propinfo, proptags = opsimdb.fetchPropInfo()
print 'hello',proptags,propinfo

#field='DD'
#numID=744
#field='WFD'


metric=Fake_Rolling(m5Col='fiveSigmaDepth',fieldname=opts.fieldname,fieldID=opts.fieldid)
#slicer = slicers.HealpixSlicer(nside=256)


slicer=slicers.OpsimFieldSlicer()
sqlconstraint = utils.createSQLWhere(opts.fieldname, proptags)
mb = metricBundles.MetricBundle(metric, slicer, sqlconstraint)

mbD = {0:mb}


mbg =  metricBundles.MetricBundleGroup(mbD, opsimdb,
                                      outDir=outDir, resultsDb=resultsDb)

mbg.runAll()
