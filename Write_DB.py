import sqlite3

class Write_DB:
    def __init__(self, dbname='test_create.db'):
        
        self.conn = sqlite3.connect(dbname)
        self.table_name="Summary"
        cursor = self.conn.cursor()
        
        print 'Creating summary table %s' %(self.table_name)
        sql = 'create table %s (obsHistID int unsigned not null, sessionID int unsigned not null, ' %(self.table_name)
        sql += 'propID int, fieldID int unsigned not null, fieldRA double, fieldDec double, '
        sql += 'filter varchar(8), expDate int unsigned, expMJD double, night int unsigned, '
        sql += 'visitTime double, visitExpTime double, finRank double, transparency double, '
        sql += 'airmass double, vSkyBright double, filtSkyBrightness double, rotSkyPos double, rotTelPos double, lst double, '
        sql += 'altitude double, azimuth double, dist2Moon double, solarElong double, moonRA double, moonDec double, '
        sql += 'moonAlt double, moonAZ double, moonPhase double, sunAlt double, sunAz double, phaseAngle double, '
        sql += 'rScatter double, mieScatter double, moonIllum double, moonBright double, darkBright double, '
        sql += 'FWHMeff double, ditheredRA double, FWHMgeom double , ditheredDec double, '
        sql += 'rawSeeing double, wind double, humidity double, slewDist double, slewTime double, fiveSigmaDepth double);'
        print 'sql',sql
        ret = self.getDbData(cursor, sql)

    def getDbData(self,cursor, sql):
        cursor.execute(sql)
        ret = cursor.fetchall()
        return ret

    def insertDbData(self,cursor, sql):
        cursor.execute(sql)

    def insert(self,dataSlice):
        cursor = self.conn.cursor()
        sql = 'insert into %s (obsHistID, sessionID, propID, fieldID, fieldRA, fieldDec, filter, ' %(self.table_name)
        sql += 'expDate, expMJD, night, visitTime, visitExpTime, finRank, transparency, airmass, vSkyBright, '
        sql += 'filtSkyBrightness, rotSkyPos, rotTelPos, lst, altitude, azimuth, dist2Moon, solarElong, moonRA, moonDec, '
        sql += 'moonAlt, moonAZ, moonPhase, sunAlt, sunAz, phaseAngle, rScatter, mieScatter, moonIllum, '
        sql += 'moonBright, darkBright, rawSeeing, wind, humidity, slewDist, slewTime, fiveSigmaDepth, '
        sql += 'FWHMeff, ditheredRA, FWHMgeom, ditheredDec) values '
        sql += '(%d, %d, %d, %d, %f, %f, "%s", %d, %f, %d, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f)' % (dataSlice['obsHistID'],dataSlice['sessionID'],dataSlice['propID'],dataSlice['fieldID'],dataSlice['fieldRA'],dataSlice['fieldDec'],dataSlice['filter'],dataSlice['expDate'],dataSlice['expMJD'],dataSlice['night'],dataSlice['visitTime'],dataSlice['visitExpTime'],dataSlice['finRank'],dataSlice['transparency'],dataSlice['airmass'],dataSlice['vSkyBright'],dataSlice['filtSkyBrightness'],dataSlice['rotSkyPos'],dataSlice['rotTelPos'],dataSlice['lst'],dataSlice['altitude'],dataSlice['azimuth'],dataSlice['dist2Moon'],dataSlice['solarElong'],dataSlice['moonRA'],dataSlice['moonDec'],dataSlice['moonAlt'],dataSlice['moonAZ'],dataSlice['moonPhase'],dataSlice['sunAlt'],dataSlice['sunAz'],dataSlice['phaseAngle'],dataSlice['rScatter'],dataSlice['mieScatter'],dataSlice['moonIllum'],dataSlice['moonBright'],dataSlice['darkBright'],dataSlice['rawSeeing'],dataSlice['wind'],dataSlice['humidity'],dataSlice['slewDist'],dataSlice['slewTime'],dataSlice['fiveSigmaDepth'],dataSlice['FWHMeff'],dataSlice['ditheredRA'],dataSlice['FWHMgeom'],dataSlice['ditheredDec'])
        self.insertDbData(cursor, sql)
        self.conn.commit()
