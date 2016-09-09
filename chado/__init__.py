from sqlalchemy import create_engine, MetaData, Table
from sqlalchemy.orm import mapper, sessionmaker

class Analysis(object):
    pass

class AnalysisProperty(object):
    pass

class Organism(object):
    pass

class Dbxref(object):
    pass

class Cvterm(object):
    pass

class Db(object):
    pass

class Cv(object):
    pass

class Feature(object):
    pass

class FeatureLocation(object):
    pass

class FeatureProperties(object):
    pass

class FeatureRelationship(object):
    pass

def ChadoAuth(parser):
    parser.add_argument('-o', '--dbhost', required=True, help='Database Host')
    parser.add_argument('-n', '--dbname', required=True, help='Database Name')
    parser.add_argument('-u', '--dbuser', help='Database Username')
    parser.add_argument('-w', '--dbpass', help='Database Password')
    parser.add_argument('-p', '--dbport', type=int, help='Database Port', default=5432)
    parser.add_argument('--dbschema', help='Database Schema (default: public)', default="public")
    parser.add_argument("-d", "--debug", help="Print debug information", action="store_true")

class ChadoInstance(object):

    def __init__(self, dbhost, dbname, dbuser, dbpass, dbschema, dbport, debug, **kwargs):
        self.dbhost = dbhost
        self.dbname = dbname
        self.dbuser = dbuser
        self.dbpass = dbpass
        self.dbport = dbport
        self.dbschema = dbschema

        self.debug = debug

    def __str__(self):
        return '<ChadoInstance at %s>' % self.dbhost

    def connect(self):

        self._engine = create_engine('postgresql://%s:%s@%s:%s/%s' % (self.dbuser, self.dbpass, self.dbhost, self.dbport, self.dbname), echo=self.debug)
        self._metadata = MetaData(self._engine, schema=self.dbschema)

        Session = sessionmaker(bind=self._engine)
        self.session = Session()

        self._test_db_access()

        self._reflect_tables()
        self._cv_id_cache = {}

    def _test_db_access(self):
        tables = self._engine.table_names(schema=self.dbschema)
        if ('analysis' not in tables or 'feature' not in tables):
            raise Exception("Could not find Chado tables in db %s" % (self._engine.url))

    def _reflect_tables(self):

        # Try to reflect tables
        analysis = Table('analysis', self._metadata, autoload=True)
        mapper(Analysis, analysis)
        analysisprop = Table('analysisprop', self._metadata, autoload=True)
        mapper(AnalysisProperty, analysisprop)
        organism = Table('organism', self._metadata, autoload=True)
        mapper(Organism, organism)
        dbxref = Table('dbxref', self._metadata, autoload=True)
        mapper(Dbxref, dbxref)
        cvterm = Table('cvterm', self._metadata, autoload=True)
        mapper(Cvterm, cvterm)
        db = Table('db', self._metadata, autoload=True)
        mapper(Db, db)
        cv = Table('cv', self._metadata, autoload=True)
        mapper(Cv, cv)
        feature = Table('feature',    self._metadata, autoload=True)
        mapper(Feature, feature)
        featureloc = Table('featureloc', self._metadata, autoload=True)
        mapper(FeatureLocation, featureloc)
        featureprop = Table('featureprop', self._metadata, autoload=True)
        mapper(FeatureProperties, featureprop)
        feature_relationship = Table('feature_relationship', self._metadata, autoload=True)
        mapper(FeatureRelationship, feature_relationship)

    def get_cvterm_id(self, type_name, cv_name):
        res = self.session.query(Cvterm, Cv).filter(Cvterm.name == type_name, Cv.name == cv_name)
        if not res.count():
            raise Exception("Could not find a cvterm with name '%s' from cv ''%s' in the database %s" % (type_name, cv_name, self._engine.url))

        return res.one().Cvterm.cvterm_id

    def get_cvterm_name(self, cv_id):
        """
        get_cvterm_name allows lookup of CV terms by their ID.

        This method caches the result in order to not hit the DB for every
        query. Maybe should investigate pre-loading popular terms? (E.g. gene,
        mRNA, etc)
        """
        if cv_id in self._cv_id_cache:
            if self._cv_id_cache[cv_id] is not None:
                return self._cv_id_cache[cv_id]
            else:
                raise Exception("Could not find a cvterm with id '%s' in the database %s" % (cv_id, self._engine.url))
        else:
            res = self.session.query(Cvterm).filter(Cvterm.cvterm_id == cv_id)
            if not res.count():
                self._cv_id_cache[cv_id] = None
            else:
                self._cv_id_cache[cv_id] = res.one().name

            return self.get_cvterm_name(cv_id)
