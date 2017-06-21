from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from sqlalchemy import create_engine, MetaData, Table
import warnings
from sqlalchemy import exc as sa_exc
from sqlalchemy.orm import mapper, sessionmaker
from chado.organism import OrganismClient
from chado.export import ExportClient
from chado.util import UtilClient
from chado.analysis import AnalysisClient
from chado.models import *


class ChadoInstance(object):

    def __init__(self, dbhost="localhost", dbname="chado", dbuser="chado", dbpass="chado", dbschema="public", dbport=5432, debug=False, **kwargs):
        self.dbhost = dbhost
        self.dbname = dbname
        self.dbuser = dbuser
        self.dbpass = dbpass
        self.dbport = dbport
        self.dbschema = dbschema
        self.debug = debug

        self._engine = create_engine('postgresql://%s:%s@%s:%s/%s' % (self.dbuser, self.dbpass, self.dbhost, self.dbport, self.dbname))
        self._metadata = MetaData(self._engine, schema=self.dbschema)
        Session = sessionmaker(bind=self._engine)
        self.session = Session()
        self._test_db_access()

        self._cv_id_cache = {}
        self._mapped = False

        if not self._mapped:
            with warnings.catch_warnings():
                # https://stackoverflow.com/a/5225951
                warnings.simplefilter("ignore", category=sa_exc.SAWarning)
                self._reflect_tables()
                self._mapped = True

        # Initialize Clients
        args = (self._engine, self._metadata, self.session, self)
        self.organism = OrganismClient(*args)
        self.export = ExportClient(*args)
        self.util = UtilClient(*args)
        self.analysis = AnalysisClient(*args)

    def __str__(self):
        return '<ChadoInstance at %s>' % self.dbhost

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

    def _test_db_access(self):
        tables = self._engine.table_names(schema=self.dbschema)
        if ('analysis' not in tables or 'feature' not in tables):
            raise Exception("Could not find Chado tables in db %s" % (self._engine.url))


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
