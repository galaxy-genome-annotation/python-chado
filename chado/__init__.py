from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import warnings

from chado.analysis import AnalysisClient
from chado.export import ExportClient
from chado.feature import FeatureClient
from chado.organism import OrganismClient
from chado.phylogeny import PhylogenyClient
from chado.util import UtilClient

from future import standard_library

from sqlalchemy import MetaData, create_engine
from sqlalchemy import exc as sa_exc
from sqlalchemy.ext.automap import automap_base
from sqlalchemy.orm import sessionmaker

standard_library.install_aliases()


class ChadoInstance(object):

    def __init__(self, dbhost="localhost", dbname="chado", dbuser="chado", dbpass="chado", dbschema="public", dbport=5432, offline=False, **kwargs):
        self.dbhost = dbhost
        self.dbname = dbname
        self.dbuser = dbuser
        self.dbpass = dbpass
        self.dbport = dbport
        self.dbschema = dbschema

        self._engine = create_engine('postgresql://%s:%s@%s:%s/%s' % (self.dbuser, self.dbpass, self.dbhost, self.dbport, self.dbname))
        self._metadata = MetaData(self._engine, schema=self.dbschema)
        Session = sessionmaker(bind=self._engine)
        self.session = Session()
        if not offline:
            self._test_db_access()

        self._cv_id_cache = {}
        self._cv_name_cache = {}
        self._mapped = False
        self.model = None

        if not self._mapped and not offline:
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
        self.feature = FeatureClient(*args)
        self.phylogeny = PhylogenyClient(*args)

    def __str__(self):
        return '<ChadoInstance at %s>' % self.dbhost

    def _reflect_tables(self):
        Base = automap_base()
        Base.prepare(self._engine, reflect=True, schema=self.dbschema)
        self.model = Base.classes

    def _test_db_access(self):
        tables = self._engine.table_names(schema=self.dbschema)
        if ('analysis' not in tables or 'feature' not in tables):
            raise Exception("Could not find Chado tables in db %s" % (self._engine.url))

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
            res = self.session.query(self.model.cvterm).filter(self.model.cvterm.cvterm_id == cv_id)
            if not res.count():
                self._cv_id_cache[cv_id] = None
            else:
                self._cv_id_cache[cv_id] = res.one().name

            return self.get_cvterm_name(cv_id)

    def get_cvterm_id(self, name, cv):
        """
        get_cvterm_id allows lookup of CV terms by their name.
        This method caches the result in order to not hit the DB for every
        query. Maybe should investigate pre-loading popular terms? (E.g. gene,
        mRNA, etc)
        """
        cvhash = cv + '____' + name
        if cvhash in self._cv_name_cache:
            if self._cv_name_cache[cvhash] is not None:
                return self._cv_name_cache[cvhash]
            else:
                raise Exception("Could not find a cvterm with name '%s' from cv '%s' in the database %s" % (name, cv, self._engine.url))
        else:
            res = self.session.query(self.model.cvterm) \
                .filter(self.model.cvterm.name == name) \
                .join(self.model.cv, self.model.cv.cv_id == self.model.cvterm.cv_id) \
                .filter_by(name=cv)

            if not res.count():
                self._cv_name_cache[cvhash] = None
            else:
                self._cv_name_cache[cvhash] = res.one().cvterm_id

            return self.get_cvterm_id(name, cv)
