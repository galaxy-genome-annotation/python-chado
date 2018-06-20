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

from sqlalchemy import BigInteger, Column, MetaData, String, TIMESTAMP, Table, Text, create_engine
from sqlalchemy import exc as sa_exc
from sqlalchemy.ext.automap import automap_base
from sqlalchemy.orm import relationship, sessionmaker

standard_library.install_aliases()


class RecordNotFoundError(Exception):
    """Raised when a db select failed."""


class ChadoModel(object):
    pass


class ChadoInstance(object):

    def __init__(self, dbhost="localhost", dbname="chado", dbuser="chado", dbpass="chado", dbschema="public", dbport=5432, offline=False, no_reflect=False, **kwargs):
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
        self._pub_id_cache = {}
        self._mapped = False
        self.model = None

        if not self._mapped and not offline:

            if no_reflect:
                # No need to do a full reflection of all tables for simple operations
                self.model = ChadoModel()

                self.model.analysis = Table('analysis', self._metadata,
                                            Column('analysis_id', BigInteger(), primary_key=True, nullable=False),
                                            Column('name', String()),
                                            Column('description', Text()),
                                            Column('program', String()),
                                            Column('programversion', String()),
                                            Column('algorithm', String()),
                                            Column('sourcename', String()),
                                            Column('sourceversion', String()),
                                            Column('sourceuri', Text()),
                                            Column('timeexecuted', TIMESTAMP()),
                                            )

                self.model.organism = Table('organism', self._metadata,
                                            Column('organism_id', BigInteger(), primary_key=True, nullable=False),
                                            Column('abbreviation', String()),
                                            Column('genus', Text()),
                                            Column('species', String()),
                                            Column('common_name', String()),
                                            Column('infraspecific_name', String()),
                                            Column('type_id', BigInteger()),
                                            Column('comment', Text()),
                                            )
            else:
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

        # ambiguous relationships to same table
        self.model.feature_relationship.subject = relationship("feature", foreign_keys=[self.model.feature_relationship.subject_id], back_populates="subject_in_relationships")
        self.model.feature.subject_in_relationships = relationship("feature_relationship", foreign_keys=[self.model.feature_relationship.subject_id])
        self.model.feature_relationship.object = relationship("feature", foreign_keys=[self.model.feature_relationship.object_id], back_populates="object_in_relationships")
        self.model.feature.object_in_relationships = relationship("feature_relationship", foreign_keys=[self.model.feature_relationship.object_id])

        self.model.featureloc.feature = relationship("feature", foreign_keys=[self.model.featureloc.feature_id], back_populates="featureloc_collection")
        self.model.feature.featureloc_collection = relationship("featureloc", foreign_keys=[self.model.featureloc.feature_id], back_populates="feature")
        self.model.featureloc.srcfeature = relationship("feature", foreign_keys=[self.model.featureloc.srcfeature_id])

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
                raise RecordNotFoundError("Could not find a cvterm with id '%s' in the database %s" % (cv_id, self._engine.url))
        else:
            res = self.session.query(self.model.cvterm.name).filter(self.model.cvterm.cvterm_id == cv_id)
            if not res.count():
                self._cv_id_cache[cv_id] = None
            else:
                self._cv_id_cache[cv_id] = res.one().name

            return self.get_cvterm_name(cv_id)

    def get_cvterm_id(self, name, cv, allow_synonyms=False):
        """
        get_cvterm_id allows lookup of CV terms by their name.
        This method caches the result in order to not hit the DB for every
        query. Maybe should investigate pre-loading popular terms? (E.g. gene,
        mRNA, etc)
        """
        cvhash = cv + '____' + name
        if allow_synonyms:
            cvhash += '____' + 'syn'

        if cvhash in self._cv_name_cache:
            if self._cv_name_cache[cvhash] is not None:
                return self._cv_name_cache[cvhash]
            else:
                raise RecordNotFoundError("Could not find a cvterm with name '%s' from cv '%s' in the database %s" % (name, cv, self._engine.url))
        else:
            res = self.session.query(self.model.cvterm.cvterm_id) \
                .join(self.model.cv, self.model.cv.cv_id == self.model.cvterm.cv_id)

            if allow_synonyms:
                res = res.join(self.model.cvtermsynonym, self.model.cvtermsynonym.cvterm_id == self.model.cvterm.cvterm_id, isouter=True) \
                    .filter((self.model.cvterm.name == name) | (self.model.cvtermsynonym.synonym == name))
            else:
                res = res.filter(self.model.cvterm.name == name)

            res = res.filter(self.model.cv.name == cv)

            if not res.count():
                self._cv_name_cache[cvhash] = None
            else:
                self._cv_name_cache[cvhash] = res[0].cvterm_id

            return self.get_cvterm_id(name, cv, allow_synonyms)

    def get_pub_id(self, name):
        """
        Allows lookup of publication by their uniquename.
        This method caches the result in order to not hit the DB for every
        query.
        """
        if name in self._pub_id_cache:
            if self._pub_id_cache[name] is not None:
                return self._pub_id_cache[name]
            else:
                raise Exception("Could not find a pub with uniquename '%s' in the database %s" % (name, self._engine.url))
        else:
            res = self.session.query(self.model.pub.pub_id) \
                .filter(self.model.pub.uniquename == name)

            if not res.count():
                self._pub_id_cache[name] = None
            else:
                self._pub_id_cache[name] = res[0].pub_id

            return self.get_pub_id(name)

    def create_cvterm(self, term, cv_name, db_name, term_definition="", cv_definition="", db_definition=""):

        try:
            cvterm = self.get_cvterm_id(term, cv_name, True)
        except RecordNotFoundError:

            # Not found, we need to create it
            # check if the db exists
            res = self.session.query(self.model.db).filter_by(name=db_name)

            if res.count() > 0:
                db = res.one()
            else:
                db = self.model.db()
                db.name = db_name
                db.definition = db_definition

                self.session.add(db)

            # check if the cv exists
            res = self.session.query(self.model.cv).filter_by(name=cv_name)

            if res.count() > 0:
                cv = res.one()
            else:
                cv = self.model.cv()
                cv.name = cv_name
                cv.definition = cv_definition

                self.session.add(cv)

            # Cvterm not found, create it
            res = self.session.query(self.model.dbxref).filter_by(accession=term, db=db)
            if res.count() > 0:
                dbxref = res.one()
            else:
                dbxref = self.model.dbxref()
                dbxref.accession = term
                dbxref.db = db

                self.session.add(dbxref)

            cvterm = self.model.cvterm()
            cvterm.name = term
            cvterm.cv = cv
            cvterm.definition = term_definition
            cvterm.dbxref = dbxref

            self.session.add(cvterm)

            self.session.flush()
            self.session.refresh(cvterm)

            cvhash = cv.name + '____' + term
            self._cv_name_cache[cvhash] = cvterm.cvterm_id

        return cvterm.cvterm_id
