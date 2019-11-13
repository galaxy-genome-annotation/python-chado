from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import re
import warnings


from chado.analysis import AnalysisClient
from chado.exceptions import RecordNotFoundError
from chado.export import ExportClient
from chado.expression import ExpressionClient
from chado.feature import FeatureClient
from chado.load import LoadClient
from chado.organism import OrganismClient
from chado.phylogeny import PhylogenyClient
from chado.util import UtilClient

from future import standard_library

from sqlalchemy import BigInteger, Column, Float, ForeignKey, Integer, MetaData, Numeric, String, TIMESTAMP, Text, create_engine
from sqlalchemy import event, exc as sa_exc
from sqlalchemy.ext.automap import automap_base
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, sessionmaker
from sqlalchemy.pool import NullPool

standard_library.install_aliases()


class ChadoModel(object):
    pass


class ChadoInstance(object):

    def __init__(self, dbhost="localhost", dbname="chado", dbuser="chado", dbpass="chado", dbschema="public", dbport=5432, dburl=None, offline=False,
                 no_reflect=False, reflect_tripal_tables=False, pool_connections=True, **kwargs):
        self.dbhost = dbhost
        self.dbname = dbname
        self.dbuser = dbuser
        self.dbpass = dbpass
        self.dbport = dbport
        self.dbschema = dbschema
        self.dburl = dburl

        if self.dburl is not None:
            if '@[' in self.dburl:
                # the url given by pglite is not in the correct format for sqlalchmemy
                self.dburl = re.sub(r"postgres://(.+)@\[(.+)\]/(.+)", r"postgres://\1@/\3?host=\2", self.dburl)
            engine_url = self.dburl
        else:
            engine_url = 'postgresql://%s:%s@%s:%s/%s' % (self.dbuser, self.dbpass, self.dbhost, self.dbport, self.dbname)

        if pool_connections:
            self._engine = create_engine(engine_url)
        else:
            # Prevent SQLAlchemy to make a connection pool.
            # Useful for galaxy dynamic options as otherwise it triggers "sorry, too many clients already" errors after a while
            self._engine = create_engine(engine_url, poolclass=NullPool)

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
                self._reflect_tables_subset()
            else:
                with warnings.catch_warnings():
                    # https://stackoverflow.com/a/5225951
                    warnings.simplefilter("ignore", category=sa_exc.SAWarning)
                    self._reflect_tables(reflect_tripal_tables=reflect_tripal_tables)
            self._mapped = True

        # Initialize Clients
        args = (self._engine, self._metadata, self.session, self)
        self.organism = OrganismClient(*args)
        self.export = ExportClient(*args)
        self.util = UtilClient(*args)
        self.analysis = AnalysisClient(*args)
        self.feature = FeatureClient(*args)
        self.load = LoadClient(*args)
        self.phylogeny = PhylogenyClient(*args)
        self.expression = ExpressionClient(*args)

    def __str__(self):
        return '<ChadoInstance at %s>' % self.dbhost

    def _reflect_tables(self, reflect_tripal_tables=False):
        Base = automap_base()

        class Blast_hit_data(Base):
            __tablename__ = 'blast_hit_data'
            __table_args__ = {'extend_existing': True}

            analysisfeature_id = Column(Integer(), primary_key=True, nullable=False)
            analysis_id = Column(Integer(), primary_key=True, nullable=False)
            feature_id = Column(Integer(), primary_key=True, nullable=False)
            db_id = Column(Integer(), primary_key=True, nullable=False)
            hit_num = Column(Integer(), primary_key=True, nullable=False)
            hit_description = Column(Text())
            hit_name = Column(Text())
            hit_url = Column(Text())
            hit_description = Column(Text())
            hit_organism = Column(Text())
            blast_org_id = Column(Integer())
            hit_accession = Column(Text())
            hit_best_eval = Column(Float())
            hit_best_score = Column(Float())
            hit_pid = Column(Float())

        class Tripal_analysis_blast(Base):
            __tablename__ = 'tripal_analysis_blast'
            __table_args__ = {'schema': 'public'}

            db_id = Column(Integer(), primary_key=True, nullable=False, default=0, index=True)
            regex_hit_id = Column(Text())
            regex_hit_def = Column(Text())
            regex_hit_accession = Column(Text())
            regex_hit_organism = Column(Text())
            hit_organism = Column(Text())
            genbank_style = Column(Integer())

        Base.prepare(self._engine, reflect=True, schema=self.dbschema)
        if reflect_tripal_tables and self.dbschema != "public":
            Base.prepare(self._engine, reflect=True, schema='public')

        # Check for schema name instead of hardcoding?
        self.model = Base.classes
        self.model.blast_hit_data = Blast_hit_data
        self.model.tripal_analysis_blast = Tripal_analysis_blast

        # ambiguous relationships to same table
        self.model.feature_relationship.subject = relationship(self.model.feature, foreign_keys=[self.model.feature_relationship.subject_id], back_populates="subject_in_relationships")
        self.model.feature.subject_in_relationships = relationship(self.model.feature_relationship, foreign_keys=[self.model.feature_relationship.subject_id])
        self.model.feature_relationship.object = relationship(self.model.feature, foreign_keys=[self.model.feature_relationship.object_id], back_populates="object_in_relationships")
        self.model.feature.object_in_relationships = relationship(self.model.feature_relationship, foreign_keys=[self.model.feature_relationship.object_id])

        self.model.featureloc.feature = relationship(self.model.feature, foreign_keys=[self.model.featureloc.feature_id], back_populates="featureloc_collection")
        self.model.feature.featureloc_collection = relationship(self.model.featureloc, foreign_keys=[self.model.featureloc.feature_id], back_populates="feature")
        self.model.featureloc.srcfeature = relationship(self.model.feature, foreign_keys=[self.model.featureloc.srcfeature_id])

    def _reflect_tables_subset(self):

        self.model = ChadoModel()

        Base = declarative_base(metadata=self._metadata)

        class Analysis(Base):
            __tablename__ = "analysis"

            analysis_id = Column(BigInteger(), primary_key=True, nullable=False)
            name = Column(String())
            description = Column(Text())
            program = Column(String())
            programversion = Column(String())
            algorithm = Column(String())
            sourcename = Column(String())
            sourceversion = Column(String())
            sourceuri = Column(Text())
            timeexecuted = Column(TIMESTAMP())

        class Organism(Base):
            __tablename__ = "organism"

            organism_id = Column(BigInteger(), primary_key=True, nullable=False)
            abbreviation = Column(String())
            genus = Column(Text())
            species = Column(String())
            common_name = Column(String())
            infraspecific_name = Column(String())
            type_id = Column(BigInteger())
            comment = Column(Text())

        class Feature(Base):
            __tablename__ = "feature"

            feature_id = Column(BigInteger(), primary_key=True, nullable=False)

        class AnalysisFeature(Base):
            __tablename__ = "analysisfeature"

            analysisfeature_id = Column(BigInteger(), primary_key=True, nullable=False)
            feature_id = Column(BigInteger(), ForeignKey('feature.feature_id'), nullable=False)
            analysis_id = Column(BigInteger(), ForeignKey('analysis.analysis_id'), nullable=False)
            rawscore = Column(Numeric())
            normscore = Column(Numeric())
            significance = Column(Numeric())
            identity = Column(Numeric())

        class Db(Base):
            __tablename__ = "db"

            db_id = Column(BigInteger(), primary_key=True, nullable=False)
            name = Column(String(), nullable=False)
            description = Column(String())
            urlprefix = Column(String())
            url = Column(String())

        self.model.analysis = Analysis
        self.model.analysisfeature = AnalysisFeature
        self.model.organism = Organism
        self.model.feature = Feature
        self.model.db = Db

        # Inspired from https://bitbucket.org/zzzeek/sqlalchemy/wiki/UsageRecipes/ManyToManyOrphan
        # In chado, a feature can be part of multiple analysis, using the analysisfeature table
        # With the following code, when we delete an analysis, it also deletes all features that were only part of this single analysis
        # FIXME This is slow for analysis with many features, it runs 1 delete query per feature
        @event.listens_for(self.session, 'after_flush')
        def delete_feature_orphans(session, ctx):
            # List all features that were attached to a deleted analysis
            if any(isinstance(i, Analysis) for i in session.deleted):
                deleted_ids = []
                for i in session.deleted:
                    if isinstance(i, AnalysisFeature):
                        deleted_ids.append(i.feature_id)
                # Look for features that were part of the deleted analysis and are not part of any other analysis
                query = session.query(Feature).filter(~Feature.analysisfeature.any()).filter(Feature.feature_id.in_(deleted_ids))
                orphans = query.all()
                for orphan in orphans:
                    session.delete(orphan)

        self.model.analysisfeature.analysis = relationship("Analysis")
        self.model.analysisfeature.feature = relationship("Feature")
        self.model.analysis.analysisfeature = relationship("AnalysisFeature", back_populates="analysis", cascade="save-update, merge, delete")
        self.model.feature.analysisfeature = relationship("AnalysisFeature", back_populates="feature")

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
            cvhash_ini = cvhash
            cvhash += '____' + 'syn'

        if cvhash in self._cv_name_cache:
            if self._cv_name_cache[cvhash] is not None:
                return self._cv_name_cache[cvhash]
            else:
                # If we allows synonyms, check the cache for the normal cvhash before throwing.
                if allow_synonyms and cvhash_ini in self._cv_name_cache and self._cv_name_cache[cvhash_ini] is not None:
                    return self._cv_name_cache[cvhash_ini]
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

    def _preload_dbxref2cvterms(self, cv):
        """
        This method caches all cvterms from a db in order to not hit the DB for every
        query.
        """

        res = self.session.query(self.model.cvterm.cvterm_id, self.model.dbxref.accession) \
                          .join(self.model.dbxref, self.model.dbxref.dbxref_id == self.model.cvterm.dbxref_id) \
                          .join(self.model.db, self.model.dbxref.db_id == self.model.db.db_id) \
                          .filter(self.model.db.name == cv) \
                          .all()

        for term in res:
            cvhash = cv + '____' + term.accession
            self._cv_name_cache[cvhash] = term.cvterm_id

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

    def create_cvterm(self, term, cv_name, db_name, term_definition="", cv_definition="", db_definition="", accession=""):

        try:
            cvterm_id = self.get_cvterm_id(term, cv_name, True)
            return cvterm_id
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
            if not accession:
                accession = term

            res = self.session.query(self.model.dbxref).filter_by(accession=accession, db=db)
            if res.count() > 0:
                dbxref = res.one()
            else:
                dbxref = self.model.dbxref()
                dbxref.accession = accession
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
