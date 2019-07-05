"""Base chado client
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

from future import standard_library

standard_library.install_aliases()


class Client(object):
    """
    Base client class implementing methods to make queries to the server
    """

    def __init__(self, engine, metadata, session, ci):
        self.engine = engine
        self.metadata = metadata
        self.session = session
        self.ci = ci
        self.model = ci.model

    def _reset_cache(self):
        """
        Reset various caches used by data loaders
        """

        self._db_cache = None
        self._xref_cache = None
        self._featxref_cache = None
        self._feature_cache = None
        self._featureloc_cache = None
        self._synonym_cache = None
        self._featsyn_cache = None
        self._featureprop_cache = None
        self._featrel_cache = None
        self._featcvterm_cache = None
        self._featured_dirty_rels = None

    def _init_db_cache(self, force=False):

        if self._db_cache is not None and force:
            self._db_cache = None

        if self._db_cache is None:
            self._db_cache = {}
            res = self.session.query(self.model.db.db_id, self.model.db.name)

            self._db_cache = {x.name: x.db_id for x in res}

    def _init_xref_cache(self, force=False):

        if self._xref_cache is not None and force:
            self._xref_cache = None

        if self._xref_cache is None:
            self._xref_cache = {}
            res = self.session.query(self.model.dbxref.accession, self.model.dbxref.dbxref_id, self.model.db.name) \
                .join(self.model.db, self.model.db.db_id == self.model.dbxref.db_id)

            self._xref_cache = {(x.name, x.accession): x.dbxref_id for x in res}

    def _init_featxref_cache(self, force=False):

        if self._featxref_cache is not None and force:
            self._featxref_cache = None

        if self._featxref_cache is None:
            self._featxref_cache = {}
            if self.cache_everything:
                res = self.session.query(self.model.feature_dbxref.feature_id, self.model.feature_dbxref.dbxref_id)

                for x in res:
                    if x.feature_id not in self._featxref_cache:
                        self._featxref_cache[x.feature_id] = []
                    self._featxref_cache[x.feature_id].append(x.dbxref_id)

    def _init_feature_cache(self, organism_id, force=False):

        if self._feature_cache is not None and force:
            self._feature_cache = None

        if self._feature_cache is None:
            self._feature_cache = {}

            if self.cache_everything:
                # We fill it even with --add_only as we typically already loaded scaffolds before
                res = self.session.query(self.model.feature.feature_id, self.model.feature.name, self.model.feature.uniquename, self.model.feature.type_id, self.model.feature.organism_id) \
                    .filter(self.model.feature.organism_id == organism_id)

                self._feature_cache = {(x.uniquename, x.organism_id, x.type_id): {'feature_id': x.feature_id, 'name': x.name, 'uniquename': x.uniquename} for x in res}

    def _init_featureloc_cache(self, organism_id, force=False):

        if self._featureloc_cache is not None and force:
            self._featureloc_cache = None

        if self._featureloc_cache is None:
            self._featureloc_cache = {}
            if self.cache_everything:
                res = self.session.query(self.model.feature.feature_id, self.model.featureloc.srcfeature_id, self.model.featureloc.fmin, self.model.featureloc.fmax, self.model.featureloc.strand) \
                    .filter(self.model.feature.organism_id == organism_id) \
                    .join(self.model.featureloc, self.model.featureloc.feature_id == self.model.feature.feature_id)

                for x in res:
                    if x.feature_id not in self._featureloc_cache:
                        self._featureloc_cache[x.feature_id] = []

                    self._featureloc_cache[x.feature_id].append((x.srcfeature_id, x.fmin, x.fmax, x.strand))

    def _init_synonym_cache(self, force=False):

        if self._synonym_cache is not None and force:
            self._synonym_cache = None

        if self._synonym_cache is None:
            self._synonym_cache = {}
            exactterm = self.ci.get_cvterm_id('exact', 'synonym_type')
            res = self.session.query(self.model.synonym.name, self.model.synonym.synonym_id) \
                .filter(self.model.synonym.type_id == exactterm)

            self._synonym_cache = {x.name: x.synonym_id for x in res}

    def _init_featsyn_cache(self, force=False):

        if self._featsyn_cache is not None and force:
            self._featsyn_cache = None

        if self._featsyn_cache is None:
            self._featsyn_cache = {}
            if self.cache_everything:
                res = self.session.query(self.model.feature_synonym.feature_id, self.model.feature_synonym.synonym_id)

                for x in res:
                    if x.feature_id not in self._featsyn_cache:
                        self._featsyn_cache[x.feature_id] = []
                    self._featsyn_cache[x.feature_id].append(x.synonym_id)

    def _init_featureprop_cache(self, organism_id, force=False):

        if self._featureprop_cache is not None and force:
            self._featureprop_cache = None

        if self._featureprop_cache is None:
            self._featureprop_cache = {}
            if self.cache_everything:
                res = self.session.query(self.model.feature.feature_id, self.model.featureprop.type_id, self.model.featureprop.value) \
                    .filter(self.model.feature.organism_id == organism_id) \
                    .join(self.model.featureprop, self.model.featureprop.feature_id == self.model.feature.feature_id)

                for x in res:
                    prehash = (x.feature_id, x.type_id)
                    if prehash not in self._featureprop_cache:
                        self._featureprop_cache[prehash] = []

                    self._featureprop_cache[prehash].append(x.value)

    def _init_featrel_cache(self, force=False):

        if self._featrel_cache is not None and force:
            self._featrel_cache = None

        if self._featrel_cache is None:
            self._featrel_cache = {}
            if self.cache_everything:
                # object= parent, subject=child
                res = self.session.query(self.model.feature_relationship.subject_id, self.model.feature_relationship.object_id, self.model.feature_relationship.type_id)
                for x in res:
                    if x.object_id not in self._featrel_cache:
                        self._featrel_cache[x.object_id] = []
                    self._featrel_cache[x.object_id].append((x.subject_id, x.type_id))

    def _init_featcvterm_cache(self, force=False):

        if self._featcvterm_cache is not None and force:
            self._featcvterm_cache = None

        if self._featcvterm_cache is None:
            self._featcvterm_cache = {}
            if self.cache_everything:
                res = self.session.query(self.model.feature_cvterm.feature_id, self.model.feature_cvterm.cvterm_id)

                for x in res:
                    if x.feature_id not in self._featcvterm_cache:
                        self._featcvterm_cache[x.feature_id] = []
                    self._featcvterm_cache[x.feature_id].append(x.cvterm_id)
