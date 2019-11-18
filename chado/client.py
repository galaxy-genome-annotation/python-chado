"""Base chado client
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import re

from chado.exceptions import RecordNotFoundError

from chakin.io import warn

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
        self._analysisfeature_cache = None
        self._analysisprop_cache = None
        self._interpro_cache = None

        self.cache_existing = True

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
            if self.cache_existing:
                res = self.session.query(self.model.feature_dbxref.feature_id, self.model.feature_dbxref.dbxref_id)

                for x in res:
                    if x.feature_id not in self._featxref_cache:
                        self._featxref_cache[x.feature_id] = []
                    self._featxref_cache[x.feature_id].append(x.dbxref_id)

    def _init_feature_cache(self, organism_id, type_id=None, match_on_name=False, force=False):

        if self._feature_cache is not None and force:
            self._feature_cache = None

        if self._feature_cache is None:
            self._feature_cache = {}

            if self.cache_existing:
                # We fill it even with --add_only as we typically already loaded scaffolds before
                res = self.session.query(self.model.feature.feature_id, self.model.feature.name, self.model.feature.uniquename, self.model.feature.type_id, self.model.feature.organism_id) \
                    .filter(self.model.feature.organism_id == organism_id)
                if type_id:
                    res = res.filter(self.model.feature.type_id == type_id)

                if match_on_name:
                    self._feature_cache = {(x.name, x.organism_id, x.type_id): {'feature_id': x.feature_id, 'name': x.name, 'uniquename': x.uniquename} for x in res}
                else:
                    self._feature_cache = {(x.uniquename, x.organism_id, x.type_id): {'feature_id': x.feature_id, 'name': x.name, 'uniquename': x.uniquename} for x in res}

    def _match_feature(self, feature_id, re_name, query_type, organism_id, skip_missing=False):

        seqterm = self.ci.get_cvterm_id(query_type, 'sequence')

        if re_name:
            re_res = re.search(re_name, feature_id)
            if re_res:
                feature_id = re_res.group(1)

        cache_id = (feature_id, organism_id, seqterm)

        if cache_id not in self._feature_cache:
            if skip_missing:
                warn('Could not find feature with name "%s", skipping it', feature_id)
                return None
            else:
                raise RecordNotFoundError('Could not find feature with name "%s"' % feature_id)

        return self._feature_cache[cache_id]['feature_id']

    def _init_analysisfeature_cache(self, analysis_id, force=False):

        if self._analysisfeature_cache is not None and force:
            self._analysisfeature_cache = None

        if self._analysisfeature_cache is None:
            self._analysisfeature_cache = {}
            res = self.session.query(self.model.analysisfeature.analysisfeature_id, self.model.analysisfeature.feature_id, self.model.analysisfeatureprop.type_id, self.model.analysisfeatureprop.rank) \
                              .filter(self.model.analysisfeature.analysis_id == analysis_id) \
                              .outerjoin(self.model.analysisfeatureprop, self.model.analysisfeature.analysisfeature_id == self.model.analysisfeatureprop.analysisfeature_id)
            for x in res:
                if x.feature_id not in self._analysisfeature_cache:
                    self._analysisfeature_cache[x.feature_id] = {'analysisfeature_id': x.analysisfeature_id, 'props': {}}
                # Add analysisfeatureprop in cache too if there are some
                if x.type_id is not None and x.rank is not None:
                    if x.type_id not in self._analysisfeature_cache[x.feature_id]['props']:
                        self._analysisfeature_cache[x.feature_id]['props'][x.type_id] = 0
                    if x.rank > self._analysisfeature_cache[x.feature_id]['props'][x.type_id]:
                        self._analysisfeature_cache[x.feature_id]['props'][x.type_id] = x.rank

    def _init_analysisprop_cache(self, force=False):

        if self._analysisprop_cache is not None and force:
            self._analysisprop_cache = None

        if self._analysisprop_cache is None:
            self._analysisprop_cache = {}
            res = self.session.query(self.model.analysisprop.analysis_id, self.model.analysisprop.type_id, self.model.analysisprop.rank)
            for x in res:
                if x.analysis_id not in self._analysisprop_cache:
                    self._analysisprop_cache[x.analysis_id] = {}
                if x.type_id not in self._analysisprop_cache[x.analysis_id]:
                    self._analysisprop_cache[x.analysis_id][x.type_id] = 0
                if x.rank > self._analysisprop_cache[x.analysis_id][x.type_id]:
                    self._analysisprop_cache[x.analysis_id][x.type_id] = x.rank

    def _add_analysis_feature(self, feature_id, analysis_id, type_id=None, value=None):

        if feature_id not in self._analysisfeature_cache:
            analysis_feature = self.model.analysisfeature()
            analysis_feature.feature_id = feature_id
            analysis_feature.analysis_id = analysis_id
            self.session.add(analysis_feature)
            self.session.flush()
            self.session.refresh(analysis_feature)
            analysisfeature_id = analysis_feature.analysisfeature_id
            self._analysisfeature_cache[feature_id] = {'analysisfeature_id': analysisfeature_id, 'props': {}}
        else:
            analysisfeature_id = self._analysisfeature_cache[feature_id]['analysisfeature_id']

        # Add analysisfeatureprop to
        if type_id and value:
            rank = 0
            if type_id in self._analysisfeature_cache[feature_id]['props']:
                rank = self._analysisfeature_cache[feature_id]['props'][type_id] + 1

            analysis_feature_prop = self.model.analysisfeatureprop()
            analysis_feature_prop.analysisfeature_id = analysisfeature_id
            analysis_feature_prop.type_id = type_id
            # Only works for a specific indent.. might be a way to do it better maybe
            analysis_feature_prop.value = value
            analysis_feature_prop.rank = rank
            self.session.add(analysis_feature_prop)
            self.session.flush()

            self._analysisfeature_cache[feature_id]['props'][type_id] = rank

        return analysisfeature_id

    def _add_analysisprop(self, analysis_id, type_id, value):
        rank = 0
        if analysis_id not in self._analysisprop_cache:
            self._analysisprop_cache[analysis_id] = {}

        if type_id in self._analysisprop_cache[analysis_id]:
            rank = self._analysisprop_cache[analysis_id][type_id] + 1

        analysis_prop = self.model.analysisprop()
        analysis_prop.analysis_id = analysis_id
        analysis_prop.type_id = type_id
        analysis_prop.value = value
        analysis_prop.rank = rank
        self.session.add(analysis_prop)
        self.session.flush()

        self._analysisprop_cache[analysis_id][type_id] = rank

        return analysis_prop

    def _init_featureloc_cache(self, organism_id, force=False):

        if self._featureloc_cache is not None and force:
            self._featureloc_cache = None

        if self._featureloc_cache is None:
            self._featureloc_cache = {}
            if self.cache_existing:
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
            if self.cache_existing:
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
            if self.cache_existing:
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
            if self.cache_existing:
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
            if self.cache_existing:
                res = self.session.query(self.model.feature_cvterm.feature_id, self.model.feature_cvterm.cvterm_id)

                for x in res:
                    if x.feature_id not in self._featcvterm_cache:
                        self._featcvterm_cache[x.feature_id] = []
                    self._featcvterm_cache[x.feature_id].append(x.cvterm_id)

    def _add_feat_cvterm(self, feat, term):
        xref = term.split(':')
        if len(xref) != 2:
            return
        xref_db = xref[0]
        xref_acc = xref[1]
        try:
            term = self.ci.get_cvterm_id(xref_acc, xref_db)
        except RecordNotFoundError:
            term = self.ci.create_cvterm(xref_acc, xref_db, xref_db)

        return self._add_feat_cvterm_with_id(feat, term)

    def _add_feat_cvterm_with_id(self, feat, cvterm_id, pub_id=None):
        if pub_id is None:
            pub_id = self.ci.get_pub_id('null')
        if feat not in self._featcvterm_cache or cvterm_id not in self._featcvterm_cache[feat]:
            cvt2feat = self.model.feature_cvterm()
            cvt2feat.cvterm_id = cvterm_id
            cvt2feat.feature_id = feat
            cvt2feat.pub_id = pub_id
            self.session.add(cvt2feat)
            if feat not in self._featcvterm_cache:
                self._featcvterm_cache[feat] = []
            self._featcvterm_cache[feat].append(cvterm_id)

    def _init_interpro_cache(self, force=False):

        if self._interpro_cache is not None and force:
            self._interpro_cache = None

        if self._interpro_cache is None:
            self._interpro_cache = {}
            if self.cache_existing:
                res = self.session.query(self.model.dbxref.accession, self.model.cvterm.cvterm_id) \
                    .join(self.model.db, self.model.db.db_id == self.model.dbxref.db_id) \
                    .filter(self.model.db.name == "INTERPRO") \
                    .join(self.model.cvterm, self.model.dbxref.dbxref_id == self.model.cvterm.dbxref_id)

                self._interpro_cache = {x.accession: x.cvterm_id for x in res}
