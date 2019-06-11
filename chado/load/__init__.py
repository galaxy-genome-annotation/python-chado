"""
Contains loader methods
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import csv
import os
import re
import tempfile
import xml.etree.ElementTree as ET

import chado
from chado.client import Client

from chakin.io import warn

from future import standard_library

standard_library.install_aliases()


class LoadClient(Client):

    def load_blast(self, analysis_id, blast_output,
                   blast_ext=None, blastdb=None, blastdb_id=None,
                   blast_parameters=None, query_re=None, query_type=None,
                   query_uniquename=False, is_concat=False, search_keywords=False,
                   no_parsed="all"):
        """
        Load a blast analysis

        :type analysis_id: int
        :param analysis_id: Analysis ID

        :type blast_output: str
        :param blast_output: Path to the Blast file to load (single XML file, or directory containing multiple XML files)

        :type blast_ext: str
        :param blast_ext: If looking for files in a directory, extension of the blast result files

        :type blastdb: str
        :param blastdb: Name of the database blasted against (must be in the Chado db table)

        :type blastdb_id: str
        :param blastdb_id: ID of the database blasted against (must be in the Chado db table)

        :type blast_parameters: str
        :param blast_parameters: Blast parameters used to produce these results

        :type query_re: str
        :param query_re: The regular expression that can uniquely identify the query name. This parameters is required if the feature name is not the first word in the blast query name.

        :type query_type: str
        :param query_type: The feature type (e.g. \'gene\', \'mRNA\', \'contig\') of the query. It must be a valid Sequence Ontology term.

        :type query_uniquename: bool
        :param query_uniquename: Use this if the --query-re regular expression matches unique names instead of names in the database.

        :type is_concat: bool
        :param is_concat: If the blast result file is simply a list of concatenated blast results.

        :type search_keywords: bool
        :param search_keywords: Extract keywords for Tripal search

        :type no_parsed: str
        :param no_parsed: Maximum number of hits to parse per feature. Default=all

        :rtype: dict
        :return: Number of processed hits

        """

        if blastdb_id:
            found_db = self.session.query(self.model.db).filter_by(db_id=blastdb_id)
            if not found_db:
                raise Exception("Invalid db ID")
        elif blastdb:
            found_db = self.session.query(self.model.db).filter_by(name=blastdb)
            if not found_db:
                raise Exception("Invalid db name")
            blastdb_id = found_db.one().db_id

        if not blastdb_id:
            raise Exception("Either blastdb or blastdb_id is required")

        res = self.session.query(self.model.analysis).filter_by(analysis_id=analysis_id)
        if not res.count():
            raise Exception("Analysis with the id {} was not found".format(analysis_id))

        if os.path.isfile(blast_output):
            count_ins = self._parse_xml(analysis_id, blastdb, blast_output, no_parsed, blast_ext, query_re, query_type, query_uniquename, is_concat, search_keywords)
            return {'inserted': count_ins}

    def load_go(self, input, organism_id, analysis_id, query_type='polypeptide', match_on_name=False,
                name_column=2, go_column=5, re_name=None, skip_missing=False):
        """
        Load GO annotation from a tabular file

        :type input: str
        :param input: Path to the input tabular file to load

        :type organism_id: int
        :param organism_id: Organism ID

        :type analysis_id: int
        :param analysis_id: Analysis ID

        :type query_type: str
        :param query_type: The feature type (e.g. \'gene\', \'mRNA\', 'polypeptide', \'contig\') of the query. It must be a valid Sequence Ontology term.

        :type match_on_name: bool
        :param match_on_name: Match features using their name instead of their uniquename

        :type name_column: int
        :param name_column: Column containing the feature identifiers (2, 3, 10 or 11; default=2).

        :type go_column: int
        :param go_column: Column containing the GO id (default=5).

        :type re_name: str
        :param re_name: Regular expression to extract the feature name from the input file (first capturing group will be used).

        :type skip_missing: bool
        :param skip_missing: Skip lines with unknown features or GO id instead of aborting everything.

        :rtype: dict
        :return: Number of inserted GO terms
        """

        if analysis_id and len(self.ci.analysis.get_analyses(analysis_id=analysis_id)) != 1:
            raise Exception("Could not find analysis with id '{}'".format(analysis_id))

        if len(self.ci.organism.get_organisms(organism_id=organism_id)) != 1:
            raise Exception("Could not find organism with id '{}'".format(organism_id))

        self.cache_everything = True
        seqterm = self.ci.get_cvterm_id(query_type, 'sequence')

        # Cache all possibly existing features
        existing = self.session.query(self.model.feature.feature_id, self.model.feature.name, self.model.feature.uniquename) \
            .filter_by(organism_id=organism_id, type_id=seqterm) \
            .all()
        if match_on_name:
            existing = {ex.name: ex.feature_id for ex in existing}
        else:
            existing = {ex.uniquename: ex.feature_id for ex in existing}

        # Cache all existing cvterms from GO cv
        db = 'GO'
        self.ci._preload_dbxref2cvterms(db)

        count_ins = 0

        # Cache anaysisfeature content for given analysis_id
        _analysisfeature_cache = []
        res = self.session.query(self.model.analysisfeature.feature_id) \
                          .filter(self.model.analysisfeature.analysis_id == analysis_id)
        for x in res:
            if x.feature_id not in _analysisfeature_cache:
                _analysisfeature_cache.append(x.feature_id)

        # Parse the tab file
        with open(input) as in_gaf:
            rd = csv.reader(in_gaf, delimiter=str("\t"))
            for row in rd:
                if row[0] and row[0][0] in ('!', '#'):
                    # skip header
                    continue

                term = row[go_column - 1]
                term_sp = term.split(':')
                if len(term_sp) != 2:
                    return
                term_db = term_sp[0]
                term_acc = term_sp[1]

                feat_id = row[name_column - 1]
                if re_name:
                    re_res = re.search(re_name, feat_id)
                    if re_res:
                        feat_id = re_res.group(1)

                if feat_id not in existing:
                    if skip_missing:
                        print('Could not find feature with name "%s", skipping it' % feat_id)
                        continue
                    else:
                        raise Exception('Could not find feature with name "%s"' % feat_id)

                try:
                    term_id = self.ci.get_cvterm_id(term_acc, term_db)
                except chado.RecordNotFoundError:
                    term_id = None

                if not term_id:
                    if skip_missing:
                        print('Could not find term with name "%s", skipping it' % term_acc)
                        continue
                    else:
                        raise Exception('Could not find term with name "%s"' % term_acc)

                # Add feature<->cvterm association
                self._populate_featcvterm_cache()
                self._add_feat_cvterm(existing[feat_id], term)

                # Associate the feature to the analysis
                if existing[feat_id] not in _analysisfeature_cache:
                    afeat = self.model.analysisfeature()
                    afeat.feature_id = existing[feat_id]
                    afeat.analysis_id = analysis_id
                    self.session.add(afeat)
                    _analysisfeature_cache.append(existing[feat_id])

                    # Add to analysisfeatureprop too (we're sure it doesn't already exist as we just created the analysisfeature)
                    afeatp = self.model.analysisfeatureprop()
                    afeatp.analysisfeature = afeat
                    afeatp.type_id = term_id
                    afeatp.value = term
                    afeatp.rank = 0
                    self.session.add(afeatp)

                count_ins += 1

        self.session.commit()

        self._reset_cache()

        return {'inserted': count_ins}

    def _parse_xml(self, an_id, blastdb, blast_output, no_parsed, blast_ext, query_re, query_type, query_uniquename, is_concat, search_keywords):

        cv_term_id = self.ci.get_cvterm_id('analysis_blast_output_iteration_hits', 'tripal')
        num_iter = 0
        if(is_concat):
            # File is concatenated, need to break it appart
            try:
                fd, path = tempfile.mkstemp()
                with os.open(blast_output) as in_fh:
                    for line in in_fh:
                        line = line.strip()
                        if not line:
                            continue
                        fd.write(line + "\n")
                        # If we have a full part, process it and delete/recreate temp file
                        if(re.search('</BlastOutput>', line)):
                            fd.close()
                            num_iter += self._parse_xml(an_id, blastdb, path, no_parsed, blast_ext, query_re, query_type, query_uniquename, False, search_keywords)
                            os.remove(path)
                            fd, path = tempfile.mkstemp()
            finally:
                fd.close()
                os.remove(path)
                return

        tree = ET.ElementTree(file=blast_output)

        for iteration in tree.iter(tag="Iteration"):
            self._manage_iteration(iteration, an_id, blastdb, blast_output, no_parsed, blast_ext, query_re, query_type, query_uniquename, is_concat, search_keywords, cv_term_id)
            num_iter += 1
        return num_iter

    def _manage_iteration(self, iteration, an_id, blastdb, blast_output, no_parsed, blast_ext, query_re, query_type, query_uniquename, is_concat, search_keywords, cv_term_id):
        feature_id = 0
        analysis_feature_id = 0
        iteration_tags_xml = ''
        num_hits = 1
        feature = None

        for child in iteration:
            if child.tag == 'Iteration_query-def':
                iteration_tags_xml += "  <{}>{}</{}>\n".format(child.tag, child.text, child.tag)
                if query_re and re.search(query_re, child.text):
                    feature = re.search(query_re, child.text).group(1)
                elif re.search(r'^(.*?)\s.*$', child.text):
                    feature = re.search(r'^(.*?)\s.*$', child.text).group(1)
                else:
                    feature = child.text
                if not feature and query_re:
                    raise Exception("Cannot find feature in {} using the regular expression: {}".format(child.text, query_re))

                # Need to find the feature in chado
                entity_cv_term_id = self.ci.get_cvterm_id(query_type, 'sequence')

                if not entity_cv_term_id:
                    raise Exception("Cannot find cvterm id for query type {}".format(query_type))

                res = self.session.query(self.model.feature).filter_by(type_id=entity_cv_term_id)
                if query_uniquename:
                    res = res.filter_by(uniquename=query_uniquename)
                else:
                    res = res.filter_by(name=feature)

                if not res:
                    raise Exception("Database query failed when searching for feature {}".format(feature))
                if not res.count():
                    raise Exception("Failed: {} cannot find a matching feature in the database".format(feature))
                if res.count() > 1:
                    # Ambiguous : Skip feature. Need to log the result.
                    warn("Ambiguous: %s matches more than one feature and is being skipped", feature)
                    continue
                feature_id = res.one().feature_id

            elif child.tag == 'Iteration_hits':
                if not feature_id:
                    # Skip line
                    warn("Cannot add blast results as feature_id is missing.")
                    continue
                xml_content = "<Iteration>\n{}    <{}>\n".format(iteration_tags_xml, child.tag)
                for hit in child:
                    if hit.tag == "Hit":
                        if (no_parsed == "all" or num_hits <= no_parsed):
                            xml_content += "        <Hit>"
                            xml_content += hit.text + ''.join(ET.tostring(e, encoding="unicode") for e in hit)
                            xml_content += "</Hit>\n"
                    num_hits += 1
                xml_content += "\n  </{}>\n</Iteration>".format(child.tag)

                analysis_feature = self.session.query(self.model.analysisfeature).filter_by(feature_id=feature_id, analysis_id=an_id)
                # Create if not existing
                if not analysis_feature.count():
                    analysis_feature = self.model.analysisfeature()
                    analysis_feature.feature_id = feature_id
                    analysis_feature.analysis_id = an_id
                    self.session.add(analysis_feature)
                    self.session.flush()
                    self.session.refresh(analysis_feature)
                    analysis_feature_id = analysis_feature.analysisfeature_id
                else:
                    analysis_feature_id = analysis_feature.one().analysisfeature_id

                analysis_feature_prop = self.session.query(self.model.analysisfeatureprop) \
                    .join(self.model.analysisfeature, self.model.analysisfeature.analysisfeature_id == self.model.analysisfeatureprop.analysisfeature_id) \
                    .filter(self.model.analysisfeature.feature_id == feature_id, self.model.analysisfeature.analysis_id == an_id, self.model.analysisfeatureprop.type_id == cv_term_id)

                if analysis_feature_prop.count():
                    an_feature_prop_id = analysis_feature_prop.first().analysisfeatureprop_id
                    self.session.query(self.model.analysisfeatureprop).filter_by(analysisfeatureprop_id=an_feature_prop_id).update({'value': xml_content})
                else:
                    analysis_feature_prop = self.model.analysisfeatureprop()
                    analysis_feature_prop.analysisfeature_id = analysis_feature_id
                    analysis_feature_prop.type_id = cv_term_id
                    analysis_feature_prop.value = xml_content
                    analysis_feature_prop.rank = 0
                    self.session.add(analysis_feature_prop)
                    self.session.flush()
                    self.session.refresh(analysis_feature_prop)

                if search_keywords:
                    # remove any existing entries. we'll replace them
                    res = self.session.query(self.model.blast_hit_data).filter_by(analysisfeature_id=analysis_feature_id)
                    res.delete(synchronize_session=False)
                    self.session.commit()

                    db = self.session.query(self.model.db).filter_by(db_id=blastdb)
                    analysis = self.session.query(self.model.analysis).filter_by(analysis_id=an_id)

                    blast_obj = self._get_blast_obj(xml_content, db.one(), feature_id, analysis.one())

                    # iterate through the hits and add the records to the blast_hit_data table
                    index = 1

                    for hit in blast_obj["hits_array"]:
                        blast_org_name = hit["hit_organism"]
                        blast_org_id = None

                        if blast_org_name:
                            res = self.session.query(self.model.blast_organisms).filter_by(blast_org_name=blast_org_name)
                            if not res.count():
                                blast_organism = self.model.blast_organisms()
                                blast_organism.blast_org_name = blast_org_name
                                self.session.add(blast_organism)
                                self.session.flush()
                                self.session.refresh(blast_organism)
                                blast_org_id = blast_organism.blast_org_id
                            else:
                                blast_org_id = res.one().blast_org_id

                        blast_hit_data = self.model.blast_hit_data()
                        blast_hit_data.analysisfeature_id = analysis_feature_id
                        blast_hit_data.analysis_id = an_id
                        blast_hit_data.feature_id = feature_id
                        blast_hit_data.db_id = blastdb
                        blast_hit_data.hit_num = index
                        blast_hit_data.hit_name = hit['hit_name']
                        blast_hit_data.hit_url = hit['hit_url']
                        blast_hit_data.hit_description = hit['description']
                        blast_hit_data.hit_organism = blast_org_name
                        blast_hit_data.blast_org_id = blast_org_id
                        blast_hit_data.hit_accession = hit['accession']
                        blast_hit_data.hit_best_eval = hit['best_evalue']
                        blast_hit_data.hit_best_score = hit['best_score']
                        blast_hit_data.hit_pid = hit['percent_identity']
                        self.session.add(blast_hit_data)
                        self.session.flush()
                        index += 1
            else:
                iteration_tags_xml += "  <{}>{}</{}>\n".format(child.tag, child.text, child.tag)

    def _get_blast_obj(self, xml_content, blast_db, feature_id, blast_analysis):
        blast_object = {}

        blast_object["xml"] = xml_content

        db_name = ""
        is_genbank = ""
        regex_hit_id = ""
        regex_hit_def = ""
        regex_hit_organism = ""
        regex_hit_accession = ""
        db_organism = ""

        parser_query = self.session.query(self.model.tripal_analysis_blast).filter_by(db_id=blast_db.db_id)
        if parser_query.count():
            parser = parser_query.one()
            db_name = parser.displayname
            is_genbank = parser.genbank_style
            regex_hit_id = parser.regex_hit_id
            regex_hit_def = parser.regex_hit_def
            regex_hit_organism = parser.regex_hit_organism
            regex_hit_accession = parser.regex_hit_accession
            db_organism = parser.hit_organism

        # Regex in Python do not take "/"
        if not regex_hit_id:
            regex_hit_id = r'^(.*?)\s.*$'

        if not regex_hit_def:
            regex_hit_def = r'^.*?\s(.*)$'

        if not regex_hit_accession:
            regex_hit_accession = r'^(.*?)\s.*$'

        # get analysis information

        blast_object["analysis"] = blast_analysis
        blast_db.displayname = db_name
        blast_object["db"] = blast_db

        if not db_name:
            blast_object["title"] = blast_analysis.name
        else:
            blast_object["title"] = db_name

        if hasattr(self.model, 'tripal_analysis'):
            res = self.session.query(self.model.tripal_analysis).filter_by(analysis_id=blast_analysis.analysis_id)
            if res.count():
                blast_analysis.nid = res.one().nid

        root = ET.fromstring(xml_content)

        for sub_iteration in root:
            if sub_iteration.tag == 'Iteration_query-def':
                blast_object["xml_tag"] = sub_iteration.text

            if sub_iteration.tag == 'Iteration_hits':
                blast_object["xml_tag"] = sub_iteration.text
                blast_object["feature_id"] = feature_id
                # Initialize hit variable
                # Replaced php array with python list (no int key)
                hits_array = []
                number_hits = 0
                accession = ''
                hit_name = ''
                description = ''
                hit_organism = 'Unknown'
                # Initialize hsp variable

                for hit in sub_iteration:
                    best_evalue = 0
                    best_score = 0
                    best_identity = 0
                    best_len = 0
                    accession = ''
                    hit_name = ''
                    description = ''
                    hit_organism = 'Unknown'

                    for hit_value in hit:
                        # This part was commented out in php version. (https://github.com/tripal/tripal_analysis_blast/blob/7.x-3.x/includes/TripalImporter/BlastImporter.inc#L943)
                        # if hit_value.tag == "Hit_id":
                        # if is_genbank:
                        if hit_value.tag == "Hit_def":
                            if is_genbank:
                                description = hit_value.text
                                regex = re.search(r'^.*\[(.*?)\].*$', description)
                                if regex:
                                    hit_organism = regex.group(1)
                            else:
                                accession = re.search(regex_hit_accession, hit_value.text).group(1) if re.search(regex_hit_accession, hit_value.text) else ''
                                hit_name = re.search(regex_hit_id, hit_value.text).group(1) if re.search(regex_hit_id, hit_value.text) else ''
                                description = re.search(regex_hit_def, hit_value.text).group(1) if re.search(regex_hit_def, hit_value.text) else ''
                                if regex_hit_organism:
                                    hit_organism = re.search(regex_hit_organism, hit_value.text).group(1) if re.search(regex_hit_organism, hit_value.text) else 'Unknown'
                                elif db_organism:
                                    hit_organism = db_organism
                        elif hit_value.tag == "Hit_accession":
                            if is_genbank:
                                accession = hit_value.text
                                hit_name = hit_value.text
                        elif hit_value.tag == "Hit_hsps":
                            # Should be only one child, but iterate anyway
                            for hsp in hit_value:
                                for hsp_type in hsp:
                                    if (best_evalue and best_score and best_len and best_identity):
                                        # No need to go past that, we don't use it
                                        break
                                    if hsp_type.tag == "Hsp_score":
                                        if not best_score:
                                            best_score = hsp_type.text
                                    elif hsp_type.tag == "Hsp_evalue":
                                        if not best_evalue:
                                            best_evalue = hsp_type.text
                                    elif hsp_type.tag == "Hsp_identity":
                                        if not best_identity:
                                            best_identity = hsp_type.text
                                    elif hsp_type.tag == "Hsp_align-len":
                                        if not best_len:
                                            best_len = hsp_type.text
                    # Finished a Hit, saving value
                    number_hits += 1
                    hit_dict = {
                        'accession': accession,
                        'hit_organism': hit_organism,
                        'hit_name': hit_name,
                        'best_evalue': best_evalue,
                        'best_score': best_score,
                        'description': description
                    }

                    if (accession and blast_db.urlprefix):
                        hit_dict['hit_url'] = blast_db.urlprefix + accession
                    else:
                        query = self.session.query(self.model.feature).filter_by(uniquename=hit_name)
                        if query:
                            hit_dict['hit_url'] = "ID" + query.one().feature_id
                        else:
                            hit_dict['hit_url'] = None

                    if best_len:
                        percent_identity = "{0:.2f}".format((float(best_identity) / float(best_len)) * 100)
                        hit_dict['percent_identity'] = percent_identity

                    hits_array.append(hit_dict)

        blast_object['number_hits'] = number_hits
        blast_object['hits_array'] = hits_array
        return blast_object
