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
import warnings
import xml.etree.ElementTree as ET

import chado
from chado.client import Client
from chado.exceptions import RecordNotFoundError

from chakin.io import warn

from future import standard_library

from sqlalchemy import Column, Float, ForeignKey, Integer, String, Table, inspect
from sqlalchemy import exc as sa_exc

standard_library.install_aliases()


class LoadClient(Client):

    def __init__(self, engine, metadata, session, ci):

        self._reset_cache()

        Client.__init__(self, engine, metadata, session, ci)

    def blast(self, analysis_id, organism_id, input, blastdb=None, blastdb_id=None,
              re_name=None, query_type="polypeptide", match_on_name=False, skip_missing=False):
        """
        Load a blast analysis, in the same way as does the tripal_analysis_blast module

        :type analysis_id: int
        :param analysis_id: Analysis ID

        :type organism_id: int
        :param organism_id: Organism ID

        :type input: str
        :param input: Path to the Blast XML file to load

        :type blastdb: str
        :param blastdb: Name of the database blasted against (must be in the Chado db table)

        :type blastdb_id: int
        :param blastdb_id: ID of the database blasted against (must be in the Chado db table)

        :type query_type: str
        :param query_type: The feature type (e.g. \'gene\', \'mRNA\', 'polypeptide', \'contig\') of the query. It must be a valid Sequence Ontology term.

        :type match_on_name: bool
        :param match_on_name: Match features using their name instead of their uniquename

        :type re_name: str
        :param re_name: Regular expression to extract the feature name from the input file (first capturing group will be used).

        :type skip_missing: bool
        :param skip_missing: Skip lines with unknown features or GO id instead of aborting everything.

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

        # Cache many things to speed up loading
        self._reset_cache()
        seqterm = self.ci.get_cvterm_id(query_type, 'sequence')
        self._init_feature_cache(organism_id, seqterm, match_on_name)

        self._init_analysisfeature_cache(analysis_id)

        self._init_analysisprop_cache()

        self._hit_details_cache = None

        if not os.path.exists(input):
            raise Exception("{} was not found".format(input))

        self._setup_tables("blast")

        count_ins = self._parse_blast_xml(analysis_id, blastdb_id, input, re_name, query_type, True, organism_id, skip_missing)

        blastdb_ap = self.ci.get_cvterm_id('analysis_blast_blastdb', 'tripal')
        self._add_analysisprop(analysis_id, type_id=blastdb_ap, value=blastdb_id)

        self.session.commit()

        self._reset_cache()

        return {'inserted': count_ins}

    def go(self, input, organism_id, analysis_id, query_type='polypeptide', match_on_name=False,
           name_column=2, go_column=5, re_name=None, skip_missing=False):
        """
        Load GO annotation from a tabular file, in the same way as does the tripal_analysis_go module

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

        seqterm = self.ci.get_cvterm_id(query_type, 'sequence')

        # Cache all possibly existing features
        self._reset_cache()
        self._init_feature_cache(organism_id, seqterm, match_on_name)

        # Cache analysisfeature content for given analysis_id
        self._init_analysisfeature_cache(analysis_id)

        self._init_featcvterm_cache()

        # Cache all existing cvterms from GO cv
        db = 'GO'
        self.ci._preload_dbxref2cvterms(db)

        count_ins = 0

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
                    raise Exception('Malformed term "%s"' % term)
                term_db = term_sp[0]
                term_acc = term_sp[1]

                feat_id = row[name_column - 1]

                feat_id = self._match_feature(feat_id, re_name, query_type, organism_id, skip_missing)
                if skip_missing and feat_id is None:
                    continue

                try:
                    term_id = self.ci.get_cvterm_id(term_acc, term_db)
                except chado.RecordNotFoundError:
                    term_id = None

                if not term_id:
                    if skip_missing:
                        warn('Could not find term with name "%s", skipping it', term_acc)
                        continue
                    else:
                        raise Exception('Could not find term with name "%s"' % term_acc)

                # Add feature<->cvterm association
                self._add_feat_cvterm_with_id(feat_id, term_id)

                # Associate the feature to the analysis
                self._add_analysis_feature(feat_id, analysis_id, term_id, term)

                count_ins += 1

        self.session.commit()

        self._reset_cache()

        return {'inserted': count_ins}

    def interpro(self, analysis_id, organism_id, input, parse_go=False, re_name=None, query_type="polypeptide",
                 match_on_name=False, skip_missing=False):
        """
        Load an InterProScan analysis, in the same way as does the tripal_analysis_intepro module

        :type analysis_id: int
        :param analysis_id: Analysis ID

        :type organism_id: int
        :param organism_id: Organism ID

        :type input: str
        :param input: Path to the InterProScan file to load

        :type parse_go: bool
        :param parse_go: Load GO annotation to the database

        :type query_type: str
        :param query_type: The feature type (e.g. \'gene\', \'mRNA\', \'polypeptide\', \'contig\') of the query. It must be a valid Sequence Ontology term.

        :type match_on_name: bool
        :param match_on_name: Match features using their name instead of their uniquename

        :type re_name: str
        :param re_name: Regular expression to extract the feature name from the input file (first capturing group will be used).

        :type skip_missing: bool
        :param skip_missing: Skip lines with unknown features or GO id instead of aborting everything.

        :rtype: dict
        :return: Number of processed hits

        """

        res = self.session.query(self.model.analysis).filter_by(analysis_id=analysis_id)
        if not res.count():
            raise Exception("Analysis with the id {} was not found".format(analysis_id))

        if len(self.ci.organism.get_organisms(organism_id=organism_id)) != 1:
            raise Exception("Could not find organism with id '{}'".format(organism_id))

        res = self.session.query(self.model.analysisfeature).filter_by(analysis_id=analysis_id)
        if res.count():
            res.delete(synchronize_session=False)

        # Cache all possibly existing features
        self._reset_cache()
        seqterm = self.ci.get_cvterm_id(query_type, 'sequence')
        self._init_feature_cache(organism_id, seqterm, match_on_name)

        # Cache analysisfeature content for given analysis_id
        self._init_analysisfeature_cache(analysis_id)
        self._init_featcvterm_cache()
        self._init_interpro_cache()

        # Cache all existing cvterms from GO cv
        db = 'GO'
        self.ci._preload_dbxref2cvterms(db)

        count_ins = 0
        self._setup_tables("interpro")
        if not os.path.exists(input):
            self.session.rollback()
            raise Exception("{} was not found".format(input))

        count_ins += self._parse_interpro_xml(analysis_id, organism_id, input, parse_go, re_name, query_type, skip_missing)

        self.session.commit()

        self._reset_cache()

        return {'inserted': count_ins}

    def _parse_interpro_xml(self, analysis_id, organism_id, interpro_output, parse_go, re_name, query_type, skip_missing):
        tree = ET.iterparse(interpro_output)
        # If it starts with 'protein-matches' or 'nucleotide-sequence-matches' then this is InterPro v5 XML
        # Need to strip namespace
        for _, elem in tree:
            if '}' in elem.tag:
                elem.tag = elem.tag.split('}', 1)[1]
        root = tree.root
        if re.search("^protein-matches", root.tag) or re.search("^nucleotide-sequence-matches", root.tag):
            counts = self._parse_interpro_xml5(analysis_id, organism_id, root, parse_go, re_name, query_type, skip_missing)
        elif re.search("^EBIInterProScanResults", root.tag) or re.search("^interpro_matches", root.tag):
            counts = self._parse_interpro_xml4(analysis_id, organism_id, root, interpro_output, parse_go, re_name, query_type, skip_missing)
        else:
            raise Exception("Xml format was not recognised")
        return counts

    def _parse_interpro_xml5(self, analysis_id, organism_id, xml, parse_go, re_name, query_type, skip_missing):
        res = self.session.query(self.model.db).filter_by(name="GO")
        if res.count():
            go_db_id = res.one().db_id
        else:
            warn("Goterm loading was requested but the GO schema is not installed in chado, skipping")
            go_db_id = False

        total_count = 0
        for entity in xml:
            total_count += 1
            for child in entity:
                child_name = child.tag
                if child_name == "xref":
                    seq_id = child.get('id')
                    try:
                        feature_id = self._match_feature(seq_id, re_name, query_type, organism_id, skip_missing=False)  # we need to have an exception if it fails
                    except RecordNotFoundError:
                        seq_name = child.get('name', "")
                        feature_id = self._match_feature(seq_name, re_name, query_type, organism_id, skip_missing)
                    if skip_missing and feature_id is None:
                        continue
                    analysisfeature_id = self._add_analysis_feature_ipr(feature_id, analysis_id, entity)
                    if not analysisfeature_id:
                        continue
                    ipr_array = self._parse_feature_xml(entity, feature_id)
                    ipr_terms = ipr_array["iprterms"]
                    self._load_ipr_terms(ipr_terms, feature_id, analysis_id, skip_missing)

                    if parse_go and go_db_id:
                        self._load_go_terms(ipr_array["goterms"], feature_id, analysis_id, go_db_id, skip_missing)
        return total_count

    def _parse_interpro_xml4(self, analysis_id, organism_id, xml, interpro_file, parse_go, re_name, query_type, skip_missing):
        # If there is an EBI header then we need to skip that
        # and set our proteins array to be the second element of the array. This
        # occurs if results were generated with the online InterProScan tool.
        # if the XML starts in with the results then this happens when InterProScan
        # is used command-line and we can just use the object as is
        res = self.session.query(self.model.db).filter_by(name="GO")
        if res.count():
            go_db_id = res.one().db_id
        else:
            warn("Goterm loading was requested but the GO schema is not installed in chado, skipping")
            go_db_id = False

        total_count = 0
        if re.search("^EBIInterProScanResults", xml.tag):
            proteins = xml[1]
        elif re.search("^interpro_matches", xml.tag):
            proteins = xml

        for protein in proteins:
            total_count += 1
            # match the protein id with the feature name
            feature_id = 0
            seqid = protein.get('id')
            # Remove _ORF from the sequence name
            seqid = re.search(r'^(.+)_\d+_ORF\d+.*', seqid).group(1)
            # match the name of the feature in the XML file to a feature in Chado
            feature_id = self._match_feature(seqid, re_name, query_type, organism_id, skip_missing)
            if not feature_id:
                continue
            # Create an entry in the analysisfeature table and add the XML for this feature
            # to the analysisfeatureprop table
            analysisfeature_id = self._add_analysis_feature_ipr(feature_id, analysis_id, protein)
            if not analysisfeature_id:
                continue

            # parse the xml
            ipr_array = self._parse_feature_xml(protein, feature_id)
            ipr_terms = ipr_array['iprterms']
            # Add IPR terms
            self._load_ipr_terms(ipr_terms, feature_id, analysis_id, skip_missing)

            if parse_go and go_db_id:
                self._load_go_terms(ipr_array["goterms"], feature_id, analysis_id, go_db_id, skip_missing)
        return total_count

    def _add_analysis_feature_ipr(self, feature_id, analysis_id, xml):

        type_id = self.ci.get_cvterm_id('analysis_interpro_xmloutput_hit', 'tripal')
        # Only works for a specific indentation, might be a way to do it better maybe
        return self._add_analysis_feature(feature_id, analysis_id, type_id, "   " + ET.tostring(xml).decode())

    def _parse_feature_xml(self, xml, feature_id):
        name = xml.tag
        attrib = xml.attrib

        if name == 'nucleotide-sequence':
            return self._parse_feature_xml5_nucleotide(xml, feature_id)
        if name == 'protein':
            # XML 5 protein key has no attributes, XML4 does
            if len(attrib) == 0:
                return self._parse_feature_xml5_protein(xml, feature_id)
            else:
                return self._parse_feature_xml4(xml, feature_id)

    def _parse_feature_xml5_nucleotide(self, xml, feature_id):
        results = {
            "format": "XML5",
            "iprterms": {},
            "goterms": []
        }
        for child in xml:
            if child.tag == "orf":
                for sub_element in child:
                    if sub_element.tag == "protein":
                        terms = self._parse_feature_xml5_protein(sub_element, feature_id)
                        for ipr_id, iprterm in terms['iprterms'].items():
                            if ipr_id not in results['iprterms']:
                                results['iprterms'][ipr_id] = {}
                            results['iprterms'][ipr_id]['ipr_desc'] = iprterm['ipr_desc']
                            results['iprterms'][ipr_id]['ipr_name'] = iprterm['ipr_name']
                        for goterm_id in terms['goterms']:
                            if goterm_id not in results['goterms']:
                                results['goterms'].append(goterm_id)

        return results

    def _parse_feature_xml5_protein(self, xml, feature_id):
        terms = {
            'format': 'XML5',
            "iprterms": {},
            "goterms": []
        }
        # iterate through each element of the 'protein' children
        for child in xml:
            if child.tag == 'matches':
                for match_element in child:
                    match = {}
                    match["locations"] = {}
                    # sometimes an alignment is made but there is no corresponding IPR term
                    # so we default the match IPR term to 'noIPR'
                    match_ipr_id = 'noIPR'
                    match_ipr_desc = ''
                    match_ipr_name = ''
                    for match_detail in match_element:
                        # the <signature> tag contains information about the match in the
                        # member database (e.g. GENE3D, PFAM, etc).
                        if match_detail.tag == 'signature':
                            # find the IPR term and GO Terms associated with this match
                            for sig_element in match_detail:
                                # Yeah, more loops !
                                # the <entry> tag contains the IPR term entry that corresponds to this match
                                if sig_element.tag == 'entry':
                                    match_ipr_id = sig_element.get('ac')
                                    match_ipr_desc = sig_element.get('desc')
                                    match_ipr_name = sig_element.get('name')
                                    # get the GO terms which are children of the <entry> element
                                    for entry_element in sig_element:
                                        if entry_element.tag == 'go-xref':
                                            go_id = entry_element.get('id')
                                            terms['goterms'].append(go_id)
                    # add this match to the IPR term key to which it is associated
                    if match_ipr_id not in terms['iprterms']:
                        terms['iprterms'][match_ipr_id] = {}

                    terms['iprterms'][match_ipr_id]['ipr_name'] = match_ipr_name
                    terms['iprterms'][match_ipr_id]['ipr_desc'] = match_ipr_desc
        return terms

    def _parse_feature_xml4(self, xml, feature_id):
        terms = {
            'format': 'XML4',
            'iprterms': {},
            'goterms': []
        }
        # iterate through each interpro results for this protein
        for interpro in xml:
            # get the interpro term for this match
            ipr_id = interpro.get('id')
            terms['iprterms'][ipr_id] = {
                'ipr_name': interpro.get('name'),
                'ipr_desc': interpro.get('name'),
            }
            # iterate through the elements of the interpro result
            for level1 in interpro:
                if level1.tag == 'classification':
                    if level1.get('class_type') == "GO":
                        go_id = level1.get('id')
                        terms['goterms'].append(go_id)
        return terms

    def _load_ipr_terms(self, ipr_terms, feature_id, analysis_id, skip_missing):
        for ipr_id, ipr_term in ipr_terms.items():
            if (ipr_term["ipr_name"] and ipr_term["ipr_name"] != 'noIPR'):
                # currently there is no InterPro Ontology OBO file so we can't
                # load the IPR terms that way, we need to just add them
                # as we encounter them. If the term already exists
                # we do not want to update it.

                # Check using IPRnumber (in case ipr_name changed at some point in time)
                if ipr_id in self._interpro_cache:
                    cvterm_id = self._interpro_cache[ipr_id]
                else:
                    cvterm_id = self.ci.create_cvterm(ipr_term['ipr_name'], 'INTERPRO', 'INTERPRO', term_definition=ipr_term['ipr_desc'], accession=ipr_id)
                    if not cvterm_id:
                        if skip_missing:
                            warn('Could not find cvterm %s %s, skipping it', ipr_id, ipr_term['ipr_name'])
                            continue
                        else:
                            raise Exception('Could not find cvterm %s %s' % ipr_id, ipr_term['ipr_name'])
                    self._interpro_cache[ipr_id] = cvterm_id

                # Insert IPR terms into the feature_cvterm table
                # the default pub_id of 1 (NULL) is used. if the cvterm already exists then just skip adding it
                self._add_feat_cvterm_with_id(feature_id, cvterm_id)

                # Insert IPR terms into the analysisfeatureprop table but only if it
                # doesn't already exist
                self._add_analysis_feature(feature_id, analysis_id, cvterm_id, ipr_id)

    def _load_go_terms(self, go_terms, feature_id, analysis_id, go_db_id, skip_missing):
        for go_id in go_terms:
            term = go_id
            term_sp = term.split(':')
            if len(term_sp) != 2:
                self.session.rollback()
                raise Exception("Cannot parse GO term {}".format(go_id))
            term_db = term_sp[0]
            term_acc = term_sp[1]

            try:
                goterm_id = self.ci.get_cvterm_id(term_acc, term_db)
            except chado.RecordNotFoundError:
                goterm_id = None

            if not goterm_id:
                if skip_missing:
                    warn('Could not find term with name "%s", skipping it', term_acc)
                    continue
                else:
                    raise Exception('Could not find term with name "%s"' % term_acc)

            # Insert GO terms into feature_cvterm table. Default pub_id = 1 (NULL) was used. But
            # only insert if not already there
            self._add_feat_cvterm_with_id(feature_id, goterm_id)

            # Insert Go terms into the analysisfeatureprop table but only if it
            # doesn't already exist
            self._add_analysis_feature(feature_id, analysis_id, goterm_id, term_acc)

    def _parse_blast_xml(self, an_id, blastdb_id, blast_output, re_name, query_type, check_concat, organism_id, skip_missing):

        num_iter = 0
        error = False
        if (check_concat):
            # File is concatenated, need to break it appart
            try:
                fd, path = tempfile.mkstemp()
                with open(blast_output) as in_fh:
                    file = open(path, 'a')
                    for line in in_fh:
                        line = line.strip()
                        if not line:
                            continue
                        file.write(line + "\n")
                        # If we have a full part, process it and delete/recreate temp file
                        if (re.search('</BlastOutput>', line)):
                            file.close()
                            num_iter += self._parse_blast_xml(an_id, blastdb_id, path, re_name, query_type, False, organism_id, skip_missing)
                            os.remove(path)
                            fd, path = tempfile.mkstemp()
                            file = open(path, 'a')
            except Exception as e:
                error = e
            finally:
                file.close()
                os.remove(path)
                if error:
                    raise error
                return num_iter
        else:
            tree = ET.ElementTree(file=blast_output)

            for iteration in tree.iter(tag="Iteration"):
                self._manage_iteration(iteration, an_id, blastdb_id, blast_output, re_name, query_type, organism_id, skip_missing)
                num_iter += 1

        return num_iter

    def _manage_iteration(self, iteration, an_id, blastdb_id, blast_output, re_name, query_type, organism_id, skip_missing):
        iteration_tags_xml = ''
        num_hits = 1

        # Need to find the feature in chado
        entity_cv_term_id = self.ci.get_cvterm_id(query_type, 'sequence')
        cv_term_id = self.ci.get_cvterm_id('analysis_blast_output_iteration_hits', 'tripal')

        if not entity_cv_term_id:
            raise Exception("Cannot find cvterm id for query type {}".format(query_type))

        feature_id = None
        feature_txt_id = None

        for child in iteration:
            if child.tag == 'Iteration_query-def':
                iteration_tags_xml += "  <{}>{}</{}>\n".format(child.tag, child.text, child.tag)
                feature_txt_id = child.text
                try:
                    feature_id = self._match_feature(child.text, re_name, query_type, organism_id, skip_missing=False)  # we need to have an exception if it fails
                except RecordNotFoundError:
                    first_word = re.search(r'^(.*?)\s.*$', child.text)
                    if first_word:
                        feature_id = self._match_feature(first_word.group(1), re_name, query_type, organism_id, skip_missing)

            elif child.tag == 'Iteration_hits':
                if feature_id is None:
                    if feature_txt_id is None:
                        raise Exception("No <Iteration_query-def> tag found before <Iteration_hits>, malformed xml while reading: %s" % child)
                    elif skip_missing:
                        continue
                    else:
                        raise Exception("Could not find a feature with id %s" % feature_txt_id)

                xml_content = "<Iteration>\n{}    <{}>\n".format(iteration_tags_xml, child.tag)
                for hit in child:
                    if hit.tag == "Hit":
                        xml_content += "        <Hit>"
                        xml_content += hit.text + ''.join(ET.tostring(e).decode() for e in hit)
                        xml_content += "</Hit>\n"
                    num_hits += 1
                xml_content += "\n  </{}>\n</Iteration>".format(child.tag)

                analysis_feature_id = self._add_analysis_feature(feature_id, an_id, cv_term_id, xml_content)

                # Load hit details in blast_hit_data table
                # remove any existing entries for current feature, we'll replace them
                res = self.session.query(self.model.blast_hit_data).filter_by(analysisfeature_id=analysis_feature_id)
                res.delete(synchronize_session=False)
                self.session.expire_all()

                db = self.session.query(self.model.db).filter_by(db_id=blastdb_id)
                analysis = self.session.query(self.model.analysis).filter_by(analysis_id=an_id)

                hits_details = self._get_hits_details(xml_content, db.one(), feature_id, analysis.one())

                # iterate through the hits and add the records to the blast_hit_data table
                hit_num = 1
                for hit in hits_details:
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
                    blast_hit_data.db_id = blastdb_id
                    blast_hit_data.hit_num = hit_num
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
                    hit_num += 1
            else:
                iteration_tags_xml += "  <{}>{}</{}>\n".format(child.tag, child.text, child.tag)

    def _get_hits_details(self, xml_content, blast_db, feature_id, blast_analysis):
        hits_details = []

        if self._hit_details_cache is None:
            self._hit_details_cache = {}
            self._hit_details_cache[blast_db.db_id] = {}
            self._hit_details_cache[blast_db.db_id]['is_genbank'] = ""
            self._hit_details_cache[blast_db.db_id]['regex_hit_id'] = ""
            self._hit_details_cache[blast_db.db_id]['regex_hit_def'] = ""
            self._hit_details_cache[blast_db.db_id]['regex_hit_organism'] = ""
            self._hit_details_cache[blast_db.db_id]['regex_hit_accession'] = ""
            self._hit_details_cache[blast_db.db_id]['db_organism'] = ""
            parser_query = self.session.query(self.model.tripal_analysis_blast).filter_by(db_id=blast_db.db_id)
            if parser_query.count():
                parser = parser_query.one()
                self._hit_details_cache[blast_db.db_id]['is_genbank'] = parser.genbank_style
                self._hit_details_cache[blast_db.db_id]['regex_hit_id'] = parser.regex_hit_id
                self._hit_details_cache[blast_db.db_id]['regex_hit_def'] = parser.regex_hit_def
                self._hit_details_cache[blast_db.db_id]['regex_hit_organism'] = parser.regex_hit_organism
                self._hit_details_cache[blast_db.db_id]['regex_hit_accession'] = parser.regex_hit_accession
                self._hit_details_cache[blast_db.db_id]['db_organism'] = parser.hit_organism

        is_genbank = self._hit_details_cache[blast_db.db_id]['is_genbank']
        regex_hit_id = self._hit_details_cache[blast_db.db_id]['regex_hit_id']
        regex_hit_def = self._hit_details_cache[blast_db.db_id]['regex_hit_def']
        regex_hit_organism = self._hit_details_cache[blast_db.db_id]['regex_hit_organism']
        regex_hit_accession = self._hit_details_cache[blast_db.db_id]['regex_hit_accession']
        db_organism = self._hit_details_cache[blast_db.db_id]['db_organism']

        # Regex in Python do not take "/"
        if not regex_hit_id:
            regex_hit_id = r'^(.*?)\s.*$'

        if not regex_hit_def:
            regex_hit_def = r'^.*?\s(.*)$'

        if not regex_hit_accession:
            regex_hit_accession = r'^(.*?)\s.*$'

        if hasattr(self.model, 'tripal_analysis'):
            res = self.session.query(self.model.tripal_analysis).filter_by(analysis_id=blast_analysis.analysis_id)
            if res.count():
                blast_analysis.nid = res.one().nid

        root = ET.fromstring(xml_content)

        for sub_iteration in root:
            if sub_iteration.tag == 'Iteration_hits':
                accession = ''
                hit_name = ''
                description = ''
                hit_organism = 'Unknown'

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
                        if query.count():
                            hit_dict['hit_url'] = "ID" + query.one().feature_id
                        else:
                            hit_dict['hit_url'] = None

                    if best_len:
                        percent_identity = "{0:.2f}".format((float(best_identity) / float(best_len)) * 100)
                        hit_dict['percent_identity'] = percent_identity

                    hits_details.append(hit_dict)
        return hits_details

    def _setup_tables(self, module):
        if module == "interpro":
            self.ci.create_cvterm(term='analysis_interpro_xmloutput_hit', term_definition='Hit in the interpro XML output. Each hit belongs to a chado feature. This cvterm represents a hit in the output', cv_name='tripal', db_name='tripal')

        # Term for blast
        if module == "blast":
            self.ci.create_cvterm(term='analysis_blast_output_iteration_hits', term_definition='Hits of a blast', cv_name='tripal', db_name='tripal')
            # Tables for blast
            added_table = False
            if not inspect(self.engine).has_table('tripal_analysis_blast', schema='public'):
                tripal_analysis_blast_table = Table(
                    'tripal_analysis_blast', self.metadata,
                    Column('db_id', Integer, primary_key=True, nullable=False, default=0, index=True),
                    Column('regex_hit_id', String, nullable=True),
                    Column('regex_hit_def', String, nullable=True),
                    Column('regex_hit_accession', String, nullable=True),
                    Column('regex_hit_organism', String, nullable=True),
                    Column('hit_organism', String, nullable=True),
                    Column('genbank_style', Integer, default=0),
                    schema="public"
                )
                tripal_analysis_blast_table.create(self.engine)
                added_table = True

            if not inspect(self.engine).has_table('blast_organisms', schema=self.ci.dbschema):
                blast_organisms_table = Table(
                    'blast_organisms', self.metadata,
                    Column('blast_org_id', Integer, primary_key=True, nullable=False),
                    Column('blast_org_name', String, index=True, unique=True),
                    schema=self.ci.dbschema
                )

                blast_organisms_table.create(self.engine)
                # Needed here for foreign key later
                with warnings.catch_warnings():
                    # https://stackoverflow.com/a/5225951
                    warnings.simplefilter("ignore", category=sa_exc.SAWarning)
                    self.ci._reflect_tables()
                    self.model = self.ci.model

            if not inspect(self.engine).has_table('blast_hit_data', schema=self.ci.dbschema):
                blast_hit_data_table = Table(
                    'blast_hit_data', self.metadata,
                    Column('analysisfeature_id', Integer, ForeignKey(self.model.analysisfeature.analysisfeature_id, ondelete="CASCADE"), nullable=False, index=True, primary_key=True),
                    Column('analysis_id', Integer, ForeignKey(self.model.analysis.analysis_id, ondelete="CASCADE"), nullable=False, index=True, primary_key=True),
                    Column('feature_id', Integer, ForeignKey(self.model.feature.feature_id, ondelete="CASCADE"), nullable=False, index=True, primary_key=True),
                    Column('db_id', Integer, ForeignKey(self.model.db.db_id, ondelete="CASCADE"), nullable=False, index=True, primary_key=True),
                    Column('hit_num', Integer, nullable=False, primary_key=True),
                    Column('hit_name', String, index=True),
                    Column('hit_url', String),
                    Column('hit_description', String),
                    Column('hit_organism', String, index=True),
                    Column('blast_org_id', Integer, ForeignKey(self.model.blast_organisms.blast_org_id, ondelete="CASCADE"), index=True),
                    Column('hit_accession', String, index=True),
                    Column('hit_best_eval', Float, index=True),
                    Column('hit_best_score', Float),
                    Column('hit_pid', Float),
                    schema=self.ci.dbschema
                )

                blast_hit_data_table.create(self.engine)

            with warnings.catch_warnings():
                # https://stackoverflow.com/a/5225951
                warnings.simplefilter("ignore", category=sa_exc.SAWarning)
                self.ci._reflect_tables()
                self.model = self.ci.model

            if added_table:
                # add swissprot, trembl and co
                res = self.session.query(self.model.db.db_id) \
                                  .filter(self.model.db.name.ilike('%swissprot%'))
                for x in res:
                    blast_db_record = self.model.tripal_analysis_blast()
                    blast_db_record.db_id = x.db_id
                    blast_db_record.displayname = 'ExPASy Swissprot'
                    blast_db_record.regex_hit_id = '^sp\\|.*?\\|(.*?)\\s.*?$'
                    blast_db_record.regex_hit_def = '^sp\\|.*?\\|.*?\\s(.*)$'
                    blast_db_record.regex_hit_accession = 'sp\\|(.*?)\\|.*?\\s.*?$'
                    blast_db_record.genbank_style = 0
                    blast_db_record.regex_hit_organism = '^.*?OS=(.*?)\\s\\w\\w=.*$'
                    self.session.add(blast_db_record)
                    self.session.flush()

                res = self.session.query(self.model.db.db_id) \
                                  .filter(self.model.db.name.ilike('%trembl%'))
                for x in res:
                    blast_db_record = self.model.tripal_analysis_blast()
                    blast_db_record.db_id = x.db_id
                    blast_db_record.displayname = 'ExPASy TrEMBL'
                    blast_db_record.regex_hit_id = '^.*?\\|(.*?)\\s.*?$'
                    blast_db_record.regex_hit_def = '^.*?\\|.*?\\s(.*)$'
                    blast_db_record.regex_hit_accession = '^(.*?)\\|.*?\\s.*?$'
                    blast_db_record.genbank_style = 0
                    blast_db_record.regex_hit_organism = '^.*?OS=(.*?)\\s\\w\\w=.*$'
                    self.session.add(blast_db_record)
                    self.session.flush()

                res = self.session.query(self.model.db.db_id) \
                                  .filter(self.model.db.name.ilike('%genbank%'))
                for x in res:
                    blast_db_record = self.model.tripal_analysis_blast()
                    blast_db_record.db_id = x.db_id
                    blast_db_record.displayname = 'Genbank'
                    blast_db_record.regex_hit_id = ''
                    blast_db_record.regex_hit_def = ''
                    blast_db_record.regex_hit_accession = ''
                    blast_db_record.genbank_style = 1
                    blast_db_record.regex_hit_organism = ''
                    self.session.add(blast_db_record)
                    self.session.flush()

            try:
                self.ci.get_cvterm_id('analysis_blast_blastdb', 'tripal')
            except chado.RecordNotFoundError:
                self.ci.create_cvterm('analysis_blast_blastdb', 'tripal', 'tripal')

    def _get_dbs(self):
        """
        Get all dbs

        :rtype: list of dict
        :return: Db information
        """

        res = self.session.query(self.model.db)

        data = []
        for db in res:
            data.append({
                'db_id': db.db_id,
                'name': db.name,
                'description': db.description,
                'urlprefix': db.urlprefix,
            })
        return data
