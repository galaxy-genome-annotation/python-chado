"""
Contains possible interactions with the Chado Features
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import csv
import hashlib
import operator
import re
import time
from functools import reduce

from BCBio import GFF

from Bio import Seq, SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature

import chado
from chado.client import Client

from future import standard_library

standard_library.install_aliases()


class FeatureClient(Client):
    """
    Access to the chado features
    """

    def __init__(self, engine, metadata, session, ci):

        self._reset_cache()

        Client.__init__(self, engine, metadata, session, ci)

    def _reset_cache(self):

        self._synonym_cache = None
        self._db_cache = None
        self._xref_cache = None
        self._feature_cache = None
        self._featureloc_cache = None
        self._featrel_cache = None
        self._featured_dirty_rels = None
        self._featxref_cache = None
        self._featcvterm_cache = None
        self._featsyn_cache = None
        self._featureprop_cache = None

    def get_features(self, organism_id=None, analysis_id=None, name=None, uniquename=None):
        """
        Get all or some features

        :type organism_id: str
        :param organism_id: organism_id filter

        :type analysis_id: str
        :param analysis_id: analysis_id filter

        :type name: str
        :param name: name filter

        :type uniquename: str
        :param uniquename: uniquename filter

        :rtype: list of dict
        :return: Features information
        """

        # check if the organism exists
        res = self.session.query(self.model.feature, self.model.analysisfeature.analysis_id)
        if organism_id:
            res = res.filter_by(organism_id=organism_id)
        if name:
            res = res.filter_by(name=name)
        if uniquename:
            res = res.filter_by(uniquename=uniquename)
        res = res.join(self.model.analysisfeature, self.model.analysisfeature.feature_id == self.model.feature.feature_id)
        if analysis_id:
            res = res.filter(self.model.analysisfeature.analysis_id == analysis_id)

        data = []
        for feat in res:
            data.append({
                'feature_id': feat.feature.feature_id,
                'dbxref_id': feat.feature.dbxref_id,
                'organism_id': feat.feature.organism_id,
                'analysis_id': feat.analysis_id,
                'name': feat.feature.name,
                'uniquename': feat.feature.uniquename,
                'residues': feat.feature.residues,
                'seqlen': feat.feature.seqlen,
                'md5checksum': feat.feature.md5checksum,
                'type_id': feat.feature.type_id,
                'is_analysis': feat.feature.is_analysis,
                'is_obsolete': feat.feature.is_obsolete,
                'timeaccessioned': feat.feature.timeaccessioned,
                'timelastmodified': feat.feature.timelastmodified,
            })
        return data

    def delete_features(self, organism_id=None, analysis_id=None, name=None, uniquename=None):
        """
        Get all or some features

        :type organism_id: str
        :param organism_id: organism_id filter

        :type analysis_id: str
        :param analysis_id: analysis_id filter

        :type name: str
        :param name: name filter

        :type uniquename: str
        :param uniquename: uniquename filter

        :rtype: list of dict
        :return: Features information
        """

        # check if the organism exists
        res = self.session.query(self.model.feature)
        if organism_id:
            res = res.filter_by(organism_id=organism_id)
        if name:
            res = res.filter_by(name=name)
        if uniquename:
            res = res.filter_by(uniquename=uniquename)
        if analysis_id:
            res = res.filter(
                self.model.analysisfeature.feature_id == self.model.feature.feature_id,
                self.model.analysisfeature.analysis_id == analysis_id
            )

        res = res.delete(synchronize_session=False)

        self.session.commit()

        self._reset_cache()

        return res

    def load_fasta(self, fasta, organism_id, sequence_type='contig', analysis_id=None,
                   re_name=None, re_uniquename=None, match_on_name=False, update=False, db=None,
                   re_db_accession=None, rel_type=None, re_parent=None, parent_type=None):
        """
        Load features from a fasta file

        :type fasta: str
        :param fasta: Path to the Fasta file to load

        :type organism_id: int
        :param organism_id: Organism ID

        :type sequence_type: str
        :param sequence_type: Sequence type

        :type analysis_id: int
        :param analysis_id: Analysis ID

        :type re_name: str
        :param re_name: Regular expression to extract the feature name from the fasta sequence id (first capturing group will be used).

        :type re_uniquename: str
        :param re_uniquename: Regular expression to extract the feature name from the fasta sequence id (first capturing group will be used).

        :type match_on_name: bool
        :param match_on_name: Match existing features using their name instead of their uniquename

        :type update: bool
        :param update: Update existing feature with new sequence instead of throwing an error

        :type db: int
        :param db: External database to cross reference to.

        :type re_db_accession: str
        :param re_db_accession: Regular expression to extract an external database accession from the fasta sequence id (first capturing group will be used).

        :type rel_type: str
        :param rel_type: Relation type to parent feature ('part_of' or 'derives_from').

        :type re_parent: str
        :param re_parent: Regular expression to extract parent uniquename from the fasta sequence id (first capturing group will be used).

        :type parent_type: str
        :param parent_type: Sequence type of the parent feature

        :rtype: dict
        :return: Number of inserted sequences
        """

        if rel_type and rel_type not in ['part_of', 'derives_from']:
            raise Exception('Unsupported parent relation type (--rel_type)')

        if (db or re_db_accession) and not (db and re_db_accession):
            raise Exception('--db and --re_db_accession should both be specified')

        if (rel_type or re_parent or parent_type) and not (rel_type and re_parent and parent_type):
            raise Exception('--rel_type, --re_parent and --rel_type should all be specified')

        if analysis_id and len(self.ci.analysis.get_analyses(analysis_id=analysis_id)) != 1:
            raise Exception("Could not find analysis with id '{}'".format(analysis_id))

        if len(self.ci.organism.get_organisms(organism_id=organism_id)) != 1:
            raise Exception("Could not find organism with id '{}'".format(organism_id))

        seqterm = self.ci.get_cvterm_id(sequence_type, 'sequence')

        # Cache all possibly existing features
        existing = self.session.query(self.model.feature.feature_id, self.model.feature.name, self.model.feature.uniquename) \
            .filter_by(organism_id=organism_id, type_id=seqterm) \
            .all()
        if match_on_name:
            existing = {ex.name: ex.feature_id for ex in existing}
        else:
            existing = {ex.uniquename: ex.feature_id for ex in existing}

        # Cache all possible parent features
        if parent_type:
            parentterm = self.ci.get_cvterm_id(parent_type, 'sequence')
            relterm = self.ci.get_cvterm_id(rel_type, 'sequence')
            existing_parent = self.session.query(self.model.feature.feature_id, self.model.feature.uniquename) \
                .filter_by(organism_id=organism_id, type_id=parentterm) \
                .all()
            existing_parent = {ex.uniquename: ex.feature_id for ex in existing_parent}

        # Cache all existing dbxref
        existing_dbxref = self.session.query(self.model.dbxref.accession, self.model.dbxref.dbxref_id) \
            .filter_by(db_id=db) \
            .all()
        existing_dbxref = {ex.accession: ex.dbxref_id for ex in existing_dbxref}

        count_ins = 0
        count_up = 0
        for seq_record in SeqIO.parse(fasta, "fasta"):

            # Prepare the dbxref stuff if needed
            dbxref_id = None
            if db and re_db_accession:

                accession = seq_record.id
                re_res = re.search(re_db_accession, accession)
                if re_res:
                    accession = re_res.group(1)

                if accession in existing_dbxref:
                    dbxref_id = existing_dbxref[accession]
                else:
                    dbx = self.model.dbxref()
                    dbx.db_id = db
                    dbx.accession = accession
                    self.session.add(dbx)

                    self.session.flush()
                    self.session.refresh(dbx)

                    dbxref_id = dbx.dbxref_id
                    existing_dbxref[accession] = dbx.dbxref_id

            # Compute md5 checksum
            md5 = hashlib.md5()
            md5.update(str(seq_record.seq).encode('utf-8'))
            md5 = md5.hexdigest()

            # Determine identifiers
            name_ok = seq_record.id
            uname_ok = seq_record.id

            if re_name:
                re_res = re.search(re_name, name_ok)
                if re_res:
                    name_ok = re_res.group(1)

            if re_uniquename:
                re_res = re.search(re_uniquename, uname_ok)
                if re_res:
                    uname_ok = re_res.group(1)

            identifier = name_ok if match_on_name else uname_ok

            # Insert or update
            if identifier not in existing:
                feat = self.model.feature()
                if dbxref_id:
                    feat.dbxref_id = dbxref_id
                feat.organism_id = organism_id
                feat.name = name_ok
                feat.uniquename = uname_ok
                feat.residues = str(seq_record.seq)
                feat.seqlen = len(seq_record)
                feat.md5checksum = md5
                feat.type_id = seqterm
                self.session.add(feat)

                # Add link to analysis
                if analysis_id:
                    afeat = self.model.analysisfeature()
                    afeat.feature = feat
                    afeat.analysis_id = analysis_id
                    self.session.add(afeat)

                existing[identifier] = feat.feature_id

                # Create relationship if needed
                if parent_type:
                    # We won't touch relationship for updated seqs (cause I'm lazy + I don't want to break an already loaded feature)
                    parent_ok = seq_record.id

                    if re_name:
                        re_res = re.search(re_parent, parent_ok)
                        if re_res:
                            parent_ok = re_res.group(1)

                    if parent_ok in existing_parent:
                        featr = self.model.feature_relationship()
                        featr.subject = feat
                        featr.object_id = existing_parent[parent_ok]
                        featr.type_id = relterm
                        self.session.add(featr)
                    else:
                        raise Exception("Could not find a parent feature with uniquename '{}' ('{}').".format(parent_ok, seq_record.id))

                count_ins += 1
            elif update:
                self.session.query(self.model.feature). \
                    filter_by(feature_id=existing[identifier]). \
                    update({
                        'residues': str(seq_record.seq),
                        'seqlen': len(seq_record),
                        'md5checksum': md5
                    })
                count_up += 1
            else:
                raise Exception("Found an existing feature with '{}': '{}' ('{}'). Use --update to update its sequence.".format('name' if match_on_name else 'uniquename', identifier, seq_record.id))

        self.session.commit()

        self._reset_cache()

        return {'inserted': count_ins, 'updated': count_up}

    def load_featureprops(self, tab_file, analysis_id, organism_id, prop_type, feature_type=None, match_on_name=False):
        """
        Load feature properties from a tabular file (Column1: feature name or uniquename, Column2: property value)

        :type tab_file: str
        :param tab_file: Path to the tabular file file to load

        :type analysis_id: int
        :param analysis_id: Analysis ID

        :type organism_id: int
        :param organism_id: Organism ID

        :type prop_type: str
        :param prop_type: Type of the feature property (cvterm will be created if it doesn't exist)

        :type feature_type: str
        :param feature_type: Type of the target features in sequence ontology (will speed up loading if specified)

        :type match_on_name: bool
        :param match_on_name: Match features using their name instead of their uniquename

        :rtype: dict
        :return: Number of inserted featureprop
        """

        if len(self.ci.analysis.get_analyses(analysis_id=analysis_id)) != 1:
            raise Exception("Could not find analysis with id '{}'".format(analysis_id))

        if len(self.ci.organism.get_organisms(organism_id=organism_id)) != 1:
            raise Exception("Could not find organism with id '{}'".format(organism_id))

        # Cache all existing features
        existing = self.session.query(self.model.feature.feature_id, self.model.feature.name, self.model.feature.uniquename) \
            .filter_by(organism_id=organism_id)

        if feature_type:
            seqterm = self.ci.get_cvterm_id(feature_type, 'sequence')
            existing = existing.filter(self.model.feature.type_id == seqterm)

        existing = existing.all()

        if match_on_name:
            existing = {ex.name: ex.feature_id for ex in existing}
        else:
            existing = {ex.uniquename: ex.feature_id for ex in existing}

        count_ins = 0
        with open(tab_file, 'r') as tsvin:
            tsvin = csv.reader(tsvin, delimiter=str("\t"))
            for row in tsvin:
                if len(row) != 2:
                    raise Exception("Malformed input tabular file '{}' (should be a 2-column tab delimited file)".format(tab_file))

                # Insert or update
                if row[0] not in existing:
                    matchon = 'name' if match_on_name else 'uniquename'
                    raise Exception("Could not find a feature with {} '{}'".format(matchon, row[0]))

                if row[1]:
                    self._add_featureprop(organism_id, existing[row[0]], prop_type, row[1])
                    count_ins += 1

        self.session.commit()

        return {'inserted': count_ins}

    def load_gff(self, gff, analysis_id, organism_id, landmark_type=None, re_protein=None, re_protein_capture="^(.*?)$", fasta=None, no_seq_compute=False, quiet=False, add_only=False):
        """
        Load features from a gff file

        :type gff: str
        :param gff: Path to the Fasta file to load

        :type analysis_id: int
        :param analysis_id: Analysis ID

        :type organism_id: int
        :param organism_id: Organism ID

        :type landmark_type: str
        :param landmark_type: Type of the landmarks (will speed up loading if provided, e.g. contig, should be a term of the Sequence ontology)

        :type re_protein: str
        :param re_protein: Replacement string for the protein name using capturing groups defined by --re_protein_capture

        :type re_protein_capture: str
        :param re_protein_capture: Regular expression to capture groups in protein name to use in --re_protein (e.g. "^(.*?)-R([A-Z]+)$", default="^(.*?)$")

        :type fasta: str
        :param fasta: Path to a Fasta containing sequences for some features. When creating a feature, if its sequence is in this fasta file it will be loaded. Otherwise for mRNA and polypeptides it will be computed from the genome sequence (if available), otherwise it will be left empty.

        :type no_seq_compute: bool
        :param no_seq_compute: Disable the computation of mRNA and polypeptides sequences based on genome sequence and positions.

        :type quiet: bool
        :param quiet: Hide progress information

        :type add_only: bool
        :param add_only: Use this flag if you're not updating existing features, but just adding new features to the selected analysis and organism. It will speedup loading, and reduce memory usage, but might produce errors in case of already existing feature.

        :rtype: None
        :return: None
        """

        if len(self.ci.analysis.get_analyses(analysis_id=analysis_id)) != 1:
            raise Exception("Could not find analysis with id '{}'".format(analysis_id))

        if len(self.ci.organism.get_organisms(organism_id=organism_id)) != 1:
            raise Exception("Could not find organism with id '{}'".format(organism_id))

        self.cache_everything = not add_only

        # Get possible landmarks
        landmarks = self.session.query(self.model.feature.name, self.model.feature.uniquename, self.model.feature.feature_id, self.model.feature.type_id, self.model.feature.organism_id) \
            .filter_by(organism_id=organism_id)
        if landmark_type:
            # Filter by landmark type if provided (else we look for all features)
            landmark_type_id = self.ci.get_cvterm_id(landmark_type, 'sequence')
            landmarks = landmarks.filter(self.model.feature.type_id == landmark_type_id)

        self._landmark_cache = {}
        for lm in landmarks:
            if lm.name not in self._landmark_cache:
                self._landmark_cache[lm.name] = []
            if lm.feature_id not in self._landmark_cache[lm.name]:
                self._landmark_cache[lm.name].append(lm.feature_id)  # There may be multiple landmarks with the same name

            # Also look for uniquename
            if lm.uniquename not in self._landmark_cache:
                self._landmark_cache[lm.uniquename] = []
            if lm.feature_id not in self._landmark_cache[lm.uniquename]:
                self._landmark_cache[lm.uniquename].append(lm.feature_id)

        examiner = GFF.GFFExaminer()
        gff_handle = open(gff)
        gff_limits = examiner.available_limits(gff_handle)
        gff_handle.close()

        # Check that we have all the cvterms in the db
        for feat_type in gff_limits['gff_type']:
            # Will raise an exception if not present + keep value in cache
            self.ci.get_cvterm_id(feat_type[0], 'sequence', True)

        # Read optional fasta file
        self._fasta_sequence_cache = {}
        if fasta:
            for record in SeqIO.parse(fasta, "fasta"):
                self._fasta_sequence_cache[record.id] = str(record.seq)

        # Check that all landmarks are there
        for seq_id in gff_limits['gff_id']:
            seq_id = seq_id[0]
            if seq_id not in self._landmark_cache:
                if landmark_type:
                    # Landmark does not exist yet, but we know how to create it
                    lm = SeqFeature(FeatureLocation(0, 1), type=landmark_type, qualifiers={'ID': [seq_id], 'Name': [seq_id]})
                    if seq_id in self._fasta_sequence_cache:
                        added_feat = self._add_feature_with_attr(None, lm, analysis_id, organism_id, have_loc=False, residues=self._fasta_sequence_cache[seq_id])
                    else:
                        added_feat = self._add_feature_with_attr(None, lm, analysis_id, organism_id, have_loc=False)
                    self._landmark_cache[seq_id] = [added_feat['feature_id']]
                else:
                    raise Exception("Could not find landmark named '{}', add --landmark_type to create it".format(seq_id))
            elif len(self._landmark_cache[seq_id]) > 1:
                raise Exception("Found {} landmarks with same name '{}'".format(len(self._landmark_cache[seq_id]), seq_id))

        count_ins = 0

        for rec in GFF.parse(gff):

            # Preload landmark seq to compute some seqs on it
            # We compare to ????... as the gff parser will populate rec.seq with a fake sequence based on the size from "sequence-region" header
            if not no_seq_compute:
                if rec.id in self._fasta_sequence_cache:
                    rec.seq = Seq.Seq(self._fasta_sequence_cache[rec.id])
                    del self._fasta_sequence_cache[rec.id]  # Save a little memory
                elif len(rec.seq) > 0 and str(rec.seq)[0:10] != "??????????":
                    seq_res = self.session.query(self.model.feature.residues) \
                        .filter(self.model.feature.uniquename == rec.id)

                    if landmark_type:
                        seq_res = seq_res.filter(self.model.feature.type_id == landmark_type_id)

                    seq_res = seq_res.all()

                    if len(seq_res) == 1 and seq_res[0].residues:
                        rec.seq = Seq.Seq(seq_res[0].residues)

            # Set a custom attr to store the chado feature_id
            rec._chado_feature_id = self._landmark_cache[rec.id][0]
            if not quiet:
                print("Loading features on {}".format(rec.id))

            for f in rec.features:

                self._load_gff_feature_with_children(rec, f, analysis_id, organism_id, re_protein_capture, re_protein, no_seq_compute=no_seq_compute)
                count_ins += 1

                if not quiet:
                    print("Inserted feature #{}".format(count_ins))

        self._update_rel_ranks()

        self.session.commit()

        self._reset_cache()

        return {'inserted': count_ins}

    def _load_gff_feature_with_children(self, rec, f, analysis_id, organism_id, re_protein_capture, re_protein, parent=None, no_seq_compute=False):

        full_transcript_seq = None
        if f.type == 'mRNA':
            seq_exons = []
            seq_cds = []
            min_cds = None
            max_cds = None

            # To compute mRNA and polypeptide
            for subrna in f.sub_features:
                if subrna.type == 'CDS':
                    seq_cds.append(rec.seq[subrna.location.nofuzzy_start:subrna.location.nofuzzy_end])

                    if min_cds is None or subrna.location.start < min_cds:
                        min_cds = subrna.location.start
                    if max_cds is None or subrna.location.end > max_cds:
                        max_cds = subrna.location.end
                if subrna.type == 'exon':
                    seq_exons.append(rec.seq[subrna.location.nofuzzy_start:subrna.location.nofuzzy_end])

            if not no_seq_compute and len(rec.seq) > 0 and str(rec.seq)[0:10] != "??????????":
                if seq_exons:
                    full_transcript_seq = reduce(operator.add, seq_exons)
                elif seq_cds:
                    full_transcript_seq = reduce(operator.add, seq_cds)
                if f.strand == -1:
                    full_transcript_seq = full_transcript_seq.reverse_complement()

        if full_transcript_seq is not None:
            added_feat = self._add_feature_with_attr(rec, f, analysis_id, organism_id, residues=str(full_transcript_seq), parent=parent)
        else:
            added_feat = self._add_feature_with_attr(rec, f, analysis_id, organism_id, parent=parent)

        mrna_has_polypeptide = False
        for subf in f.sub_features:

            self._load_gff_feature_with_children(rec, subf, analysis_id, organism_id, re_protein_capture, re_protein, parent=added_feat['feature_id'], no_seq_compute=no_seq_compute)

            if f.type == 'mRNA':
                mrna_has_polypeptide = mrna_has_polypeptide or (subf.type == 'polypeptide')

        # Create a polypeptide feature
        if f.type == 'mRNA' and not mrna_has_polypeptide:

            if re_protein:
                pep_uname = re.sub(re_protein_capture, re_protein, added_feat['uniquename'])
            else:
                pep_uname = added_feat['uniquename'] + '-protein'
            polypeptide = SeqFeature(FeatureLocation(min_cds, max_cds), type="polypeptide", strand=f.location.strand, qualifiers={"source": subrna.qualifiers['source'], 'ID': [pep_uname], 'Name': [added_feat['name']]})

            protein_seq = None
            if not no_seq_compute and len(rec.seq) > 0 and str(rec.seq)[0:10] != "??????????":
                full_cds_seq = reduce(operator.add, seq_cds)
                if f.strand == -1:
                    full_cds_seq = full_cds_seq.reverse_complement()
                protein_seq = str(full_cds_seq.translate())

            self._add_feature_with_attr(rec, polypeptide, analysis_id, organism_id, residues=protein_seq, parent=added_feat['feature_id'], parent_rel='derives_from')

    def _add_feature_with_attr(self, rec, f, analysis_id, organism_id, residues=None, dbxref_id=None, is_analysis=None, is_obsolete=None, parent=None, parent_rel='part_of', have_loc=True):

        # Prepare name and uniquename
        if 'ID' in f.qualifiers:
            f_uname = f.qualifiers['ID'][0]
        elif 'Name' in f.qualifiers:
            f_uname = f.qualifiers['Name'][0]
        elif have_loc:
            f_uname = "{}-{}-{}:{}..{}".format(time.time(), f.type, rec.id, f.location.start, f.location.end)
        else:
            f_uname = "{}-{}".format(time.time(), f.type)

        if 'Name' in f.qualifiers:
            f_name = f.qualifiers['Name'][0]
        else:
            f_name = f_uname

        feat_term = self.ci.get_cvterm_id(f.type, 'sequence', True)

        # Fill the existing feature cache if not already done
        if self._feature_cache is None:
            self._feature_cache = {}

            if self.cache_everything:
                # We fill it even with --add_only as we typically already loaded scaffolds before
                res = self.session.query(self.model.feature.feature_id, self.model.feature.name, self.model.feature.uniquename, self.model.feature.type_id, self.model.feature.organism_id) \
                    .filter(self.model.feature.organism_id == organism_id)

                self._feature_cache = {(x.uniquename, x.organism_id, x.type_id): {'feature_id': x.feature_id, 'name': x.name, 'uniquename': x.uniquename} for x in res}

        if self._featureloc_cache is None:
            self._featureloc_cache = {}
            if self.cache_everything:
                res = self.session.query(self.model.feature.feature_id, self.model.featureloc.srcfeature_id, self.model.featureloc.fmin, self.model.featureloc.fmax, self.model.featureloc.strand) \
                    .filter(self.model.feature.organism_id == organism_id) \
                    .join(self.ci.model.featureloc, self.ci.model.featureloc.feature_id == self.ci.model.feature.feature_id)

                for x in res:
                    if x.feature_id not in self._featureloc_cache:
                        self._featureloc_cache[x.feature_id] = []

                    self._featureloc_cache[x.feature_id].append((x.srcfeature_id, x.fmin, x.fmax, x.strand))

        # See if we have a sequence to load from fasta file
        if f_uname in self._fasta_sequence_cache:
            residues = self._fasta_sequence_cache[f_uname]

        feat_uid = (f_uname, organism_id, feat_term)

        if feat_uid in self._feature_cache:
            feat_id = self._feature_cache[feat_uid]['feature_id']
            rank = 0
            if feat_id in self._featureloc_cache:
                rank = len(self._featureloc_cache[feat_id])
        else:
            md5checksum = None
            seqlen = None
            if residues is not None:
                seqlen = len(residues)

                # Compute md5 checksum
                md5checksum = hashlib.md5()
                md5checksum.update(str(residues).encode('utf-8'))
                md5checksum = md5checksum.hexdigest()
            else:
                # We need an md5 even if empty...
                md5checksum = hashlib.md5()
                md5checksum.update("".encode('utf-8'))
                md5checksum = md5checksum.hexdigest()

            score = None
            if 'score' in f.qualifiers and f.qualifiers['score']:
                score = f.qualifiers['score'][0]

            feat = self._add_feature(analysis_id, organism_id, f_uname, feat_term, name=f_name, residues=residues, seqlen=seqlen, md5checksum=md5checksum, dbxref_id=dbxref_id, is_analysis=is_analysis, is_obsolete=is_obsolete, score=score)
            feat_id = feat.feature_id
            rank = 0

            self._feature_cache[feat_uid] = {'feature_id': feat_id, 'name': feat.name, 'uniquename': feat.uniquename}

        if have_loc:
            self._add_featureloc(rec._chado_feature_id, f, feat_id, rank)

        self._load_feat_alias(f, feat_id)

        self._load_feat_dbxref(f, feat_id)

        if 'Gap' in f.qualifiers:
            for gap in f.qualifiers['Gap']:
                self._add_featureprop(organism_id, feat_id, 'Gap', gap)

        if 'Note' in f.qualifiers:
            for gap in f.qualifiers['Note']:
                self._add_featureprop(organism_id, feat_id, 'Note', gap)

        special_qualifiers = [
            'ID',
            'Name',
            'Alias',
            'Parent',
            'Target',
            'Derives_from',
            'Dbxref',
            'Ontology_term',
            'Is_circular',
            'target_organism',
            'target_type',
            'organism',
            'Gap',
            'Note',
            'source',
            'phase',
            'score'
        ]

        for qual in f.qualifiers:
            if qual not in special_qualifiers:
                for gap in f.qualifiers[qual]:
                    self._add_featureprop(organism_id, feat_id, qual, gap)

        if parent:
            self._set_feature_parent(feat_id, parent, parent_rel)

        if 'Target' in f.qualifiers:
            for target in f.qualifiers['Target']:
                self._add_target(feat_id, target)

        if 'Derives_from' in f.qualifiers:
            for parent in f.qualifiers['Derives_from']:
                for x in self._feature_cache:
                    if x[0] == parent:
                        parent = self._feature_cache[x]['feature_id']
                        self._set_feature_parent(feat_id, parent, 'derives_from')
                        break

        return {'feature_id': feat_id, 'name': f_name, 'uniquename': f_uname}

    def _add_featureloc(self, src, f, feat, rank=0):
        phase = None
        if 'phase' in f.qualifiers:
            phase = f.qualifiers['phase'][0]

        self._do_add_featureloc(src, feat, rank, f.location.start, f.location.end, f.location.strand, phase)

    def _do_add_featureloc(self, src, feat, rank, start, end, strand, phase=None):
        loc_hash = (src, start, end, strand)

        if feat not in self._featureloc_cache or loc_hash not in self._featureloc_cache[feat]:
            feat_loc = self.model.featureloc()
            feat_loc.feature_id = feat
            feat_loc.srcfeature_id = src
            feat_loc.fmin = start
            feat_loc.fmax = end
            feat_loc.strand = strand
            feat_loc.rank = rank
            if phase is not None:
                feat_loc.phase = phase

            self.session.add(feat_loc)
            self.session.flush()
            self.session.refresh(feat_loc)

            if feat not in self._featureloc_cache:
                self._featureloc_cache[feat] = []

            self._featureloc_cache[feat].append(loc_hash)

    def _load_feat_alias(self, f, feat):

        if 'Alias' in f.qualifiers:

            exactterm = self.ci.get_cvterm_id('exact', 'synonym_type')
            pub_id = self.ci.get_pub_id('null')

            if self._synonym_cache is None:
                self._synonym_cache = {}
                res = self.session.query(self.model.synonym.name, self.model.synonym.synonym_id) \
                    .filter(self.model.synonym.type_id == exactterm)

                self._synonym_cache = {x.name: x.synonym_id for x in res}

            if self._featsyn_cache is None:
                self._featsyn_cache = {}
                if self.cache_everything:
                    res = self.session.query(self.model.feature_synonym.feature_id, self.model.feature_synonym.synonym_id)

                    for x in res:
                        if x.feature_id not in self._featsyn_cache:
                            self._featsyn_cache[x.feature_id] = []
                        self._featsyn_cache[x.feature_id].append(x.synonym_id)

            for alias in f.qualifiers['Alias']:
                if alias not in self._synonym_cache:
                    syn = self.model.synonym()
                    syn.name = alias
                    syn.type_id = exactterm
                    syn.synonym_sgml = ''
                    self.session.add(syn)

                    self.session.flush()
                    self.session.refresh(syn)

                    self._synonym_cache[alias] = syn.synonym_id

                if feat not in self._featsyn_cache or self._synonym_cache[alias] not in self._featsyn_cache[feat]:
                    syn2feat = self.model.feature_synonym()
                    syn2feat.synonym_id = self._synonym_cache[alias]
                    syn2feat.feature_id = feat
                    syn2feat.pub_id = pub_id
                    self.session.add(syn2feat)

                if feat not in self._featsyn_cache:
                    self._featsyn_cache[feat] = []

                self._featsyn_cache[feat].append(self._synonym_cache[alias])

    def _add_featureprop(self, organism_id, feat, prop, value):

        try:
            propterm = self.ci.get_cvterm_id(prop, 'feature_property')
        except chado.RecordNotFoundError:
            propterm = self.ci.create_cvterm(prop, 'feature_property', 'internal')

        cache_hash = (feat, propterm)
        rank = 0

        if self._featureprop_cache is None:
            self._featureprop_cache = {}
            if self.cache_everything:
                res = self.session.query(self.model.feature.feature_id, self.model.featureprop.type_id, self.model.featureprop.value) \
                    .filter(self.model.feature.organism_id == organism_id) \
                    .join(self.ci.model.featureprop, self.ci.model.featureprop.feature_id == self.ci.model.feature.feature_id)

                for x in res:
                    prehash = (x.feature_id, x.type_id)
                    if prehash not in self._featureprop_cache:
                        self._featureprop_cache[prehash] = []

                    self._featureprop_cache[prehash].append(x.value)

        if cache_hash in self._featureprop_cache:
            if value in self._featureprop_cache[cache_hash]:
                # Don't add two times the same featureprop
                return

            rank = len(self._featureprop_cache[cache_hash])

        prop = self.model.featureprop()
        prop.type_id = propterm
        prop.feature_id = feat
        prop.value = value
        prop.rank = rank
        self.session.add(prop)

        if cache_hash not in self._featureprop_cache:
            self._featureprop_cache[cache_hash] = []
        self._featureprop_cache[cache_hash].append(value)

    def _load_feat_dbxref(self, f, feat):

        if self._db_cache is None:
            self._db_cache = {}
            res = self.session.query(self.model.db.db_id, self.model.db.name)

            self._db_cache = {x.name: x.db_id for x in res}

        if self._xref_cache is None:
            self._xref_cache = {}
            res = self.session.query(self.model.dbxref.accession, self.model.dbxref.dbxref_id, self.model.db.name) \
                .join(self.model.db, self.model.db.db_id == self.model.dbxref.db_id)

            self._xref_cache = {(x.name, x.accession): x.dbxref_id for x in res}

        if self._featxref_cache is None:
            self._featxref_cache = {}
            if self.cache_everything:
                res = self.session.query(self.model.feature_dbxref.feature_id, self.model.feature_dbxref.dbxref_id)

                for x in res:
                    if x.feature_id not in self._featxref_cache:
                        self._featxref_cache[x.feature_id] = []
                    self._featxref_cache[x.feature_id].append(x.dbxref_id)

        if self._featcvterm_cache is None:
            self._featcvterm_cache = {}
            if self.cache_everything:
                res = self.session.query(self.model.feature_cvterm.feature_id, self.model.feature_cvterm.cvterm_id)

                for x in res:
                    if x.feature_id not in self._featcvterm_cache:
                        self._featcvterm_cache[x.feature_id] = []
                    self._featcvterm_cache[x.feature_id].append(x.cvterm_id)

        if 'Dbxref' in f.qualifiers:

            for xref in f.qualifiers['Dbxref']:

                self._add_feat_dbxref(feat, xref)

        if 'source' in f.qualifiers:

            for source in f.qualifiers['source']:

                self._add_feat_dbxref(feat, 'GFF_source:{}'.format(source))

        if 'Ontology_term' in f.qualifiers:

            for term in f.qualifiers['Ontology_term']:

                self._add_feat_cvterm(feat, term)

    def _add_feat_dbxref(self, feat, xref):

        xref = xref.split(':')
        if len(xref) != 2:
            return
        xref_db = xref[0]
        xref_acc = xref[1]

        if xref_db not in self._db_cache:
            db = self.model.db()
            db.name = xref_db
            db.description = 'Added automatically by the GFF loader'
            self.session.add(db)

            self.session.flush()
            self.session.refresh(db)

            self._db_cache[xref_db] = db.db_id

        if (xref_db, xref_acc) not in self._xref_cache:
            dbxref = self.model.dbxref()
            dbxref.db_id = self._db_cache[xref_db]
            dbxref.accession = xref_acc
            dbxref.version = ''
            self.session.add(dbxref)

            self.session.flush()
            self.session.refresh(dbxref)

            self._xref_cache[(xref_db, xref_acc)] = dbxref.dbxref_id

        if feat not in self._featxref_cache or self._xref_cache[(xref_db, xref_acc)] not in self._featxref_cache[feat]:
            dbx2feat = self.model.feature_dbxref()
            dbx2feat.dbxref_id = self._xref_cache[(xref_db, xref_acc)]
            dbx2feat.feature_id = feat
            self.session.add(dbx2feat)

            if feat not in self._featxref_cache:
                self._featxref_cache[feat] = []

            self._featxref_cache[feat].append(self._xref_cache[(xref_db, xref_acc)])

    def _add_feat_cvterm(self, feat, term):

        xref = term.split(':')
        if len(xref) != 2:
            return
        xref_db = xref[0]
        xref_acc = xref[1]

        try:
            term = self.ci.get_cvterm_id(xref_acc, xref_db)
        except chado.RecordNotFoundError:
            term = self.ci.create_cvterm(xref_acc, xref_db, xref_db)

        pub_id = self.ci.get_pub_id('null')

        if feat not in self._featcvterm_cache or term not in self._featcvterm_cache[feat]:
            cvt2feat = self.model.feature_cvterm()
            cvt2feat.cvterm_id = term
            cvt2feat.feature_id = feat
            cvt2feat.pub_id = pub_id
            self.session.add(cvt2feat)

            if feat not in self._featcvterm_cache:
                self._featcvterm_cache[feat] = []

            self._featcvterm_cache[feat].append(term)

    def _add_target(self, feat, target_str):

        target = target_str.split(' ')
        if len(target) != 3 and len(target) != 4:
            raise Exception('Malformed Target value: {}'.format(target_str))

        strand = 1
        if len(target) == 4:
            if target[3] == '+':
                strand = 1
            elif target[3] == '-':
                strand = -1
            else:
                raise Exception('Malformed Target value (bad strand): {}'.format(target_str))

        landmark_str = target[0]
        landmark = None
        start = int(target[1])
        end = int(target[2])
        rank = 0
        if feat in self._featureloc_cache:
            rank = len(self._featureloc_cache[feat])

        for x in self._feature_cache:
            if x[0] == landmark_str:
                landmark = self._feature_cache[x]['feature_id']
                break

        if landmark is None:
            raise Exception('Malformed Target value (unknown target): {}'.format(target_str))

        self._do_add_featureloc(landmark, feat, rank, start, end, strand)

    def _set_feature_parent(self, feat, parent, parent_rel='part_of'):

        partofterm = self.ci.get_cvterm_id('part_of', 'sequence', True)
        reltypeterm = self.ci.get_cvterm_id(parent_rel, 'sequence', True)

        if self._featrel_cache is None:
            self._featrel_cache = {}
            if self.cache_everything:
                # object= parent, subject=child
                res = self.session.query(self.model.feature_relationship.subject_id, self.model.feature_relationship.object_id, self.model.feature_relationship.type_id)
                for x in res:
                    if x.object_id not in self._featrel_cache:
                        self._featrel_cache[x.object_id] = []
                    self._featrel_cache[x.object_id].append((x.subject_id, x.type_id))

        if parent not in self._featrel_cache or (feat, reltypeterm) not in self._featrel_cache[parent]:

            rel = self.model.feature_relationship()
            rel.subject_id = feat
            rel.object_id = parent
            rel.type_id = reltypeterm

            self.session.add(rel)

            if parent not in self._featrel_cache:
                self._featrel_cache[parent] = []
            self._featrel_cache[parent].append((feat, reltypeterm))

            if reltypeterm == partofterm:
                if self._featured_dirty_rels is None:
                    self._featured_dirty_rels = []
                if parent not in self._featured_dirty_rels:
                    self._featured_dirty_rels.append(parent)

    def _update_rel_ranks(self):
        """
        Updates the rank columns in feature_relationship table based on the order of child features
        Only do this for part_of relationships
        """

        if self._featured_dirty_rels:

            partofterm = self.ci.get_cvterm_id('part_of', 'sequence', True)

            for parent in self._featured_dirty_rels:
                if parent in self._featureloc_cache:
                    parent_src = self._featureloc_cache[parent][0][0]
                    children = self._featrel_cache[parent]
                    children_locs = []
                    for x in children:
                        if x[1] == partofterm and x[0] in self._featureloc_cache:
                            for y in self._featureloc_cache[x[0]]:
                                if parent_src == y[0]:
                                    children_locs.append((x[0], y[1]))
                    children_locs = sorted(children_locs, key=lambda x: x[1])

                for rank, child in enumerate(children_locs):
                    self.session.query(self.model.feature_relationship) \
                        .filter(
                            self.model.feature_relationship.object_id == parent,
                            self.model.feature_relationship.subject_id == child[0],
                            self.model.feature_relationship.type_id == partofterm) \
                        .update({"rank": rank})

    def _add_feature(self, analysis_id, organism_id, uniquename, type_id, name=None, residues=None, seqlen=None, md5checksum=None, dbxref_id=None, is_analysis=None, is_obsolete=None, score=None):

        feat = self.model.feature()
        feat.organism_id = organism_id
        feat.uniquename = uniquename
        feat.type_id = type_id
        if name is not None:
            feat.name = name
        if residues is not None:
            feat.residues = residues
        if seqlen is not None:
            feat.seqlen = seqlen
        if md5checksum is not None:
            feat.md5checksum = md5checksum
        if is_analysis is not None:
            feat.is_analysis = is_analysis
        if is_obsolete is not None:
            feat.is_obsolete = is_obsolete

        self.session.add(feat)

        self.session.flush()
        self.session.refresh(feat)

        afeat = self.model.analysisfeature()
        afeat.feature = feat
        afeat.analysis_id = analysis_id
        if score:
            afeat.significance = score
        self.session.add(afeat)

        return feat
