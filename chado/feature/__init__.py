"""
Contains possible interactions with the Chado Features
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import hashlib
import re

from Bio import SeqIO

from chado.client import Client

from future import standard_library

standard_library.install_aliases()


class FeatureClient(Client):
    """
    Access to the chado features
    """

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

        :rtype: None
        :return: None
        """

        if rel_type and rel_type not in ['part_of', 'derives_from']:
            raise Exception('Unsupported parent relation type (--rel_type)')

        if (db or re_db_accession) and not (db and re_db_accession):
            raise Exception('--db and --re_db_accession should both be specified')

        if (rel_type or re_parent or parent_type) and not (rel_type and re_parent and parent_type):
            raise Exception('--rel_type, --re_parent and --rel_type should all be specified')

        seqterm = self.ci.get_cvterm_id(sequence_type, 'sequence')

        # Cache all possibly existing features
        existing = self.session.query(self.model.feature) \
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
            existing_parent = self.session.query(self.model.feature) \
                .filter_by(organism_id=organism_id, type_id=parentterm) \
                .all()
            existing_parent = {ex.uniquename: ex.feature_id for ex in existing_parent}

        # Cache all existing dbxref
        existing_dbxref = self.session.query(self.model.dbxref) \
            .filter_by(db_id=db) \
            .all()
        existing_dbxref = {ex.accession: ex.dbxref_id for ex in existing_dbxref}

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

                    dbxref_id = dbx
                    existing_dbxref[accession] = dbx

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
                    if isinstance(dbxref_id, int):  # I'm not proud of this
                        feat.dbxref_id = dbxref_id
                    else:
                        feat.dbxref = dbxref_id
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
            elif update:
                self.session.query(self.model.feature). \
                    filter_by(feature_id=existing[identifier]). \
                    update({
                        'residues': str(seq_record.seq),
                        'seqlen': len(seq_record),
                        'md5checksum': md5
                    })
            else:
                raise Exception("Found an existing feature with '{}': '{}' ('{}'). Use --update to update its sequence.".format('name' if match_on_name else 'uniquename', identifier, seq_record.id))

        self.session.commit()

    def load_gff(self, gff, analysis_id, organism_id, sequence_type):
        """
        Load features from a gff file

        :type gff: str
        :param gff: Path to the Fasta file to load

        :type analysis_id: int
        :param analysis_id: Analysis ID

        :type organism_id: int
        :param organism_id: Organism ID

        :type sequence_type: str
        :param sequence_type: Sequence type

        :rtype: None
        :return: None
        """
        raise NotImplementedError()
