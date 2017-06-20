"""
Export data from chado
"""
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from builtins import open
from future import standard_library
standard_library.install_aliases()
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation as BioFeatureLocation
from chado.client import Client
from chado.models import *


class ExportClient(Client):

    def export_fasta(self, organism_id, file=False):
        """
        Export reference sequences as fasta.

        :type organism_id: int
        :param organism_id: Organism ID

        :type file: bool
        :param file: If true, write to files in CWD

        :rtype: None
        :return: None
        """

        # check if the organism exists
        res = self.session.query(Organism) \
            .filter(Organism.organism_id.in_([organism_id]))

        for org in res:
            sequence_features = self.session.query(Feature) \
                .filter_by(organism_id=org.organism_id) \
                .filter(Feature.seqlen > 0)

            if file:
                output = open('{0.organism_id}.{0.genus}.{0.species}-{0.common_name}.fa'.format(org), 'w')
            else:
                output = sys.stdout

            for seq in sequence_features:
                output.write('>{0.uniquename} [ID={0.feature_id}]'.format(seq))
                output.write('\n')
                output.write(seq.residues)
                output.write('\n')

            if file:
                output.close()

    def export_gff3(self, organism_id, file=False):
        """
        Export organism features as GFF3

        :type organism_id: int
        :param organism_id: Organism ID

        :rtype: None
        :return: None
        """
        # check if the organism exists
        res = self.session.query(Organism).filter(Organism.organism_id.in_([organism_id]))
        sys.stderr.write("Processing %s sequences\n" % res.count())

        for org in res:
            # TODO: can we do this properly?
            seq = Seq("A" * 1, IUPAC.unambiguous_dna)

            # Annotation features
            features = self.session.query(Feature, FeatureLocation) \
                .filter_by(organism_id=org.organism_id) \
                .filter(Feature.seqlen==None) \
                .join(FeatureLocation, Feature.feature_id == FeatureLocation.feature_id, isouter=True)
            sys.stderr.write("\tProcessing %s features\n" % features.count())

            biopy_features = {}
            for idx, (feature, featureloc) in enumerate(features):
                if idx % 5000 == 0:
                    sys.stderr.write("\t%s / %s\n" % (idx, features.count()))

                #[u'dbxref_id', u'feature_id', u'is_analysis', u'is_obsolete',
                # u'md5checksum', u'name', u'organism_id', u'residues', u'seqlen',
                # u'timeaccessioned', u'timelastmodified', u'type_id',
                # u'uniquename']
                #[u'feature_id', u'featureloc_id', u'fmax', u'fmin',
                # u'is_fmax_partial', u'is_fmin_partial', u'locgroup', u'phase',
                # u'rank', u'residue_info', u'srcfeature_id', u'strand']
                qualifiers = {
                    self.ci.get_cvterm_name(prop.type_id): prop.value for prop in
                    self.session.query(FeatureProperties).filter_by(feature_id=feature.feature_id).all()
                }

                qualifiers['ID'] = feature.uniquename

                biopy_features[feature.feature_id] = SeqFeature(
                    BioFeatureLocation(featureloc.fmin, featureloc.fmax),
                    id=feature.uniquename,
                    type=self.ci.get_cvterm_name(feature.type_id),
                    strand=featureloc.strand,
                    qualifiers=qualifiers
                )

        #res = self.session.query(Organism).filter(Organism.organism_id.in_(organism_id))
            relationships = self.session.query(FeatureRelationship) \
                .filter(FeatureRelationship.subject_id.in_(biopy_features.keys()))
            sys.stderr.write("\tProcessing %s relationships\n" % relationships.count())

     #feature_relationship_id | subject_id | object_id | type_id | value | rank
    #-------------------------+------------+-----------+---------+-------+------
                           #1 |          4 |         3 |      37 |       |    0
                           #2 |          5 |         4 |      37 |       |    0
                           #3 |          6 |         4 |      37 |       |    0
                           #4 |          7 |         4 |      37 |       |    0

            features = []

            def findById(feature_list, id):
                for feature in feature_list:
                    if feature.id == id:
                        yield feature

                    if hasattr(feature, 'sub_features'):
                        for x in findById(feature.sub_features, id):
                            yield x

            # Now to re-parent things properly
            # This is BioPython 1.67 ONLY since they broke compatability with bcbio-gff
            # https://github.com/biopython/biopython/issues/928
            for idx, rel in enumerate(relationships):
                if idx % 5000 == 0:
                    sys.stderr.write("\t%s / %s\n" % (idx, relationship.count()))

                term = self.ci.get_cvterm_name(rel.type_id)
                if term != 'part_of':
                    sys.stderr.write("\tCannot handle non-part_of relationships (%s %s %s)\n" % (rel.subject_id, term, rel.object_id))
                    continue

                # Try and find the features in features.
                child = list(findById(features, biopy_features[rel.subject_id].id))
                parent = list(findById(features, biopy_features[rel.object_id].id))

                assert len(child) <= 1
                assert len(parent) <= 1
                alreadyProcessedParent = False
                alreadyProcessedChild = False

                # If they aren't there, pull them from the complete set.
                if len(child) == 0:
                    child = biopy_features[rel.subject_id]
                else:
                    child = child[0]
                    alreadyProcessedChild = True

                if len(parent) == 0:
                    parent = biopy_features[rel.object_id]
                else:
                    parent = parent[0]
                    alreadyProcessedParent = True

                if not hasattr(parent, 'sub_features'):
                    parent.sub_features = []

                parent.sub_features.append(child)
                if alreadyProcessedChild and alreadyProcessedParent:
                    # Here we've seen both (they're BOTH in the list), so we need to remove
                    # child and not touch parent since we added to parent already
                    features.remove(child)
                elif alreadyProcessedChild and not alreadyProcessedParent:
                    # Here our child is already in features, so we need to remove it from
                    # the feature set, add to the parent (done) and re-place in features.
                    features.remove(child)
                    features.append(parent)
                elif not alreadyProcessedChild and alreadyProcessedParent:
                    # In this case we've seen the parent before, already in list, no need to do anything
                    # features.append(parent)
                    pass
                else:
                    # Otherwise, completely new feature.
                    features.append(parent)

            n = org.common_name if org.common_name else 'org_%s' % org.organism_id
            record = SeqRecord(
                seq, id=n, name=n,
                description="%s %s" % (org.genus, org.species),
            )
            record.features = sorted(features, key=lambda f: f.location.start)

            GFF.write([record], sys.stdout)

    def export_gbk(self, organism_id):
        """
        Export organism features as genbank

        :type organism_id: int
        :param organism_id: Organism ID

        :rtype: None
        :return: None
        """

        # check if the organism exists
        res = self.ci.session.query(Organism) \
            .filter(Organism.organism_id.in_([organism_id]))

        sys.stderr.write("Processing %s sequences\n" % res.count())
        for org in res:
            seq = None

            record_features = []
            features = self.ci.session.query(Feature, FeatureLocation) \
                .filter_by(organism_id=org.organism_id) \
                .join(FeatureLocation, Feature.feature_id == FeatureLocation.feature_id, isouter=True)

            sys.stderr.write("\tProcessing %s features\n" % features.count())
            for idx, (feature, featureloc) in enumerate(features):
                if idx % 5000 == 0:
                    sys.stderr.write("\t%s / %s\n" % (idx, features.count()))
                # Sequence containing feature
                if feature.residues:
                    # This seems bad? What if multiple things have seqs?
                    seq = Seq(feature.residues, IUPAC.unambiguous_dna)
                else:
                    qualifiers = {
                        self.ci.get_cvterm_name(prop.type_id): prop.value for prop in
                        self.ci.session.query(FeatureProperties).filter_by(feature_id=feature.feature_id).all()
                    }
                    record_features.append(
                        SeqFeature(
                            BioFeatureLocation(featureloc.fmin, featureloc.fmax),
                            id=feature.uniquename,
                            type=self.ci.get_cvterm_name(feature.type_id),
                            strand=featureloc.strand,
                            qualifiers=qualifiers
                        )
                    )

            record = SeqRecord(
                seq, id=org.common_name,
                name=org.common_name,
                description="%s %s" % (org.genus, org.species),
            )
            record.features = record_features

            SeqIO.write([record], sys.stdout, 'genbank')
