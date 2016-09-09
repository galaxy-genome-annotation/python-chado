#!/usr/bin/env python
import argparse
from BCBio import GFF
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation as BioFeatureLocation
from chado import ChadoAuth, ChadoInstance, Organism, Feature, FeatureLocation, FeatureProperties, FeatureRelationship
import logging
logging.basicConfig()
log = logging.getLogger(name='export_gff3')
try:
    import tqdm
    HAS_PROGRESS_BAR = True
except ImportError:
    # Progress bar library, not required to run.
    HAS_PROGRESS_BAR = False


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Export a GFF3 formatted dataset for an organism')

    parser.add_argument("orgId", nargs='+', help="The id of the organism")
    ChadoAuth(parser)
    args = parser.parse_args()

    ci = ChadoInstance(**vars(args))
    ci.connect()

    # check if the organism exists
    res = ci.session.query(Organism).filter(Organism.organism_id.in_(args.orgId))

    if HAS_PROGRESS_BAR:
        pbar = tqdm.tqdm(total=res.count(), desc='organism')

    for org in res:
        if HAS_PROGRESS_BAR:
            pbar.update(1)

        # TODO: can we do this properly?
        seq = Seq("A" * 1, IUPAC.unambiguous_dna)

        record_features = []
        # Annotation features
        features = ci.session.query(Feature, FeatureLocation) \
            .filter_by(organism_id=org.organism_id) \
            .filter(Feature.seqlen==None) \
            .join(FeatureLocation, Feature.feature_id == FeatureLocation.feature_id, isouter=True)

        if HAS_PROGRESS_BAR:
            fbar = tqdm.tqdm(total=features.count(), desc='features')

        biopy_features = {}
        for feature, featureloc  in features:
            if HAS_PROGRESS_BAR:
                fbar.update(1)
            #[u'dbxref_id', u'feature_id', u'is_analysis', u'is_obsolete',
            # u'md5checksum', u'name', u'organism_id', u'residues', u'seqlen',
            # u'timeaccessioned', u'timelastmodified', u'type_id',
            # u'uniquename']
            #[u'feature_id', u'featureloc_id', u'fmax', u'fmin',
            # u'is_fmax_partial', u'is_fmin_partial', u'locgroup', u'phase',
            # u'rank', u'residue_info', u'srcfeature_id', u'strand']
            qualifiers = {
                ci.get_cvterm_name(prop.type_id): prop.value for prop in
                ci.session.query(FeatureProperties).filter_by(feature_id=feature.feature_id).all()
            }

            qualifiers['ID'] = feature.uniquename

            biopy_features[feature.feature_id] = SeqFeature(
                BioFeatureLocation(featureloc.fmin, featureloc.fmax),
                id=feature.uniquename,
                type=ci.get_cvterm_name(feature.type_id),
                strand=featureloc.strand,
                qualifiers=qualifiers
            )

        if HAS_PROGRESS_BAR:
            fbar.close()

    #res = ci.session.query(Organism).filter(Organism.organism_id.in_(args.orgId))
        relationships = ci.session.query(FeatureRelationship) \
            .filter(FeatureRelationship.subject_id.in_(biopy_features.keys()))

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
        for rel in relationships:
            term = ci.get_cvterm_name(rel.type_id)
            if term != 'part_of':
                log.error("Cannot handle non-part_of relationships (%s %s %s)", rel.subject_id, term, rel.object_id)
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

    if HAS_PROGRESS_BAR:
        pbar.close()
