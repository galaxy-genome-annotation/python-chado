#!/usr/bin/env python
import argparse
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation as BioFeatureLocation
try:
    import tqdm
    HAS_PROGRESS_BAR = True
except ImportError:
    # Progress bar library, not required to run.
    HAS_PROGRESS_BAR = False


from chado import ChadoAuth, ChadoInstance, Organism, Feature, FeatureLocation, FeatureProperties

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Export a GenBank formatted dataset for an organism')

    parser.add_argument("orgId", nargs='+', help="The id of the organism")
    ChadoAuth(parser)
    args = parser.parse_args()

    ci = ChadoInstance(**vars(args))
    ci.connect()

    # check if the organism exists
    res = ci.session.query(Organism) \
        .filter(Organism.organism_id.in_(args.orgId))

    if HAS_PROGRESS_BAR:
        pbar = tqdm.tqdm(total=res.count(), desc='organism')

    for org in res:
        seq = None

        if HAS_PROGRESS_BAR:
            pbar.update(1)

        record_features = []
        features = ci.session.query(Feature, FeatureLocation) \
            .filter_by(organism_id=org.organism_id) \
            .join(FeatureLocation, Feature.feature_id == FeatureLocation.feature_id, isouter=True) \
            .all()

        for feature, featureloc  in features:
            # Sequence containing feature
            if feature.residues:
                # This seems bad? What if multiple things have seqs?
                seq = Seq(feature.residues, IUPAC.unambiguous_dna)
            else:
                qualifiers = {
                    ci.get_cvterm_name(prop.type_id): prop.value for prop in
                    ci.session.query(FeatureProperties).filter_by(feature_id=feature.feature_id).all()
                }
                record_features.append(
                    SeqFeature(
                        BioFeatureLocation(featureloc.fmin, featureloc.fmax),
                        id=feature.uniquename,
                        type=ci.get_cvterm_name(feature.type_id),
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

    if HAS_PROGRESS_BAR:
        pbar.close()
