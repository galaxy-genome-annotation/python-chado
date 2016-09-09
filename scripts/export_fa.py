#!/usr/bin/env python
import argparse
import sys
try:
    import tqdm
    HAS_PROGRESS_BAR = True
except ImportError:
    # Progress bar library, not required to run.
    HAS_PROGRESS_BAR = False

from chado import ChadoAuth, ChadoInstance, Organism, Feature

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Export any sequence records associated with an organism')

    parser.add_argument("orgId", nargs='+', help="The id of the organism")
    parser.add_argument("--file", action='store_true', help="Store sequences in files split by organism")
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
        if HAS_PROGRESS_BAR:
            pbar.update(1)

        sequence_features = ci.session.query(Feature) \
            .filter_by(organism_id=org.organism_id) \
            .filter(Feature.seqlen > 0)

        if args.file:
            output = open('{0.organism_id}.{0.genus}.{0.species}-{0.common_name}.fa'.format(org), 'w')
        else:
            output = sys.stdout

        if HAS_PROGRESS_BAR:
            sbar = tqdm.tqdm(total=sequence_features.count(), desc='sequences')

        for seq in sequence_features:
            output.write('>{0.uniquename} [ID={0.feature_id}]'.format(seq))
            output.write('\n')
            output.write(seq.residues)
            output.write('\n')

            if HAS_PROGRESS_BAR:
                sbar.update(1)

        if HAS_PROGRESS_BAR:
            sbar.close()

        if args.file:
            output.close()

    if HAS_PROGRESS_BAR:
        pbar.close()
