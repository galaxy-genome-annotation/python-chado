#!/usr/bin/env python
import argparse
from chado import ChadoAuth, ChadoInstance, Organism

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create a new organism')

    parser.add_argument("--genus", help="The genus of the organism")
    parser.add_argument("--species", help="The species of the organism")
    parser.add_argument("--common", help="The common name of the organism")
    parser.add_argument("--abbr", help="The abbreviation of the organism")
    parser.add_argument("--comment", help="The comment of the organism")
    parser.add_argument("-q", "--quiet", action='store_true', help="Only print out IDs")

    ChadoAuth(parser)
    args = parser.parse_args()

    ci = ChadoInstance(**vars(args))

    ci.connect()

    # check if the organism exists
    res = ci.session.query(Organism)
    if args.genus:
        res = res.filter_by(genus=args.genus)
    if args.species:
        res = res.filter_by(species=args.species)
    if args.common:
        res = res.filter_by(common_name=args.common)
    if args.abbr:
        res = res.filter_by(abbreviation=args.abbr)
    if args.comment:
        res = res.filter_by(comment=args.comment)

    for org in res:
        if args.quiet:
            print org.organism_id
        else:
            print '\t'.join(map(str, (
                org.organism_id,
                org.genus,
                org.species,
                org.abbreviation,
                org.common_name,
                org.comment,
            )))
