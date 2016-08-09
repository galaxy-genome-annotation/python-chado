#!/usr/bin/env python
import json
import argparse
from chado import ChadoAuth, ChadoInstance, Organism

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create a new organism')

    parser.add_argument("--genus", required=True, help="The genus of the organism")
    parser.add_argument("--species", help="The species of the organism")
    parser.add_argument("--common", required=True, help="The common name of the organism")
    parser.add_argument("--abbr", required=True, help="The abbreviation of the organism")
    parser.add_argument("--description", help="The abbreviation of the organism")

    ChadoAuth(parser)
    args = parser.parse_args()

    ci = ChadoInstance(args.dbhost, args.dbname, args.dbuser, args.dbpass, args.dbschema, args.debug)

    ci.connect()

    # check if the organism exists
    res = ci.session.query(Organism).filter_by(common_name = args.common)

    if (res.count() > 0):
        raise Exception("Found a preexisting organism with the same attributes in the database %s" % (ci._engine.url))

    org = Organism()
    org.abbreviation = args.abbr
    org.genus = args.genus
    org.species = args.species
    org.common_name = args.common
    org.comment = args.description
    ci.session.add(org)
    ci.session.commit()

    print "New organism created with ID: %s" % org.organism_id
