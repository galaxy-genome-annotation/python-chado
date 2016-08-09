#!/usr/bin/env python
import json
import argparse
from chado import ChadoAuth, ChadoInstance, Organism

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create a new analysis')

    parser.add_argument("--name", required=True, help="Analysis name")
    parser.add_argument("--program", required=True, help="Program name")
    parser.add_argument("--program-version", required=True, help="Program version")
    parser.add_argument("--algorithm", help="Algorithm name")
    parser.add_argument("--source-name", required=True, help="Source name")
    parser.add_argument("--source-version", help="Source version")
    parser.add_argument("--source-uri", help="Source URI")

    ChadoAuth(parser)
    args = parser.parse_args()

    ci = ChadoInstance(args.dbhost, args.dbname, args.dbuser, args.dbpass, args.dbschema, args.debug)

    ci.connect()

    # check if the analysis exists
    res = ci.session.query(Analysis).filter_by(name = args.name)
    found_res = (res.count() > 0)

    if not found_res:
        newa = Analysis()
        newa.name = args.name
        newa.program = args.program
        newa.programversion = args.program_version
        newa.algorithm = args.algorithm
        newa.sourcename = args.source_name
        newa.sourceversion = args.source_version
        newa.sourceuri = args.source_uri
        newa.timeexecuted = args.timeexecuted # This should be a timestamp
        session.add(newa)
        session.commit()
    else:
        raise Exception("Found a preexisting analysis with the same attributes in the database %s" % (self._engine.url))
