#!/usr/bin/env python
import json
import argparse
from datetime import datetime
from chado import ChadoAuth, ChadoInstance, Analysis

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create a new analysis')

    parser.add_argument("--name", required=True, help="Analysis name")
    parser.add_argument("--program", required=True, help="Program name")
    parser.add_argument("--program-version", required=True, help="Program version")
    parser.add_argument("--algorithm", help="Algorithm name")
    parser.add_argument("--source-name", required=True, help="Source name")
    parser.add_argument("--source-version", help="Source version")
    parser.add_argument("--source-uri", help="Source URI")
    parser.add_argument("--date-executed", help="Date of execution of the analysis (format=YYYY-MM-DD, default=today)")

    ChadoAuth(parser)
    args = parser.parse_args()

    ci = ChadoInstance(args.dbhost, args.dbname, args.dbuser, args.dbpass, args.dbschema, args.debug)

    ci.connect()

    # check if the analysis exists
    res = ci.session.query(Analysis).filter_by(name = args.name)

    if res.count() > 0:
        raise Exception("Found a preexisting analysis with the same attributes in the database %s" % (ci._engine.url))

    date = datetime.today()
    if args.date_executed:
        date = datetime.strptime(args.date_executed, '%Y-%m-%d')

    newa = Analysis()
    newa.name = args.name
    newa.program = args.program
    newa.programversion = args.program_version
    newa.algorithm = args.algorithm
    newa.sourcename = args.source_name
    newa.sourceversion = args.source_version
    newa.sourceuri = args.source_uri
    newa.timeexecuted = date
    ci.session.add(newa)
    ci.session.commit()

    print "New analysis created with ID: %s" % newa.analysis_id
