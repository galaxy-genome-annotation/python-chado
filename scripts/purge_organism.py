#!/usr/bin/env python
import json
import argparse
from chado import ChadoAuth, ChadoInstance, Organism

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create a new organism')

    ChadoAuth(parser)
    args = parser.parse_args()

    ci = ChadoInstance(**vars(args))

    ci.connect()

    # check if the organism exists
    res = ci.session.query(Organism).delete()
    ci.session.commit()

    print "Removed %i organism" % res
