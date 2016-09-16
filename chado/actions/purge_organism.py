#!/usr/bin/env python
import json
import argparse
from chado import ChadoAuth, ChadoInstance, Organism

class purge_organism(object):

    def run(self, args):
        parser = argparse.ArgumentParser(prog=('chado %s' % self.__class__.__name__), description='Remove all organisms from the database')

        ChadoAuth(parser)
        args = parser.parse_args(args)

        ci = ChadoInstance(**vars(args))

        ci.connect()

        # check if the organism exists
        res = ci.session.query(Organism).delete()
        ci.session.commit()

        print "Removed %i organism" % res
