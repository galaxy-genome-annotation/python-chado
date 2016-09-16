#!/usr/bin/env python
import os
import subprocess
import argparse
from chado import ChadoAuth

class dbshell(object):

    def run(self, args):
        parser = argparse.ArgumentParser(prog=('chado %s' % self.__class__.__name__), description='Invoke a database shell session')

        ChadoAuth(parser)
        args = parser.parse_args(args)
        subprocess.call(
            [
                'psql',
                '-h', args.dbhost,
                '-U', args.dbuser,
                '-p', str(args.dbport),
                args.dbname
            ],
            env=dict(os.environ, PGPASSWORD=args.dbpass)
        )
