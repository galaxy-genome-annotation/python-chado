#!/usr/bin/env python
import os
import subprocess
import argparse
from chado import ChadoAuth

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create a new organism')

    ChadoAuth(parser)
    args = parser.parse_args()
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
