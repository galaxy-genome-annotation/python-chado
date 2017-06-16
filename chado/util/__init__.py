#!/usr/bin/env python
import os
import subprocess
from chado.client import Client


class UtilClient(Client):
    """
    Class containing some chado utilities
    """

    def dbshell(self):
        """
        Open a psql session to the database

        :rtype: None
        :return: None
        """
        env = dict(
            os.environ,
            PGPASSWORD=self.ci.dbpass,
            PGHOST=self.ci.dbhost,
            PGUSER=self.ci.dbuser,
            PGDATABASE=self.ci.dbname,
            PGPORT=str(self.ci.dbport),
            PGSCHEMA=self.ci.dbschema,
        )

        subprocess.call(
            ['psql'],
            env=env
        )
