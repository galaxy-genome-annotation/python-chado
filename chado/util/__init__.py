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

    def launch_docker_image(self, background=False):
        """
        Launch a chado docker image.

        :type background: bool
        :param background: Launch the image in the background

        :rtype: None
        :return: None
        """

        subprocess.call([
            'docker',
            'run', '-it',
            '-d' if background else '',
            '-e', 'INSTALL_YEAST_DATA=1',
            '-p', '5432:5432',
            'erasche/chado:1.31-jenkins97-pg9.5'
        ])
