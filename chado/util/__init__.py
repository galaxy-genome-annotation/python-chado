#!/usr/bin/env python
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from builtins import dict
from builtins import str
from future import standard_library
standard_library.install_aliases()
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
