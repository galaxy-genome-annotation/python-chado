#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import subprocess
from builtins import dict
from builtins import str

from chado.client import Client

from future import standard_library

standard_library.install_aliases()


class UtilClient(Client):
    """
    Some chado utilities
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

    def launch_docker_image(self, background=False, no_yeast=False):
        """
        Launch a chado docker image.

        :type background: bool
        :param background: Launch the image in the background

        :type no_yeast: bool
        :param no_yeast: Disable loading of example yeast data

        :rtype: None
        :return: None
        """

        cmd = [
            'docker',
            'run', '-it',
        ]

        if background:
            cmd.append('-d')

        if not no_yeast:
            cmd += ['-e', 'INSTALL_YEAST_DATA=1']

        cmd += [
            '-p', '5432:5432',
            'erasche/chado:1.31-jenkins110.2-pg9.5'
        ]

        subprocess.call(cmd)
