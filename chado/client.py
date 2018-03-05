"""Base chado client
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

from future import standard_library

standard_library.install_aliases()


class Client(object):
    """
    Base client class implementing methods to make queries to the server
    """

    def __init__(self, engine, metadata, session, ci):
        self.engine = engine
        self.metadata = metadata
        self.session = session
        self.ci = ci
        self.model = ci.model
