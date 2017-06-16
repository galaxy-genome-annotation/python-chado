"""Base chado client
"""
import json


class Client(object):
    """
    Base client class implementing methods to make queries to the server
    """

    def __init__(self, engine, metadata, session):
        self.engine = engine
        self.metadata = metadata
        self.session = session
