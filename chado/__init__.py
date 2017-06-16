from sqlalchemy import create_engine, MetaData, Table
from sqlalchemy.orm import mapper, sessionmaker
from chado.organism import OrganismClient


class ChadoInstance(object):

    def __init__(self, dbhost="localhost", dbname="chado", dbuser="chado", dbpass="chado", dbschema="public", dbport=5432, debug=False, **kwargs):
        self.dbhost = dbhost
        self.dbname = dbname
        self.dbuser = dbuser
        self.dbpass = dbpass
        self.dbport = dbport
        self.dbschema = dbschema
        self.debug = debug

        self._engine = create_engine('postgresql://%s:%s@%s:%s/%s' % (self.dbuser, self.dbpass, self.dbhost, self.dbport, self.dbname))
        self._metadata = MetaData(self._engine, schema=self.dbschema)
        Session = sessionmaker(bind=self._engine)
        self.session = Session()

        self.organism = OrganismClient(self._engine, self._metadata, self.session)

    def __str__(self):
        return '<ChadoInstance at %s>' % self.dbhost

    # def connect(self):
        # # self._test_db_access()
        # self._cv_id_cache = {}

    # def _test_db_access(self):
        # tables = self._engine.table_names(schema=self.dbschema)
        # if ('analysis' not in tables or 'feature' not in tables):
            # raise Exception("Could not find Chado tables in db %s" % (self._engine.url))

    def _reflect_tables(self):
        organism = Table('organism', self._metadata, autoload=True)
        mapper(Organism, organism)
