from chado import *

ci = ChadoInstance(dbuser="postgres", dbpass="postgres", dbname="postgres")

def setup_package():
    global ci
