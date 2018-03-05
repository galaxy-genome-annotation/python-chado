from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

from future import standard_library

standard_library.install_aliases()


def ChadoAuth(parser):
    parser.add_argument('-o', '--dbhost', required=True, help='Database Host')
    parser.add_argument('-n', '--dbname', required=True, help='Database Name')
    parser.add_argument('-u', '--dbuser', help='Database Username')
    parser.add_argument('-w', '--dbpass', help='Database Password')
    parser.add_argument('-p', '--dbport', type=int, help='Database Port', default=5432)
    parser.add_argument('--dbschema', help='Database Schema (default: public)', default="public")
