from getpass import getpass
import MySQLdb

__author__ = 'eran'

__initlized = False

USER = ''
PASSWD = ''
DB_NAME = ''

def init_connection(user=None, passwd=None, db=None):
    """
    Initializes connection details to database
    :param user: database user
    :param passwd: database password
    :param db: database name
    """
    global USER, PASSWD, DB_NAME, __initlized
    if user is None:
        user = raw_input('Enter db user: ')
    if passwd is None:
        passwd = getpass('password: ')
    if db is None:
        db = raw_input('Enter db name: ')
    USER = user
    PASSWD = passwd
    DB_NAME = db
    __initlized = True


def get_connection():
    """
    Get connection object to local database
    """
    global USER, PASSWD, DB_NAME, __initlized
    if not __initlized:
        init_connection()
    return MySQLdb.connect(host='localhost', user=USER, passwd=PASSWD, db=DB_NAME, local_infile=1)
