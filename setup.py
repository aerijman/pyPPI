from distutils.core import setup
from distutils.command.install import INSTALL_SCHEMES


for scheme in INSTALL_SCHEMES.values():
    scheme['data'] = scheme['purelib']

setup(
    name='pyPPI',
    version='0.1',
    packages=['pyPPI', 'pyPPI.alignment', 'pyPPI.surfaceComplementarity'],
    scripts=['bin/setupPpiDb.py'],
    url=['https://github.com/eranroz/pyPPI','https://github.com/aerijman/pyPPI'],
    license='MIT',
    authors=['eranroz','aerijman'],
    author_email=['eranroz@cs.huji.ac.il','aerijman@fredhutch.org'],
    description='Protein-Protein interactions calculations', requires=['requests', 'pymysql'],
    data_files=[('pyPPI/molprobity', ['pyPPI/molprobity/molprobity.pl',
                                      'pyPPI/molprobity/remediator.pl',
                                      'pyPPI/molprobity/reduce']),
                ('pyPPI/sqls', ['pyPPI/sqls/createDB.sql',
                                'pyPPI/sqls/createInterface.sql',
                                'pyPPI/sqls/donors2.sql',
                                'pyPPI/sqls/proteinComplex.sql']),
                ('pyPPI', ['pyPPI/DonAcc2.txt']),
                ]
)
