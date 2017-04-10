"""
Setups a protein database in MySQL: a database of interesting properties of the proteins based on scripts of this library.
Create a folder and place there some file with list of PDBs to analyze.
The program will create the following directory structure in the same directory:
    ./pdbs/ - list of pdbs downloaded
    ./results/ - some results of the analysis scripts
    ./debug/ - scripts results for specific structure
"""

import string
import subprocess
import urllib2
import argparse
import os
import sys

if __name__ == "__main__":
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))

from kdtree2 import KDTree
import math
from ASA import ASA
from pdbReader import PDBReader
from hbonds import hbonds
import DBConfig
import surfaceComplementarity.VDW as VDW
import surfaceComplementarity.interfaceDepth as Periphery
__author__ = 'eran'

"""
Distance in angtroms between the chains that is relevant for defining the interface
"""
INTERFACE_DISTANCE = 4
WORKING_DIRECTORY = './'
PDBS_DIR = "./pdbs/"
RESULTS_DIR = "./results/"
MOLPROBITY_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, 'molprobity')) + os.path.sep
#os.path.relpath('../molprobity/', __file__)


def downloadPDB(pdb):
    """
    Downloads a PDB from protein data base
    :param pdb: pdb identifier
    """
    url = 'http://www.rcsb.org/pdb/files/{0}.pdb'.format(pdb)
    print 'downloading ', pdb
    print url
    req = urllib2.Request(url)
    f = urllib2.urlopen(req)
    res = f.read()
    print 'end downloading'

    newPDB = getFile(pdb)
    print>> newPDB, res
    newPDB.close()


def getFile(name):
    global PDBS_DIR
    return file(os.path.join(PDBS_DIR, name + ".pdb"), "w")


def downloadDB(pdbList):
    print "Downloading pdbs according to list"
    for pdb in pdbList:
        #dont download twice the same PDB
        if os.path.exists(os.path.join(PDBS_DIR, pdb + "_FH.pdb")): continue
        downloadPDB(pdb)
        molprobity(pdb)
    print "Finished downloading pdbs"


def molprobity(pdbName):
    global MOLPROBITY_DIR, PDBS_DIR
    if os.path.exists(MOLPROBITY_DIR + pdbName + '_FH.pdb'):
        return True#already exist
    print 'Starting molprobity ', pdbName
    os.path.join(MOLPROBITY_DIR, 'remediator.pl')
    subprocess.check_output('perl ' + os.path.join(MOLPROBITY_DIR, 'remediator.pl') + ' ' + os.path.join(PDBS_DIR,
                                                                                                         pdbName + ".pdb") + ' > a',
                            shell=True)
    try:
        subprocess.check_output(MOLPROBITY_DIR + 'reduce a > b', shell=True)
    except:
        print 'error prasing PDB ', pdbName
        pass#yakky kaky, but reduce returns 1 exit
    subprocess.check_output(
        'perl ' + MOLPROBITY_DIR + 'remediator.pl b -oldout> ' + os.path.join(PDBS_DIR, pdbName + "_FH.pdb"),
        shell=True)

    os.remove(os.path.join(PDBS_DIR, pdbName + ".pdb"))
    return True


def buildASAperAtomForComplex(pdb, result):
    asaCalc = ASA(pdb)
    asaCalc.execute()
    for atom, asa in asaCalc.interPerAtom.iteritems():
        #complex inter
        res = [pdb.name, atom.chain, atom.residue, atom.resId, atom.symbol, atom.atomType, asa, atom.tempFactor, 0]
        toPrint = string.join([str(a) for a in res], ',')
        print>> result, toPrint
        #complex intra (seperatd)
        asa = asaCalc.diffASAperAtom[atom] + asa
        res = [pdb.name, atom.chain, atom.residue, atom.resId, atom.symbol, atom.atomType, asa, atom.tempFactor, 1]
        toPrint = string.join([str(a) for a in res], ',')
        print>> result, toPrint


def calcInterfaceDist(pdb, result):
    """
    Defines interface by distance
    """
    global INTERFACE_DISTANCE
    partA = [a for a in pdb.atoms if a.chain in pdb.interfaceParts[0]]
    partB = [a for a in pdb.atoms if a.chain in pdb.interfaceParts[1]]
    if len(partA) == 0 or len(partB) == 0:
        print 'WARNING: ', pdb.name, ' doesnt have atoms in one its chains'
        return
    aTree = KDTree.construct_from_data(partA[:])
    bTree = KDTree.construct_from_data(partB[:])
    complexChains = string.join(pdb.interfaceParts, ':')
    for part, tree in [(partA, bTree), (partB, aTree)]:
        for atom in part:
            near = tree.findNearest(query_point=atom.coord, num=1)[0]
            dist = math.sqrt(atom.distance(near))
            if dist < INTERFACE_DISTANCE:
                print>> result, string.join(
                    [pdb.name, complexChains, atom.chain, str(atom.resId), atom.symbol, atom.atomType, str(dist)], ',')


def createInterfaceCSV(pdbsToAnalyze):
    """
    interface can be defined by either ASA or distance
    we use both of them
    """
    global PDBS_DIR, RESULTS_DIR
    if all(os.path.exists(os.path.join(RESULTS_DIR, resFile)) for resFile in ['PerAtomASA.csv', 'PerAtomASA.csv']):
        print 'Data already exist in result directory.'
        return

    asaPerAtom = file(os.path.join(RESULTS_DIR, 'PerAtomASA.csv'), 'w')
    distancePerAtom = file(os.path.join(RESULTS_DIR, 'PerAtomDistance.csv'), 'w')
    pdbs = os.listdir(PDBS_DIR)
    print>> distancePerAtom, 'PDB,Chains,Chain,ResId,Symbol,Atom,MinDistance'
    print>> asaPerAtom, 'PDB,Chain,Residue,ResId,Symbol,AtomType,ASA,tempFactor,Seperated'
    failedPDBs = []
    pdbsNamesToChains = dict((p[0], p[1].split(':') if len(p) > 1 else None) for p in pdbsToAnalyze)
    for pdbName in pdbs:
        if pdbName[0:4] not in pdbsNamesToChains: continue
        pdb = PDBReader.readFile(os.path.join(PDBS_DIR, pdbName), pdbsNamesToChains[pdbName[0:4]])
        try:
            print 'Writing ASA for ', pdb.name
            buildASAperAtomForComplex(pdb, asaPerAtom)
            print 'Writing distance for ', pdb.name
            calcInterfaceDist(pdb, distancePerAtom)
        except IndexError:
            failedPDBs.append(pdb.name)

    asaPerAtom.close()
    distancePerAtom.close()
    print 'Finished'
    if len(failedPDBs) > 0:
        print 'Failed to process:', ','.join(failedPDBs)


def createDataBase():
    print 'Creating DB:', DBConfig.DB_NAME
    installDB = os.path.abspath(os.path.join(os.path.dirname(__file__), 'createDB.sql'))
    metadataDB = os.path.abspath(os.path.join(os.path.dirname(__file__), 'donors2.sql'))
    createInterfaceSql = os.path.abspath(os.path.join(os.path.dirname(__file__), 'createInterface.sql'))

    subprocess.call(
        "mysql -u %s -p%s -e 'create database if not exists %s'" % (DBConfig.USER, DBConfig.PASSWD, DBConfig.DB_NAME),
        shell=True)
    #create schema
    subprocess.call('mysql %s -u%s -p%s < %s ' % (DBConfig.DB_NAME, DBConfig.USER, DBConfig.PASSWD, installDB),
                    shell=True)
    #insert metadata
    subprocess.call('mysql %s -u%s -p%s < %s ' % (DBConfig.DB_NAME, DBConfig.USER, DBConfig.PASSWD, metadataDB),
                    shell=True)
    conn = DBConfig.get_connection()
    cursor = conn.cursor()
    cursor.execute('''
    load data local infile '%s' into table interfaceDist fields terminated by ',' optionally enclosed by '"'  lines terminated by '\n' ignore 1 lines (PDB,Chains,Chain,ResId,Symbol,Atom,MinDist);
    ''' % (os.path.join(RESULTS_DIR, 'PerAtomDistance.csv')))
    cursor.execute('''
    load data local infile '%s' into table perAtomASA fields terminated by ',' optionally enclosed by '"'  lines terminated by '\n' ignore 1 lines (PDB,Chain,Residue,ResId,Symbol,Atom,ASA,Bfactor,Seperated);
    ''' % (os.path.join(RESULTS_DIR, 'PerAtomASA.csv')))
    conn.commit()

    #create interface table
    subprocess.call('mysql %s -u%s -p%s < %s ' % (DBConfig.DB_NAME, DBConfig.USER, DBConfig.PASSWD, createInterfaceSql),
                    shell=True)

    #add metadata table with complexs in the database
    pdbsNamesToChains = dict((p[0], p[1].split(':') if len(p) > 1 else None) for p in pdbsToAnalyzeWithChains)
    dataToInsert=[]
    for pdbName, chains in pdbsNamesToChains.iteritems():
        pdb = PDBReader.readFile(os.path.join(PDBS_DIR, '%s_FH.pdb'%pdbName), pdbsNamesToChains[pdbName[0:4]])
        if chains is None:
            compunds = pdb.compunds.split(' - ')
            dataToInsert.append((pdbName, pdb.interfaceParts[0], compunds[0] if len(compunds)>1 else compunds, pdb.interfaceParts[1], compunds[1] if len(compunds)>1 else '' ))
        else:
            dataToInsert.append((pdbName, pdb.interfaceParts[0], '', pdb.interfaceParts[1], ''))

    cursor = conn.cursor()
    cursor.executemany('''
    INSERT INTO proteinComplex (PDB,UnboundChainA,NameA,UnboundChainB,NameB)
    values (%s,%s,%s,%s,%s)
    ''', dataToInsert)
    conn.commit()
    conn.close()
    print 'database created!'


def getInterfaceAtoms(cur, pdb):
    """
    Gets interface atoms from database
    :param cur: cursor to database
    :param pdb: pdb object to get atoms from
    :return: list of interface atoms
    """
    cur.execute('''
    select Chain,ResId,Symbol from NinterfaceAtoms
    where PDB='%s'
    ''' % pdb.name)
    interfaceAtoms = []
    for chain, resid, symbol in cur.fetchall():
        interfaceAtoms.append(
            (a for a in pdb.atoms if a.chain == chain and a.resId == resid and a.symbol == symbol).next())
    return interfaceAtoms

def fillInterfacePeriphrial(pdbsToAnalyze):
    global PDBS_DIR, RESULTS_DIR

    if os.path.exists(os.path.join(RESULTS_DIR, 'interfacePeriphrial.csv')):
        print 'Data already exist in result directory for interface periphery.'
        return

    pdbsNamesToChains = dict((p[0], p[1].split(':') if len(p) > 1 else None) for p in pdbsToAnalyze)

    Periphery.PDBS_DIR = PDBS_DIR
    with file(os.path.join(RESULTS_DIR, 'interfacePeriphrial.csv'), 'w') as interfacePeriphrial:
        print>>interfacePeriphrial,'PDB,Chain,ResId,Symbol,Peripherial,PropPeri'
        for pdbName, chains in pdbsNamesToChains.iteritems():
            print 'Calculating peripheral table for %s '%pdbName
            depthL, peripherialL = Periphery.calcPeriphrialForPDB(pdbName, chains)
            for atom, peri, propPeri in peripherialL:
                print>> interfacePeriphrial, ','.join([pdbName, atom.chain, str(atom.resId), atom.symbol, str(peri), str(propPeri)])

    conn = DBConfig.get_connection()
    cursor = conn.cursor()
    cursor.execute('''
    load data local infile '%s' into table interfacePeriphrial
    fields terminated by ',' optionally enclosed by '"'  lines terminated by '\n'
    ignore 1 lines (PDB,Chain,ResId,Symbol,Peri,PropPeri);
    ''' % (os.path.join(RESULTS_DIR, 'interfacePeriphrial.csv')))
    conn.commit()
    conn.close()

def calcEnergyTerms(pdbsToAnalyze):
    """
    Finds hydrogen bonds near interface atoms and calculates their energy,
    and calculates VDW and electrostatic energy for PDB
    """
    global PDBS_DIR, RESULTS_DIR

    if all(os.path.exists(os.path.join(RESULTS_DIR, resFile)) for resFile in ['Ndrieding.csv', 'interfaceVDW.csv']):
        print 'Data already exist in result directory for energy terms.'
        return

    conn = DBConfig.get_connection()
    cursor = conn.cursor()
    pdbsNamesToChains = dict((p[0], p[1].split(':') if len(p) > 1 else None) for p in pdbsToAnalyze)
    with file(os.path.join(RESULTS_DIR, 'Ndrieding.csv'), 'w') as driedingResult:
        print>> driedingResult, 'PDB,DonorChain,DonorResId,DonorSymbol,AccChain,AccResId,AccSymbol,Energy'

        pdbs = os.listdir(PDBS_DIR)
        for pdbName in pdbs:
            if pdbName[0:4] not in pdbsNamesToChains: continue
            pdb = PDBReader.readFile(os.path.join(PDBS_DIR, pdbName), pdbsNamesToChains[pdbName[0:4]])
            interfaceAtoms = getInterfaceAtoms(cursor, pdb)
            bonds = hbonds(pdb)
            bonds.HDPlusDefinition = False
            cBondList = bonds.hbonds(interfaceAtoms)
            print('Calcing Hbonds for %s' % pdb.name)
            for donor, acceptor, eng in cBondList:
                toPrint = [pdb.name, donor.chain, donor.resId, donor.symbol, acceptor.chain, acceptor.resId,
                           acceptor.symbol, eng]
                print>> driedingResult, ','.join([str(a) for a in toPrint])
    cursor.execute('''
    load data local infile '%s' into table Ndrieding
    fields terminated by ',' optionally enclosed by '"'  lines terminated by '\n'
    ignore 1 lines (PDB,DonorChain,DonorResId,DonorSymbol,AccChain,AccResId,AccSymbol,Energy);
    ''' % (os.path.join(RESULTS_DIR, 'Ndrieding.csv')))
    conn.commit()

    print 'Calculating VDW energy between interfaces'
    with file(os.path.join(RESULTS_DIR, 'interfaceVDW.csv'), 'w') as vdw_result:
        print>> vdw_result, 'PDB,VDV,VDVx,clashV,clashS'
        VDW.PDBS_DIR = PDBS_DIR
        for pdb, chains in pdbsNamesToChains.iteritems():
            print 'Calcing VDW for %s'%pdb
            sumVDW, sumVDWx, clashV, clashS = VDW.calcCompl(pdb, chains)
            print>> vdw_result, ','.join([pdb, str(sumVDW),str(sumVDWx), str(clashV), str(clashS)])
    cursor.execute('''
    load data local infile '%s' into table interfaceVDW
    fields terminated by ',' optionally enclosed by '"'  lines terminated by '\n'
    ignore 1 lines (PDB,VDV,VDVx6,ClashV,ClashS);
    ''' % (os.path.join(RESULTS_DIR, 'interfaceVDW.csv')))
    conn.commit()
    cursor.close()
    conn.close()



if __name__ == "__main__":
    __package__ = "scripts.setupDB.setupDB"
    parser = argparse.ArgumentParser(description="Setup/download protein database based on PDB")
    #parser.add_argument("--downloadList",required=True,help="A file with a list of PDB to download")
    parser.add_argument("pdbList", help="A file with a list of PDB to download")
    parser.add_argument("--folder", help="Name of the folder to contain downloaded files")
    parser.add_argument("--dbName", help="Name of the database to create.")
    args = parser.parse_args()
    if args.pdbList is None:
        sys.exit("Please provide a file with list of PDBs to anaylze")

    WORKING_DIRECTORY = args.folder if args.folder is not None else os.path.dirname(os.path.abspath(args.pdbList))
    print 'WORKING DIR:', WORKING_DIRECTORY

    PDBS_DIR = os.path.join(WORKING_DIRECTORY, 'pdbs')
    RESULTS_DIR = os.path.join(WORKING_DIRECTORY, 'results')
    for dir in [PDBS_DIR, RESULTS_DIR]:
        if not os.path.exists(dir):
            os.mkdir(dir)

    pdbsToAnalyzeWithChains = [pdb.strip().upper().split("_") for pdb in file(args.pdbList, 'r') if
                               pdb[0:1] != '#']#todo: add treatment for chains specificatin instad of [0:4]
    pdbsToAnalyze = [pdb[0] for pdb in pdbsToAnalyzeWithChains]
    downloadDB(pdbsToAnalyze)#download from PDB bank and add hydrogens
    createInterfaceCSV(pdbsToAnalyzeWithChains)#define interface by distance and by asa
    print '''The script will now create DB. DB is required for extra calculations
    including VDW and hydrogen bonds
     '''
    try:
        if args.dbName:
            DBConfig.DB_NAME=args.dbName
        DBConfig.init_connection()
        createDataBase(pdbsToAnalyzeWithChains)

        #post database creation scripts
        fillInterfacePeriphrial(pdbsToAnalyzeWithChains)
        calcEnergyTerms(pdbsToAnalyzeWithChains)
    except KeyboardInterrupt:
        print 'DB will not be created. Use ./results table to see the results'


