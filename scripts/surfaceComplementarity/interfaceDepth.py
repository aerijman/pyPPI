"""
Finds depth and periphriality for atoms
"""

import sys, os, string, math, MySQLdb, logging, datetime

sys.path.append('../')
sys.path.append('../algorithms')
import itertools
from ASA import ASA, radio_atom, spiral
from pdbReader import PDBReader
from kdtree2 import KDTree
import DBConfig

PDBS_DIR = '/home/eran/Eran/pdbs/'


def getNamesAndChainFromDB():
    conn = DBConfig.get_connection()
    cursor = conn.cursor()

    #create interface table
    cursor.execute("""select Complex from origTable2""")
    res = []
    for c in cursor.fetchall():
        c = c[0]
        res.append((c[0:4], c[5:].split(':')))
    return res


def getInterfaceAtoms(pdb, includingDistance=False):
    conn = DBConfig.get_connection()
    cursor = conn.cursor()

    #create interface table
    if includingDistance:
        cursor.execute("""select Chain,ResId,Symbol from
                             NinterfaceAtoms
                             where PDB='%s'""" % ( pdb.name ))
    else:
        cursor.execute("""select Chain,ResId,Symbol from
                             interfaceAtoms
                             where PDB='%s'""" % ( pdb.name ))
        #find unbounds relating to
    interfaceAtoms = []
    for chain, resId, symbol in cursor.fetchall():
        interfaceAtoms += [a for a in pdb.atoms if a.chain == chain and a.resId == resId and a.symbol[0:3] == symbol]
    return interfaceAtoms


def oldCode():
    aTree = KDTree.construct_from_data(components[0][:])
    interDist = dict()
    for atom in components[0]:
        if interDist.has_key(atom):
            continue
        far = aTree.findFarest(query_point=atom.coord, num=1)
        far = far[0]
        dist = far.distance(atom)
        interDist[far] = dist
        interDist[atom] = dist
    ll = [(k, v) for k, v in interDist.iteritems()]
    ll = sorted(ll, key=lambda x: -x[1])
    maxDist = ll[0][1]
    print 'Peripherial'
    peripherial, nperipherial = set(), set()
    for atom, dist in ll:
        if dist > 2.0 / 3.0 * maxDist:
            print atom.resId, atom.symbol
            if atom.chain == 'C':
                peripherial.add(str(atom.resId))

    print 'non peripherial'
    for atom, dist in ll:
        if dist <= 2.0 / 3.0 * maxDist:
            nperipherial.add(str(atom.resId))
            if atom.chain == 'C':
                print atom.resId, atom.symbol

    print 'sel peripherial, res ', string.join(peripherial, '+')
    print 'sel nonperipherial, res ', string.join(nperipherial, '+')


def assignDepth(interface):
    """
     Finds the distance of farest atom in interface for each atom
     """
    aTree = KDTree.construct_from_data(interface)
    depth = []
    res = []
    maxDepth = 0.0
    for atom in interface:
        far = aTree.findFarest(query_point=atom.coord, num=1)[0]
        distance2 = atom.distance(far)
        maxDepth = max(distance2, maxDepth)
        depth.append((atom, distance2))

    maxDepth = math.sqrt(maxDepth)
    for atom, distance2 in depth:
        dist = math.sqrt(distance2)
        res.append((atom, dist, dist / maxDepth))
    return res


def assignPeriphrial(interface, surfaceAtoms):
    """
     Finds the nearest surface atom which is non interface to each atom
     """
    peripherial = []
    res = []
    aTree = KDTree.construct_from_data([atom for atom in surfaceAtoms if atom not in interface])
    maxPeripherial = 0
    for atom in interface:
        near = aTree.findNearest(query_point=atom.coord, num=1)[0]
        distance2 = atom.distance(near)
        maxPeripherial = max(maxPeripherial, distance2)
        peripherial.append((atom, distance2))
    maxPeripherial = math.sqrt(maxPeripherial)
    for atom, distance2 in peripherial:
        dist = math.sqrt(distance2)
        res.append((atom, dist, dist / maxPeripherial))
    return res


def calcPeriphrialForPDB(pdbName, chains):
    global PDBS_DIR
    pdbPath = os.path.join(PDBS_DIR, '%s_FH.pdb' % pdbName)
    pdb = PDBReader.readFile(pdbPath, chains)
    asa = ASA(pdb)
    asa.execute()
    interface = getInterfaceAtoms(pdb, True)
    components = []
    surfaceComponents = []
    for part in pdb.interfaceParts:
        partAtoms = [a for a in interface if a.chain in part]
        surfaceAtoms = [atom for atom in pdb.atoms if atom.chain in part and not asa.isBuried(atom)]
        components.append(partAtoms)
        surfaceComponents.append(surfaceAtoms)
    depth = []
    peripherial = []
    for i in range(0, 2):
        depth += assignDepth(components[i])
        peripherial += assignPeriphrial(components[i], surfaceComponents[i])
    return depth, peripherial


def main():
    depthRes = file('depthRes2.csv', 'w')
    periRes = file('periRes2.csv', 'w')
    print>> depthRes, string.join(['PDB', 'Chain', 'ResId', 'Symbol', 'Depth', 'PropDepth'], ',')
    print>> periRes, string.join(['PDB', 'Chain', 'ResId', 'Symbol', 'Peripherial', 'PropPeri'], ',')
    for pdbName, chain in getNamesAndChainFromDB():
        print 'Proccessing %s' % pdbName
        depthL, peripherialL = calcPeriphrialForPDB(pdbName, chain)
        for atom, depth, propDepth in depthL:
            print>> depthRes, string.join(
                [pdbName, atom.chain, str(atom.resId), atom.symbol, str(depth), str(propDepth)], ',')
        for atom, peri, propPeri in peripherialL:
            print>> periRes, string.join([pdbName, atom.chain, str(atom.resId), atom.symbol, str(peri), str(propPeri)],
                                         ',')
        depthRes.flush()
        periRes.flush()
    depthRes.close()
    periRes.close()
    print 'Finished'


if __name__ == "__main__":
    main()
"""
pdb = PDBReader.readFile('/home/eran/Eran/pdbs/1A2K_FH.pdb',['C','AB'])
asa = ASA(pdb)
asa.execute()
interface = getInterfaceAtoms(pdb)
components=[]
surfaceComponents=[]
for part in pdb.interfaceParts:
     partAtoms = [a for a in interface if a.chain in part]
     surfaceAtoms = [atom for atom in pdb.atoms if atom.chain in part and not asa.isBuried(atom)]
     components.append(partAtoms)
     surfaceComponents.append(surfaceAtoms)

depth = assignDepth(components[0])
peripherial = assignPeriphrial(components[0], surfaceComponents[0])
"""


