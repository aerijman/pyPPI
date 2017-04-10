"""

H2O:
    from 2 atoms - less than 3.5

References:
    *pl files
    *http://www.wwpdb.org/documentation/format33/sect9.html#ATOM -pdb RFC
"""

import re, os, math, string
from math import sqrt
from atom import atom
from water import water
from kdtree import KDTree
import time
import logging

NEIGHBOR_DISTANCE = 16

def angle(DH,distance,HA,radians=True):
    ang= math.acos((DH**2 - distance**2 + HA**2)/(2*DH*HA))
    if not radians:
        ang=ang*(180/math.pi)
    return ang

def radianToAngle(radian):
    return radian*(180/math.pi)


"""
TODO:
*use chains as parameters and not 'a' 'b'. started as pseudoChain

improvments:
* AVAL and BVAL smart selection
"""
class PDBReader:

    @staticmethod
    def readFile(path, interfaceParts=None):
        name = os.path.basename(path)[0:4]
        pdbFile = file(path)
        atoms = []
        waters = []
        pH = None
        cmpnds=[]
        cmpndChains=[]
        hetAtms=[]
        logging.info('reading pdb file (atoms and HETATM of HOH) %s', path)
        #read file for chain atoms and water
        for line in pdbFile.readlines():
            if line[0:6] == 'ATOM  ' or line[0:6] == 'HETATM':
                atomNum, atomSymbol = line[6:11].strip(),line[12:16].strip()
                residue, chain,resId, x, y, z, atomType = line[17:20].strip(),line[21:22],line[22:26].strip(),line[30:38],line[38:46],line[46:54],line[77:78]
                tempFactor = float(line[60:66])

                #skip: when two rotamers appear, such as AVAL and BVAL we ignore the BVAL
                if line[16:17]!=' ' and line[16:17]!='A':
                    continue
                #skip on multimodel. see 1EAW Res 60 chain A
                if line[26:27]!=' ':
                        continue
                coord = (float(x), float(y), float(z))
                if line[0:6]=='ATOM  ':

                    atoms.append( atom(atomNum, atomSymbol, residue, chain,resId, atomType, coord, tempFactor=tempFactor ) )
                elif line[0:6]=='HETATM':# and residue.strip()=='HOH':
                    if residue=='HOH':
                        waters.append( water( atomNum, atomSymbol, residue, chain,resId, atomType, coord) )
                    else:
                        hetAtms.append( atom(atomNum, atomSymbol, residue, chain,resId, atomType, coord, tempFactor=tempFactor ) )
            elif line[0:6] == "COMPND" and 'MOLECULE:' in line:
                cmpnds.append(line.split('MOLECULE:')[1].strip('; \n').replace(',',' '))
            elif line[0:6] == "COMPND" and 'CHAIN:' in line:
                cmpndChains.append(line.split('CHAIN: ')[1].strip('; \n').replace(', ',''))
            elif line[0:14] == "REMARK 200  PH":
                pH = line[45:50]
            elif line[0:3] == "END" and len(line)==4:
                break
            elif line[0:14] == "MODEL        2":
                print 'model2'
                logging.info('using model 1 (more models are ignored)')
                break
        logging.info('finished reading file. waters:' + str(len(waters)) + ' atoms: ' +str(len(atoms)))
        if interfaceParts is None:
            interfaceParts = cmpndChains
        return PDBReader(name, atoms, waters, cmpnds, hetAtms, interfaceParts=interfaceParts)

    """ class that handles PDB files and utilities """
    PDB_REGEX = re.compile('(?:ATOM|HETATM)\s*([0-9]+)\s+([\w\d]{,4})\s*(\w{2,4})\s+(\w)\s+(\d+)\s+(-?[0-9\\.]+)\s*(-?[0-9\\.]+)\s*(-?[0-9\\.]+)\s*.*?([A-Z])')
    def __init__(self, name, atoms, waters, compunds,hetAtms, interfaceParts=['A','B']):
        self.compunds = string.join(compunds,' - ')
        self.interfaceParts = interfaceParts
        self.name = name
        self.atoms = atoms
        self.waters = waters
        self.hetatms = hetAtms
        self.chains = set([a.chain for a in atoms])
        self.interfaceCache=None
        self.cacheDistance = 0
        lastAtom = None
        self.__buildIndex()

    def __buildIndex(self):
        """ Init the internal indexs for the atoms and their pseudoChain, and k-tree """
        logging.debug('building k-tree')
        self.ktree = KDTree.construct_from_data(self.atoms[:])
        logging.debug('end building k-tree')
        logging.debug('building indexs')
        for i,a in enumerate(self.atoms):
            a.atomIndex=i
            #assign pseudo chain to residue which symbols the interface chains (example: A:B)
            for j, interPart in enumerate(self.interfaceParts):
                if a.chain in interPart:
                    a.pseudoChain=chr(ord('A')+j)

        logging.debug('end building index')

    def getFile(self, name):
        if not os.path.exists("./debug/"):
            os.makedirs("./debug/")

        return file("./debug/"+self.name+name,"w")

    def getPseudoChains(self):
        for j, interPart in enumerate(self.interfaceParts):
            yield chr(ord('A')+j)

    def getInterface(self, maxDistance=NEIGHBOR_DISTANCE):
        if self.interfaceCache!=None and self.cacheDistance==maxDistance:
            return self.interfaceCache

        """
        Get atoms not from same chain, having distance less than maxDistance ignores H
        """
        i= time.clock()
        self.cacheDistance=maxDistance
        self.interfaceCache=set()
        interfacesT=self.interfaceCache
        #we assume dimers

        for atom in [atom for atom in self.atoms if atom.pseudoChain=='A']:
            for atom2 in self.ktree.findByDistance(query_point=atom.coord, distance=maxDistance):
                if atom2.pseudoChain==atom.pseudoChain or atom2.pseudoChain==None:
                    continue
                interfacesT.add(atom)
                interfacesT.add(atom2)

        #extend with H
        print 'interface atoms:',len(interfacesT)
        print 'finished interfaces', (time.clock()-i)
        return self.interfaceCache
        """
        i= time.clock()
        interfaces=dict()
        interfaces['A']=[]
        interfaces['B']=[]
        print 'start interfaces'
        for atom, atom2 in ([(atomA,atomB) for atomA,atomB in combinations(self.atoms,2) if atomFilter(atomA,atomB)]):
            if atom.distance(atom2)<maxDistance:
                interfaces[atom.chain].append(atom.resId)
                interfaces[atom2.chain].append(atom2.resId)
                atom.interface.append(atom2)
                #print atom
                atom2.interface.append(atom)
        print 'interface atoms:',len(interfaces['A'])
        i= time.clock()
        print 'finished interfaces', (time.clock()-i)
        return interfaces
        """


    def atoms(self):
            return self.atoms

    def getNextAtoms(self,atom,i):
        if atom.atomIndex+i>=len(self.atoms):
            logging.debug('error getting atom %s',atom)
            return atom
        return self.atoms[atom.atomIndex+i]

def main():
    #a=PDBFile(u"1AY7_FH.pdb")
    #a=PDBFile(u"./pdbs/1RV6_FH.pdb",interfaceParts=['VW','X'])
    #a.getInterface()
    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)
    #a=PDBReader.readFile(u"./unboundPdbs/1EGL_FH.pdb")

    #a=PDBReader.readFile(u"./unboundPdbs/1BU6_FH.pdb")
    #a=PDBReader.readFile(u"./pdbs/1EAW_FH.pdb")
    #PDBReader.readFile(u"./pdbs/3SGB_FH.pdb")
    print 'READ'
    #print len(a.atoms)
    #print len([at for at in a.atoms if at.chain=='A' and at.resId==60])

    #a=PDBFile(u"3SZJ.pdb")
    #a=PDBFile(u"2VDB_FH.pdb")
    #a=PDBFile(u"1R0R_FH.pdb", interfaceParts=['E','I'] )
    #a.h2OStuff()
    #a.hbondsStuff()

if __name__ == "__main__":
    main()
