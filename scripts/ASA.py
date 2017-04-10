import math
import logging
import string

from pdbReader import PDBReader
from resReduce import mapUnbound


R_WATER = 1.4
#signtificant distance
KNOWN_RADIUS = {
    'C': 1.700,
    'N': 1.625,
    'O': 1.480,
    'S': 1.782,
    'H': 1.000,
    'X': 1.700 #for undetermined atom assume Carbon
}
SPIRAL_POINTS = 20

"""
Original script:
*find atoms on interface having distance less than 25
*foreach atom:
        initate spiral according to atom
        foreach atomB that is 25 from a
                if distance to one of the spiral points is less than ATOMRADIUS+1.4**2
                        clash point (v=0)
        foreach point on spiral:
                if not clash
                        qpoint++
                        print point
                        calc area
TODO:
        preformence improve - when using kdtree use the radio_atom+R_WATER+1.782 instead of 5
        clashs - removed. just to check for sure
"""


class ASA(object):
    def __init__(self, pdb, unboundA=None, unboundB=None):
        """
                Calculates asa for pdb. you can specify pdbA and pdbB as unbound version
                otherwise it just extract it from the pdb
                """
        self.pdb = pdb
        if unboundA != None:
            self.pdbA = unboundA[0]
            self.pdbB = unboundB[0]
            self.chainA = unboundA[1]
            self.chainB = unboundB[1]
        else:
            self.pdbA = None
            self.pdbB = None
            self.chainA = None
            self.chainB = None
        self.interASA = dict()            #calculated from complex
        self.intraASA = dict()            #calculated from complex
        self.diffASAperAtom = dict()      #calculated from complex
        self.interPerAtom = dict()        #calculated from complex

        #not used now
        self.unboundASA = dict()  #unbound ASA per chain
        self.intraASAUnB = dict() #calculated from unbound
        self.diffASAperAtomUnbound = dict()
        #results
        self.interfaceASA = 0            #ASA in interface bound against unbound
        self.combinedASA = 0            #TOTAL diff ASA (intra-inter from complex)
        self.boundVsUnboundASA = 0      #NOT USED
        self.unboundIntraTotal = 0      #NOT USED
        self.unboundMisses = 0
        #for outputing
        self.perAtom = self.pdb.getFile('per-atom.pdb')
        self.q_points = 0
        self.interfaceASAComp = {'A': 0, 'B': 0}

    def execute(self):
        #calculate asa only for interface chains
        for chain in string.join(self.pdb.interfaceParts, ''):
            logging.info('calculating ASA for chain %s', chain)
            self.interASA[chain], self.intraASA[chain] = self.mesh(chain)

        intraTotal = sum([sum([self.intraASA[chain] for chain in part]) for part in self.pdb.interfaceParts])
        interTotal = sum([sum([self.interASA[chain] for chain in part]) for part in self.pdb.interfaceParts])
        self.combinedASA = intraTotal - interTotal

        if self.pdbA != None and self.pdbB != None:
            self.calcInterfaceASA()
        """
                # global ASA bound vs unbound
                if self.pdbA!=None and self.pdbB!=None:
                        for i, pseudo in enumerate(self.pdb.getPseudoChains()):
                                interTransformed=dict()
                                area = self.meshUnbound(pseudo,interTransformed)
                                print 'pseudo chain',pseudo,' Area: ',area
                                self.unboundASA[pseudo]=area
                                self.unboundIntraTotal+=area
                        self.boundVsUnboundASA=self.unboundIntraTotal-interTotal
                """
        self.output()
        self.perAtom.close()

    @staticmethod
    def calcForUnboundPDB(pdb):
        asaPerAtom = dict()
        interstringChains = string.join(pdb.interfaceParts, '')
        calcASAforAtom = ASA.calcASAforAtom
        for chain in pdb.interfaceParts[0]:
            for atom, neigbors in ASA.nearAtomsOnChain(chain, pdb.atoms, pdb.ktree):
                intraNeigbors = [a for a in neigbors if a.pseudoChain == atom.pseudoChain]
                intraArea = calcASAforAtom(atom, intraNeigbors, False)
                asaPerAtom[atom] = intraArea

        return asaPerAtom

    def calcInterfaceASA(self):
        diffASAA, missesA = self.calcASAUnboundMissingReduced('A', self.pdbA, self.chainA)
        diffASAB, missesB = self.calcASAUnboundMissingReduced('B', self.pdbB, self.chainB)

        self.interfaceASAComp['A'] = (diffASAA, missesA)
        self.interfaceASAComp['B'] = (diffASAB, missesB)
        self.interfaceASA = diffASAA + diffASAB
        self.unboundMisses = missesA + missesB
        logging.info('INTERFACE ASA: %s', str(self.interfaceASA))
        return


    def calcASAUnboundMissingReduced(self, pseudoChain, unboundPdb, unboundChain):
        """
                calc interface ASA for bound-unbound
                interface: ASA=0 in inter
                """

        ktree = unboundPdb.ktree
        interPerAtom = self.interPerAtom

        def unboundNeigbors(atom):
            return [atom2 for atom2 in ktree.findByDistance(query_point=atom.coord, distance=(radio_atom(
                atom.atomType) + R_WATER * 2 + 1.8) ** 2) if atom2 != atom
                    and atom2.chain in unboundChain]

        mapping = mapUnbound(self.pdb, unboundPdb, self.pdb.interfaceParts[0 if pseudoChain == 'A' else 1],
                             unboundChain)
        missInterface = self.pdb.getFile("missingASAInterface{0}.txt".format(unboundPdb.name))
        diffASA, misses = 0, 0
        calcASAforAtom = ASA.calcASAforAtom
        for aBound, asa in self.diffASAperAtom.iteritems():
            if asa > 0 and aBound.pseudoChain == pseudoChain:
                if aBound not in mapping:
                    print>> missInterface, '{0:^4}{1:^6}{2:^4}'.format(aBound.chain, aBound.resId, aBound.symbol)
                    logging.warn('interface atom isnt mapped %s', aBound)
                    misses += 1
                    continue
                aUnbound = mapping[aBound]
                neigbors = unboundNeigbors(aUnbound)
                intraArea = calcASAforAtom(aUnbound, neigbors, False)
                interArea = interPerAtom[aBound]
                diffASA += (intraArea - interArea)
        missInterface.close()
        logging.info('INTERFACE ASA: %s', str(self.interfaceASA))
        return diffASA, misses

    def completeASA(self, missingAtoms):
        """
            Completes ASA for pdbUnbound based on missingAtoms

            Parameters:
            pdb - pdb for the complex
            pdbUnbound - pdb class for unbound
            missingAtoms - list of atoms, from complex on the same chain

            Returns:
                ASA to add the unbound
                
            How it works:
            On a specific residue, calcuates the ASA from the first missing atom to
            the last atom in this residue
            The calculation is done on a residue alone (not dependent on neigbors residues)
            """

        #first group missing atoms by res
        grouped = dict()
        for a in missingAtoms:
            if not grouped.has_key(a.resId):
                grouped[a.resId] = []
            grouped[a.resId].append(a)

        asaToAdd = 0
        pdb = self.pdb
        calcASAforAtom = ASA.calcASAforAtom
        for resId, missAtoms in grouped.iteritems():
            missAtoms.sort(key=lambda x: x.index)
            #first atom that appear iin both complex and unbound on the same residue
            firstNonMissing = pdb[missAtoms[0].index - 1]

            tail = [a for a in pdb.atoms if
                    a.resId == resId and a.index > firstNonMissing and a.chain == firstNonMissing.chain]
            neigb = ASA.nearAtomsOnChain(firstNonMissing.chain, pdb.atoms, pdb.kdtree)
            neigb = [a for a in neigb if a.resId == resId and a.index < firstNonMissing.index]
            #remove asa for the last atom in the residue
            asaToAdd -= calcASAforAtom(firstNonMissing, neigb, False)
            #add to asa area of the remaining tail
            for atom in tail:
                neigb = ASA.nearAtomsOnChain(atom.chain, pdb.atoms, pdb.kdtree)
                neigb = [a for a in neigb if a.resId == resId]
                asaToAdd += calcASAforAtom(atom, neigb, False)
        return asaToAdd

    def output(self):
        for chain in string.join(self.pdb.interfaceParts, ''):
            logging.info('-------------')
            logging.info('ASA for chain %s', chain)
            logging.info('--inter: %s', self.interASA[chain])
            logging.info('--intra: %s', self.intraASA[chain])
            logging.info('diff ASA: %s', self.intraASA[chain] - self.interASA[chain])

        interfaceDistance, interfaceASA = dict(), dict()#for debugging output the interface with asa and with distance
        interfaces = self.pdb.getInterface()
        for atom in interfaces:
            if atom.chain not in interfaceDistance:
                interfaceDistance[atom.chain] = set()
            interfaceDistance[atom.chain].add(str(atom.resId))

        for atom in [a for a in self.pdb.atoms if (a in self.diffASAperAtom) and self.diffASAperAtom[a] > 0]:
            if atom.chain not in interfaceASA:
                interfaceASA[atom.chain] = set()
            interfaceASA[atom.chain].add(str(atom.resId))

        selDist, selASA = [], []
        for chain, reses in interfaceDistance.iteritems():
            selDist.append(' res ' + string.join(reses, '+') + ' and chain ' + chain)
        for chain, reses in interfaceASA.iteritems():
            selASA.append(' res ' + string.join(reses, '+') + ' and chain ' + chain)
        interfaceOutput = self.pdb.getFile('.interface.pml')
        print>> interfaceOutput, 'sel ' + string.join(selDist, ' or ')
        print>> interfaceOutput, 'show lines,sele'
        print>> interfaceOutput, 'set_name sele,Distance-int'
        print>> interfaceOutput, 'sel ' + string.join(selASA, ' or ')
        print>> interfaceOutput, 'show lines,sele'
        print>> interfaceOutput, 'set_name sele,ASA-int'
        interfaceOutput.close()

    def meshUnbound(self, pseudoChain, interTransformed):
        pdb = None
        pChain = None
        intraASA = 0
        if pseudoChain == 'A':
            pdb = self.pdbA
            pChain = self.chainA
        else:
            pdb = self.pdbB
            pChain = self.chainB

        calcASAforAtom = ASA.calcASAforAtom
        for atom in [atom for atom in pdb.atoms if atom.chain in pChain]:
            neigbors = [atom2 for atom2 in pdb.ktree.findByDistance(query_point=atom.coord, distance=25)
                        if atom2 != atom and atom2.chain in pChain]

            intraArea = calcASAforAtom(atom, neigbors, False)
            self.intraASAUnB[(atom.chain, atom.resId, atom.symbol)] = intraArea
            #TODO calc something like diffASAperAtom
            """
                        origChain=chainMapping[atom.chain]
                        if (origChain,atom.resId,atom.symbol) not in interTransformed:
                                print 'not found', origChain,atom.resId,atom.symbol
                        else:
                                print 'found'
                                inter=interTransformed[(origChain,atom.resId,atom.symbol)]
                                self.diffASAperAtomUnbound[(origChain,atom.resId,atom.symbol)]=intraArea-inter
                        """
            intraASA += intraArea
        return intraASA

    def mesh(self, chain):
        interASA, intraASA = 0, 0
        interstringChains = string.join(self.pdb.interfaceParts, '')
        pdb = self.pdb
        calcASAforAtom = ASA.calcASAforAtom
        for atom, neigbors in ASA.nearAtomsOnChain(chain, pdb.atoms, pdb.ktree):
            #pseudo chain or chain?
            intraNeigbors = [a for a in neigbors if a.pseudoChain == atom.pseudoChain]
            interOptimized = [a for a in neigbors if a.chain in interstringChains]
            intraArea = calcASAforAtom(atom, intraNeigbors, True, printPoints=self.print_point)
            interArea = calcASAforAtom(atom, interOptimized, False)

            self.interPerAtom[atom] = interArea
            self.diffASAperAtom[atom] = intraArea - interArea
            interASA += interArea
            intraASA += intraArea
        return interASA, intraASA

    def getASA(self, atom):
        return self.diffASAperAtom[atom]

    def isBuried(self, atom):
        return self.interPerAtom[atom] == 0

    @staticmethod
    def calcASAforAtom(atom, neigbors, toPrint, printPoints=None):
        r = radio_atom(atom.atomType) + R_WATER
        spiralPoints = list(spiral(r, atom))
        #clashs=[]
        for atom2 in neigbors:
            for x, y, z in spiralPoints[:]:
                r2 = radio_atom(atom2.atomType)
                dist = atom2.distanceFromXYZ((x, y, z))
                if dist < (r2 + R_WATER) ** 2:#water is added to two sides of the equation?
                    #clashs.append((x,y,z))
                    """
                                        if atom.x==14.235 and abs(x-15.1311577535)<0.1:
                                                print (x,y,z),atom2.coord,'|',(r2 + R_WATER)**2
                                                print r2
                                                print '{0:.7f},{1:.7f}'.format(dist,(r2 + R_WATER)**2)
                                        """
                    spiralPoints.remove((x, y, z))
        if toPrint:
            for x, y, z in spiralPoints[:]:
                printPoints(x, y, z)
                #self.print_point( x, y, z )

        area = area_calc(r, len(spiralPoints), SPIRAL_POINTS)
        return area

    @staticmethod
    def nearAtomsOnChain(chain, atoms, ktree):
        extraRad = R_WATER * 2 + 1.8
        for atom in [atom for atom in atoms if atom.chain == chain]:
            neigbors = []
            for atom2 in ktree.findByDistance(query_point=atom.coord,
                                              distance=(radio_atom(atom.atomType) + extraRad) ** 2):#25
                if atom2 != atom:
                    neigbors.append(atom2)
            yield atom, neigbors

    def print_point(self, x, y, z):
        self.q_points += 1
        toPrint = 'HETATM {0:^4}  O   HOH I {1:^4} {2:^11.3f} {3:^4.3f} {4:^8.3f}  0.83 56.58           O '.format(
            self.q_points, self.q_points, x, y, z)
        print>> self.perAtom, toPrint


def radio_atom(atom):
    if atom == 'X':
        logging.error('undetermined atom while calcing ASA. Assuimg C')
    return KNOWN_RADIUS[atom]


def area_calc(x, y, z):
    return (4 * math.pi * (x) ** 2) * y / (z);


def spiral(r, atom, nPoints=SPIRAL_POINTS):
    move = 0;
    z = 1 - (1.0 / nPoints)
    zMove = (float(2) / nPoints)
    mvPlus = ( math.pi * ( 3 - math.sqrt(5) ) )
    for k in range(0, nPoints):
        rUnit = math.sqrt(1 - z * z);
        xCoord = ( math.cos(move) * rUnit ) * r + atom.x
        yCoord = ( math.sin(move) * rUnit ) * r + atom.y
        zCoord = z * r + atom.z
        move = move + mvPlus
        z -= zMove
        yield ( xCoord, yCoord, zCoord )


res = None


def main():
    import resReduce

    resReduce.PRINT_ALIGNED = True
    global res
    print 'ASA script'
    import time

    i = time.clock()
    pdb = PDBReader.readFile(u"./pdbs/1FQJ_FH.pdb", interfaceParts=['A', 'B'])
    pdbA = PDBReader.readFile(u"./unboundPdbs/1TND_FH.pdb")
    pdbB = PDBReader.readFile(u"./unboundPdbs/1FQI_FH.pdb")
    asa = ASA(pdb, unboundA=(pdbA, "C"), unboundB=(pdbB, "A"))

    asa.execute()
    print 'TOOK ', (time.clock() - i)

    print asa.combinedASA
    print asa.boundVsUnboundASA
    res = asa
    """
        How ASA can change dramitcly in no known interface (B54 in 1A2K)
        for atom, diff in res.diffASAperAtom.iteritems():
	if (atom.chain,atom.resId,atom.symbol) in res.diffASAperAtomUnbound:
		if abs(diff-res.diffASAperAtomUnbound[(atom.chain,atom.resId,atom.symbol)])>10:
			print atom
	"""


if __name__ == "__main__":
    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)
    main()
