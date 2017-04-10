import string
import sys
import logging
from algorithms.NWAlignment import NWAlignment

from pdbReader import PDBReader
from atom import ResiduesCodes

#sys.path.append( 'algorithms' )
#from algorithms.NWAlignment import NWAlignment
from databaseReader import readDatabase

PRINT_ALIGNED=False

def getResidues(pdbBound, pdbUnbound, chainBound, chainUnbound):
    resBound, resUnbound = [], []
    for chain in chainBound:
        resBound += [ResiduesCodes[a.residue if a.residue!='GLX' else 'GLN'] for a in pdbBound.atoms if a.chain == chain and a.symbol=='CA']
    for chain in chainUnbound:
        resUnbound += [ResiduesCodes[a.residue if a.residue!='GLX' else 'GLN'] for a in pdbUnbound.atoms if a.chain == chain and a.symbol=='CA']
    return resBound, resUnbound

def getAtoms(pdbBound, pdbUnbound, chainBound, chainUnbound):
    resBound, resUnbound = [], []
    for chain in chainBound:
        resBound += [a for a in pdbBound.atoms if a.chain == chain]
    for chain in chainUnbound:
        resUnbound += [a for a in pdbUnbound.atoms if a.chain == chain]
    return resBound, resUnbound




def alignResidues(bound, unbound):
    global PRINT_ALIGNED
    seq1 = string.join(bound,'')
    seq2 = string.join(unbound,'')
    align = NWAlignment(seq1,seq2)
    align.fillIn()
    nSeq1, nSeq2 = align.getTraceback()
    if PRINT_ALIGNED:
        printAligned( string.join(nSeq1,''), string.join(nSeq2,'') )
    isStart=False
    resIdA, resIdB = -1, -1

    missesA, missesB = [], []
    for a,b in zip(nSeq1,nSeq2):
        if a!='-':
            resIdA+=1
        if b!='-':
            resIdB+=1

        if a != '-' and b!='-':
            isStart = True
        
        if a == '-':
             missesB.append(resIdB)
        if b == '-':
            missesA.append(resIdA)
            
    return missesA, missesB

def printAligned(seqA, seqB):
    linesLength = 100
    for i in range(0,min(len(seqA), len(seqB)),linesLength):
        print 'A  '+seqA[i:i+linesLength]
        print 'B  '+seqB[i:i+linesLength]
        
def redcueAtoms(pdbBound, pdbUnbound, chainBound, chainUnbound):
    """
    Returns atoms of bound align to unbound
    When there is a miss within the range it is removed
    """
    resBound, resUnbound = getResidues(pdbBound, pdbUnbound, chainBound, chainUnbound)
    missesBound, missesUnbound = alignResidues(resBound, resUnbound)
    boundAtoms, unboundAtoms = getAtoms(pdbBound, pdbUnbound, chainBound, chainUnbound)

    mappingB= dict(((a.chain,a.resId),i) for i,a in enumerate([a for a in boundAtoms if a.symbol=='CA']))
    boundReduce = [a for a in boundAtoms if (a.chain,a.resId) in mappingB and mappingB[(a.chain,a.resId)] not in missesBound]
    for a in boundAtoms:
        if (a.chain,a.resId) not in mappingB:
            logging.warn('WARNING missing CA in ', a)
    mappingUn = dict(((a.chain,a.resId),i) for i,a in enumerate([a for a in unboundAtoms if a.symbol=='CA']))
    unboundReduce = [a for a in unboundAtoms if (a.chain,a.resId) in mappingUn and mappingUn[(a.chain,a.resId)] not in missesUnbound]
    for a in unboundAtoms:
        if (a.chain,a.resId) not in mappingUn:
            logging.warn('WARNING missing CA in ', a)
    return boundReduce, unboundReduce

def mapUnbound(pdbBound, pdbUnbound, chainBound, chainUnbound):
    """
        return a dictionary mapping between atom bound to atom unbound
    """
    boundReduce, unboundReduce = redcueAtoms(pdbBound, pdbUnbound, chainBound, chainUnbound) 
    mapping=dict()
    unboundReduce.append(None)
    unboundIter=iter(unboundReduce)
    boundRes, unBoundRes = -1, -1
    residueUnbound = []
    y = unboundIter.next()
    for x in boundReduce:
        if boundRes!=x.resId and y!=None:
            boundRes = x.resId 
            residueUnbound=[y]
            unBoundRes = y.resId
            y = unboundIter.next()
            while y!=None and y.resId==unBoundRes:
                residueUnbound.append(y)
                y = unboundIter.next()

        for aUnbound in residueUnbound:
            if aUnbound.symbol==x.symbol:
                mapping[x]=aUnbound
    return mapping


def main():
    calls = readDatabase('Database.csv')
    #1FC2_FH
    for path, chains, unboundA, unboundB, kd in calls:
        pdb = PDBReader.readFile(path,interfaceParts = chains.split(':'))
        pdbB = PDBReader.readFile('./unboundPdbs/{0}_FH.pdb'.format(unboundB.split('_')[0]))
        chainBound = chains.split(':')[1]
        chainUnbound = unboundB.split('_')[1]
        print unboundB.split('_')[0]
        mapUnbound( pdb, pdbB, chainBound, chainUnbound )
    print 'finished'

TEST_MAPPING = None
def test():
    global TEST_MAPPING, PRINT_ALIGNED
    from pdbReader import PDBReader
    PRINT_ALIGNED=True
    pdb=PDBReader.readFile(u"./pdbs/1FQJ_FH.pdb",interfaceParts=['A','B'])
    pdbA=PDBReader.readFile(u"./unboundPdbs/1TND_FH.pdb")
    TEST_MAPPING = mapUnbound(pdb, pdbA, "A", "C")
    
    print 'done'
    
if __name__ == '__main__':
    #main()
    test()
