import os
import string

from pdbReader import PDBReader
from atom import ResiduesCodes
from algorithms.NWAlignment import NWAlignment


def readDatabase( csvName ):
    """
    for each record returns (path,chains,unboundA,unboundB,kD)

    Reads database csv file, checks for existing PDB and unbound pdbs
    Assumes:
        ./unbondPdbs/ID_FH.pdb - for each ID unbound pdb
        ./pdbs/ID_FH.pdb - for each ID bound pdb

    csvName - database filename

    """
    calls = []
    command = file( csvName )
    missingFiles, missedUnbound = [], []
    for line in command.readlines():
            protinInput=line.split(',')[0]
            path='./pdbs/{0}_FH.pdb'.format(protinInput.split('_')[0])
            if not os.path.isfile(path):
                missingFiles.append(protinInput.split('_')[0])
                continue
            unboundA=line.split(',')[2]
            unboundB=line.split(',')[4]
            pathA='./unboundPdbs/{0}_FH.pdb'.format(unboundA.split('_')[0])
            pathB='./unboundPdbs/{0}_FH.pdb'.format(unboundB.split('_')[0])
            if not os.path.isfile(pathA):
                missedUnbound.append(pathA.split('_')[0])
                continue
            if not os.path.isfile(pathB):
                missedUnbound.append(pathB.split('_')[0])
                continue
            chains=protinInput.split('_')[1]
            kD = line.split(',')[7]
            calls.append((path,chains,unboundA,unboundB,kD))

    if len(missingFiles)>0 or len(missedUnbound)>0:
        print 'The following files are missing:'
        print string.join(missingFiles,'\n')
        print 'The following unbound files are missing:'
        print string.join(missedUnbound,'\n')
        if raw_input('Continue (y to continue)? \n')!="y":
            print 'Terminated'
            raise Exception("Missing files")
    return calls



test=None
   
def getResidues(pdb, pdbA, pdbB, chains, chainA, chainB, bone = False):
    residuesA = []
    for chain in chains.split(':')[0]:
        residuesA += [a.residue for a in pdb.atoms if a.chain == chain and (not bone or a.symbol=='CA')]
    residuesB = []
    for chain in chains.split(':')[1]:
        residuesB += [a.residue for a in pdb.atoms if a.chain == chain and (not bone or a.symbol=='CA')]
    
    unboundA, unboundB = [], []
    for chain in chainA:
        unboundA += [a.residue for a in pdbA.atoms if a.chain == chain and (not bone or a.symbol=='CA')]
    for chain in chainB:
        unboundB += [a.residue for a in pdbB.atoms if a.chain == chain and (not bone or a.symbol=='CA')]
    return residuesA, residuesB, unboundA, unboundB

def getAtoms(pdb, pdbA, pdbB, chains, chainA, chainB):
    residuesA = []
    for chain in chains.split(':')[0]:
        residuesA += [a for a in pdb.atoms if a.chain == chain]
    residuesB = []
    for chain in chains.split(':')[1]:
        residuesB += [a for a in pdb.atoms if a.chain == chain]
    
    unboundA, unboundB = [], []
    for chain in chainA:
        unboundA += [a for a in pdbA.atoms if a.chain == chain]
    for chain in chainB:
        unboundB += [a for a in pdbB.atoms if a.chain == chain]
    return residuesA, residuesB, unboundA, unboundB

def findDiffForSameLength(residuesA, residuesB, unboundA, unboundB):
    isAEq, isBEq = True, True
    
    for a,b in zip(residuesA, unboundA):
        if a.residue!=b.residue:
            print 'a:',a,'|',b
            isAEq=False
    
    for a,b in zip(residuesB, unboundB):
        if a.residue!=b.residue:
            print 'b:', a,'|',b
            isBEq=False
    
    return isAEq and isBEq

def reduceChain(atoms, start, end):
    cas = [a for a in atoms if a.symbol=='CA']
    mapping = dict((a.resId,i) for i,a in enumerate(cas))
   
    return [a for a in atoms if mapping[a.resId]>=start and mapping[a.resId]<end]

def check2( path, chains, unboundA, unboundB ):
    pdb = PDBReader.readFile(path,interfaceParts = chains.split(':'))
    pdbA = PDBReader.readFile('./unboundPdbs/{0}_FH.pdb'.format(unboundA.split('_')[0]))
    pdbB = PDBReader.readFile('./unboundPdbs/{0}_FH.pdb'.format(unboundB.split('_')[0]))
    chainA = unboundA.split('_')[1]
    chainB = unboundB.split('_')[1]
    boundA, boundB, unboundA, unboundB = getResidues(pdb, pdbA, pdbB, chains, chainA, chainB, bone=True)
    #to 1 char code
    for r in boundA+boundB:
        if r not in ResiduesCodes:#GLX
            print 'ERROR'
            return pdb.name,[], [], [], []
    boundA = [] + [ResiduesCodes[r] for r in boundA]
    boundB = [] + [ResiduesCodes[r] for r in boundB]
    unboundA = [] + [ResiduesCodes[r] for r in unboundA]
    unboundB = [] + [ResiduesCodes[r] for r in unboundB]
    print 'for ' , pdb.name
    mutationA, missesA, startA, endA = alignResidues(boundA, unboundA)
    mutationB, missesB, startB, endB = alignResidues(boundB, unboundB)
    print missesB
    boundAtoms, boundBatoms, unboundAatoms, unboundBatoms = getAtoms(pdb, pdbA, pdbB, chains, chainA, chainB)
    print string.join([ResiduesCodes[a.residue] for a in reduceChain(unboundAatoms, startA, endA) if a.symbol=='CA'],'')
    
    return pdb.name,mutationA, missesA, mutationB, missesB

def printAligned(seqA, seqB):
    linesLength = 100
    for i in range(0,min(len(seqA), len(seqB)),linesLength):
        print 'A  '+seqA[i:i+linesLength]
        print 'B  '+seqB[i:i+linesLength]

def alignResidues(bound, unbound):
    global DEBUG
    seq1 = string.join(bound,'')
    seq2 = string.join(unbound,'')

    #print 'before alignment:'
    #printAligned( seq1, seq2 )
    align = NWAlignment(seq1,seq2)
    align.fillIn()
    nSeq1, nSeq2 = align.getTraceback()
    #print 'After alignment:'
    if DEBUG:
        printAligned( string.join(nSeq1,''), string.join(nSeq2,'') )

    isStart=False
    tmpMiss = []
    misses, mutations = [], []
    resIdA = 0
    resIdB = 0
    startIndex, endIndex = -1, 0
    for a,b in zip(nSeq1,nSeq2):
        endIndex += 1
        if isStart==False:
                startIndex+=1
        resIdA+=1 if a!='-' else 0
        resIdB+=1 if b!='-' else 0
        if a == '-' or b == '-':
            tmpMiss.append(1)
            if isStart:
                if a == '-':
                    tmpMiss.append(resIdB)
                else:
                    #print 'BOUND HAS MISSING!!!', resIdA
                    tmpMiss.append(resIdA)
        elif a!=b:
            mutations.append((resIdA, resIdB))
        else:
            isStart = True
            if len(tmpMiss)>0:
                misses += tmpMiss
                tmpMiss = []
    endIndex = endIndex - len(tmpMiss)
    if DEBUG==True:
        print startIndex,endIndex
    return mutations, misses, startIndex, endIndex

DEBUG=False
def main():
    calls = readDatabase('Database.csv')
    qa = file('databaseQA.csv','w')
    print>>qa,'pdb,mutationA, missesA, mutationB, missesB'
    for path, chains, unboundA, unboundB, kd in calls:
        pdb,mutationA, missesA, mutationB, missesB = check2( path, chains, unboundA, unboundB )
        toPrint=string.join([pdb,str(len(mutationA)), str(len(missesA)), str(len(mutationB)), str(len(missesB))],',')
        print>>qa,toPrint
    qa.close()
    print 'finished'

    """
    different order
    if pdb.name not in ['1BVK', '1EAW']:
         findDiffForSameLength(pdb, pdbA, pdbB, chains, chainA, chainB)
        return
    """
        
    """
    Checks the PDB repository whether the chains are complete
    and have the same residues
    """
    
if __name__ == '__main__':
    main()




"""
def check( path, chains, unboundA, unboundB ):
    pdb = PDBReader.readFile(path,interfaceParts = chains.split(':'))
   
    pdbA = PDBReader.readFile('./unboundPdbs/{0}_FH.pdb'.format(unboundA.split('_')[0]))
    pdbB = PDBReader.readFile('./unboundPdbs/{0}_FH.pdb'.format(unboundB.split('_')[0]))
    chainA = unboundA.split('_')[1]
    chainB = unboundB.split('_')[1]
    residuesA, residuesB, unboundA, unboundB = getResidues(pdb, pdbA, pdbB, chains, chainA, chainB)
    print 'for ', pdb.name
    if len(residuesA) == len(unboundA) and len(residuesB) == len(unboundB):
        if findDiffForSameLength(residuesA, residuesB, unboundA, unboundB):
            print 'OK! ', path
        else:
            print 'Same length, but different residues'
    elif len(chains.split(':')[0])==len(chainA) and len(chains.split(':')[1])==len(chainB):
        findDiffForSameChainNums(pdb, pdbA, pdbB, chains, chainA, chainB)
    else:
        print pdb.name
        print 'Not same number of chains'

def findDiffForSameChainNums(pdb, pdbA, pdbB, chains, chainA, chainB):
    #Finds diffs for chain by chain 
    def reduceChains( chainMap, oriAtoms, newAtoms ):
        for oriChain, newChain in chainMap.iteritems():
            oriResIds= [] + [a.resId for a in oriAtoms if a.chain==oriChain]
            newResIds= [] + [a.resId for a in newAtoms if a.chain==newChain]
            minOri, maxOri = min(oriResIds), max(oriResIds)
            minNew, maxNew = min(newResIds), max(newResIds)
            minRes=max(minNew, minOri)
            maxRes=min(maxNew, maxOri)
            print 'reducing chain',oriChain, '({0})'.format(newChain), minRes, ':', maxRes
            oriReduced =[] + [a for a in oriAtoms if a.chain==oriChain and a.resId>minRes and a.resId<maxRes  and a.atomType!='H']
            newReduced =[] + [a for a in newAtoms if a.chain==newChain and a.resId>minRes and a.resId<maxRes  and a.atomType!='H']
            if len(newReduced)==0:
                global test
                test=newAtoms
                print len(newAtoms)
                raw_input("Error?")
            yield oriReduced, newReduced

    def reduceBackbonesChains( chainMap, oriAtoms, newAtoms ):
        for oriChain, newChain in chainMap.iteritems():
            oriResIds= [] + [a.resId for a in oriAtoms if a.chain==oriChain]
            newResIds= [] + [a.resId for a in newAtoms if a.chain==newChain]
            minOri, maxOri = min(oriResIds), max(oriResIds)
            minNew, maxNew = min(newResIds), max(newResIds)
            minRes=max(minNew, minOri)
            maxRes=min(maxNew, maxOri)
            print 'reducing chain',oriChain, '({0})'.format(newChain), minRes, ':', maxRes
            oriReduced =[] + [a for a in oriAtoms if a.chain==oriChain and a.resId>minRes and a.resId<maxRes  and a.symbol == 'CA']
            newReduced =[] + [a for a in newAtoms if a.chain==newChain and a.resId>minRes and a.resId<maxRes  and a.symbol == 'CA']
            if len(newReduced)==0:
                global test
                test=newAtoms
                print len(newAtoms)
                raw_input("Error?")
            yield oriReduced, newReduced
        
    toOrigA=dict(zip(chains.split(':')[0], chainA))
    toOrigB=dict(zip(chains.split(':')[1], chainB))
    atomsA = []
    for chain in chains.split(':')[0]:
        atomsA += [a for a in pdb.atoms if a.chain == chain]
    atomsB = []
    for chain in chains.split(':')[1]:
        atomsB += [a for a in pdb.atoms if a.chain == chain]
    unboundA, unboundB = [], []
    for chain in chainA:
        unboundA += [a for a in pdbA.atoms if a.chain == chain]
    for chain in chainB:
        unboundB += [a for a in pdbB.atoms if a.chain == chain]

    hasMutatin = False

    
    atomsARed, atomsBRed =[], []
    atomsUnARed, atomsUnBRed =[], []
    missing=False
    for redBound, redUnbound in reduceBackbonesChains(toOrigA, atomsA, unboundA):
        if len(redBound)!=len(redUnbound):
            resRedBound=dict([(a.resId, a.residue) for a in redBound])
            resRedUnbound=dict([(a.resId, a.residue) for a in redUnbound])
            mutations=[str(k)+':'+ v+'-'+resRedUnbound[k] for k, v in resRedBound.iteritems() if v!=resRedUnbound[k]]
            if len(mutations)>0:
                print 'mutant',string.join(mutations, '\n')
            else:
                print 'some atom is missing. chain A (after reduce): ', pdb.name,'    ', len(redBound), ':', len(redUnbound)
            raw_input("ok?")
            missing=True
        else:
           atomsARed+=redBound
           atomsUnARed+=redUnbound
    
    for redBound, redUnbound in reduceBackbonesChains(toOrigB, atomsB, unboundB):
        if len(redBound)!=len(redUnbound):
            resRedBound=dict([(a.resId, a.residue) for a in redBound])
            resRedUnbound=dict([(a.resId, a.residue) for a in redUnbound])
            mutations=[str(k)+':'+ v+'-'+resRedUnbound[k] for k, v in resRedBound.iteritems() if v!=resRedUnbound[k]]
            if len(mutations)>0:
                print 'mutant',string.join(mutations, '\n')
            else:
                print 'some atom is missing. chain B (after reduce): ', pdb.name,'    ', len(redBound), ':', len(redUnbound)
            raw_input("ok?")
            missing=True
        else:
           atomsBRed+=redBound
           atomsUnBRed+=redUnbound
    if not missing:
        print 'find diff for reduced ', pdb.name
        findDiffForSameLength(atomsARed, atomsBRed, atomsUnARed, atomsUnBRed)
    
 
"""
