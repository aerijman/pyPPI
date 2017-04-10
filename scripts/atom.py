"""
Class representing atom
"""
ResiduesCodes = {
    'ALA': 'A',
    'ARG': 'R',
    'ASN': 'N',
    'ASP': 'D',
    'CYS': 'C',
    'GLU': 'E',
    'GLN': 'Q',
    'GLY': 'G',
    'HIS': 'H',
    'ILE': 'I',
    'LEU': 'L',
    'LYS': 'K',
    'MET': 'M',
    'PHE': 'F',
    'PRO': 'P',
    'SER': 'S',
    'THR': 'T',
    'TRP': 'W',
    'TYR': 'Y',
    'VAL': 'V',
    'GLX': 'Q',#FAKE
    'ASX': 'N',#FAKE
    'UNK': 'T'#SUPER-FAKE
}
class atom(object):
    
    
    #def __init__(self, atomNum, atomSymbol, residue, chain,resId, x, y, z,atomType, atomIndex=0 ):
    
    def __init__(self, atomNum, atomSymbol, residue, chain,resId,atomType, coord, atomIndex=0, tempFactor=0 ):
        self.neighbors = []
        self.symbol = atomSymbol
        self.residue = residue
        self.chain = chain
        self.resId = int(resId)
        self.x = coord[0]
        self.y = coord[1]
        self.z = coord[2]
        self.coord = coord
        self.interface = []
        self.atomType=atomType
        self.atomNum=atomNum
        #used by pdb:
        self.atomIndex=atomIndex
        self.tempFactor=tempFactor
        
        self.pseudoChain=None#instead of chain

    def distance(self, atom2):
        return (self.x-atom2.x)**2 +(self.y-atom2.y)**2 +(self.z-atom2.z)**2

    def distanceFromXYZ(self,XYZ):
        return (self.x-XYZ[0])**2 +(self.y-XYZ[1])**2 +(self.z-XYZ[2])**2

    def isInterface(self):
        return len(self.interface)>0

    def resCode(self):
        global ResiduesCodes
        return ResiduesCodes[self.residue]

    def __str__(self):
        return  self.chain+str(self.resId)+"|"+self.atomNum
