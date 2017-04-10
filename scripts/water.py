import atom

class water(atom.atom):
    def __init__(self, atomNum, atomSymbol, residue, chain,resId, atomType, coord):
        super(water, self).__init__(atomNum, atomSymbol, residue, chain,resId, atomType, coord)
