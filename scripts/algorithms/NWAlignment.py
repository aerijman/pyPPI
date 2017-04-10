import os
from DynamicProgramming import DynamicProgramming

class NWAlignment(DynamicProgramming):
    """ implementation of Needleman-Wunsch"""
    
    space = -2 #TODO: define score for space
    MATCH = 1
    MISS_MATCH = -1
    def __init__(self, seq1, seq2):
        self.seq1 = seq1
        self.seq2 = seq2
        super(NWAlignment, self).__init__(len(seq2)+1, len(seq1)+1)   
    
    def getInitialPointer(self, cell):
        if cell.row == 0 and cell.col != 0:
            return self.matrix[cell.row][cell.col-1]
        elif cell.row != 0 and cell.col == 0:
            return self.matrix[cell.row-1][cell.col]
        else:
            return None

    def getInitialScore(self, cell):
        if cell.row == 0 and cell.col != 0:
            return cell.col*NWAlignment.space
        elif cell.row != 0 and cell.col == 0:
            return cell.row*NWAlignment.space
        else:
            return 0

    def fillInCell(self, currentCell, cellAbove, cellToLeft,  cellAboveLeft):
        rowSpaceScore = cellAbove.score+NWAlignment.space
        colSpaceScore = cellToLeft.score+NWAlignment.space
        matchMissmatchScore = cellAboveLeft.score
        #give match/missmatch score
        matchMissmatchScore +=BLOSUM62.getValues(self.seq2[currentCell.row-1], self.seq1[currentCell.col-1])
        
        if rowSpaceScore >= colSpaceScore:
            if matchMissmatchScore >= rowSpaceScore:
                currentCell.score = matchMissmatchScore
                currentCell.prevCell = cellAboveLeft
            else:
                currentCell.score = rowSpaceScore
                currentCell.prevCell = cellAbove
        else:
            if matchMissmatchScore >= colSpaceScore:
                currentCell.score = matchMissmatchScore
                currentCell.prevCell = cellAboveLeft
            else:
                currentCell.score = colSpaceScore
                currentCell.prevCell = cellToLeft        

    def getTraceback(self):
        align1 = []
        align2 = []
        cell = self.getTracebackStartingCell()
        while not self.traceBackIsDone(cell):
            prevCell = cell.prevCell
            if cell.row-prevCell.row == 1:
                align2.insert(0, self.seq2[cell.row-1])
            else:
                align2.insert(0, '-')

            if cell.col-prevCell.col == 1:
                align1.insert(0, self.seq1[cell.col-1])
            else:
                align1.insert(0, '-')
            cell = prevCell
        return align1, align2

    def getTracebackStartingCell(self):
        matrix = self.matrix
        return matrix[len(matrix) - 1][len(matrix[0]) - 1];
    
    def traceBackIsDone(self, cell):
        return cell.prevCell == None

class BLOSUM62:
    __instance = None

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
        'VAL': 'V'
    }

    def __init__(self):
        f = file('BLOSUM62.csv')
        headers = f.readline().strip().split(',')[1:]
        self.subValues = dict()
        for l in f.readlines():
            values = l.split(',')
            resDest = values[0]
            for resSrc, v in zip(headers, values[1:]):
                self.subValues[(resSrc,resDest)]=int(v)
        f.close()

    @staticmethod
    def getValues(resA, resB):
        return BLOSUM62.instance().subValues[(resA, resB)]
    
    @staticmethod
    def instance():
        if not BLOSUM62.__instance:
            BLOSUM62.__instance = BLOSUM62()
        return BLOSUM62.__instance
            
def main():
    import string
    seq1 = 'VAVAAA'
    seq2 = 'VWAAL'
    align = NWAlignment(seq1,seq2)
    align.fillIn()

    nSeq1, nSeq2 = align.getTraceback()
    print nSeq1
    print nSeq2
    
    
if __name__ == '__main__':
    main()
