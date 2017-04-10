"""Sequence alignment"""

class Cell(object):
    def __init__(self, r, c):
        self.prevCell = None
        self.score = 0
        self.row = r
        self.col = c

    def setScore(self, score):
        self.score = score

    def setPrevCell(self, prev):
        self.prevCell = prev

    def getScore(self):
        return self.score

class DynamicProgramming(object):
    def __init__(self, rows, cols):
        self.matrix=[]
        for r in range(0, rows):
            row = list([Cell(r, c) for c in range(0, cols)])
            self.matrix.append(row)

        self.initializeScores();
        self.initializePointers();

    def initializeScores(self):
        for row in self.matrix:
            for col in row:
                col.setScore( self.getInitialScore(col) )

    def initializePointers(self):
        for row in self.matrix:
            for col in row:
                col.setPrevCell( self.getInitialPointer(col) )

    def getInitialScore(self, cell):
        print 'TODO: implement getInitialScore'

    def getInitialPointer(self, cell):
        print 'TODO: implement getInitialPointer'

    def fillIn(self):
        for row in self.matrix[1:]:
            for cell in row[1:]:
             cellAbove = self.matrix[cell.row-1][cell.col]
             cellToLeft =  self.matrix[cell.row][cell.col-1]
             cellAboveLeft =  self.matrix[cell.row-1][cell.col-1]
             self.fillInCell(cell, cellAbove, cellToLeft, cellAboveLeft)


    def fillInCell(self, currentCell, cellAbove, cellToLeft,  cellAboveLeft):
        print 'TODO: implement fillInCell'

    def getTraceback(self):
        print 'TODO: implement getTraceback'


