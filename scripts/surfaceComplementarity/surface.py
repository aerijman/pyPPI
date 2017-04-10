import sys, os, string,math, MySQLdb, logging, datetime
import scipy

sys.path.append( '../' )
sys.path.append( '../algorithms' )
from pdbReader import PDBReader
from ASA import radio_atom, spiral
from ASA import KNOWN_RADIUS, R_WATER
import time
from kdtree2 import KDTree

"""
PDB coordinates can be negative
"""
#1j = sqrt -1

MAX_WW_RAD = max(KNOWN_RADIUS.values())
MOL_SURFACE = 1
OUTSIDE_MOL = 0
RO_INSIDE_A = -1#large negative - their 15
DELTA_INSIDE_B = 1.1#small positive
STEP_SIZE=0.3#they used 0.7-0.8

def getInterfaceAtoms(pdb):
     conn=MySQLdb.connect(host='localhost', user='root', passwd='Admin', db='ppi')
     cursor=conn.cursor()
     
     #create interface table
     cursor.execute("""select Chain,ResId,Symbol from
                        interfaceAtoms
                        where PDB='%s'"""%( pdb.name ))
     #find unbounds relating to
     interfaceAtoms=[]
     for chain,resId,symbol in cursor.fetchall():
          interfaceAtoms+=[a for a in pdb.atoms if a.chain==chain and a.resId==resId and a.symbol[0:3]==symbol]
     return interfaceAtoms

def getNamesAndChainFromDB():
     conn=MySQLdb.connect(host='localhost', user='root', passwd='Admin', db='ppi')
     cursor=conn.cursor()
     
     #create interface table
     cursor.execute("""select Complex, niceLog from origTable2""")
     res=[]
     for c,lg in cursor.fetchall():
          res.append((c[0:4], c[5:].split(':'),lg))
     return res

def getGridBound(atoms,stepSize):
    #todo:change to interfce bound
    xs=min(a.x for a in atoms)-MAX_WW_RAD, max(a.x for a in atoms)+MAX_WW_RAD
    xs=xs[0]-xs[0]%stepSize,xs[1]-xs[1]%stepSize
    ys=min(a.y for a in atoms)-MAX_WW_RAD, max(a.y for a in atoms)+MAX_WW_RAD
    ys=ys[0]-ys[0]%stepSize,ys[1]-ys[1]%stepSize
    zs=min(a.z for a in atoms)-MAX_WW_RAD, max(a.z for a in atoms)+MAX_WW_RAD
    zs=zs[0]-zs[0]%stepSize,zs[1]-zs[1]%stepSize
    return xs,ys,zs


def getSampleGrid(pdb, compChains):
    samplingPoint=1000
    compA=[a for a in pdb.atoms if a.chain in compChains]
    compABound=getGridBound(compA)
    for j in range(0,samplingPoint):
        randPoint=[random.uniform(compABound[d][0],compABound[d][1]) for d in range(0,3)]
        atomsAround=[a for a in pdb.ktree.findByDistance(query_point=randPoint, distance=MAX_WW_RAD) if atom.chain in compChains]
        if len(atomsAround)==0:
            v=OUTSIDE_MOL
        elif len(atomsAround)==1:
            #thin layer around described by only 1 atom touching it
            #we may think of making it broder (adding more point around)
            v=RO_INSIDE_A
        else:
            v=MOL_SURFACE
        yield (randPoint,v)

def frange(x, y, jump):
  #float point suprise...
  while x < y:
    yield x
    x += jump

def getGrid(pdb,compA,gridSize,insideVal):
     global STEP_SIZE
     
     SURFACE_WIDTH=STEP_SIZE/2.0
     MAX_WW_RAD = max(KNOWN_RADIUS.values())     
     compABound=getGridBound(compA,STEP_SIZE)#only on interface instead of compA
     
     #grid=scipy.zeros([gridSize[i] for i in range(0,3)])
     grid=scipy.zeros([(compABound[i][1]-compABound[i][0])/STEP_SIZE for i in range(0,3)])
     #compABound=zip([0,0,0],[i/STEP_SIZE for i in grid.shape])
     logging.debug('getting grid of %s, insideval %s'%(grid.shape,insideVal))
     ktree=KDTree.construct_from_data(pdb.atoms[:])
     #ktree=pdb.ktree
     #i= time.clock()
     ti=datetime.datetime.now()
     numInside=0
     interAChain=set(a.chain for a in compA)
     """
     Method a
     """
     
     logging.debug('grid coords:'+string.join(['(%s, %s)'%(compABound[i][0],compABound[i][1]) for i in range(0,3)],' | '))
     for gx,x in enumerate(frange(compABound[0][0],compABound[0][1], STEP_SIZE)):
          if gx>=grid.shape[0]:
               continue
          for gy,y in enumerate(frange(compABound[1][0],compABound[1][1], STEP_SIZE)):
               if gy>=grid.shape[1]:
                    continue
               for gz,z in enumerate(frange(compABound[2][0],compABound[2][1], STEP_SIZE)):
                    if gz>=grid.shape[2]:
                         continue
                    #for atom in pdb.ktree.findByDistance(query_point=(x,y,z), distance=MAX_WW_RAD**2):
                    for atom in ktree.findByDistance(query_point=(x,y,z), distance=MAX_WW_RAD**2):
                         if grid[gx,gy,gz]==insideVal:#dont override inside
                              break
                         if atom.chain in interAChain:
                              atomRad2=radio_atom(atom.atomType)
                              distanceXYZ=math.sqrt(atom.distanceFromXYZ((x,y,z)))
                              if distanceXYZ<=atomRad2+SURFACE_WIDTH:
                                   grid[gx,gy,gz]=insideVal if distanceXYZ<atomRad2 else 1
                                   if distanceXYZ<atomRad2:
                                        numInside+=1
                                        
                    
     #logging.debug('building grid took %s. inside grid %s'%(time.clock()-i,numInside))
     print 'building grid took %s. inside grid %s'%((datetime.datetime.now()-ti).total_seconds(),numInside)
     return grid
     """Method b
     
     grid=scipy.zeros([gridSize[i] for i in range(0,3)])
     i= time.clock()
     numInside=0
     for atom in [a for a in pdb.atoms if a.chain in interAChain]:
          atomRad=radio_atom(atom.atomType)
          atomRad+=R_WATER*2
          coordInRange=[i for i in range(0,3) if atom.coord[i]-atomRad<compABound[i][1] and atom.coord[i]+atomRad>compABound[i][0]]
          if len(coordInRange)!=3:
               continue
          atomRad2=atomRad**2
          ranges=[(max(compABound[i][0],atom.coord[i]-atomRad),min(compABound[i][1],atom.coord[i]+atomRad)) for i in range(0,3)]
          ranges=[(r[0]-r[0]%STEP_SIZE, r[1]-r[1]%STEP_SIZE) for r in ranges]
          for gx,x in enumerate(frange(ranges[0][0],ranges[0][1],STEP_SIZE)):
               for gy,y in enumerate(frange(ranges[1][0],ranges[1][1],STEP_SIZE)):
                    for gz,z in enumerate(frange(ranges[2][0],ranges[2][1],STEP_SIZE)):
                         point=(x,y,z)
                         #gridP=[(point[i]-compABound[i][0])/STEP_SIZE for i in range(0,3)]
                         if grid[gx,gy,gz]!=insideVal:
                              distanceXYZ=math.sqrt(atom.distanceFromXYZ((x,y,z)))
                              if distanceXYZ<=atomRad+SURFACE_WIDTH:
                                   grid[gx,gy,gz]=insideVal if distanceXYZ<atomRad2 else 1
                                   if distanceXYZ<atomRad2:
                                        numInside+=1

     print 'building grid took', (time.clock()-i), " inside grid ",numInside
     raise Exception("check")
     return grid
     """


def getShift(pdb,comp):
    """
    Calculates shift between proteins.
    Is needed for complexs: i think no
    Now it uses averagem maybe should use minimum
    """
    global STEP_SIZE
    interA=comp[0]
    interB=comp[1]
    interAChain=set(a.chain for a in interA)
    distances=[[],[],[]]
    for iA in interB:
        atomRad=radio_atom(iA.atomType)
        distanceAround=(atomRad+R_WATER*2+MAX_WW_RAD)**2
        atomsAround=[a for a in pdb.ktree.findByDistance(query_point=iA.coord, distance=distanceAround) if a.chain in interAChain]
        if len(atomsAround)>0:
             distances[0]+=[iA.x-a.x for a in atomsAround]
             distances[1]+=[iA.y-a.y for a in atomsAround]
             distances[2]+=[iA.z-a.z for a in atomsAround]
    #shift definition:
    #article - number of grid steps by which molecule B is shift with respect to A
    #here - AVG or minimum of distance of surface coords

    #it may be not optimal for long proteins
    #shift=[sum(d/abs(d)*math.sqrt(abs(d)) if d!=0 else 0 for d in distances[i])/len(distances[i]) for i in range(0,3)]
    #shift=[sum(math.sqrt(d))/len(distances[i]) for i in range(0,3)  for d in distances[i]]
    #minimum
    shift=[math.sqrt(min(abs(d) for d in distances[i] if d!=0)) for i in range(0,3) ]
    shift=[int((shift[i]-shift[i]%STEP_SIZE)/STEP_SIZE) for i in range(0,3) ]
    return shift
    

global GRID1, GRID2

def showPlot(g1,g2,shift):
     import matplotlib.pyplot as plt
     from mpl_toolkits.mplot3d import Axes3D
     a,b,y=shift
     fig = plt.figure()
     ax = fig.add_subplot(111, projection='3d')
     #g1
     xs,ys,zs=[],[],[]
     for l in range(0,g1.shape[0]):
          for m in range(0,g1.shape[1]):
               for n in range(0,g1.shape[2]):
                    if g1[l,m,n]!=1:
                         xs.append(l)
                         ys.append(m)
                         zs.append(n)
     ax.scatter(xs, ys, zs, c='r')
     #g2
     xs,ys,zs=[],[],[]
     for l in range(0,g2.shape[0])[::2]:
          for m in range(0,g2.shape[1])[::2]:
               for n in range(0,g2.shape[2])[::2]:
                    if g2.shape[0]>l-a and g2.shape[1]>m-b and g2.shape[1]>n-y and g2[l-a,m-b,n-y]!=1:
                         xs.append(l-a)
                         ys.append(m-b)
                         zs.append(n-y)
     ax.scatter(xs, ys, zs, c='b')
     #ax.
     
     ax.set_xlabel('X Label')
     ax.set_ylabel('Y Label')
     ax.set_zlabel('Z Label')

     plt.show()


def calcCompl(pdbName, chains):
     global GRID1, GRID2
     global RO_INSIDE_A, DELTA_INSIDE_B, STEP_SIZE
     pdb = PDBReader.readFile('/home/eran/Eran/pdbs/%s_FH.pdb'%pdbName,chains)
     
     interface = getInterfaceAtoms(pdb)
     components=[]
     for part in pdb.interfaceParts:
          components.append([a for a in interface if a.chain in part])
     #A is the larger component and b is the smaller
     components.sort(key=lambda y:-len(y))
     bounds=[]
     normalBound=[1,1,1]
     for comp in components:
          b=getGridBound(comp,STEP_SIZE)
          #print b
          bounds.append(b)
          for i in range(0,3):
               normalBound[i]=max((b[i][1]-b[i][0])*(1/STEP_SIZE),normalBound[i])
     normalBound=[max(normalBound[i] for i in range(0,3))]*3
     logging.debug('using grid of size %s', normalBound)

     #shift=getShift(pdb,components)
     shift=[int((bounds[1][i][0]-bounds[0][i][0])*(1.0/STEP_SIZE)) for i in range(0,3)]
     #print [((bounds[1][i][0]-bounds[0][i][0])) for i in range(0,3)]
     g1=getGrid(pdb,components[0],normalBound,RO_INSIDE_A)
     g2=getGrid(pdb,components[1],normalBound,DELTA_INSIDE_B)
     
     a,b,y=shift
     GRID1, GRID2=g1,g2
     c=0.0
     """
     for a in frange(-g2.shape[0],g2.shape[0],STEP_SIZE):
          for b in frange(-g2.shape[1],g2.shape[1],STEP_SIZE):
               for y in frange(-g2.shape[2],g2.shape[2],STEP_SIZE):
                    shift=a,b,y
                    print 'using shift %s,%s,%s'%(a,b,y)
                    g1F=scipy.fft(g1)
                    #g2Translate=scipy.zeros([normalBound[i] for i in range(0,3)])
                    #rngs=[range(shift[i], g2.shape[i]+shift[i]) for i in range(0,3)]
                    #for l in range(0,g2.shape[0]):
                    #     for m in range(0,g2.shape[1]):
                    #          for n in range(0,g2.shape[2]):
                    #               if g2Translate.shpae[0]>l+a and g2Translate.shpae[1]>m+b and g2Translate.shpae[2]>n+y:
                    #                    g2Translate[l+a,m+b,n+y]=g2[l,m,n]
                    g2F=scipy.fft(g2)#g2Translate
                    c=scipy.ifft(g1F*g2F)
                    for k in c:
                         print 'k corr: ',k
     """
     i=time.clock()
     numOfClashes=0
     contact=0
     #clashCoords=[]
     for l in range(0,g1.shape[0]):
          for m in range(0,g1.shape[1]):
               for n in range(0,g1.shape[2]):
                    if g2.shape[0]>l-a and g2.shape[1]>m-b and g2.shape[2]>n-y and l-a>0 and m-b>0 and n-y>0:
                         if g1[l,m,n]*g2[l-a,m-b,n-y]<0:
                              numOfClashes+=1
                              #clashCoords.append((bounds[0][0][0]+l*STEP_SIZE,bounds[0][1][0]+m*STEP_SIZE,bounds[0][2][0]+n*STEP_SIZE))
                         elif g1[l,m,n]*g2[l-a,m-b,n-y]>0:
                              contact+=1
                         c+=g1[l,m,n]*g2[l-a,m-b,n-y]

     logging.debug('calcing corr took %s'%(time.clock()-i))
     #print 'Clash coords'
     #for cla in clashCoords:
     #     print cla, '  :', string.join([str(a.chain)+''+str(a.resId) for a in pdb.ktree.findByDistance(query_point=cla, distance=3)])
     print 'calcing corr took ',(time.clock()-i)
     #print 'clashes: ',numOfClashes
     #print 'using shift %s,%s,%s'%(a,b,y)
     print 'Corr: %s   | Shift: %s      Clashes: %surface       Conactas: %s'%(c,(a,b,y), numOfClashes,contact)
     #showPlot(g1,g2,shift)
     return c,numOfClashes,contact
     
     """
     i=time.clock()
     g1F=scipy.fft(g1,int(normalBound[0]))
     g2Translate=scipy.zeros([g1.shape[i] for i in range(0,3)])
     #rngs=[range(shift[i], g2.shape[i]+shift[i]) for i in range(0,3)]
     for l in range(0,g2.shape[0]):
          for m in range(0,g2.shape[1]):
               for n in range(0,g2.shape[2]):
                    if g2Translate.shape[0]>l+a and g2Translate.shape[1]>m+b and g2Translate.shape[2]>n+y:
                         g2Translate[l+a,m+b,n+y]=g2[l,m,n]
     g2F=scipy.fft(g2,int(normalBound[0]))
     c=sum(scipy.ifft(g1F*g2F))
     
     print 'calcing corr took', (time.clock()-i)
     
     print 'using shift %s,%s,%s'%(a,b,y)
     print 'Corr: ',c
     """

tRES=dict()

def check():
     #c1=calcCompl('1A2K',['C','AB'])
     #calcCompl('1EZU',['C','AB'])
     #calcCompl('1AKJ',['DE','AB'])
     calcCompl('1DFJ',['E','I'])

def main():
     global tRES
     
     #c1=calcCompl('1A2K',['C','AB'])#medium
     #c2=calcCompl('1S1Q',['A','B'])#low
     #c3=calcCompl('1DFJ',['E','I'])#high
     res=dict()
     result=file('complementarity.csv','w')
     print>>result,string.join(['PDB,complementarity,clash,contact,niceLog'],',')
     for pdb, chains,niceLog in getNamesAndChainFromDB():
          logging.info('Checking %s'%pdb)
          corr,clash,contact=calcCompl(pdb,chains)
          res[pdb]=corr
          #print>>result,string.join([pdb,str(res[pdb]),str(niceLog)],',')
          print>>result,string.join([pdb,str(corr),str(clash),str(contact),str(niceLog)],',')
          result.flush()

     tRES=res
     result.close()
     print 'finished'


if __name__ == "__main__":
    logging.basicConfig(filename=None,format='%(levelname)s: %(message)s', level=logging.INFO)
    main()
    #check()
