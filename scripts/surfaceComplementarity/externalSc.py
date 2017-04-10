import subprocess, re, MySQLdb, string
def useExternalSC(pdbName, chains):
    print 'Calcing Sc for '+pdbName
    pdbPath='/home/eran/Eran/pdbs/%s_FH.pdb'%pdbName
    res=subprocess.check_output('./sc/sc %s -1:%s -2:%s'%(pdbPath,chains[0],chains[1]),shell=True)
    return float(res.split('\t')[1])

def getNamesAndChainFromDB():
     conn=MySQLdb.connect(host='localhost', user='root', passwd='Admin', db='ppi')
     cursor=conn.cursor()
     
     #create interface table
     cursor.execute("""select Complex, niceLog from origTable2""")
     res=[]
     for c,lg in cursor.fetchall():
          res.append((c[0:4], c[5:].split(':'),lg))
     return res

def main():
    result=file('scResults.csv','w')
    print>>result,string.join(['complex','sc','log'],',')
    for pdb, chains,niceLog in getNamesAndChainFromDB():
        sc=useExternalSC(pdb,chains)
        print>>result,string.join([pdb,str(sc),str(niceLog)],',')
        result.flush()
    result.close()
    print 'Finished'


if __name__ == "__main__":
    main()
    

