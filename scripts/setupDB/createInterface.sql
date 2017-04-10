/* Insert to interface table all atoms that have diffASA>0*/
insert into NinterfaceAtoms (PDB,Chain,Residue,ResId,Symbol,atom,diffASA)
select PDB,Chain,Residue,ResId,Symbol,Atom,max(ASA)-min(ASA) from perAtomASA
group by PDB,Chain,ResId,Symbol
having stddev(ASA)>0;
/* Insert to interface table all atoms that have enough distance */
insert ignore into NinterfaceAtoms (PDB,Chain,Residue,ResId,Symbol,atom)
select asa.PDB,asa.Chain,asa.Residue,asa.ResId,asa.Symbol,dist.Atom from interfaceDist dist
inner join
	perAtomASA asa
on
	dist.PDB=asa.PDB and
	dist.Chain=asa.Chain and
	dist.ResId=asa.ResId and
	dist.Symbol=asa.Symbol and
	Seperated=0
