set @engCutoff:=-0.6;
set @unsetisfiedMax:=-0.5;/* mb to lower to -0.5 */
set @unsatifiedHASA:=5;/* dont include exposed hbonds that may form HB with water */
select 
	concat(orig.PDB,'_',UnboundChainA,':',UnboundChainB), orig.NameA, orig.NameB,
	count(distinct(if(((LOCATE(hb.DonorChain,UnboundChainA)>0 and 
	LOCATE(hb.AccChain,UnboundChainA)=0) OR (LOCATE(hb.AccChain,UnboundChainA)>0 and 
	LOCATE(hb.DonorChain,UnboundChainB)=0)),hb.fullDrieding_id,0)))-1 as inter,
	count(distinct(if(((LOCATE(hb.DonorChain,UnboundChainA)>0 and 
	LOCATE(hb.AccChain,UnboundChainA)>0)
	OR
	(LOCATE(hb.AccChain,UnboundChainB)>0 and 
	LOCATE(hb.DonorChain,UnboundChainB)>0)),hb.fullDrieding_id,0)))-1 as intra,
	interfaceVDW.VDV as VDW,
	fASA.diffASA
from 
	proteinComplex orig
inner join
	NinterfaceAtoms atms
on
	atms.PDB=left(orig.PDB,4)
inner join
     interfacePeriphrial peri
on
     peri.PDB=atms.PDB and
     peri.Chain=atms.Chain and
     peri.Resid=atms.Resid and
     peri.Symbol=atms.Symbol

left join
	Ndrieding hb
on 
	hb.Energy<@engCutoff and
	atms.PDB=hb.PDB and
	((atms.Chain=hb.DonorChain and
	atms.ResId=hb.DonorResId and
	atms.Symbol=hb.DonorSymbol) or 
	(atms.Chain=hb.AccChain and
	atms.ResId=hb.AccResId and
	atms.Symbol=hb.AccSymbol))
inner join
	(select orig.PDB as Comp,sum(asa.diffASA) as diffASA
	from proteinComplex orig
	inner join NinterfaceAtoms asa
	on asa.PDB=left(orig.PDB,4)
group by orig.PDB
) fASA
on fASA.Comp=orig.PDB
left join
	interfaceVDW
on
	atms.PDB=interfaceVDW.PDB
where peri.Peri>=3
group by orig.PDB
