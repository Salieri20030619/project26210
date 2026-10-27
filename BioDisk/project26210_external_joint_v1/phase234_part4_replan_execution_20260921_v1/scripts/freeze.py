from common import *
import h5py, scipy, traceback, shutil, math
from anndata.io import read_elem
from scipy.sparse.csgraph import dijkstra
from scipy.spatial import cKDTree
def main():
 save('STATUS.json',dict(status='A0_FREEZING',time=stamp()))
 old=json.loads((F/'audit/RDS_BEFORE.json').read_text());save('audit/RDS_BEFORE.json',old)
 plan=R/'phase233_part4_complete_replan_20260921_v1/PART4_COMPLETE_ANALYSIS_PLAN_CN.md'
 sources=[meta(plan)]
 defs=json.loads((R/'phase218_spatial_identity_alignment_20260920_v1/PROGRAM_DEFINITIONS.json').read_text())
 author=json.loads((F/'PROGRAM_DEFINITIONS.json').read_text())['programs']
 programs={k:defs[k] for k in ['spp1_state_no_spp1','c1qc_state','polyamine_mta_context','cytotoxicity','tcell_dysfunction','immunosuppression_core','invasion_emt','vascular_remodeling']}
 programs.update(author)
 programs['spp1_no_LGALS3']=[g for g in programs['spp1_state_no_spp1'] if g!='LGALS3']
 programs['vascular_no_LGALS3']=[g for g in programs['vascular_remodeling'] if g!='LGALS3']
 custompath=Path('/mnt/BioDisk/ludingyu/project26210/results/part4_spatial/customer_style/run_metadata.json')
 custom=json.loads(custompath.read_text())['gene_sets']
 names=['Pyrimidine metabolism','Cysteine and methionine metabolism','Arginine and proline metabolism','Glycerophospholipid metabolism']
 short=['pyrimidine','cysteine_methionine','arginine_proline','glycerophospholipid']
 gmt=R/'phase199_fig3_marker_GSEA_20260918_v5/inputs/KEGG_2021_Human.gmt'
 kg={l.split('\t')[0]:l.strip().split('\t')[2:] for l in gmt.read_text().splitlines()}
 for name,k in zip(names,short):
  matches=[v for key,v in kg.items() if key.lower().startswith(name.lower())];assert len(matches)==1,(name,len(matches))
  programs['KEGG_'+k]=sorted(set(matches[0])-{'TYMP'})
  programs['custom_'+k]=sorted(set(custom[name])-{'TYMP'})
 for path in [custompath,gmt,R/'phase218_spatial_identity_alignment_20260920_v1/PROGRAM_DEFINITIONS.json',F/'PROGRAM_DEFINITIONS.json']:
  sources.append(meta(path));sources[-1]['sha256']=sha(path)
 lrpath=R/'phase213_immune_restart_execution_20260919_v1/I4_directed_communication/FROZEN_RECEPTOR_CANDIDATES.tsv'
 lr=pd.read_csv(lrpath,sep='\t');lr['family']='SPP1_LR';lr['database']='CellChatDB v1; frozen Phase213 extraction'
 dbpath=R/'phase35_part1_c01_c04_cell_communication_v1/resource_lock_v1/cellchat_interactions_v1.tsv'
 db=pd.read_csv(dbpath,sep='\t');postn=db[db.ligand=='POSTN']
 for rr in postn.itertuples():
  lr.loc[len(lr)]=[rr.interaction_name,rr.ligand,rr.receptor,rr.receptor.replace('_',';'),rr.evidence,'POSTN_LR_comparison',str(rr.version)]
 tab('tables/FROZEN_LR_CANDIDATES.tsv',lr)
 for path in [lrpath,dbpath]:sources.append(dict(**meta(path),sha256=sha(path)))
 save('PROGRAM_DEFINITIONS.json',dict(programs=programs,score='phase234 equal mean within-section population gene Z of log1p(counts/total*10000); not GSVA',eligibility='at least half listed genes present and at least 3 present (2 for two-gene set), at least20 spots detect >=min(2,n_present); nonzero variance. Exact present genes retained; no association-based selection',all_metabolism_exclude_TYMP=True))
 ep=[]
 def add(fam,x,y,role='primary'):
  ep.append(dict(family=fam,endpoint=x+'__'+y,exposure=x,target=y,role=role))
 for k in short:
  add('KEGG_four_axes','TYMP_RNA','KEGG_'+k)
  add('custom_four_axes','TYMP_RNA','custom_'+k)
 vlist=['V1_SPP1_q05','V2_SPP1_fraction','V3_SPP1_macro_share','V4_SPP1_program']
 for v in vlist:
  add('MTA_context_SPP1','polyamine_mta_context',v)
  for t in ['endothelial_fraction','pericyte_fraction']:add('SPP1_vascular',v,t)
  for t in ['POSTNCAF','ECM4']:add('SPP1_stroma_program',v,t)
  for t in ['cytotoxicity','tcell_dysfunction','immunosuppression_core','invasion_emt']:add('SPP1_function',v,t)
 for x in ['total_macro_fraction','C1QC_fraction','APOEMph']:
  for t in ['endothelial_fraction','pericyte_fraction','POSTNCAF','ECM4']:add('comparison_partners',x,t,'comparison')
 add('shared_gene_program','V4_SPP1_program','vascular_remodeling')
 add('shared_gene_removed','spp1_no_LGALS3','vascular_no_LGALS3','sensitivity')
 tab('tables/FROZEN_ENDPOINTS.tsv',pd.DataFrame(ep))
 rules=dict(frozen_time=stamp(),plan_sha256=sha(plan),known_prior_results=True,blind_preregistration=False,
 broad_scales=['W0','W1_4x4','W1_5x5'],broad='all QC finite; >=50 spots; block medians on original coordinate rank(method=first) qcut; >=8 occupied blocks; no new spot/block P',
 radius_um=[800,400,200],radius_distance='weighted shortest path on frozen first-ring graph; excludes self; edge length pixels*55/spot_diameter_fullres; median edge must 90..110um; no cross-component edges',
 source_gate='Q0 finite; Q1 macro fraction>=0.05; Q2>=0.10; at least3 neighbors and50 sources',
 scales_core='TYMP->8metabolic; MTA->V1/V3; V1/V3->8vascular/stroma/function targets; Q1 A0',
 qualification='400um V3->8 targets Q0 Q1 Q2; own masks and Q2 common-source mask',
 adjustment='400um Q1 V1/V3->8 targets; A0 A1(depth source/neighbor) A2(+macro fraction source/neighbor and malignant fraction) A3(+source target). V1 uses other-macro fraction excluding SPP1 to avoid source as covariate; V3 uses totalmacro; A4 coordinates for two stroma targets only',
 boundary='frozen malignant fraction>=0.20; signed nearest graph distance to expression interface. tumor bands<=800/400/200 vs tumor core>800; exterior bands<=800/400/200 vs exterior>800; >=20 per region. Continuous200um bins and slope when >=3 bins with20 spots. No histology claims',
 LR='9 SPP1 and all frozen POSTN DB interactions; log-normalized ligand vs neighborhood min-subunit receptor expression, requiring all subunits present and jointly detected>=20 spots. At least20 ligand positive; 800/400/200/ring, Q0 A0; compatibility not cell-type communication',
 patient_summary='Fisher-z equal preparation mean for rho, arithmetic for differences/slopes; source-separated; n>=10 exact two-sided sign P + Wilcoxon companion,10000 patient bootstrap median CI; smaller n descriptive',
 FDR='dataset x family x layer; parallel variables and all scales/models in same layer included; unestimated registered rows P=1 only inside correction, displayed P/FDR NA',
 physical_reference='https://cdn.10xgenomics.com/image/upload/v1758134784/support-documents/CG000686_VisiumHDSpatialApplicationsPerformance_TechnicalNote_RevB.pdf',
 new_figures=0)
 save('IMPLEMENTATION_FREEZE_BEFORE_ASSOCIATION.json',rules)
 records=json.loads((A/'ASSET_ROLE_MANIFEST.json').read_text())
 pat=pd.read_csv(S/'FROZEN_SECTION_PATIENT_MAP.tsv',sep='\t').fillna('').set_index('section')
 for sid,v in pat.iterrows():
  records.append(dict(sample=sid,dataset='GSE175540',patient=v.patient_id,input_path=str(R.parent/'project26210_external_joint_v1_phase54_gse175540_spatial_v1/run_v1/sections_v1'/sid/'processed_section_v1.h5ad'),abundance_path=str(R.parent/'project26210_phase78_fine_state_spatial_v1/run_v1/sections_v1'/sid/'fine_state_abundance_v1.tsv.gz'),preparation=v.preparation,cohort=v.cohort,treatment=v.prior_treatment or 'UNKNOWN'))
 images=pd.read_csv(F/'tables/IMAGE_COORDINATE_SUPPORT.tsv',sep='\t').set_index('sample')
 qcs=[];sections=[];scales=[]
 for rec in records:
  sid=rec['sample'];fp=F/'spots'/(sid+'_F1_SPATIAL_INPUTS.tsv.gz');d=pd.read_csv(fp,sep='\t')
  gp=(S/'frozen_inputs' if rec['dataset']=='GSE175540' else A/'inputs')/(sid+'_adjacency.npz')
  G=sparse.load_npz(gp).tocsr();ab=pd.read_csv(rec['abundance_path'],sep='\t',index_col=0)
  if rec['dataset']=='GSE175540':ab=ab[[c for c in ab if c.startswith('q05__')]];ab.columns=[c.removeprefix('q05__') for c in ab]
  else:ab.columns=[c.removeprefix('q05cell_abundance_w_sf_') for c in ab]
  assert len(ab.columns)==13 and np.array_equal(ab.index,d.spot_id)
  assert np.isfinite(ab.to_numpy()).all() and (ab.to_numpy()>=0).all()
  # Current continuous-TYMP reference names verified against existing frozen fractions.
  spp=[c for c in ab if 'SPP1' in c];c1=[c for c in ab if 'C1QC' in c];oth=[c for c in ab if c=='Other Macrophage']
  assert len(spp)==len(c1)==len(oth)==1,ab.columns.tolist()
  total=ab.sum(axis=1).to_numpy();macro=ab[spp+c1+oth].sum(axis=1).to_numpy();sv=ab[spp[0]].to_numpy()
  d['V1_SPP1_q05']=sv;d['V2_SPP1_fraction']=sv/total;d['V3_SPP1_macro_share']=sv/macro
  d['C1QC_fraction']=ab[c1[0]].to_numpy()/total
  d['other_macro_fraction']=(macro-sv)/total;d['total_q05']=total
  np.testing.assert_allclose(d.V3_SPP1_macro_share,d.SPP1_share,atol=1e-12)
  np.testing.assert_allclose(macro/total,d.total_macro_fraction,atol=1e-12)
  for key,col in [('endothelial_fraction','Endothelial'),('pericyte_fraction','SMC/Pericyte')]:
   np.testing.assert_allclose(ab[col].to_numpy()/total,d[key],atol=1e-12)
  for key in ['T cell','T/NK','T Cell','T/NK cell']:
   if key in ab:d['Tcell_q05']=ab[key].to_numpy();d['Tcell_fraction']=ab[key].to_numpy()/total
  with h5py.File(rec['input_path']) as h:
   obs=read_elem(h['obs']);var=read_elem(h['var']);assert np.array_equal(obs.index,d.spot_id)
   X=read_elem(h['layers/counts'] if 'counts' in h.get('layers',{}) else h['X']).astype(float)
   assert np.isfinite(X.data).all() and (X.data>=0).all() and np.allclose(X.data,np.round(X.data))
   totals=np.asarray(X.sum(axis=1)).ravel();assert (totals>0).all()
   req=sorted(set(sum(programs.values(),[]))|{'TYMP','SPP1','POSTN'}|set(';'.join(lr.subunits).split(';')))
   present=[g for g in req if g in var.index];raw=X[:,var.index.get_indexer(present)].toarray()
   log=np.log1p(raw/totals[:,None]*10000);lookup={g:i for i,g in enumerate(present)}
   for gene in ['TYMP','SPP1','POSTN']:
    d[gene+'_RNA']=log[:,lookup[gene]] if gene in lookup else np.nan
   for key,genes in programs.items():
    used=[g for g in genes if g in lookup];idx=[lookup[g] for g in used]
    if idx:
     zz=log[:,idx];sd=zz.std(axis=0);score=((zz-zz.mean(axis=0))/np.where(sd==0,1,sd)).mean(axis=1)
     nd=int(((raw[:,idx]>0).sum(axis=1)>=min(2,len(idx))).sum())
    else:score=np.full(len(d),np.nan);nd=0
    gate=len(used)>=max(min(3,len(genes)),math.ceil(len(genes)*.5)) and nd>=20 and np.nanstd(score)>0
    d[key]=score if gate else np.nan
    qcs.append(dict(sample=sid,program=key,listed=len(genes),present=len(used),detected_spots=nd,eligible=gate,used_genes=';'.join(used),missing_genes=';'.join(sorted(set(genes)-set(used)))))
    if key in author and gate and len(used)==len(genes):np.testing.assert_allclose(score,d['program_'+key],atol=1e-10)
   d['V4_SPP1_program']=d.spp1_state_no_spp1
   for rr in lr.itertuples():
    genes=rr.subunits.split(';');ix=[lookup[g] for g in genes if g in lookup]
    nd=int((raw[:,ix]>0).all(axis=1).sum()) if len(ix)==len(genes) else 0
    ligand_nd=int((raw[:,lookup[rr.ligand]]>0).sum()) if rr.ligand in lookup else 0
    gate=len(ix)==len(genes) and nd>=20 and ligand_nd>=20
    d['LR_'+rr.receptor]=np.min(log[:,ix],axis=1) if gate else np.nan
    qcs.append(dict(sample=sid,program='LR_'+rr.interaction,listed=len(genes),present=len(ix),detected_spots=nd,ligand_detected_spots=ligand_nd,eligible=gate,used_genes=';'.join(g for g in genes if g in lookup),missing_genes=';'.join(g for g in genes if g not in lookup)))
  d['log1p_UMI']=np.log1p(totals)
  xy=d[['x','y']].to_numpy();rr,cc=G.nonzero();factor=55/float(images.loc[sid,'spot_diameter_fullres'])
  distances=np.linalg.norm(xy[rr]-xy[cc],axis=1)*factor;med=float(np.median(distances))
  physical=90<=med<=110 and np.isfinite(distances).all()
  WG=sparse.csr_matrix((distances,(rr,cc)),shape=G.shape);sparse.save_npz(P/'graphs'/(sid+'_ring.npz'),G)
  d['expression_distance_um']=np.nan
  if physical:
   rows=[];cols=[];vals=[]
   for start in range(0,len(d),256):
    ix=np.arange(start,min(start+256,len(d)));dst=dijkstra(WG,directed=False,indices=ix,limit=800)
    r,c=np.where(np.isfinite(dst)&(dst>0));rows.extend((r+start).tolist());cols.extend(c.tolist());vals.extend(dst[r,c].tolist())
   D=sparse.csr_matrix((vals,(rows,cols)),shape=G.shape)
   assert (D-D.T).max()<1e-8
   sparse.save_npz(P/'graphs'/(sid+'_distance800.npz'),D)
   ids=np.flatnonzero(d.expression_interface.to_numpy(bool))
   if len(ids):
    dd=dijkstra(WG,directed=False,indices=ids,min_only=True)
    d['expression_distance_um']=np.where(np.isfinite(dd),dd,np.nan)
   for radius in [800,400,200]:
    H=D.copy();H.data=(H.data<=radius).astype(float);H.eliminate_zeros()
    count=np.asarray(H.sum(axis=1)).ravel()
    scales.append(dict(sample=sid,scale=str(radius),physical_eligible=True,sources_ge3=int((count>=3).sum()),median_neighbors=float(np.median(count)),min_neighbors=int(count.min()),max_neighbors=int(count.max()),median_section_fraction=float(np.median(count/max(len(d)-1,1))),median_edge_um=med))
  else:
   for radius in [800,400,200]:scales.append(dict(sample=sid,scale=str(radius),physical_eligible=False,reason='IMAGE_DIAMETER_SPACING_INCONSISTENT',median_edge_um=med))
  deg=np.asarray(G.sum(axis=1)).ravel()
  scales.append(dict(sample=sid,scale='ring',physical_eligible=True,sources_ge3=int((deg>=3).sum()),median_neighbors=float(np.median(deg)),median_edge_um=med))
  for grid in [4,5]:
   bins=[pd.qcut(pd.Series(xy[:,axis]).rank(method='first'),q=grid,labels=False).to_numpy() for axis in range(2)]
   d['block_'+str(grid)]=bins[0]*grid+bins[1]
  tab('inputs/'+sid+'_INPUTS.tsv.gz',d)
  section=dict(spots=len(d),physical_eligible=physical,pixel_to_um=factor,median_edge_um=med,reference_columns=';'.join(ab.columns),treatment=rec.get('treatment','UNKNOWN'),preparation=rec.get('preparation','UNKNOWN'),cohort=rec.get('cohort',rec['dataset']))
  section={**rec,**section}
  sections.append(section)
  for path in [fp,gp,rec['input_path'],rec['abundance_path']]:
   m=meta(path)
   if m['bytes']<50000000:m['sha256']=sha(path)
   sources.append(m)
  print('A0',sid,len(d),'median edge um',round(med,3),'physical',physical,flush=True)
 tab('tables/INPUT_SECTION_REGISTER.tsv',pd.DataFrame(sections));tab('tables/PROGRAM_LR_ELIGIBILITY.tsv',pd.DataFrame(qcs));tab('tables/SCALE_ELIGIBILITY.tsv',pd.DataFrame(scales))
 save('audit/SOURCE_MANIFEST.json',sources)
 assert len(sections)==38 and sum(x['spots'] for x in sections)==114640
 save('A0_COMPLETE.json',dict(time=stamp(),sections=38,spots=114640,physical_sections=sum(x['physical_eligible'] for x in sections),endpoint_count=len(ep),association_not_started=True))
 save('SOFTWARE.json',dict(python=sys.version,numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,h5py=h5py.__version__))
 save('STATUS.json',dict(status='A0_COMPLETE',time=stamp()))
if __name__=='__main__':
 try:main()
 except Exception:
  (P/'ERROR_A0.txt').write_text(traceback.format_exc());raise
