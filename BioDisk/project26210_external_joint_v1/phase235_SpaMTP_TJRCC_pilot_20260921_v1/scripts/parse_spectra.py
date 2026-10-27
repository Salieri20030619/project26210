from pathlib import Path
import json,hashlib,xml.etree.ElementTree as ET,uuid,re,gc,sys,datetime
import numpy as np,pandas as pd
P=Path(__file__).resolve().parents[1]
ns={'m':'http://psi.hupo.org/ms/mzml'}
def cvs(node):return {c.get('accession'):c.get('value','') for c in node.findall('.//m:cvParam',ns)}
def tab(n,d):d.to_csv(P/'tables'/n,sep='\t',index=False)
def sha(f):return hashlib.sha256(Path(f).read_bytes()).hexdigest()
def main():
 # Frozen targeted mass eligibility only, no molecule identity from mass alone.
 formulas={'MTA':'C11H15N5O3S','SAM_neutral_zwitterion':'C15H22N6O5S','SAH':'C14H20N6O5S','Methionine':'C5H11NO2S','Spermidine':'C7H19N3','Spermine':'C10H26N4','Putrescine':'C4H12N2','Ornithine':'C5H12N2O2','Thymidine':'C10H14N2O5','Thymine':'C5H6N2O2','Uridine':'C9H12N2O6','Uracil':'C4H4N2O2','Adenosine':'C10H13N5O4'}
 atom={'C':12.,'H':1.00782503223,'N':14.00307400443,'O':15.99491461957,'S':31.9720711744}
 mass=lambda f:sum(atom[a]*int(n or 1) for a,n in re.findall(r'([A-Z][a-z]?)(\d*)',f))
 proton=1.007276466621
 adducts={'pos':{'M+H':(proton,1),'M+Na':(22.989218,1),'M+K':(38.963158,1),'M+NH4':(18.033823,1),'M+2H':(2*proton,2)},'neg':{'M-H':(-proton,1),'M+Cl':(34.969402,1),'M+FA-H':(44.998201,1),'M-H2O-H':(-18.010565-proton,1)}}
 lock=dict(time=datetime.datetime.now().astimezone().isoformat(),samples=['R114_T','Y12_T'],selection='two smallest compressed complete tumor neg/pos pairs, before spectra inspection',mass_tolerance_ppm=5,formulas=formulas,adducts=adducts,interpretation='mass-compatible candidate, not confirmed metabolite; unknown unit grid coords not micron',no_spot_P=True,new_figures=False)
 (P/'metadata/TECHNICAL_PILOT_FREEZE.json').write_text(json.dumps(lock,indent=2))
 qcs=[];candidates=[];axes=[]
 (P/'results/matrices').mkdir(exist_ok=True)
 for xmlfile in sorted((P/'extracted').glob('*.imzML')):
  label=xmlfile.stem;sample=label.rsplit('-',1)[0].replace('-','_');pol=label.rsplit('-',1)[1]
  ibd=xmlfile.with_suffix('.ibd');blob=ibd.read_bytes()
  root=ET.parse(xmlfile).getroot()
  content=cvs(root.find('m:fileDescription/m:fileContent',ns))
  expectedsha=content['IMS:1000091'];actualsha=hashlib.sha1(blob).hexdigest()
  assert actualsha==expectedsha,(label,actualsha,expectedsha)
  actualuuid=str(uuid.UUID(bytes=blob[:16]));assert actualuuid==content['IMS:1000080']
  groups={x.get('id'):cvs(x) for x in root.findall('m:referenceableParamGroupList/m:referenceableParamGroup',ns)}
  settings=cvs(root.find('m:scanSettingsList',ns))
  sl=root.find('m:run/m:spectrumList',ns);n=int(sl.get('count'));coords=[];values=[];axis=None;offsets=[]
  for s in sl:
   sc=cvs(s.find('m:scanList',ns));coords.append((int(sc['IMS:1000050']),int(sc['IMS:1000051'])))
   arrays={}
   for b in s.findall('m:binaryDataArrayList/m:binaryDataArray',ns):
    ref=b.find('m:referenceableParamGroupRef',ns).get('ref');a={**groups[ref],**cvs(b)}
    off=int(a['IMS:1000102']);count=int(a['IMS:1000103']);length=int(a['IMS:1000104'])
    dt='<f8' if 'MS:1000523' in a else '<f4';assert length==count*np.dtype(dt).itemsize and off+length<=len(blob)
    arr=np.frombuffer(blob,dtype=dt,count=count,offset=off)
    arrays[ref]=arr;offsets.append((ref,off,count,length))
   if axis is None:axis=arrays['mzArray'].copy()
   else:assert np.array_equal(axis,arrays['mzArray'])
   values.append(arrays['intensityArray'].copy())
  X=np.vstack(values);xy=np.array(coords)
  assert X.shape==(n,len(axis)) and len(np.unique(xy,axis=0))==n
  assert np.isfinite(X).all() and (X>=0).all() and np.isfinite(axis).all() and (np.diff(axis)>0).all()
  tic=X.sum(axis=1,dtype=np.float64);active=tic>0
  qcs.append(dict(sample=sample,polarity=pol,pixels=n,features=len(axis),nonzero_pixels=int(active.sum()),zero_pixels=int((~active).sum()),nonzero_features=int((X.sum(axis=0)>0).sum()),mz_min=float(axis.min()),mz_max=float(axis.max()),x_min=int(xy[:,0].min()),x_max=int(xy[:,0].max()),y_min=int(xy[:,1].min()),y_max=int(xy[:,1].max()),pixel_size_x=settings.get('IMS:1000046'),pixel_size_y=settings.get('IMS:1000047'),units='array grid; physical units unspecified in deposited file',uuid_verified=True,ibd_sha1_verified=True,matrix_nonnegative=True))
  np.savez_compressed(P/'results/matrices'/(label+'.npz'),intensity=X,mz=axis,xy=xy,tic=tic)
  rows=[]
  for j,mz in enumerate(axis):
   rows.append(dict(sample=sample,polarity=pol,feature_index=j,mz=float(mz),nonzero_pixels=int((X[:,j]>0).sum()),sum_intensity=float(X[:,j].sum(dtype=np.float64))))
  axes.extend(rows)
  for molecule,formula in formulas.items():
   for add,(shift,charge) in adducts[pol].items():
    theo=(mass(formula)+shift)/charge;ppm=(axis-theo)/theo*1e6;j=int(np.argmin(abs(ppm)))
    matches=np.flatnonzero(abs(ppm)<=5)
    ids=matches if len(matches) else [j]
    for k in ids:
     val=X[:,k];positive=val>0
     candidates.append(dict(sample=sample,polarity=pol,molecule=molecule,formula=formula,adduct=add,theoretical_mz=theo,observed_mz=float(axis[k]),ppm=float(ppm[k]),within_5ppm=bool(abs(ppm[k])<=5),feature_index=int(k),detected_pixels=int(positive.sum()),detected_fraction_all=float(positive.mean()),median_positive_intensity=float(np.median(val[positive])) if positive.any() else 0,identity='mass_candidate_only' if abs(ppm[k])<=5 else 'no_mass_candidate',is_nearest_nonmatch=not bool(len(matches))))
  # Preserve raw intensity exactly enough for round-trip float32 when feeding official ReadSM_mtx.
  dd=pd.DataFrame(X[active],columns=[format(m,'.12g') for m in axis]);dd.insert(0,'y',xy[active,1]);dd.insert(0,'x',xy[active,0])
  dd.to_csv(P/'results/matrices'/(label+'_SpaMTP_input.csv'),index=False,float_format='%.9g')
  tab(label+'_PIXEL_QC.tsv',pd.DataFrame(dict(x=xy[:,0],y=xy[:,1],tic=tic,nonzero_features=(X>0).sum(axis=1),eligible_nonzero=active)))
  print('SPECTRA',qcs[-1],flush=True)
  del root,sl,blob,values,X,dd;gc.collect()
 tab('SPECTRAL_TECHNICAL_QC.tsv',pd.DataFrame(qcs));tab('ALL_MZ_FEATURES.tsv',pd.DataFrame(axes));cc=pd.DataFrame(candidates)
 tab('TARGETED_MASS_CANDIDATES_ALL.tsv',cc);tab('TARGETED_MASS_CANDIDATES_PASS5PPM.tsv',cc[cc.within_5ppm])
 print(cc[cc.within_5ppm].to_string(index=False),flush=True)
 (P/'RAW_SPECTRA_VALIDATION.json').write_text(json.dumps(dict(status='PASS',spectra_files=len(qcs),all_ibd_sha1_uuid_pass=True,matrix_count=len(qcs),pixels=sum(q['pixels'] for q in qcs),mass_candidates=int(cc.within_5ppm.sum()),official_molecule_identity_claim=False),indent=2))
if __name__=='__main__':main()
