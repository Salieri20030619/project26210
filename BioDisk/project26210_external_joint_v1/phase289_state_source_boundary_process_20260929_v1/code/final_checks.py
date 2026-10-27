from common import *
import h5py
from anndata.io import read_elem
from scipy import sparse,stats
def main():
 gate('final_raw_and_topology_checks')
 reg=read(R/'phase234_part4_replan_execution_20260921_v1/tables/INPUT_SECTION_REGISTER.tsv')
 checked=0
 for sample in ['GLMF1','H46_T','GSM5924030']:
  z=np.load(P/'inputs'/f'{sample}_PANELS.npz');r=reg[reg['sample']==sample].iloc[0]
  with h5py.File(r.input_path,'r') as h:
   gs=np.asarray(read_elem(h['var/gene_symbol_original'])) if 'gene_symbol_original' in h['var'] else np.asarray(read_elem(h['var'][h['var'].attrs['_index']]))
   ids=np.asarray(read_elem(h['obs'][h['obs'].attrs['_index']]))
   node=h['layers/counts'] if 'layers/counts' in h else h['X']
   for k in np.linspace(0,len(z['spot_id'])-1,8,dtype=int):
    i=int(np.flatnonzero(ids==z['spot_id'][k])[0]);lo,hi=node['indptr'][i:i+2]
    row=sparse.csr_matrix((node['data'][lo:hi],node['indices'][lo:hi],[0,hi-lo]),shape=(1,len(gs)))
    vals=row[:,[np.flatnonzero(gs==g)[0] for g in z['genes']]].toarray().ravel()
    assert np.array_equal(vals,z['counts'][k]);checked+=1
 top=read(P/'tables/BOUNDARY_CONTINUITY_DESCRIPTIVE.tsv')
 for r in top.itertuples():
  d=read(P/'inputs'/f'{r.sample}_{r.definition}_{r.radius}_EDGES.tsv.gz')
  x=d.macro_out_g_SPP1.to_numpy();y=d.out_p_CAF_activation.to_numpy();ok=np.isfinite(x)&np.isfinite(y);high=np.zeros(len(d),bool)
  if ok.any():
   high[ok]=(x[ok]>=np.quantile(x[ok],.75))&(y[ok]>=np.quantile(y[ok],.75))&(x[ok]>0)&(y[ok]>0)
   joint=(stats.rankdata(x[ok])/(ok.sum()+1))*(stats.rankdata(y[ok])/(ok.sum()+1))
   assert np.array_equal(joint,d.loc[ok,'joint_rank'].to_numpy())
  assert np.array_equal(high,d.joint_high)
  good=ok&d.geometry_eligible.to_numpy()
  assert good.sum()==r.evaluable_edges and (high&good).sum()==r.joint_high_edges
 pt=read(P/'tables/LOCAL_PATIENT_EFFECTS.tsv');rows=[]
 keys=['dataset','radius','model','family','endpoint']
 for key,g in pt.groupby(keys):
  z=g.pivot(index='patient',columns='definition',values='effect').dropna()
  for name in z.columns:
   vals=z[name].to_numpy();rows.append(dict(zip(keys,key),definition=name,patients=len(vals),positive=int((vals>0).sum()),negative=int((vals<0).sum()),median_effect=float(np.median(vals)) if len(vals) else np.nan,P=np.nan,interpretation='common_patient_descriptive_only'))
 tab(rows,'LOCAL_COMMON_PATIENT_DESCRIPTIVE.tsv')
 # Explicit documented numerical/data fingerprints for complete delivery.
 for n in ['states','boundary','process','local_tests','validate_state_process','validate_boundary_local','descriptive_diagnostics']:assert (P/'logs'/f'{n}_exit.txt').read_text().strip()=='0'
 freeze=json.loads((P/'FREEZE.json').read_text());assert sha(P/'PROTOCOL_BEFORE_ANALYSIS.md')==freeze['protocol_sha256'] and sha(P/'PANELS_FROZEN.json')==freeze['panels_sha256']
 dump(P/'qc/FINAL_VALIDATION.json',dict(status='PASS',new_spatial_panel_source_rows=checked,source_formats=3,all_joint_rank_and_coverage_counts_match=True,common_patient_comparison_saved=True,source_program_and_boundary_validations='PASS',final_analysis_exits_zero=True,protocol_unchanged=True,no_new_figures=True))
 print('FINAL_VALIDATION_PASS',flush=True)
if __name__=='__main__':main()
