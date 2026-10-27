from source import *
from spatial import CORE
def main():
 checks=[]
 for job in ['source','spatial','validation','context','external']:
  assert (P/'logs'/f'{job}_exit.txt').read_text().strip()=='0',job
 for name in ['INTERNAL_VALIDATION.json','CONTEXT_VALIDATION.json','EXTERNAL_VALIDATION.json']:assert json.loads((P/'qc'/name).read_text())['status']=='PASS'
 d=read(P/'tables/CONTEXT_SENSITIVITY_SUMMARY.tsv');e=read(P/'tables/CONTEXT_SENSITIVITY_PATIENT_EFFECTS.tsv')
 for r in d.itertuples():
  v=e[(e.group==r.group)&(e.endpoint==r.endpoint)].effect.dropna().to_numpy()
  assert len(v)==r.n
  if len(v)>=10:assert abs(stats.binomtest(np.count_nonzero(v>0),np.count_nonzero(v),.5).pvalue-r.P)<1e-12
 assert np.allclose(stats.false_discovery_control(d.P.fillna(1)),d.FDR)
 checks.append(dict(check='context_sign_BH',endpoints=len(d)))
 inter=read(P/'tables/SPATIAL_INTERACTION_DETAILS.tsv');err=0.;nn=0
 for r in inter[inter.effect.notna()].itertuples():
  z=read(P/'inputs'/f'{r.sample}_SPLIT_INTERFACE.tsv.gz');cov=CORE+(['outside_TYMP'] if r.model=='outside_TYMP_adjusted' else [])
  a=z.loc[z.split_eligible,['inside_TYMP','SPP1_POSTN_coverage',r.endpoint]+cov].dropna().to_numpy();a=(a-a.mean(0))/a.std(0)
  D=np.column_stack([np.ones(len(a)),a[:,0],a[:,1],a[:,0]*a[:,1],a[:,3:]])
  Q,RR=np.linalg.qr(D,mode='reduced');beta=np.linalg.solve(RR,Q.T@a[:,2])
  assert abs(beta[3]-r.effect)<1e-10 and abs(beta[2]-beta[3]-r.coverage_slope_TYMP_low)<1e-10 and abs(beta[2]+beta[3]-r.coverage_slope_TYMP_high)<1e-10
  err=max(err,abs(beta[3]-r.effect));nn+=1
 checks.append(dict(check='interaction_coefficients_and_conditional_slopes_QR',n=nn,max_error=err))
 # Same-library CD8 values, all eight source estimates, no inferential P with nine patients.
 df=read(P/'tables/TYMP_CD8_FUNCTION_SUMMARY.tsv');pa=read(P/'tables/TYMP_CD8_MATCHED_PATIENTS.tsv')
 for r in df.itertuples():
  z=pa[pa.endpoint==r.endpoint];other='macro_TYMP' if r.source=='tumor_TYMP' else 'tumor_TYMP';A=z[[r.source,'y']].rank().to_numpy()
  cols=[np.ones(len(z))]+[z.study.eq(s).to_numpy(float) for s in sorted(z.study.unique())[1:]]
  if r.model=='study_other_source_adjusted':cols.append(z[other].rank().to_numpy())
  D=np.column_stack(cols);Q=np.linalg.qr(D,mode='reduced')[0];res=A-Q@(Q.T@A)
  assert abs(stats.pearsonr(res[:,0],res[:,1]).statistic-r.rho)<1e-10 and r.n==9 and pd.isna(r.P)
 checks.append(dict(check='CD8_rhos_and_no_small_n_P',n=len(df)))
 external=read(P/'tables/EXTERNAL_FIXED_GENE_SAMPLE_TPM.tsv').set_index('gene');g=read(P/'tables/EXTERNAL_FIXED_GENE_RESULTS.tsv').set_index('gene')
 for gene in g.index:
  control=external.loc[gene,['GSM5226831','GSM5226832','GSM5226833']].to_numpy(float);treated=external.loc[gene,['GSM5226837','GSM5226838','GSM5226839']].to_numpy(float)
  assert abs(np.log2(treated.mean()+1)-np.log2(control.mean()+1)-g.at[gene,'TYMP_vs_MCSF_log2_meanTPMplus1'])<1e-10
 for gene in ['SPP1','LGALS3']:assert external.loc[gene,['GSM5226837','GSM5226838','GSM5226839']].max()<external.loc[gene,['GSM5226831','GSM5226832','GSM5226833']].min()
 checks.append(dict(check='external_log_effects_and_all_control_comparison',genes=12))
 assert sha(P/'PROTOCOL_BEFORE_TESTING.md')==json.loads((P/'FREEZE.json').read_text())['protocol_sha256']
 final=dict(status='PASS',at=now(),all_analysis_exit_codes=0,checks=checks,inherited_validation_files=['INTERNAL_VALIDATION.json','CONTEXT_VALIDATION.json','EXTERNAL_VALIDATION.json'],source_inputs_unchanged=True,figures=0,manuscript_changed=False,wetlab_reanalysed=False)
 dump(P/'qc/FINAL_VALIDATION.json',final)
 print('FINAL_VALIDATION_PASS',checks,flush=True)
if __name__=='__main__':main()
