from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
P=Path(__file__).resolve().parents[1];R=P.parent
for d in ['plot_data','delivery']: (P/d).mkdir(exist_ok=True)
L=['Endothelial','Malignant','Myeloid','SMC_Pericyte']; TF=['HNF4A','CEBPD','JUNB','STAT3']
stat=pd.read_csv(P/'differential/RCC100_chr1_statistics.txt',sep='\t').set_index('Motif')
assert len(stat)==4 and (stat.Num>=30).all() and np.isfinite(stat.to_numpy()).all()
motif_plan=json.loads((P/'metadata/FOUR_MOTIFS_FROZEN.json').read_text()); motifs={r['TF']:r['matrix_id'] for r in motif_plan}
assert set(stat.index)==set(motifs.values())
union=pd.read_csv(P/'motif_matching/UNION_UNIQUE_mpbs.bed',sep='\t',header=None)
checks=[]
for tf,motif in motifs.items():
    m=union[union[3]==motif]; lens=(m[2]-m[1]).unique();assert len(lens)==1
    n=int(lens[0]); assert len(m)==int(stat.loc[motif,'Num'])
    f=P/'differential/Lineplots'/f'{motif}.txt'
    profiles=pd.read_csv(f,sep='\t');assert profiles.shape==(200,4) and np.isfinite(profiles.to_numpy()).all()
    for lin in L:
        v=profiles[lin].to_numpy();nc=v[100-n//2:100+n//2].sum()
        left=v[100-n//2-n:100-n//2].sum();right=v[100+n//2:100+n//2+n].sum()
        ps=(left+right-2*nc)/n;tc=(v.sum()-nc)/(200-n)
        error=max(abs(ps-stat.loc[motif,'Protection_Score_'+lin]),abs(tc-stat.loc[motif,'TC_'+lin]))
        assert error<5e-7,(motif,lin,error)
        checks.append(dict(motif=motif,lineage=lin,profile_formula_max_error=float(error)))
activity=pd.DataFrame({tf:[stat.loc[motifs[tf],'Protection_Score_'+l]+stat.loc[motifs[tf],'TC_'+l] for l in L] for tf in TF},index=L)
activity.index.name='lineage';activity.to_csv(P/'tables/ATAC_CHR1_TF_ACTIVITY.tsv',sep='\t')
reg=pd.read_csv(R/'phase234_part4_replan_execution_20260921_v1/tables/INPUT_SECTION_REGISTER.tsv',sep='\t')
groups=json.loads((P/'metadata/LINEAGE_TYPING_PLAN.json').read_text())['groups'];summary=[];inputs=[]
for name in ['GLMF1','GLMF2']:
    path=Path(reg[reg['sample']==name].iloc[0].abundance_path)
    ab=pd.read_csv(path,sep='\t',index_col=0);ab.columns=ab.columns.str.replace('q05cell_abundance_w_sf_','',regex=False)
    coords=pd.read_csv(R/f'phase234_part4_replan_execution_20260921_v1/inputs/{name}_INPUTS.tsv.gz',sep='\t').set_index('spot_id')
    assert ab.index.is_unique and coords.index.is_unique and set(ab.index)==set(coords.index)
    ab=ab.loc[coords.index]; assert set(ab.columns)==set(sum(groups.values(),[]))
    assert (ab.to_numpy()>=0).all() and np.isfinite(ab.to_numpy()).all()
    broad=pd.DataFrame({l:ab[groups[l]].sum(axis=1) for l in L})
    coverage=broad.sum(axis=1)/ab.sum(axis=1)
    weights=broad.div(broad.sum(axis=1),axis=0);assert np.allclose(weights.sum(axis=1),1)
    projected=weights.to_numpy()@activity.to_numpy()
    independent=np.column_stack([sum(weights[l].to_numpy()*activity.loc[l,tf] for l in L) for tf in TF])
    err=float(np.max(np.abs(projected-independent)));assert err<1e-12
    assert np.all(projected>=activity.min(axis=0).to_numpy()-1e-12) and np.all(projected<=activity.max(axis=0).to_numpy()+1e-12)
    eligible=coverage>=.8;projected[~eligible.to_numpy(),:]=np.nan
    out=coords[['x_hires_px','y_hires_px','physical_x_um','physical_y_um']].copy()
    out['represented_fraction']=coverage;out['projection_eligible']=eligible
    for l in L:out['conditional_weight_'+l]=weights[l]
    for j,tf in enumerate(TF):out[tf]=projected[:,j]
    out.to_csv(P/f'plot_data/{name}_spatial_ATAC_projection.tsv',sep='\t',index_label='spot_id')
    summary.append(dict(sample=name,spots=len(out),eligible=int(eligible.sum()),masked=int((~eligible).sum()),eligible_percent=float(eligible.mean()*100),median_represented_fraction=float(coverage.median()),projection_formula_max_error=err))
    inputs.append(dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
meta=dict(status='PASS',scope='RCC100 chromosome 1 technical pilot',formula='Sum of conditional broad-lineage fractions * (HINT Protection_Score + TC)',normalization='HINT factors estimated over the four predefined motifs only',minimum_unique_motif_sites=30,minimum_represented_fraction=.8,missing_lineages='Lymphoid, Mast, Normal_epithelial: omitted from conditional projection, never assigned zero activity',samples=summary,profile_formula_checks=checks,inputs=inputs,inferential_tests=0)
(P/'qc/ATAC_PROJECTION_VALIDATION.json').write_text(json.dumps(meta,indent=2))
print(activity.to_string());print(json.dumps(summary,indent=2));print('PROJECTION_EXIT_STATUS=0')
