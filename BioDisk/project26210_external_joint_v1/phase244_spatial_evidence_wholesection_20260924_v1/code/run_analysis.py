from pathlib import Path
import sys,json,datetime,hashlib,subprocess,traceback,ast,math
import numpy as np,pandas as pd
from scipy.stats import spearmanr,binomtest
P=Path(__file__).resolve().parents[1];ROOT=P.parent;MS=ROOT/'phase242_SpaMTP_12tumor_analysis_20260924_v1';ST=ROOT/'phase234_part4_replan_execution_20260921_v1'
def now():return datetime.datetime.now().astimezone().isoformat()
def dump(f,o):Path(f).write_text(json.dumps(o,ensure_ascii=False,indent=2,allow_nan=False))
def sha(f):return hashlib.sha256(Path(f).read_bytes()).hexdigest()
def gate(stage):
    r=subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True)
    with open(P/'logs/memory_gate.log','a') as f:f.write(stage+'\n'+r.stdout+r.stderr)
    if r.returncode:raise RuntimeError('memory gate '+stage)
def tab(d,f):d.to_csv(P/'tables'/f,sep='\t',index=False)
def snap(f):return {'path':str(f),'sha256':sha(f),'size':f.stat().st_size,'mtime_ns':f.stat().st_mtime_ns}
CFG={'statistics':{'permutations':100000,'bootstrap':10000}}
tree=ast.parse((ROOT/'phase240_bulk_immune_analysis_20260924_v1/code/pipeline.py').read_text())
exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ['rank','corr','spear','bh','evaluate']],type_ignores=[]),'frozen_numerical_functions','exec'),globals())
def main():
    gate('specification')
    cfg=json.loads((MS/'CONFIG.json').read_text())
    plan={'at':now(),'scope':'post-hoc exploratory whole-section between-patient association, distinct from within-section spatial co-localization; chosen to remove registration dependency','samples':cfg['samples'],'molecules':['Spermidine','Spermine'],'polarity':'pos','adduct':'M+H','ppm':5,'mass_identity':'candidate only','ST_endpoints':['TYMP_RNA_mean','SPP1_RNA_mean','SPP1_macrophage_ratio_of_sums'],'RNA_summary':'arithmetic mean existing log1p(counts/spot_total*10000) over all frozen QC spots','macro_summary':'sum V1_SPP1_q05 / sum total_macro_q05 across all frozen QC spots; avoids unstable mean of per-spot ratios','MSI_summary':'arithmetic mean peak_intensity/TIC*feature_count over every positive-TIC acquired pixel; zeros for absent peak intensity at measured pixels retained; missing mass candidate remains missing','normalization_sensitivity':'ratio of summed peak intensity to summed TIC times feature_count; same all-positive-TIC pixels, not selected by result','hypothesis_family':'six molecule-endpoint comparisons, one family per normalization specification; both fully reported, main is pixelwise mean','statistics':CFG['statistics'],'seed':26210244,'minimum_n':10,'no_registration_selection':True,'no_patient_selection_by_outcome':True,'no_new_figures':True,'old_local_spatial_tests':'retained unchanged, no test swapping'}
    dump(P/'PLAN_BEFORE_TESTS.json',plan)
    dump(P/'STATUS.json',{'status':'EXTRACTING_WHOLE_SECTION','at':now()})
    paths=[MS/'CONFIG.json',MS/'tables/TARGETED_MASS_CANDIDATES_PASS5PPM.tsv',MS/'tables/REGISTRATION_SUMMARY.tsv',MS/'tables/PEAK_COVERAGE.tsv',MS/'tables/PATIENT_LEVEL_PRIMARY_AND_SENSITIVITY.tsv',ST/'scripts/freeze.py']
    paths += [ST/'inputs'/f'{s}_INPUTS.tsv.gz' for s in cfg['samples']]
    paths += [MS/'results/matrices'/f'{s.replace("_","-")}-pos.npz' for s in cfg['samples']]
    before=[snap(f) for f in paths];dump(P/'INPUT_MANIFEST_BEFORE.json',before)
    peaks=pd.read_csv(paths[1],sep='\t')
    rows=[];massqc=[]
    for sample,patient in zip(cfg['samples'],cfg['patient_ids']):
        gate('sample_'+sample)
        d=pd.read_csv(ST/'inputs'/f'{sample}_INPUTS.tsv.gz',sep='\t')
        assert d.patient.nunique()==1 and d.patient.iloc[0]==patient and d.spot_id.is_unique
        assert (d.total_macro_q05>=0).all() and d.total_macro_q05.sum()>0 and (d.V1_SPP1_q05>=0).all()
        macro=float(d.V1_SPP1_q05.sum()/d.total_macro_q05.sum());assert 0<=macro<=1
        record=dict(sample=sample,patient=patient,ST_spots=len(d),TYMP_RNA_mean=float(d.TYMP_RNA.mean()),SPP1_RNA_mean=float(d.SPP1_RNA.mean()),SPP1_macrophage_ratio_of_sums=macro,mean_log1p_UMI=float(d.log1p_UMI.mean()),malignant_fraction_mean=float(d.malignant_fraction.mean()))
        z=np.load(MS/'results/matrices'/f'{sample.replace("_","-")}-pos.npz');X=z['intensity'];tic=z['tic'];valid=tic>0
        assert valid.sum()>50 and np.allclose(X.sum(1,dtype=np.float64),tic)
        for molecule in plan['molecules']:
            p=peaks[(peaks['sample']==sample)&(peaks.polarity=='pos')&(peaks.adduct=='M+H')&(peaks.molecule==molecule)]
            assert len(p)<=1
            q=dict(sample=sample,patient=patient,molecule=molecule,MSI_pixels=len(tic),positive_TIC_pixels=int(valid.sum()),features=X.shape[1],candidate_present=bool(len(p)))
            if len(p):
                j=int(p.feature_index.iloc[0]);assert abs(float(z['mz'][j])-float(p.observed_mz.iloc[0]))<1e-8
                signal=X[valid,j].astype(float);fr=signal/tic[valid]*X.shape[1]
                a=float(fr.mean());b=float(signal.sum()/tic[valid].sum()*X.shape[1])
                q.update(observed_mz=float(z['mz'][j]),feature_index=j,ppm=float(p.ppm.iloc[0]),detected_pixels=int((signal>0).sum()),detected_fraction=float((signal>0).mean()),mean_TIC_fraction=a,pooled_TIC_fraction=b)
                record[molecule+'_mean_TIC']=a;record[molecule+'_pooled_TIC']=b
            else:
                record[molecule+'_mean_TIC']=np.nan;record[molecule+'_pooled_TIC']=np.nan
                q['exclusion_reason']='no frozen mass-compatible candidate; not zero concentration'
            massqc.append(q)
        rows.append(record)
    data=pd.DataFrame(rows);assert len(data)==12 and data.patient.is_unique
    tab(data,'WHOLE_SECTION_PATIENT_INPUT.tsv');tab(pd.DataFrame(massqc),'WHOLE_SECTION_MASS_QC.tsv')
    gate('patient_associations')
    rng=np.random.default_rng(plan['seed']);out=[]
    for spec,suffix in [('main_pixelwise_mean','mean_TIC'),('normalization_sensitivity_pooled_TIC','pooled_TIC')]:
        for molecule in plan['molecules']:
            for target in plan['ST_endpoints']:
                s=data[[molecule+'_'+suffix,target]].dropna();n=len(s)
                r=evaluate(s.iloc[:,0].to_numpy(),s.iloc[:,1].to_numpy(),rng) if n>=10 else {'rho':spear(s.iloc[:,0].to_numpy(),s.iloc[:,1].to_numpy()),'p':None,'ci_low':None,'ci_high':None,'loo_min':None,'loo_max':None,'nonestimable':'fewer_than_10_independent_patients'}
                if n>=10:assert abs(r['rho']-spearmanr(s.iloc[:,0],s.iloc[:,1]).statistic)<1e-12
                out.append(dict(specification=spec,molecule=molecule,endpoint=target,n_patients=n,evidence='posthoc_exploratory_between_patient_whole_section',**r))
        rr=[r for r in out if r['specification']==spec]
        for r,q in zip(rr,bh([r['p'] if r['p'] is not None else 1. for r in rr])):r['FDR']=float(q) if r['p'] is not None else None
    result=pd.DataFrame(out);tab(result,'WHOLE_SECTION_ASSOCIATIONS_ALL.tsv')
    # Diagnose geometric ambiguity without changing or filtering original spatial inference.
    reg=pd.read_csv(MS/'tables/REGISTRATION_SUMMARY.tsv',sep='\t');pos=reg[reg.polarity=='pos']
    pq=pd.read_csv(MS/'tables/PEAK_COVERAGE.tsv',sep='\t')
    pq=pq[(pq['rank']==1)&(pq.coverage_mode=='candidate_specific')&(pq.polarity=='pos')&(pq.adduct=='M+H')&pq.molecule.isin(plan['molecules'])]
    tab(pos,'POSITIVE_MODE_REGISTRATION_DIAGNOSTIC.tsv');tab(pq,'PRIMARY_MAPPED_PEAK_COVERAGE.tsv')
    spatial={'positive_spectra':len(pos),'ambiguous_orientation_samples':int((pos.retained_candidates>1).sum()),'single_candidate_samples':int((pos.retained_candidates==1).sum()),'mask_variants_not_unanimous':int((pos.preferred_orientation_wins_mask_variants<pos.total_mask_variants).sum()),'min_Dice':float(pos.best_Dice.min()),'max_Dice':float(pos.best_Dice.max()),'anatomically_validated':int(pos.validated_registration.sum()),'main_mass_candidate_patients':data['Spermidine_mean_TIC'].notna().sum().item()}
    dump(P/'SPATIAL_TECHNICAL_DIAGNOSIS.json',spatial)
    gate('legacy_evidence_verification')
    # Keep complete original summary tables, and verify every estimable primary sign P from patient effects.
    validations=[];preserved=[]
    keys=['dataset','family','endpoint','exposure','target','layer','scale','qualification','model','effect_type']
    for prefix in ['A1','A2_A4']:
        summary=pd.read_csv(ST/'tables'/f'{prefix}_COHORT_SUMMARY.tsv',sep='\t',dtype={k:str for k in keys})
        effects=pd.read_csv(ST/'tables'/f'{prefix}_PATIENT_EFFECTS.tsv',sep='\t',dtype={k:str for k in keys})
        source=ST/'tables'/f'{prefix}_COHORT_SUMMARY.tsv'
        preserved.append(snap(source));tab(summary,f'LEGACY_{prefix}_COHORT_SUMMARY_UNCHANGED.tsv')
        for row in summary.itertuples():
            if pd.isna(row.P):continue
            mask=np.ones(len(effects),dtype=bool)
            for k in keys:mask &= effects[k].fillna('NA').eq(str(getattr(row,k)) if pd.notna(getattr(row,k)) else 'NA').to_numpy()
            g=effects[mask];v=g.effect.dropna().to_numpy()
            assert g.patient.is_unique and len(v)==row.n
            positive=int((v>0).sum());negative=int((v<0).sum())
            pv=float(binomtest(positive,positive+negative,.5,alternative='two-sided').pvalue)
            assert abs(pv-row.P)<1e-10,(prefix,row.endpoint,pv,row.P)
            validations.append(dict(source=prefix,dataset=row.dataset,endpoint=row.endpoint,scale=row.scale,model=row.model,n=row.n,P_saved=row.P,P_recomputed=pv))
        # Preserve all 400um SPP1 ligand and broad metabolism findings, including negative and unestimable.
        sel=summary[(summary.exposure.eq('SPP1_RNA')&summary.scale.eq('400'))|(summary.exposure.eq('polyamine_mta_context')&summary.scale.isin(['W0','W1_4x4','W1_5x5','400']))]
        tab(sel,f'LEGACY_{prefix}_RELEVANT_ALL_DIRECTIONS.tsv')
    tab(pd.DataFrame(validations),'LEGACY_INDEPENDENT_SIGN_CHECKS.tsv')
    for rec in before:
        f=Path(rec['path']);assert f.stat().st_size==rec['size'] and f.stat().st_mtime_ns==rec['mtime_ns'] and sha(f)==rec['sha256']
    dump(P/'PRESERVED_LEGACY_SOURCES.json',preserved)
    lines=['# 空间证据复核与免配准整片对应','','## 已有空间证据可保留，MSI局部对应仍受配准约束','Phase234的多胺/MTA转录背景—SPP1绝对丰度具有患者间一致性：GSE175540全组织23/23正，中位rho0.2256，P2.384e-7/FDR7.153e-7；TJ-RCC12/12正，rho0.1634，P0.0004883/FDR0.005859。转录背景不是实测MTA。','固定400um SPP1—ITGAV/ITGB1表达兼容性在TJ-RCC为11/12正，rho0.2082，P0.006348/FDR0.05713，属于名义显著；GSE175540同一主检验P0.09314，未达到名义显著。CD44、CAF/ECM及既有血管主要阴性完整保留。不能据表达兼容性声称直接细胞通讯。',f'既有患者层符号检验本轮独立复核{len(validations)}项一致；完整原汇总表和对应阴性已保存。',f'MSI正离子12组中{spatial["ambiguous_orientation_samples"]}组有多个保留朝向，{spatial["mask_variants_not_unanimous"]}组在技术掩模设置下朝向选择不一致；独立解剖确认0组。原六个局部端点与其P/FDR保持不变。','','## 整片患者间对应的完整结果','整片汇总不依赖ST—MSI坐标配准，回答患者间组织平均相对信号是否对应；本轮为观察既有结果后预先冻结的新探索性问题，不是对原局部空间主终点的重新验证。']
    for r in out:
        if r['specification']!='main_pixelwise_mean':continue
        lines.append(f"- {r['molecule']}—{r['endpoint']}: n={r['n_patients']}, rho={r['rho']:.6f}, P={r['p']:.6g}, FDR={r['FDR']:.6g}, bootstrap95%CI=[{r['ci_low']:.6f},{r['ci_high']:.6f}]。")
    lines+=['','## Methods','全部12例登记肿瘤患者按固定名单纳入。沿用正离子M+H、5ppm峰资格，峰缺失保留为不可估计；所有正TIC采集像素参与，未检出强度为0的像素仍保留。主量为逐像素峰/TIC相对信号平均，乘以固定特征数；总强度之比作为全部六端点的归一化对照，均报告不择优。','TYMP/SPP1采用所有已冻结QC spot的log1p每万计数平均；SPP1巨噬比例为整片SPP1 q05之和/总成熟巨噬q05之和，不作跨片绝对q05比较。患者层双侧Spearman置换100000次、bootstrap10000次，seed26210244；每个归一化口径六项共同BH。病例内区域或像素不充当独立样本。','','## Discussion与下一步','MSI为相对离子强度和质量相容候选，缺少结构确证和内标；邻近切片的组织构成、采集/归一化及离子抑制仍可影响整片比较，阳性也不能写成局部共定位或MTA中介。主分析、所有归一化对照及原局部阴性并列保留。','局部空间关系的优先改进为取得MSI自身光学图或作者变换，并由独立形态标志核对朝向；不以TYMP/SPP1表达或相关方向选配准。若无此类依据，应以已验证的ST空间主线展示，MSI保留为候选探索。可以设计新验证数据，但不继续围绕P<0.05切换半径、ROI、患者或显著性检验。','GLMF两例完整组织可保留主案例展示，但不能扩增患者数；目前没有新增图或修改既有正文。']
    (P/'REPORT_CN.md').write_text('\n'.join(lines))
    qa=dict(status='PASS',at=now(),registered_patients=12,new_exploratory_main_tests=6,normalization_sensitivity_tests=6,original_spatial_tests_unchanged=True,legacy_sign_tests_verified=len(validations),old_input_hashes_rechecked=len(before),new_figures=0,registration_independent_whole_section=True)
    dump(P/'FINAL_VALIDATION.json',qa);dump(P/'STATUS.json',dict(status='COMPLETE_VERIFIED',at=now(),validation=qa));(P/'EXIT_STATUS').write_text('0\n');(P/'RUN_FINISHED').write_text(now())
    print(result.to_string(index=False));print(json.dumps(spatial))
try:main()
except BaseException:
    (P/'logs/ERROR.txt').write_text(traceback.format_exc());dump(P/'STATUS.json',{'status':'FAILED','at':now(),'error':traceback.format_exc()});(P/'EXIT_STATUS').write_text('1\n');(P/'RUN_FINISHED').write_text(now());raise

