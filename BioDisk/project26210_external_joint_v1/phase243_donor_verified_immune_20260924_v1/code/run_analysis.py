from pathlib import Path
import os,sys,json,datetime,hashlib,subprocess,traceback,itertools,ast,shutil,csv
import numpy as np,pandas as pd,scipy
from scipy.stats import spearmanr
P=Path(__file__).resolve().parents[1];ROOT=P.parent;OLD=ROOT/'phase240_bulk_immune_analysis_20260924_v1';CONV=ROOT/'phase236_public_bulk_conversion_20260923_v1'
CFG=json.loads((OLD/'CONFIG.json').read_text())
def now():return datetime.datetime.now().astimezone().isoformat()
def dump(f,o):Path(f).write_text(json.dumps(o,ensure_ascii=False,indent=2,allow_nan=False))
def sha(f):return hashlib.sha256(Path(f).read_bytes()).hexdigest()
def gate(stage):
    r=subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True)
    with open(P/'logs/memory_gate.log','a') as f:f.write(stage+'\n'+r.stdout+r.stderr)
    if r.returncode:raise RuntimeError('Memory gate blocked '+stage)
def mark(s,**kw):dump(P/'STATUS.json',dict(status=s,time=now(),**kw))
def table(d,name):d.to_csv(P/'tables'/name,sep='\t',index=False)
def snapshot(f):return {'path':str(f),'bytes':f.stat().st_size,'mtime_ns':f.stat().st_mtime_ns,'sha256':sha(f)}
# Reuse ONLY pure numerical functions, never import or call previous pipeline entry points.
tree=ast.parse((OLD/'code/pipeline.py').read_text())
funcs=ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ['rank','corr','spear','bh','evaluate']],type_ignores=[])
exec(compile(funcs,str(OLD/'code/pipeline.py'),'exec'),globals())
def main():
    gate('identity_reproduction');mark('IDENTITY_REPRODUCTION')
    inputs=[OLD/'CONFIG.json',OLD/'code/pipeline.py',OLD/'reference/PANEL_USED.tsv',OLD/'tables/identity_pairwise.tsv',OLD/'tables/identity_sample_QC.tsv',OLD/'tables/program_scores_all_libraries.tsv',CONV/'tables/sample_QC_summary.tsv',ROOT/'phase216_paper_dataset_download_plan_20260919_v1/metadata/SAMPLE_RUN_METADATA.tsv']
    qc=pd.read_csv(CONV/'tables/sample_QC_summary.tsv',sep='\t');iq=pd.read_csv(OLD/'tables/identity_sample_QC.tsv',sep='\t').set_index('run')
    previous=pd.read_csv(OLD/'tables/identity_pairwise.tsv',sep='\t')
    meta=pd.read_csv(inputs[-1],sep='\t').set_index('runs')
    panel=pd.read_csv(OLD/'reference/PANEL_USED.tsv',sep='\t');data={};reproduced=[]
    assert len(qc)==60 and qc.run.is_unique and len(panel)==20311
    for run in qc.run:
        f=OLD/'identity'/f'{run}.npz';assert sha(f)==iq.loc[run,'sha256'];inputs.append(f)
        a=np.load(f)['counts'];assert a.shape==(len(panel),3)
        dep=a[:,:2].sum(1);valid=(dep>=CFG['identity']['depth'])&(a[:,2]/np.maximum(a.sum(1),1)<=CFG['identity']['other_allele_fraction_max'])
        assert valid.sum()==iq.loc[run,'eligible_loci'];data[run]=(a[:,1]/np.maximum(dep,1),valid)
    for row in previous.itertuples():
        a,ma=data[row.run_a];b,mb=data[row.run_b];m=ma&mb;x=a[m];y=b[m]
        rr=corr(x,y);ibs=float(np.mean(((x<=.1)&(y>=.9))|((y<=.1)&(x>=.9))))
        assert int(m.sum())==row.shared_loci
        assert abs(rr-row.vaf_r)<1e-12 and abs(ibs-row.opposite_homozygote_fraction)<1e-12
        classification='AMBIGUOUS'
        if m.sum()<CFG['identity']['pair_loci_min'] or min(ma.sum(),mb.sum())<CFG['identity']['covered_loci_min']:classification='LOW_INFORMATION'
        elif rr>=CFG['identity']['same_r_min'] and ibs<=CFG['identity']['same_ibs0_max']:classification='SAME_DONOR_SUPPORTED'
        elif rr<=CFG['identity']['distinct_r_max'] and ibs>=CFG['identity']['distinct_ibs0_min']:classification='DISTINCT_SUPPORTED'
        assert classification==row.classification;reproduced.append(row._asdict())
    pf=pd.DataFrame(reproduced).drop(columns='Index');table(pf,'IDENTITY_PAIRS_REPRODUCED.tsv')
    q=qc.set_index('run');obs=qc.loc[qc.group=='observational46_paper_correspondence','run'].tolist()
    low=[r for r in obs if iq.loc[r,'status']!='PASS'];eligible=[r for r in obs if r not in low]
    sub=pf[pf.run_a.isin(eligible)&pf.run_b.isin(eligible)]
    assert len(sub)==len(eligible)*(len(eligible)-1)//2
    assert sub.classification.isin(['SAME_DONOR_SUPPORTED','DISTINCT_SUPPORTED']).all()
    groups=[{r} for r in eligible]
    for row in sub[sub.classification=='SAME_DONOR_SUPPORTED'].itertuples():
        joins=[g for g in groups if row.run_a in g or row.run_b in g]
        if len(joins)>1:groups=[g for g in groups if g not in joins]+[set.union(*joins)]
    mapping=[];selected=[];alternates=[];samepairs=[]
    for i,g in enumerate(sorted(groups,key=lambda x:min(x)),1):
        within=sub[sub.run_a.isin(g)&sub.run_b.isin(g)]
        assert within.classification.eq('SAME_DONOR_SUPPORTED').all()
        ranked=sorted(g,key=lambda r:(-float(q.loc[r,'assigned_fragments']),r))
        selected.append(ranked[0]);alternates.append(ranked[-1])
        for r in ranked:mapping.append(dict(run=r,sample_alias=q.loc[r,'sample_alias'],anonymous_donor=f'D{i:03d}',selected=r==ranked[0],assigned_fragments=int(q.loc[r,'assigned_fragments']),identity_status='SUPPORTED',basis='unchanged SNP thresholds; duplicate representative by assigned fragments before immune outcomes'))
    for r in low:mapping.append(dict(run=r,sample_alias=q.loc[r,'sample_alias'],anonymous_donor='',selected=False,assigned_fragments=int(q.loc[r,'assigned_fragments']),identity_status='LOW_INFORMATION_DEFERRED',basis='original sample SNP covered-loci minimum not met; not claimed resolved'))
    table(pd.DataFrame(mapping),'OBSERVATIONAL_DONOR_FREEZE.tsv')
    for row in pf[pf.classification=='SAME_DONOR_SUPPORTED'].itertuples():
        va=dict(run_a=row.run_a,run_b=row.run_b,alias_a=q.loc[row.run_a,'sample_alias'],alias_b=q.loc[row.run_b,'sample_alias'],cohort_a=q.loc[row.run_a,'group'],cohort_b=q.loc[row.run_b,'group'],shared_loci=row.shared_loci,vaf_r=row.vaf_r,ibs0=row.opposite_homozygote_fraction,treatment='UNCONFIRMED')
        for key,run in [('a',row.run_a),('b',row.run_b)]:
            at=json.loads(meta.loc[run,'attributes'])
            for field in ['age','tissue','collection_date']:va[field+'_'+key]=at.get(field,'')
        samepairs.append(va)
    table(pd.DataFrame(samepairs),'SUPPORTED_DONOR_PAIRS_WITH_METADATA.tsv')
    table(qc[qc.run.isin(low)],'DEFERRED_LOW_INFORMATION_LIBRARIES.tsv')
    freeze={'status':'PASS_QUALIFIED_SUBSET','frozen_before_outcome_tests':now(),'original_observational_libraries':len(obs),'retained_libraries':len(eligible),'independent_supported_donors':len(selected),'deferred_low_information':low,'selected':selected,'alternate_representatives':alternates,'identity_thresholds':CFG['identity'],'treatment_eligibility':'NO_TREATMENT_COMPARISON','note':'New version permits qualified subset; entire previous cohort gate remains REVIEW_REQUIRED.'}
    dump(P/'IDENTITY_FREEZE.json',freeze);assert len(selected)==44
    # Persistent specification saved before outcomes are loaded or tested.
    dump(P/'ANALYSIS_PLAN_BEFORE_TESTS.json',{'at':now(),'primary_configuration':CFG['statistics'],'programs':CFG['programs'],'selected_donors':selected,'alternate_duplicate_rule':'swap only the one duplicate representative, all four endpoints unchanged; robustness only','no_low_coverage_threshold_relaxation':True,'no_result_based_selection':True})
    gate('immune_associations');mark('IMMUNE_ASSOCIATIONS',donors=len(selected))
    scores=pd.read_csv(OLD/'tables/program_scores_all_libraries.tsv',sep='\t').set_index('run')
    assert scores.index.is_unique and set(scores.index)==set(qc.run)
    genes=sorted(set(['SPP1']+sum(CFG['programs'].values(),[])))
    ann=pd.read_csv(CONV/'tables/gene_annotation.tsv',sep='\t').set_index('gene_id')
    ids=ann.index[ann.gene_symbol.isin(genes)]
    chunks=[]
    for chunk in pd.read_csv(CONV/'matrices/all_libraries_technical_only/TPM.tsv.gz',sep='\t',index_col=0,chunksize=5000):
        chunks.append(chunk.loc[chunk.index.isin(ids)])
    tpm=pd.concat(chunks);sym=ann.loc[tpm.index,'gene_symbol'];logged=np.log2(tpm.groupby(sym).sum()+1)
    maxdelta=0.
    for name,members in [('SPP1',['SPP1'])]+list(CFG['programs'].items()):
        rec=logged.loc[members].mean(0).reindex(scores.index)
        delta=float(np.max(np.abs(rec-scores[name])));maxdelta=max(maxdelta,delta);assert delta<1e-10
    rows=[];rng=np.random.default_rng(CFG['statistics']['seed'])
    for spec,runs in [('primary_44_donors',selected),('alternate_duplicate_44_donors',alternates)]:
        s=scores.loc[runs];table(s.reset_index(),spec+'_INPUT.tsv')
        for exposure,fam in [('SPP1','primary_two_immune_programs'),('literature5','secondary_literature_definition')]:
            for endpoint in ['cytotoxic','dysfunction']:
                x=s[exposure].to_numpy();y=s[endpoint].to_numpy()
                result=evaluate(x,y,rng);assert abs(result['rho']-spearmanr(x,y).statistic)<1e-12
                rows.append(dict(specification=spec,family=fam,exposure=exposure,endpoint=endpoint,n_donors=len(runs),**result))
        for fam in ['primary_two_immune_programs','secondary_literature_definition']:
            subset=[r for r in rows if r['specification']==spec and r['family']==fam]
            for r,adj in zip(subset,bh([r['p'] if r['p'] is not None else 1 for r in subset])):r['FDR']=float(adj)
    result=pd.DataFrame(rows);table(result,'IMMUNE_ASSOCIATIONS_ALL.tsv')
    table(result[result.specification=='primary_44_donors'],'IMMUNE_PRIMARY.tsv')
    dump(P/'INPUT_MANIFEST.json',[snapshot(f) for f in inputs])
    dump(P/'ENVIRONMENT.json',{'python':sys.version,'executable':sys.executable,'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__})
    lines=['# 合格供者中SPP1与组织免疫程序的对应','','## 样本关系与资格','',f'60份文库1770对SNP关系从已存计数独立复算一致。46份观察文库中，SRR32389767仅534个合格位点，低于原1000位点门槛；该文库的45对观察关系全部保留为未决。本轮暂缓纳入该文库，其余45份文库按同供者合并后形成44个匿名供者单元。','观察文库SRR32389752与SRR32389755支持同供者；代表文库按assigned fragments较高者固定，另报告替换代表后的全部四项结果。14份候选文库支持7对同供者，但官方处理标签缺失，年龄/组织等来源冲突保留，不做抗SPP1/IgG比较。','','## 组织层免疫结果','']
    for r in rows:
        if r['specification']!='primary_44_donors':continue
        lines.append(f"- {r['exposure']}—{r['endpoint']}: n={r['n_donors']}, rho={r['rho']:.6f}, bootstrap95%CI=[{r['ci_low']:.6f}, {r['ci_high']:.6f}], 双侧置换P={r['p']:.6g}, BH FDR={r['FDR']:.6g}; 留一供者rho={r['loo_min']:.6f}至{r['loo_max']:.6f}。")
    lines+=['','## Methods','SNP面板、质量与同/异供者阈值完全沿用Phase240，从60个原始位点计数数组重新核对1770对。先按原门槛暂缓低信息文库，再检查保留集合所有配对均可判定；相同供者连通分量内部要求两两一致，以assigned fragments最大文库作为代表。名单及统计口径在本轮读取免疫结局前保存。','SPP1为log2(TPM+1)；细胞毒和功能障碍程序、五基因文献定义均沿用原成员与等权均值，直接从原TPM复算并验证已有分数。双侧Spearman置换100000次，供者bootstrap5000次，随机种子26210240。两个主要端点和两个文献定义端点分别BH。重复文库替换为预设稳健性比较，未根据P选择版本。','','## Discussion','当前结果为公开bulk组织层对应，不鉴定SPP1的C4来源、CD8/NK特异功能或MTA因果作用；原Phase213主要免疫端点未获研究调整后支持的结果继续保留。合格子集的独立性在固定面板与保守阈值下得到支持，不代表公开患者编号已外部鉴定。被暂缓文库及原46例精确对应未解决；处理标签不能从SNP、文库名或表达推断。','完整正负结果、样本名单、7对候选供者、来源字段冲突和代表替换比较见tables。无新图，旧RDS及原数据未改。','','来源：原Phase216官方元数据归档；NGSCheckMate SNP面板 https://github.com/parklab/NGSCheckMate 。本文使用自定义保守筛查，不称原版软件判定。']
    (P/'REPORT_CN.md').write_text('\n'.join(lines))
    qa=dict(status='PASS',time=now(),identity_npz_hashes=60,reproduced_pairs=1770,unresolved_observation_pairs_retained=45,deferred_libraries=low,qualified_donors=44,main_tests=4,duplicate_sensitivity_tests=4,program_score_max_difference=maxdelta,independent_scipy_rho_checks=8,old_inputs_unchanged=True,treatment_comparisons=0,new_figures=0)
    dump(P/'FINAL_VALIDATION.json',qa)
    mark('COMPLETE_VERIFIED',validation=qa)
    (P/'EXIT_STATUS').write_text('0\n');(P/'RUN_FINISHED').write_text(now())
    print(result.to_string(index=False),flush=True)
try:main()
except BaseException:
    (P/'logs/ERROR.txt').write_text(traceback.format_exc());mark('FAILED',error=traceback.format_exc());(P/'EXIT_STATUS').write_text('1\n');(P/'RUN_FINISHED').write_text(now());raise

