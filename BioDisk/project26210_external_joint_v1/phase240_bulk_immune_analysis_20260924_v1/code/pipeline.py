import os,sys,json,pathlib,hashlib,datetime,subprocess,traceback,itertools,concurrent.futures,shutil,csv,time
import numpy as np
import pandas as pd
import pysam
P=pathlib.Path(__file__).resolve().parents[1]
ROOT=P.parent
CONV=ROOT/'phase236_public_bulk_conversion_20260923_v1'
CFG=json.loads((P/'CONFIG.json').read_text())
def now(): return datetime.datetime.now().astimezone().isoformat()
def dump(path,obj):
    path=pathlib.Path(path); tmp=path.with_suffix(path.suffix+'.tmp'); tmp.write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)); tmp.replace(path)
def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(2**20),b''): h.update(b)
    return h.hexdigest()
def status(stage,**kwargs):
    v={'status':stage,'updated_at':now(),**kwargs}; dump(P/'STATUS.json',v)
    with open(P/'logs/events.jsonl','a') as f:f.write(json.dumps(v)+'\n')
def rank(x):
    _,inv,c=np.unique(x,return_inverse=True,return_counts=True)
    return (np.cumsum(c)-(c-1)/2)[inv].astype(float)
def corr(x,y):
    x=x-x.mean();y=y-y.mean()
    den=np.linalg.norm(x)*np.linalg.norm(y)
    return float(x@y/den) if den else None
def spear(x,y): return corr(rank(x),rank(y))
def extract(task):
    run,bam,panel=task
    out=P/'identity'/f'{run}.npz'
    counts=np.zeros((len(panel),3),dtype=np.uint32)
    cfg=CFG['identity']
    with pysam.AlignmentFile(bam,'rb') as b:
        assert b.has_index()
        for i,(chrom,pos,ref,alt,rs) in enumerate(panel):
            for col in b.pileup(chrom,pos,pos+1,truncate=True,stepper='samtools',min_base_quality=cfg['baseq'],
                min_mapping_quality=cfg['mapq'],max_depth=cfg['max_depth'],ignore_overlaps=True,ignore_orphans=True,
                compute_baq=False,flag_filter=3844):
                if col.reference_pos!=pos: continue
                for pr in col.pileups:
                    if pr.is_del or pr.is_refskip or pr.query_position is None:continue
                    a=pr.alignment; q=pr.query_position
                    if a.query_qualities[q]<cfg['baseq']:continue
                    base=a.query_sequence[q]
                    counts[i,0 if base==ref else 1 if base==alt else 2]+=1
            if (i+1)%2000==0:
                dump(P/'identity'/f'{run}.progress.json',{'run':run,'loci_processed':i+1,'total':len(panel),'updated_at':now()})
    np.savez_compressed(out,counts=counts)
    good=(counts[:,:2].sum(1)>=cfg['depth'])&(counts[:,2]/np.maximum(counts.sum(1),1)<=cfg['other_allele_fraction_max'])
    result={'run':run,'panel_loci':len(panel),'eligible_loci':int(good.sum()),'status':'PASS' if good.sum()>=cfg['covered_loci_min'] else 'LOW_INFORMATION','sha256':sha(out),'finished_at':now()}
    dump(P/'identity'/f'{run}.QA.json',result)
    return result
def identity(qc):
    ref=list((CONV/'reference').glob('*.fa'))
    assert len(ref)==1,ref
    panel=[];excluded=[];seen=set()
    with pysam.FastaFile(str(ref[0])) as fa:
        for line in (P/'reference/SNP_GRCh38_hg38_wChr.bed').read_text().splitlines():
            c,s,e,rs,a,b=line.split()[:6];s=int(s);e=int(e)
            why=None
            if c not in ['chr'+str(i) for i in range(1,23)]:why='non_autosomal'
            elif e-s!=1 or a not in 'ACGT' or b not in 'ACGT' or a==b:why='non_biallelic_SNV'
            elif (c,s) in seen:why='duplicate_coordinate'
            else:
                base=fa.fetch(c,s,e).upper()
                if base not in [a,b]:
                    comp=str.maketrans('ACGT','TGCA')
                    a,b=a.translate(comp),b.translate(comp)
                if base not in [a,b]:why='reference_mismatch_after_strand_normalization'
                elif base==b:a,b=b,a
            if why:excluded.append([c,s,e,rs,a,b,why]);continue
            seen.add((c,s));panel.append((c,s,a,b,rs))
    panel.sort(key=lambda v:(int(v[0][3:]),v[1]))
    assert len(panel)>10000,len(panel)
    pd.DataFrame(panel,columns=['chrom','pos0','ref','alt','rsid']).to_csv(P/'reference/PANEL_USED.tsv',sep='\t',index=False)
    pd.DataFrame(excluded,columns=['chrom','start','end','rsid','ref','alt','reason']).to_csv(P/'reference/PANEL_EXCLUDED.tsv',sep='\t',index=False)
    tasks=[(row.run,row.BAM,panel) for row in qc.itertuples()]
    results=[]
    with concurrent.futures.ProcessPoolExecutor(max_workers=CFG['identity']['workers']) as ex:
        futures={ex.submit(extract,t):t[0] for t in tasks}
        for f in concurrent.futures.as_completed(futures):
            results.append(f.result())
            status('IDENTITY_RUNNING',completed_libraries=len(results),total_libraries=len(tasks),last_run=futures[f])
    pd.DataFrame(results).sort_values('run').to_csv(P/'tables/identity_sample_QC.tsv',sep='\t',index=False)
    arrays={r:np.load(P/'identity'/f'{r}.npz')['counts'] for r in qc.run}
    data={}
    for r,a in arrays.items():
        depth=a[:,:2].sum(1)
        valid=(depth>=CFG['identity']['depth'])&(a[:,2]/np.maximum(a.sum(1),1)<=CFG['identity']['other_allele_fraction_max'])
        data[r]=(a[:,1]/np.maximum(depth,1),valid)
    pairs=[]
    for a,b in itertools.combinations(qc.run,2):
        va,ma=data[a];vb,mb=data[b];mask=ma&mb;n=int(mask.sum());x=va[mask];y=vb[mask]
        r=corr(x,y) if n else None
        ibs=float(np.mean(((x<=.1)&(y>=.9))|((y<=.1)&(x>=.9)))) if n else None
        label='AMBIGUOUS'
        if n<CFG['identity']['pair_loci_min'] or min(ma.sum(),mb.sum())<CFG['identity']['covered_loci_min']:label='LOW_INFORMATION'
        elif r is not None and r>=CFG['identity']['same_r_min'] and ibs<=CFG['identity']['same_ibs0_max']:label='SAME_DONOR_SUPPORTED'
        elif r is not None and r<=CFG['identity']['distinct_r_max'] and ibs>=CFG['identity']['distinct_ibs0_min']:label='DISTINCT_SUPPORTED'
        pairs.append({'run_a':a,'run_b':b,'shared_loci':n,'vaf_r':r,'opposite_homozygote_fraction':ibs,'classification':label})
    pf=pd.DataFrame(pairs);pf.to_csv(P/'tables/identity_pairwise.tsv',sep='\t',index=False)
    obs=qc.loc[qc.group=='observational46_paper_correspondence','run'].tolist()
    sub=pf[pf.run_a.isin(obs)&pf.run_b.isin(obs)]
    unresolved=sub[~sub.classification.isin(['SAME_DONOR_SUPPORTED','DISTINCT_SUPPORTED'])]
    groups=[{x} for x in obs]
    for row in sub[sub.classification=='SAME_DONOR_SUPPORTED'].itertuples():
        joins=[g for g in groups if row.run_a in g or row.run_b in g]
        if len(joins)==2:groups=[g for g in groups if g not in joins]+[set.union(*joins)]
    inconsistent=[]
    for g in groups:
        for row in sub[sub.run_a.isin(g)&sub.run_b.isin(g)].itertuples():
            if row.classification!='SAME_DONOR_SUPPORTED':inconsistent.append([row.run_a,row.run_b])
    selected=[];mapping=[]
    q=qc.set_index('run')
    for i,g in enumerate(sorted(groups,key=lambda x:min(x)),1):
        chosen=sorted(g,key=lambda r:(-float(q.loc[r,'assigned_fragments']),r))[0]
        selected.append(chosen)
        for r in sorted(g):mapping.append({'run':r,'anonymous_unit':f'D{i:03d}','selected':r==chosen,'basis':'conservative_SNP_consistency_screen','sample_alias':q.loc[r,'sample_alias']})
    pd.DataFrame(mapping).to_csv(P/'tables/observational_units.tsv',sep='\t',index=False)
    gate={'status':'PASS' if len(unresolved)==0 and not inconsistent else 'REVIEW_REQUIRED','unresolved_observation_pairs':len(unresolved),'inconsistent_components':inconsistent,'selected_units':len(selected),'all_pair_classes':pf.classification.value_counts().to_dict(),'treatment_eligibility':'NO_TREATMENT_COMPARISON'}
    dump(P/'IDENTITY_GATE.json',gate)
    return gate,selected
def bh(v):
    v=np.array(v);order=np.argsort(v);adj=np.minimum.accumulate((v[order]*len(v)/np.arange(1,len(v)+1))[::-1])[::-1];out=np.empty(len(v));out[order]=np.minimum(adj,1);return out
def evaluate(x,y,rng):
    r=spear(x,y)
    if r is None:return {'rho':None,'p':None,'ci_low':None,'ci_high':None,'loo_min':None,'loo_max':None,'nonestimable':'constant_input'}
    a=rank(x);a=(a-a.mean())/np.linalg.norm(a-a.mean());b=rank(y);b=(b-b.mean())/np.linalg.norm(b-b.mean())
    exceed=0;total=CFG['statistics']['permutations']
    for start in range(0,total,1000):
        n=min(1000,total-start);z=np.array([rng.permutation(b) for _ in range(n)]);exceed+=int((np.abs(z@a)>=abs(r)-1e-12).sum())
    boots=[]
    for _ in range(CFG['statistics']['bootstrap']):
        ids=rng.integers(0,len(x),len(x));rr=spear(x[ids],y[ids])
        if rr is not None:boots.append(rr)
    loo=[spear(np.delete(x,i),np.delete(y,i)) for i in range(len(x))];loo=[v for v in loo if v is not None]
    lo,hi=np.quantile(boots,[.025,.975])
    return {'rho':r,'p':(exceed+1)/(total+1),'ci_low':float(lo),'ci_high':float(hi),'loo_min':min(loo),'loo_max':max(loo),'nonestimable':''}
def expression(qc,gate,selected):
    status('EXPRESSION_ANALYSIS',identity_gate=gate['status'])
    matrix=CONV/'matrices/all_libraries_technical_only/TPM.tsv.gz'
    tpm=pd.read_csv(matrix,sep='\t',index_col=0)
    annotation=pd.read_csv(CONV/'tables/gene_annotation.tsv',sep='\t').set_index('gene_id')
    assert tpm.index.is_unique and set(qc.run)==set(tpm.columns)
    required=sorted(set(['SPP1']+sum(list(CFG['programs'].values()),[])))
    symbol=annotation.reindex(tpm.index).gene_symbol.fillna('')
    membermap=annotation[annotation.gene_symbol.isin(required)].reset_index()
    membermap.to_csv(P/'tables/program_gene_row_map.tsv',sep='\t',index=False)
    grouped=tpm.loc[symbol.isin(required)].groupby(symbol[symbol.isin(required)]).sum()
    assert set(required).issubset(grouped.index),set(required)-set(grouped.index)
    logged=np.log2(grouped+1)
    sample=qc[['run','sample_alias','group']].copy().set_index('run')
    sample['SPP1']=logged.loc['SPP1']
    coverage=[]
    for name,genes in CFG['programs'].items():
        sample[name]=logged.loc[genes].mean(axis=0)
        for g in genes:
            for cohort,rs in qc.groupby('group'):
                v=grouped.loc[g,rs.run]
                coverage.append({'program':name,'gene':g,'cohort':cohort,'libraries':len(v),'zero_libraries':int((v==0).sum()),'TPM_min':float(v.min()),'TPM_median':float(v.median()),'TPM_max':float(v.max())})
    sample.to_csv(P/'tables/program_scores_all_libraries.tsv',sep='\t')
    pd.DataFrame(coverage).to_csv(P/'tables/program_coverage.tsv',sep='\t',index=False)
    if gate['status']!='PASS' or len(selected)<CFG['statistics']['independent_n_min']:
        dump(P/'ASSOCIATION_GATE.json',{'status':'NOT_STARTED_IDENTITY_REVIEW','reason':'donor-unit eligibility unresolved or fewer than 10 units','n_candidate_units':len(selected)})
        return []
    s=sample.loc[selected];s.to_csv(P/'tables/observational_analysis_input.tsv',sep='\t')
    rows=[];rng=np.random.default_rng(CFG['statistics']['seed'])
    for exposure,family in [('SPP1','primary_two_immune_programs'),('literature5','secondary_literature_definition')]:
        for endpoint in ['cytotoxic','dysfunction']:
            result=evaluate(s[exposure].to_numpy(),s[endpoint].to_numpy(),rng)
            rows.append({'family':family,'exposure':exposure,'endpoint':endpoint,'n_units':len(s),**result})
    for family in set(r['family'] for r in rows):
        subset=[r for r in rows if r['family']==family]
        vals=[r['p'] if r['p'] is not None else 1. for r in subset]
        for r,q in zip(subset,bh(vals)):r['FDR']=float(q)
    pd.DataFrame(rows).to_csv(P/'tables/immune_associations.tsv',sep='\t',index=False)
    dump(P/'ASSOCIATION_GATE.json',{'status':'COMPLETE','n_units':len(s),'tests':len(rows),'family_sizes':[2,2],'bootstrap':CFG['statistics']['bootstrap'],'permutations':CFG['statistics']['permutations']})
    return rows
def records(note):
    stamp=datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    for name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md','RESULTS_MANIFEST.tsv']:
        shutil.copy2(ROOT/name,P/'root_record_backups'/f'{stamp}_{name}')
    for name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md']:
        f=ROOT/name;f.write_text(note+'\n\n'+f.read_text())
    manifest=ROOT/'RESULTS_MANIFEST.tsv'
    with open(manifest,'r') as f: header=f.readline().rstrip('\n').split('\t')
    record={h:'' for h in header}
    record.update({k:v for k,v in {'result_id':'phase240_bulk_immune_v1','artifact_type':'report','validation_status':'PASS','validation_evidence':str(P/'FINAL_VALIDATION.json')}.items() if k in record})
    for h in header:
        low=h.lower()
        if 'path' in low:record[h]=str(P/'REPORT_CN.md')
        elif low in ['phase','stage']:record[h]='Phase240'
        elif 'status' in low:record[h]='COMPLETE'
        elif any(x in low for x in ['description','note','summary']):record[h]='Bulk identity screen and conditional immune analysis; treatment identity unresolved'
        elif 'date' in low or 'time' in low:record[h]=now()
    with open(manifest,'a') as f:csv.DictWriter(f,fieldnames=header,delimiter='\t').writerow(record)
def main():
    status('INPUT_VALIDATION')
    qa=json.loads((CONV/'MATRIX_QA.json').read_text())
    assert qa['status']=='PASS' and qa['sample_libraries']==60
    assert (CONV/'EXIT_STATUS').read_text().strip()=='0'
    checked=[]
    for o in qa['outputs']:
        assert sha(o['path'])==o['sha256'],o['path']
        checked.append(o['path'])
    qc=pd.read_csv(CONV/'tables/sample_QC_summary.tsv',sep='\t')
    assert len(qc)==60 and qc.run.is_unique and (qc.status=='PASS').all()
    dump(P/'INPUT_VALIDATION.json',{'status':'PASS','matrix_hashes_verified':len(checked),'libraries':len(qc),'time':now()})
    gate,selected=identity(qc)
    rows=expression(qc,gate,selected)
    lines=['# Phase240：公开bulk样本一致性与组织免疫对应','',
        f'完成时间：{now()}。60份文库的9个既有矩阵SHA256全部匹配。',
        f'供者一致性筛查状态：{gate["status"]}；观察队列候选独立单元数：{gate["selected_units"]}。','']
    if rows:
        lines+=['## 组织层免疫对应','']
        for r in rows:
            lines.append(f'- {r["exposure"]}—{r["endpoint"]}：n={r["n_units"]}，rho={r["rho"]}，95%患者重采样区间[{r["ci_low"]}, {r["ci_high"]}]，双侧置换P={r["p"]}，家族BH FDR={r["FDR"]}；留一供者rho范围[{r["loo_min"]}, {r["loo_max"]}]。')
    else:lines+=['供者核对尚有未决项；已计算所有文库预设程序的描述性分数与覆盖，未开展观察队列推断检验。']
    lines+=['','## 方法',
        '供者核对采用固定NGSCheckMate GRCh38常见SNP面板的自定义保守VAF/对立纯合筛查，具体质量和判定阈值在首次计算前冻结于CONFIG.json。并非原版NGSCheckMate软件判定。相同供者连接必须形成两两一致分量；未决配对阻止供者层检验。',
        '表达为同symbol行TPM先相加再log2(TPM+1)。SPP1连续暴露为主；两个组织程序按固定成员等权平均。五基因文献定义另列次要家族。两个家族分别BH；双侧置换100000次，供者bootstrap5000次，随机种子固定。','',
        '## 证据整合与边界',
        '本次组织bulk结果与Phase213原患者层细胞特异分析分别报告；Phase213主要免疫端点未获研究调整后支持仍保留。bulk程序不能指定CD8/NK来源，也不证明C4来源或MTA因果链。',
        '14份候选仍缺少官方处理及配对映射，不能按名称、SPP1表达或聚类猜anti-SPP1/IgG。供者一致性只核对样本间生物学重复，不鉴定处理。组织、年龄与取材年份冲突仍按Phase216原记录保留。',
        'SNP结果为在所检面板和质量阈值下的样本一致性证据，未作外部身份鉴定。独立性及与原论文的对应仍需结合来源记录。临床字段和可靠无共享基因的组成估计未补齐，未进行临床/组成调整。',
        '原C4、接收端、前两部分和所有旧RDS保持不变；本轮先交付文字统计，不生成新图。','',
        '来源：https://github.com/parklab/NGSCheckMate；https://academic.oup.com/nar/article/45/11/e103/3079509',
        '完整逐文库、逐配对及正负统计见tables/；参考和代码版本见reference/、CODE_SHA256.json。']
    (P/'REPORT_CN.md').write_text('\n'.join(lines))
    assert len(pd.read_csv(P/'tables/identity_pairwise.tsv',sep='\t'))==1770
    assert len(pd.read_csv(P/'tables/program_scores_all_libraries.tsv',sep='\t'))==60
    validation={'status':'PASS','input_matrix_hashes':9,'identity_pairs':1770,'sample_scores':60,'identity_gate':gate['status'],'association_tests':len(rows),'treatment_comparison_run':False,'new_figures':False,'completed_at':now()}
    dump(P/'FINAL_VALIDATION.json',validation)
    (P/'EXIT_STATUS').write_text('0\n');(P/'RUN_FINISHED').write_text(now())
    status('COMPLETE_WITH_IDENTITY_REVIEW' if not rows else 'COMPLETE_VERIFIED',validation=validation)
    note=f'## {now()} Phase240 bulk身份及免疫分析完成\n60文库输入完整性通过，1770对样本一致性检查完成；观察供者资格{gate["status"]}，免疫检验{len(rows)}项。14候选处理未知，不执行干预比较。原Phase213阴性保留；无新图、旧RDS不变。\n报告：{P}/REPORT_CN.md；FINAL_VALIDATION PASS，EXIT_STATUS 0。'
    records(note)
if __name__=='__main__':
    try:main()
    except Exception:
        (P/'logs/ERROR.txt').write_text(traceback.format_exc())
        (P/'EXIT_STATUS').write_text('1\n')
        (P/'RUN_FINISHED').write_text(now())
        status('FAILED',error=traceback.format_exc())
        raise
