from pathlib import Path
import shutil,json,hashlib
R=Path(__file__).resolve().parents[2];P=Path(__file__).resolve().parents[1];O=R/'phase329_immune_spatial_clinical_A4_20261006_v2';S=R/'phase243_donor_verified_immune_20260924_v1/tables'
for sub in ['data','assets']:
 for f in (O/sub).iterdir():
  if f.is_file():shutil.copy2(f,P/sub/f.name)
manifest=json.loads((O/'qc/SOURCE_MANIFEST.json').read_text())
for name in ['primary_44_donors_INPUT.tsv','IMMUNE_PRIMARY.tsv']:
 f=S/name;shutil.copy2(f,P/'data'/name);manifest.append({'path':str(f),'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
(P/'qc/SOURCE_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
s=(O/'code/draw.py').read_text();a=s.index('# F:');b=s.index('# G-H:')
new='''# F: frozen tissue-level correlations, one point per qualified donor.
panel('F',116,197,'SPP1 and tissue immune programs','Bulk RNA | 44 independent donors')
bi=read('primary_44_donors_INPUT.tsv');bs=read('IMMUNE_PRIMARY.tsv')
assert len(bi)==44 and bi.run.is_unique
for x,endpoint,title,col in [(127,'cytotoxic','Cytotoxicity',blue),(171,'dysfunction','Dysfunction',purple)]:
 st=bs[(bs.exposure=='SPP1')&(bs.endpoint==endpoint)].iloc[0]
 f=axmm(x,211,28,18);style(f)
 f.scatter(bi.SPP1,bi[endpoint],s=12,color=col,alpha=.78,edgecolors='white',linewidths=.3)
 f.set(xlim=(0,14),xticks=[0,6,12],xlabel='SPP1 log₂(TPM+1)')
 f.set_ylabel('Program score' if endpoint=='cytotoxic' else '')
 txt(x+14,206,title+'  ρ = '+f'{st.rho:.2f}',6.3,ha='center')
 txt(x+14,239,f'P={st.p:.1g}; FDR={st.FDR:.1g}',6.3,ha='center')
 assert bi.SPP1.between(0,14).all()
bi.to_csv(P/'data/PANEL_F_44_DONORS.tsv',sep='\\t',index=False)
'''
s=s[:a]+new+s[b:];s=s.replace('26210_immune_spatial_clinical_A4_v2','26210_immune_spatial_clinical_A4_v3');(P/'code/draw.py').write_text(s)
cap=(O/'delivery/26210_immune_spatial_clinical_A4_v2_图注.txt').read_text(encoding='utf-8-sig');a=cap.index('F｜');b=cap.index('G｜')
newcap='''F｜SPP1与组织免疫程序的共同变化。替换v2公开CRC干预面板，采用PRJNA1223824公开组织bulk既有44名合格匿名供者。按冻结SNP筛查规则排除低信息文库并合并同供者重复，每个点对应一位独立供者；全部44人均展示，不截取高/低组、不删离群点。横轴为SPP1 log₂(TPM+1)，纵轴为程序成员的log₂(TPM+1)等权均值，不按供者再标准化。细胞毒程序由NKG7、GNLY、GZMB、PRF1、CTSW组成；功能障碍相关程序由PDCD1、CTLA4、LAG3、HAVCR2、TIGIT、TOX、CXCL13、ENTPD1组成。散点不拟合新趋势线，图中统计均直接沿用Phase243冻结汇总。\nSPP1–细胞毒：Spearman ρ=0.670754，供者bootstrap 95%CI 0.455353–0.814443，双侧100000次置换P=0.000009999900001，BH FDR=0.000019999800002。\nSPP1–功能障碍：ρ=0.491896，95%CI 0.220634–0.705998，P=0.000789992100079，FDR=0.000789992100079。两项为预定义primary_two_immune_programs校正家族；5000次bootstrap，原统计种子26210240，未在本次重算检验或区间。\n两类组织程序均随SPP1升高而升高，提示组织层免疫转录共同变化；不能解释为细胞毒功能降低，也不能单独确定SPP1的巨噬来源或CD8特异耗竭。原细胞特异接收端未形成统一免疫抑制联系的结果保留在完整文章中。原v2的外部Spp1干预结果未删除，仍保留于旧图、正文及源表。\n\n'''
cap=cap[:a]+newcap+cap[b:];cap=cap.replace('（A4综合图 v2）','（A4综合图 v3）');cap='本次修改：仅替换F及对应图注；A–E、G–H保持v2。\n\n'+cap
cap='\n'.join(l for l in cap.splitlines() if not l.startswith('F / Ext Fig'))
(P/'delivery/26210_immune_spatial_clinical_A4_v3_图注.txt').write_text(cap,encoding='utf-8-sig')
print('REVISION_READY')
