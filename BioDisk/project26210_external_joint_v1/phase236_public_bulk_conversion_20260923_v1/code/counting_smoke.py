from common import *
import tempfile
T=P/'software/counting_smoke';T.mkdir(exist_ok=True)
if (T/'QA.json').exists() and json.loads((T/'QA.json').read_text()).get('status')=='PASS':
 print('Existing counting smoke PASS reused');raise SystemExit(0)
gtf=T/'genes.gtf';gtf.write_text('chrT\ttest\texon\t101\t500\t.\t+\t.\tgene_id "plus";\nchrT\ttest\texon\t1001\t1400\t.\t-\t.\tgene_id "minus";\n')
rows=['@HD\tVN:1.6\tSO:queryname','@SQ\tSN:chrT\tLN:2000']
for name,pos,flag1,flag2,mate in [('fwd',121,99,147,221),('rev',1021,99,147,1121)]:
 for flag,start,nextpos,tlen in [(flag1,pos,mate,150),(flag2,mate,pos,-150)]:rows.append('\t'.join(map(str,[name,flag,'chrT',start,60,'50M','=',nextpos,tlen,'A'*50,'I'*50,'NH:i:1'])))
sam=T/'paired.sam';sam.write_text('\n'.join(rows)+'\n');observed={}
for mode in ['no','yes','reverse']:
 out=T/(mode+'.tsv')
 with out.open('w') as f:run([str(P/'software/venv/bin/htseq-count'),'-f','sam','-r','name','-s',mode,'-a','10','-t','exon','-i','gene_id','-m','union','--nonunique','none',str(sam),str(gtf)],T/(mode+'.log'),stdout=f)
 values={l.split('\t')[0]:int(l.split('\t')[1]) for l in out.read_text().splitlines()};observed[mode]=values
assert observed['no']['plus']==observed['no']['minus']==1
assert observed['yes']['plus']==1 and observed['yes']['minus']==0
assert observed['reverse']['plus']==0 and observed['reverse']['minus']==1
save(T/'QA.json',dict(status='PASS',paired_fragments_not_reads=True,stranded_orientation_mapping_confirmed=True,observed=observed))
