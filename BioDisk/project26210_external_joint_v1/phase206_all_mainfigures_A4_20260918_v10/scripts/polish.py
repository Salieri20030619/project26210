from pathlib import Path
p=Path('phase206_all_mainfigures_A4_20260918_v10/scripts/assemble_pages.py');s=p.read_text()
s=s.replace("('A','Cell atlas',10,19,94,39)","('A','Cell atlas',10,19,94,33)").replace("('C','cNMF programs',10,69,94,39)","('C','cNMF programs',10,63,94,33)").replace("('D','Patient-program expression',10,119,94,77)","('D','Patient-program expression',10,107,94,90)")
s=s.replace("'Cell type','Malignant epithelial']","'Cell type','Malignant epithelial','query','Monocle trajectory on the original UMAP','11,950 malignant cells','224 original principal','log2 relative signal']")
s=s.replace("if any(s in txt for s in remove_patterns) or len(txt)>78:","if any(s in txt for s in remove_patterns) or len(txt)>78 or (fn==1 and letter=='D' and ('Original cNMF program' in txt or (txt.startswith('GSE') and '/' not in txt))) or (fn==1 and letter=='B' and txt.isdigit() and 1<=int(txt)<=22 and int(txt) not in [1,5,10,15,20,22]):")
s=s.replace("('GSE' in txt or 'RCC' in txt or 'SS_' in txt)","('/' in txt and ('GSE' in txt or 'RCC' in txt or 'SS_' in txt))")
start=s.index("  page.insert_text((x*MM,(y-2)*MM),letter")
end=s.index("  records.append",start)
s=s[:start]+s[end:]
s=s.replace(" out=d/'figures'/f'Main_Figure_{fn}_A4_v10.pdf'", " for letter,title,x,y,w,h in layouts[fn]:\n  page.insert_text((x*MM,(y-2)*MM),letter,fontsize=10,fontname='hebo');page.insert_text(((x+5)*MM,(y-2)*MM),title,fontsize=7.5,fontname='hebo')\n out=d/'figures'/f'Main_Figure_{fn}_A4_v10.pdf'")
p.write_text(s)
p=Path('phase206_all_mainfigures_A4_20260918_v10/scripts/redraw_panels.py');s=p.read_text().replace("[.22,.20,.65,.72]","[.20,.20,.62,.72]").replace("[.91,.40,.022,.38]","[.87,.40,.022,.38]")
s=s.replace("'Paired columns: GENIE3 / GRNBoost2'","'G3 / GB2'")
# Restore the unchanged regulon dendrogram, occupying the former spare left margin.
s=s.replace("fig=plt.figure(figsize=(90/25.4,101/25.4));ax=fig.add_axes([.025,.12,.72,.79])", "fig=plt.figure(figsize=(90/25.4,101/25.4));ax=fig.add_axes([.13,.12,.62,.79])")
s=s.replace("ca=fig.add_axes([.20,.955,.37,.018]);", "raw=(r/'phase130_section2_4_BE_correction_remote_v1/scripts/redraw_B_annotations.py').read_text();prefix=raw[:raw.index('plt.rcParams.update')];scope={'__file__':str(r/'phase130_section2_4_BE_correction_remote_v1/scripts/redraw_B_annotations.py')};exec(compile(prefix,'frozen_regulon_tree','exec'),scope);assert np.allclose(v.values,scope['Zall'][scope['leaves']]);ta=fig.add_axes([.005,.12,.105,.79]);scope['dendrogram'](scope['tree'],orientation='left',no_labels=True,color_threshold=0,above_threshold_color='#888888',ax=ta);ta.invert_yaxis();ta.axis('off')\nca=fig.add_axes([.20,.955,.37,.018]);")
p.write_text(s)
p=Path('phase206_all_mainfigures_A4_20260918_v10/scripts/redraw_density.py');s=p.read_text().replace("mark(ax);clean(ax);ax.set_title(['SPP1'", "ax.contour(np.linspace(ee[0],ee[1],v.shape[1]),np.linspace(ee[2],ee[3],v.shape[0]),v,levels=lev[1:-1],colors='#5C6566',linewidths=.25,alpha=.65);mark(ax);clean(ax);ax.set_title(['SPP1'")
p.write_text(s)
