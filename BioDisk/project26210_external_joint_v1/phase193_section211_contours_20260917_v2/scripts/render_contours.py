from pathlib import Path
import json,subprocess,hashlib
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1')
old=r/'phase192_sections28_212_figures_20260917_v1'
exec((old/'scripts/render_figures.py').read_text().split('# 01 actual')[0])
d=r/'phase193_section211_contours_20260917_v2'
for sub in ['figures','qc','tables','scripts']:(d/sub).mkdir(parents=True,exist_ok=True)
from matplotlib.lines import Line2D
cmap=LinearSegmentedColormap.from_list('contourcontrast',['#DBEAF1','#92BBC9','#F8F5E9','#E7B47C','#B95534'])
fig,axs=plt.subplots(2,2,figsize=(7.8,6.65))
title(fig,'2.11','Functional expression contours across macrophage states')
fields=['expr:SPP1','score:Cholesterol metabolism','score:Lysosome','score:Antigen processing and presentation']
audit=[];matrices={}
c4h,ce,_=density(x[m.candidate=='4'],y[m.candidate=='4'],bounds=bounds)
for ax,f,letter in zip(axs.ravel(),fields,'ABCD'):
 h,e,_=density(x,y,weights=m[f].values,bounds=bounds)
 v=np.ma.masked_where(base<base.max()*.008,np.divide(h,base,out=np.zeros_like(h),where=base>0))
 lo,hi=np.percentile(v.compressed(),[2,98]);levels=np.linspace(lo,hi,8)
 xx=np.linspace(e[0],e[1],v.shape[1]);yy=np.linspace(e[2],e[3],v.shape[0])
 im=ax.contourf(xx,yy,v,levels=levels,cmap=cmap,extend='both',antialiased=True)
 ax.contour(xx,yy,v,levels=levels[1:-1],colors='#5C6566',linewidths=.32,alpha=.65)
 ax.contour(c4h,levels=[c4h.max()*.12],extent=ce,colors=[purple],linewidths=1.5,linestyles='--',zorder=5)
 clean(ax)
 label=f.split(':')[1].replace('Antigen processing and presentation','Antigen presentation')
 ax.set_title(letter+'   '+label,loc='left',fontsize=9)
 cb=fig.colorbar(im,ax=ax,fraction=.033,pad=.025,shrink=.79,ticks=levels[::2])
 from matplotlib.ticker import FormatStrFormatter
 cb.ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'));cb.ax.tick_params(labelsize=6)
 cb.ax.set_title('Local mean\nlog1p(CP10K)',fontsize=6,pad=7)
 audit.append({'feature':f,'lower_display_p2':lo,'upper_display_p98':hi,'levels':levels.tolist(),'visible_min':float(v.min()),'visible_max':float(v.max())})
 matrices[f]=v.filled(np.nan)
fig.legend(handles=[Line2D([0],[0],color=purple,lw=1.5,ls='--',label='C4 boundary'),Line2D([0],[0],color='#5C6566',lw=.65,label='Equal-expression contour')],loc='lower center',bbox_to_anchor=(.5,.055),ncol=2,frameon=False,fontsize=8)
fig.text(.5,.027,'Cool to warm: low to high within each panel; contour intervals are equally spaced.',ha='center',fontsize=7,color='#555555')
fig.subplots_adjust(top=.88,left=.025,right=.94,bottom=.14,wspace=.19,hspace=.24)
name='06_C4_expression_contours_v2'
save(fig,name)
np.savez_compressed(d/'tables/UNCHANGED_LOCAL_MEANS.npz',**matrices)
(d/'tables/CONTOUR_LEVELS.json').write_text(json.dumps(audit,indent=2))
subprocess.run(['pdftoppm','-r','150','-singlefile','-jpeg',str(d/'figures'/f'{name}.pdf'),str(d/'qc/PDF_RENDER')],check=True)
(d/'FIGURE_LEGEND_CN.md').write_text("""# 2.11 功能表达等高线图（替换原06图）

同一冻结UMAP中展示SPP1、胆固醇代谢、溶酶体及抗原呈递的局部平均表达。灰色细线为等表达线；紫色虚线为C4密度轮廓。填色由冷色至暖色表示每个面板内部从低至高的表达，色标仍为实际log1p(CP10K)单位，而非百分位或z分数。功能程序沿用排除TYMP与SPP1后的平均基因表达评分。

本次保留原图全部22,444个成熟巨噬细胞的冻结坐标、表达值、180×180网格、标准差2网格的高斯平滑，以及低于密度峰值0.8%的遮蔽条件。仅更改显示方式：各面板以可见网格值第2至98百分位作为主显示区间，划分7个等距色阶，并叠加等表达细线；两端超出部分使用端色，色条三角表示延伸。该设置增强局部差异的可读性，不改变细胞数据、局部均值或既有统计结论。面板间色标独立，不以相同颜色直接比较不同程序的绝对强弱。

C4轮廓由该群自身密度峰值12%的等高线确定，与原图相同。图中UMAP是转录状态嵌入，不是组织空间图；SPP1在C4外也有表达。患者级功能比较仍以原报告配对结果为准。本版替换phase192的06_C4_expression_landscapes；05功能程序热图保持不变。
""",encoding='utf8')
print('CONTOUR_DONE',flush=True)
