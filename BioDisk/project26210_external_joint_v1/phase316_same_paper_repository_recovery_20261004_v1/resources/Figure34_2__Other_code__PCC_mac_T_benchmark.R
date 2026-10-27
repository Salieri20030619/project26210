cor.test.result <- cor.test(tumor_prop_avg_UCASpatial_without_ent$`C1QC+ Mac`, tumor_prop_avg_UCASpatial_without_ent$`CD8+ Tem`)
p1<-ggplot(tumor_prop_avg_UCASpatial_without_ent,aes(x=`C1QC+ Mac`,y=`CD8+ Tem`,color=subTME_category))+geom_point()+
  geom_smooth(method = 'lm', se = F, color = 'black')+annotate("text", x = -Inf, y = Inf, 
                                                               label = paste("R = ", round(cor.test.result$estimate, digits=2),
                                                                             "\np = ", format.pval(cor.test.result$p.value, digits=2)),
                                                               hjust = -0.1, vjust = 2)+
  scale_color_manual(values = brewer.pal(10,"Paired")[c(2,6,3,4)])+guides(color=F)


cor.test.result <- cor.test(tumor_prop_avg_UCASpatial_without_ent$`C1QC+ Mac`, tumor_prop_avg_UCASpatial_without_ent$`CD8+ Tex`)
p2<-ggplot(tumor_prop_avg_UCASpatial_without_ent,aes(x=`C1QC+ Mac`,y=`CD8+ Tex`,color=subTME_category))+geom_point()+
  geom_smooth(method = 'lm', se = F, color = 'black')+annotate("text", x = -Inf, y = Inf, 
                                                               label = paste("R = ", round(cor.test.result$estimate, digits=2),
                                                                             "\np = ", format.pval(cor.test.result$p.value, digits=2)),
                                                               hjust = -0.1, vjust = 2)+
  scale_color_manual(values = brewer.pal(10,"Paired")[c(2,6,3,4)])+guides(color=F)
cor.test.result <- cor.test(tumor_prop_avg_UCASpatial_without_ent$`C1QC+ Mac`, tumor_prop_avg_UCASpatial_without_ent$`CD8+ Trm`)
p3<-ggplot(tumor_prop_avg_UCASpatial_without_ent,aes(x=`C1QC+ Mac`,y=`CD8+ Trm`,color=subTME_category))+geom_point()+
  geom_smooth(method = 'lm', se = F, color = 'black')+annotate("text", x = -Inf, y = Inf, 
                                                               label = paste("R = ", round(cor.test.result$estimate, digits=2),
                                                                             "\np = ", format.pval(cor.test.result$p.value, digits=2)),
                                                               hjust = -0.1, vjust = 2)+
  scale_color_manual(values = brewer.pal(10,"Paired")[c(2,6,3,4)])+guides(color=F)


cor.test.result <- cor.test(tumor_prop_avg_UCASpatial_without_ent$`FN1+ Mac`, tumor_prop_avg_UCASpatial_without_ent$`CD8+ Tem`)
p4<-ggplot(tumor_prop_avg_UCASpatial_without_ent,aes(x=`FN1+ Mac`,y=`CD8+ Tem`,color=subTME_category))+geom_point()+
  geom_smooth(method = 'lm', se = F, color = 'black')+annotate("text", x = -Inf, y = Inf, 
                                                               label = paste("R = ", round(cor.test.result$estimate, digits=2),
                                                                             "\np = ", format.pval(cor.test.result$p.value, digits=2)),
                                                               hjust = -0.1, vjust = 2)+
  scale_color_manual(values = brewer.pal(10,"Paired")[c(2,6,3,4)])+guides(color=F)

cor.test.result <- cor.test(tumor_prop_avg_UCASpatial_without_ent$`FN1+ Mac`, tumor_prop_avg_UCASpatial_without_ent$`CD8+ Tex`)
p5<-ggplot(tumor_prop_avg_UCASpatial_without_ent,aes(x=`FN1+ Mac`,y=`CD8+ Tex`,color=subTME_category))+geom_point()+
  geom_smooth(method = 'lm', se = F, color = 'black')+annotate("text", x = -Inf, y = Inf, 
                                                               label = paste("R = ", round(cor.test.result$estimate, digits=2),
                                                                             "\np = ", format.pval(cor.test.result$p.value, digits=2)),
                                                               hjust = -0.1, vjust = 2)+
  scale_color_manual(values = brewer.pal(10,"Paired")[c(2,6,3,4)])+guides(color=F)
cor.test.result <- cor.test(tumor_prop_avg_UCASpatial_without_ent$`FN1+ Mac`, tumor_prop_avg_UCASpatial_without_ent$`CD8+ Trm`)
p6<-ggplot(tumor_prop_avg_UCASpatial_without_ent,aes(x=`FN1+ Mac`,y=`CD8+ Trm`,color=subTME_category))+geom_point()+
  geom_smooth(method = 'lm', se = F, color = 'black')+annotate("text", x = -Inf, y = Inf, 
                                                               label = paste("R = ", round(cor.test.result$estimate, digits=2),
                                                                             "\np = ", format.pval(cor.test.result$p.value, digits=2)),
                                                               hjust = -0.1, vjust = 2)+
  scale_color_manual(values = brewer.pal(10,"Paired")[c(2,6,3,4)])+guides(color=F)


cor.test.result <- cor.test(tumor_prop_avg_UCASpatial_without_ent$`FN1+ Mac`, tumor_prop_avg_UCASpatial_without_ent$`CD8+ Tem`)
p7<-ggplot(tumor_prop_avg_UCASpatial_without_ent,aes(x=`PDPN+ Mac`,y=`CD8+ Tem`,color=subTME_category))+geom_point()+
  geom_smooth(method = 'lm', se = F, color = 'black')+annotate("text", x = -Inf, y = Inf, 
                                                               label = paste("R = ", round(cor.test.result$estimate, digits=2),
                                                                             "\np = ", format.pval(cor.test.result$p.value, digits=2)),
                                                               hjust = -0.1, vjust = 2)+
  scale_color_manual(values = brewer.pal(10,"Paired")[c(2,6,3,4)])+guides(color=F)

cor.test.result <- cor.test(tumor_prop_avg_UCASpatial_without_ent$`PDPN+ Mac`, tumor_prop_avg_UCASpatial_without_ent$`CD8+ Tex`)
p8<-ggplot(tumor_prop_avg_UCASpatial_without_ent,aes(x=`PDPN+ Mac`,y=`CD8+ Tex`,color=subTME_category))+geom_point()+
  geom_smooth(method = 'lm', se = F, color = 'black')+annotate("text", x = -Inf, y = Inf, 
                                                               label = paste("R = ", round(cor.test.result$estimate, digits=2),
                                                                             "\np = ", format.pval(cor.test.result$p.value, digits=2)),
                                                               hjust = -0.1, vjust = 2)+
  scale_color_manual(values = brewer.pal(10,"Paired")[c(2,6,3,4)])+guides(color=F)
cor.test.result <- cor.test(tumor_prop_avg_UCASpatial_without_ent$`PDPN+ Mac`, tumor_prop_avg_UCASpatial_without_ent$`CD8+ Trm`)
p9<-ggplot(tumor_prop_avg_UCASpatial_without_ent,aes(x=`PDPN+ Mac`,y=`CD8+ Trm`,color=subTME_category))+geom_point()+
  geom_smooth(method = 'lm', se = F, color = 'black')+annotate("text", x = -Inf, y = Inf, 
                                                               label = paste("R = ", round(cor.test.result$estimate, digits=2),
                                                                             "\np = ", format.pval(cor.test.result$p.value, digits=2)),
                                                               hjust = -0.1, vjust = 2)+
  scale_color_manual(values = brewer.pal(10,"Paired")[c(2,6,3,4)])+guides(color=F)


ggsave(plot = p1+p2+p3+p4+p5+p6+p7+p8+p9+plot_layout(ncol=3),filename = "/data/huangzr/Spatial/Kras/version5/SI_Fig/total_TME_benchmark/UCASpatial_without_ent_PCC.pdf",width = 12,height=12)
