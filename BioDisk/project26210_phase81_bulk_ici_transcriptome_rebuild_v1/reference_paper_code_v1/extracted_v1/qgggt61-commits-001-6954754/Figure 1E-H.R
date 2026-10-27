###########Figure 1E-H######################
###########Figure 1E-H######################
library(tidyverse)
library(jjPlot)
library(ggnewscale)
dff <- df %>% rownames_to_column(var="id") %>% as_tibble() %>% 
  mutate(EGFR=as.character(EGFR),
         Recurrence=as.character(Recurrence),
         Smoking_status=as.character(Smoking_status),
         Sample_no=as.character(Sample_no)) %>% pivot_longer(-id)

dff$name <- factor(dff$name,levels =dff$name %>% unique() %>% rev())
plot <- dff %>% 
  filter(!name %in% c("Pretreatment_detectio2","Postoperative_detection2")) %>% 
  ggplot(aes(id,name))+
  geom_tile(color="white",fill="grey90",size=1)+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  geom_text(data=dff %>% filter(name=="Stage"),
            aes(id,name,label = value),inherit.aes = F)+
  geom_text(data=dff %>% filter(name=="Sample_no"),
            aes(id,name,label = value),inherit.aes = F)+
  geom_point(data=dff %>% filter(name %in% c("Smoking_status","EGFR","Recurrence")) %>% 
               filter(value !="1"),
             aes(id,name,color = value),inherit.aes = F,shape=19,size=5,show.legend = F)+
  scale_color_manual(values = "black")+
  geom_point(data=dff %>% filter(name %in% c("Treatment")) %>% 
               dplyr::rename("Treatment"="value"),
             aes(id,name,fill = Treatment),inherit.aes = F,shape=22,size=11,color="grey90")+
  scale_fill_manual(values = c("#F98400","#9986A5"))+
  new_scale_fill()+
  geom_point(data=dff %>% filter(name %in% c("Histology")) %>% 
               dplyr::rename("Histology"="value"),
             aes(id,name,fill = Histology),inherit.aes = F,shape=22,size=11,color="grey90")+
  scale_fill_manual(values = c("#E6A0C4","#6C8645","#CCC591"))+
  new_scale_fill()+
  geom_jjtriangle(data=dff %>% filter(name %in% c("Pretreatment_detectio","Pretreatment_detectio2")) %>% 
                    pivot_wider(names_from = "name") %>% 
                    mutate(name="Pretreatment_detectio") %>% 
                    filter(Pretreatment_detectio !="NA") %>% 
                    dplyr::rename("Detection"="Pretreatment_detectio"),
                  aes(id,name,fill =Detection),type = 'ul',color="white",size=0.5)+
  geom_jjtriangle(data=dff %>% filter(name %in% c("Pretreatment_detectio","Pretreatment_detectio2")) %>% 
                    pivot_wider(names_from = "name") %>% 
                    mutate(name="Pretreatment_detectio") %>% 
                    filter(Pretreatment_detectio !="NA"),
                  aes(id,name,fill =Pretreatment_detectio2),type = 'br',
                  color="white",size=0.5)+
  geom_jjtriangle(data=dff %>% filter(name %in% c("Postoperative_detection",
                                                  "Postoperative_detection2")) %>% 
                    pivot_wider(names_from = "name") %>% 
                    mutate(name="Postoperative_detection") %>% 
                    filter(Postoperative_detection !="NA"),
                  aes(id,name,fill =Postoperative_detection),type = 'ul',
                  color="white",size=0.5)+
  geom_jjtriangle(data=dff %>% filter(name %in% c("Postoperative_detection",
                                                  "Postoperative_detection2")) %>% 
                    pivot_wider(names_from = "name") %>% 
                    mutate(name="Postoperative_detection") %>% 
                    filter(Postoperative_detection !="NA"),
                  aes(id,name,fill =Postoperative_detection2),type = 'br',
                  color="white",size=0.5)+
  scale_fill_manual(values = c("#7294D4","#0A9F9D"))+
  guides(fill=guide_legend(theme=theme(legend.direction = "horizontal",
                                       legend.title.position = "top",
                                       legend.key.height = unit(0.8,"cm"),
                                       legend.key.width = unit(0.8,"cm"))
  ))+
  theme(legend.background = element_blank(),
        legend.key = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x=element_blank(),
        axis.title = element_blank(),
        axis.text.y=element_text(color="black",size=9,face="bold"),
        legend.key.spacing.y  = unit(-0.3,"cm"),
        legend.text = element_text(color="black",size=9,
                                   margin = margin(c(l=0),unit="cm")),
        plot.margin = margin(c(0.5,0.2,0.5,0.2), unit="cm"),
        legend.spacing.y = unit(-0.3,"cm"),
        legend.title = element_text(color="black",size=9,face="bold",
                                    margin = margin(c(t=0,b=0.1),unit="cm")))
ggsave(plot,file="heatmap.pdf",width=9.59,height = 3.69,unit="in",dpi=300)


