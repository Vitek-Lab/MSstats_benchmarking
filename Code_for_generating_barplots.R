####File for Graphical representation of bar plots

##########################################  WITH CONVERTER  #######################################################


plot_info3 = data.frame("Dataset" = rep("Dataset 1: Controlled Mixture - DDA - Skyline", 8),
                        "Version" = rep(c("MSstats","DeqMS","MSqRob2","pmartR","DEP","ProDA","Prolfqua","limma"),2),
                        "Stat" = c("TP", "TP", "TP", "TP","TP","TP","TP","TP", "FP", "FP", "FP", "FP","FP","FP","FP","FP"),
                        "Value" = c(-173,-167,-154,-154,-164,-152,-169,-167,92,140,101,170,250,178,179,142))


## Method Comp
plot_info3_c<-plot_info3 %>% mutate(Version = factor(Version, levels=c("MSstats","DeqMS","MSqRob2","pmartR","DEP","ProDA","Prolfqua","limma"))) %>% 
  ggplot(aes(x=Version, y=Value, fill=Version)) + 
  geom_bar(stat="identity", position="identity") +
  scale_fill_manual(values=c("violet","#33FFFF","yellow","green","purple","magenta","orange","blue"))+
  geom_hline(yintercept = 0) +ylim(-200,250)+
  labs(y = "TP                                           FP",
       x="", title="Dataset 1: Controlled Mixture - DDA - Skyline") + theme(axis.text.y=element_blank(),
                                                                            axis.title.y=element_text(size=10, face='bold'),
                                                                            axis.text.x=element_text(size=10, face='bold',angle = 90, vjust = 0.5, hjust=1))+
  annotate("text", x = c(1), y=-24.5, label = c("173"), fontface =2, size=8) +
  annotate("text", x = c(1), y=30, label = c("92"), fontface =2, size=8) +
  annotate("text", x = c(2), y=-22.5, label = c("167"), fontface =2, size=8) + 
  annotate("text", x = c(2), y=30, label = c("140"), fontface =2, size=8) +
  annotate("text", x = c(3), y=-15, label = c("154"), fontface =2, size=8) + 
  annotate("text", x = c(3), y=30, label = c("101"), fontface =2, size=8)+
  annotate("text", x = c(4), y=-24.5, label = c("154"), fontface =2, size=8) +
  annotate("text", x = c(4), y=30, label = c("170"), fontface =2, size=8) +
  annotate("text", x = c(5), y=-22.5, label = c("164"), fontface =2, size=8) + 
  annotate("text", x = c(5), y=30, label = c("304"), fontface =2, size=8) +
  annotate("text", x = c(6), y=-22.5, label = c("152"), fontface =2, size=8) + 
  annotate("text", x = c(6), y=30, label = c("178"), fontface =2, size=8) +
  annotate("text", x = c(7), y=-22.5, label = c("169"), fontface =2, size=8) + 
  annotate("text", x = c(7), y=30, label = c("179"), fontface =2, size=8) +
  annotate("text", x = c(8), y=-22.5, label = c("167"), fontface =2, size=8) + 
  annotate("text", x = c(8), y=30, label = c("142"), fontface =2, size=8) 


########################################################################################################


plot_info4 = data.frame("Dataset" = rep("Dataset 1: Controlled Mixture - DDA - Progenesis", 8),
                        "Version" = rep(c("MSstats","DeqMS","MSqRob2","pmartR","DEP","ProDA","Prolfqua","limma"),2),
                        "Stat" = c("TP", "TP", "TP", "TP","TP","TP", "TP","TP","FP", "FP", "FP", "FP","FP","FP","FP","FP"),
                        "Value" = c(-193,-199,-181,-203,-194,-191,-204,-202,50,112,131,450,278,214,188,111))


## Method Comp
plot_info4_c<-plot_info4 %>% mutate(Version = factor(Version, levels=c("MSstats","DeqMS","MSqRob2","pmartR","DEP","ProDA","Prolfqua","limma"))) %>% 
  ggplot(aes(x=Version, y=Value, fill=Version)) + 
  geom_bar(stat="identity", position="identity") +
  scale_fill_manual(values=c("violet","#33FFFF","yellow","green","purple","magenta","orange","blue"))+
  geom_hline(yintercept = 0) + ylim(-300,450)+
  labs(y = " TP                                         FP",
       x="", title="Dataset 1: Controlled Mixture - DDA - Progenesis") + theme(axis.text.y=element_blank(),
                                                                               axis.title.y=element_text(size=10, face='bold'),
                                                                               axis.text.x=element_text(size=10, face='bold',angle = 90, vjust = 0.5, hjust=1))+
  annotate("text", x = c(1), y=-30, label = c("193"), fontface =2, size=8) +
  annotate("text", x = c(1), y=55, label = c("50"), fontface =2, size=8) +
  annotate("text", x = c(2), y=-30, label = c("199"), fontface =2, size=8) + 
  annotate("text", x = c(2), y=55, label = c("112"), fontface =2, size=8) +
  annotate("text", x = c(3), y=-30, label = c("181"), fontface =2, size=8) +
  annotate("text", x = c(3), y=55, label = c("131"), fontface =2, size=8) +
  annotate("text", x = c(4), y=-30, label = c("203"), fontface =2, size=8) + 
  annotate("text", x = c(4), y=130, label = c("1048"), fontface =2, size=8)+
  annotate("text", x = c(5), y=-30, label = c("194"), fontface =2, size=8) + 
  annotate("text", x = c(5), y=55, label = c("278"), fontface =2, size=8) +
  annotate("text", x = c(6), y=-30, label = c("191"), fontface =2, size=8) + 
  annotate("text", x = c(6), y=55, label = c("214"), fontface =2, size=8)+
  annotate("text", x = c(7), y=-30, label = c("204"), fontface =2, size=8) + 
  annotate("text", x = c(7), y=55, label = c("188"), fontface =2, size=8) +
  annotate("text", x = c(8), y=-30, label = c("202"), fontface =2, size=8) + 
  annotate("text", x = c(8), y=55, label = c("111"), fontface =2, size=8)



#################################################################################################


plot_info5 = data.frame("Dataset" = rep("Dataset 1: Controlled Mixture - DDA - PD", 8),
                        "Version" = rep(c("MSstats","DeqMS","MSqRob2","pmartR","DEP","ProDA","Prolfqua","limma"),2),
                        "Stat" = c("TP", "TP", "TP", "TP","TP","TP","TP","TP", "FP", "FP", "FP", "FP","FP","FP","FP","FP"),
                        "Value" = c(-199,-186,-188,-99,-201,-165,-185,-190,205,58,127,165,124,86,56,52))


## Method Comp
plot_info5_c<-plot_info5 %>% mutate(Version = factor(Version, levels=c("MSstats","DeqMS","MSqRob2","pmartR","DEP","ProDA","Prolfqua","limma"))) %>% 
  ggplot(aes(x=Version, y=Value, fill=Version)) + 
  geom_bar(stat="identity", position="identity") +
  scale_fill_manual(values=c("violet","#33FFFF","yellow","green","purple","magenta","orange","blue"))+
  geom_hline(yintercept = 0) +ylim(-205,300)+
  labs(y = " TP                                        FP",
       x="", title="Dataset 1: Controlled Mixture - DDA - PD") + theme(axis.text.y=element_blank(),
                                                                       axis.title.y=element_text(size=10, face='bold'),
                                                                       axis.text.x=element_text(size=10, face='bold',angle = 90, vjust = 0.5, hjust=1))+
  annotate("text", x = c(1), y=-24.5, label = c("199"), fontface =2, size=8) +
  annotate("text", x = c(1), y=30, label = c("205"), fontface =2, size=8) +
  annotate("text", x = c(2), y=-24.5, label = c("186"), fontface =2, size=8) +
  annotate("text", x = c(2), y=30, label = c("58"), fontface =2, size=8) +
  annotate("text", x = c(3), y=-15, label = c("188"), fontface =2, size=8) + 
  annotate("text", x = c(3), y=30, label = c("127"), fontface =2, size=8)+
  annotate("text", x = c(4), y=-22.5, label = c("99"), fontface =2, size=8) + 
  annotate("text", x = c(4), y=30, label = c("165"), fontface =2, size=8) +
  annotate("text", x = c(5), y=-22.5, label = c("201"), fontface =2, size=8) + 
  annotate("text", x = c(5), y=30, label = c("124"), fontface =2, size=8) +
  annotate("text", x = c(6), y=-22.5, label = c("165"), fontface =2, size=8) + 
  annotate("text", x = c(6), y=30, label = c("86"), fontface =2, size=8) +
  annotate("text", x = c(7), y=-22.5, label = c("185"), fontface =2, size=8) + 
  annotate("text", x = c(7), y=30, label = c("56"), fontface =2, size=8) +
  annotate("text", x = c(8), y=-22.5, label = c("190"), fontface =2, size=8) + 
  annotate("text", x = c(8), y=30, label = c("52"), fontface =2, size=8) 

#################################################################################################


plot_info6 = data.frame("Dataset" = rep("Dataset 2: Choi2017 - DDA - MAxQuant", 6),
                        "Version" = rep(c("MSstats","DeqMS","MSqRob2","pmartR","DEP","ProDA"),2),
                        "Stat" = c("TP", "TP", "TP", "TP","TP","TP", "FP", "FP", "FP", "FP","FP","FP"),
                        "Value" = c(-21,-20,-20,-19,-21,-20,27,20,10,202,187,26))


## Method Comp
plot_info6 %>% mutate(Version = factor(Version, levels=c("MSstats","DeqMS","MSqRob2","pmartR","DEP","ProDA"))) %>% 
  ggplot(aes(x=Version, y=Value, fill=Version)) + 
  geom_bar(stat="identity", position="identity") +
  scale_fill_manual(values=c("violet","#33FFFF","yellow","green","purple","magenta"))+
  geom_hline(yintercept = 0) + ylim(-25,230)+
  labs(y = " TP                                        FP",
       x="", title="Dataset 2: Choi2017 - DDA - MAxQuant") + theme(axis.text.y=element_blank(),
                                                                   axis.title.y=element_text(size=10, face='bold'),
                                                                   axis.text.x=element_text(size=10, face='bold',angle = 90, vjust = 0.5, hjust=1))+
  annotate("text", x = c(1), y=-24.5, label = c("21"), fontface =2, size=8) +
  annotate("text", x = c(1), y=30, label = c("27"), fontface =2, size=8) +
  annotate("text", x = c(2), y=-22.5, label = c("20"), fontface =2, size=8) + 
  annotate("text", x = c(2), y=30, label = c("20"), fontface =2, size=8) +
  annotate("text", x = c(3), y=-15, label = c("20"), fontface =2, size=8) + 
  annotate("text", x = c(3), y=30, label = c("10"), fontface =2, size=8)+
  annotate("text", x = c(4), y=-22.5, label = c("19"), fontface =2, size=8) + 
  annotate("text", x = c(4), y=30, label = c("202"), fontface =2, size=8) +
  annotate("text", x = c(5), y=-22.5, label = c("21"), fontface =2, size=8) + 
  annotate("text", x = c(5), y=30, label = c("187"), fontface =2, size=8) +
  annotate("text", x = c(6), y=-24.5, label = c("20"), fontface =2, size=8) +
  annotate("text", x = c(6), y=30, label = c("26"), fontface =2, size=8) 


#################################################################################################


plot_info7 = data.frame("Dataset" = rep("Dataset 2: Choi2017 - DDA - Skyline", 6),
                        "Version" = rep(c("MSstats","DeqMS","MSqRob2","pmartR","DEP","ProDA"),2),
                        "Stat" = c("TP", "TP", "TP", "TP","TP","TP", "FP", "FP", "FP", "FP","FP","FP"),
                        "Value" = c(-23,-13,-24,-23,-19,-15,6,60,15,257,350,94))


## Method Comp
plot_info7 %>% mutate(Version = factor(Version, levels=c("MSstats","DeqMS","MSqRob2","pmartR","DEP","ProDA"))) %>% 
  ggplot(aes(x=Version, y=Value, fill=Version)) + 
  geom_bar(stat="identity", position="identity") +
  scale_fill_manual(values=c("violet","#33FFFF","yellow","green","purple","magenta"))+
  geom_hline(yintercept = 0) + ylim(-25,360)+
  labs(y = " TP                                        FP",
       x="", title="Dataset 2: Choi2017 - DDA - Skyline") + theme(axis.text.y=element_blank(),
                                                                  axis.title.y=element_text(size=10, face='bold'),
                                                                  axis.text.x=element_text(size=10, face='bold',angle = 90, vjust = 0.5, hjust=1))+
  annotate("text", x = c(1), y=-20, label = c("23"), fontface =2, size=6) +
  annotate("text", x = c(1), y=30, label = c("6"), fontface =2, size=6) +
  annotate("text", x = c(2), y=-20, label = c("13"), fontface =2, size=6) + 
  annotate("text", x = c(2), y=30, label = c("60"), fontface =2, size=6)+
  annotate("text", x = c(3), y=-20, label = c("24"), fontface =2, size=6) + 
  annotate("text", x = c(3), y=30, label = c("15"), fontface =2, size=6) +
  annotate("text", x = c(4), y=-20, label = c("23"), fontface =2, size=6) + 
  annotate("text", x = c(4), y=30, label = c("257"), fontface =2, size=6) +
  annotate("text", x = c(5), y=-20, label = c("19"), fontface =2, size=6) + 
  annotate("text", x = c(5), y=30, label = c("350"), fontface =2, size=6) +
  annotate("text", x = c(6), y=-20, label = c("15"), fontface =2, size=6) +
  annotate("text", x = c(6), y=30, label = c("94"), fontface =2, size=6) 

#################################################################################################


plot_info8 = data.frame("Dataset" = rep("Dataset 2: Choi2017 - DDA - Progenesis", 6),
                        "Version" = rep(c("MSstats","DeqMS","MSqRob2","pmartR","DEP","ProDA"),2),
                        "Stat" = c("TP", "TP", "TP", "TP","TP","TP", "FP", "FP", "FP", "FP","FP","FP"),
                        "Value" = c(-22,-24,-24,-28,-28,-25,110,204,254,413,413,219))


## Method Comp
plot_info8 %>% mutate(Version = factor(Version, levels=c("MSstats","DeqMS","MSqRob2","pmartR","DEP","ProDA"))) %>% 
  ggplot(aes(x=Version, y=Value, fill=Version)) + 
  geom_bar(stat="identity", position="identity") +
  scale_fill_manual(values=c("violet","#33FFFF","yellow","green","purple","magenta"))+
  geom_hline(yintercept = 0) + ylim(-30,600)+
  labs(y = " TP                                             FP",
       x="", title="Dataset 2: Choi2017 - DDA - Progenesis") + theme(axis.text.y=element_blank(),
                                                                     axis.title.y=element_text(size=10, face='bold'),
                                                                     axis.text.x=element_text(size=10, face='bold',angle = 90, vjust = 0.5, hjust=1))+
  annotate("text", x = c(1), y=-20, label = c("22"), fontface =2, size=6) +
  annotate("text", x = c(1), y=50, label = c("110"), fontface =2, size=6) +
  annotate("text", x = c(2), y=-20, label = c("24"), fontface =2, size=6) + 
  annotate("text", x = c(2), y=50, label = c("204"), fontface =2, size=6) +
  annotate("text", x = c(3), y=-20, label = c("24"), fontface =2, size=6) +
  annotate("text", x = c(3), y=20, label = c("254"), fontface =2, size=6) +
  annotate("text", x = c(4), y=-20, label = c("28"), fontface =2, size=6) + 
  annotate("text", x = c(4), y=20, label = c("413"), fontface =2, size=6) +
  annotate("text", x = c(5), y=-20, label = c("28"), fontface =2, size=6) + 
  annotate("text", x = c(5), y=20, label = c("413"), fontface =2, size=6) +
  annotate("text", x = c(6), y=-20, label = c("25"), fontface =2, size=6) + 
  annotate("text", x = c(6), y=50, label = c("219"), fontface =2, size=6)


#########################################################################################################


#Without converter

##########################################################################################################




plot_info2 = data.frame("Dataset" = rep("Dataset 1: Controlled Mixture - DDA - MaxQuant", 8),
                        "Version" = rep(c("MSstats","DeqMS","MSqRob2","pmartR","DEP","ProDA","Prolfqua","limma"),2),
                        "Stat" = c("TP", "TP", "TP", "TP","TP","TP","TP","TP", "FP", "FP", "FP", "FP","FP","FP","FP","FP"),
                        "Value" = c(-191,-196,-123,-117,-198,-160,-185,-196,77,80,85,197,94,72,124,83))


## Method Comp
plot_info2<-plot_info2 %>% mutate(Version = factor(Version, levels=c("MSstats","DeqMS","MSqRob2","pmartR","DEP","ProDA","Prolfqua","limma"))) %>% 
  ggplot(aes(x=Version, y=Value, fill=Version)) + 
  geom_bar(stat="identity", position="identity") +
  scale_fill_manual(values=c("violet","#33FFFF","yellow","green","purple","magenta","orange","blue"))+
  geom_hline(yintercept = 0) +
  labs(y = "TP                                               FP",
       x="", title="Dataset 1: Controlled Mixture - DDA - MaxQuant") + theme(axis.text.y=element_blank(),
                                                                             axis.title.y=element_text(size=10, face='bold'),
                                                                             axis.text.x=element_text(size=10, face='bold',angle = 90, vjust = 0.5, hjust=1))+
  annotate("text", x = c(1), y=-24.5, label = c("191"), fontface =2, size=8) +
  annotate("text", x = c(1), y=30, label = c("77"), fontface =2, size=8) +
  annotate("text", x = c(2), y=-22.5, label = c("196"), fontface =2, size=8) + 
  annotate("text", x = c(2), y=30, label = c("80"), fontface =2, size=8) +
  annotate("text", x = c(3), y=-22.5, label = c("123"), fontface =2, size=8) + 
  annotate("text", x = c(3), y=30, label = c("85"), fontface =2, size=8) +
  annotate("text", x = c(4), y=-22.5, label = c("117"), fontface =2, size=8) + 
  annotate("text", x = c(4), y=30, label = c("197"), fontface =2, size=8) +
  annotate("text", x = c(5), y=-24, label = c("198"), fontface =2, size=8) + 
  annotate("text", x = c(5), y=30, label = c("94"), fontface =2, size=8)+
  annotate("text", x = c(6), y=-24.5, label = c("160"), fontface =2, size=8) +
  annotate("text", x = c(6), y=30, label = c("72"), fontface =2, size=8) +
  annotate("text", x = c(7), y=-24, label = c("185"), fontface =2, size=8) + 
  annotate("text", x = c(7), y=30, label = c("124"), fontface =2, size=8)+
  annotate("text", x = c(8), y=-24.5, label = c("196"), fontface =2, size=8) +
  annotate("text", x = c(8), y=30, label = c("83"), fontface =2, size=8) 

#################################################################################################


plot_info3 = data.frame("Dataset" = rep("Dataset 1: Controlled Mixture - DDA - Skyline", 8),
                        "Version" = rep(c("MSstats","DeqMS","MSqRob2","pmartR","DEP","ProDA","Prolfqua","limma"),2),
                        "Stat" = c("TP", "TP", "TP", "TP","TP","TP", "TP","TP","FP", "FP", "FP", "FP","FP","FP","FP","FP"),
                        "Value" = c(-173,-136,-137,-131,-31,-8,-190,-140,92,100,51,142,151,11,159,93))


## Method Comp
plot_info3<-plot_info3 %>% mutate(Version = factor(Version, levels=c("MSstats","DeqMS","MSqRob2","pmartR","DEP","ProDA","Prolfqua","limma"))) %>% 
  ggplot(aes(x=Version, y=Value, fill=Version)) + 
  geom_bar(stat="identity", position="identity") +
  scale_fill_manual(values=c("violet","#33FFFF","yellow","green","purple","magenta","orange","blue"))+
  geom_hline(yintercept = 0) +ylim(-200,250)+
  labs(y = " TP                                         FP",
       x="", title="Dataset 1: Controlled Mixture - DDA - Skyline") + theme(axis.text.y=element_blank(),
                                                                            axis.title.y=element_text(size=10, face='bold'),
                                                                            axis.text.x=element_text(size=10, face='bold',angle = 90, vjust = 0.5, hjust=1))+
  annotate("text", x = c(1), y=-24.5, label = c("173"), fontface =2, size=8) +
  annotate("text", x = c(1), y=30, label = c("92"), fontface =2, size=8) +
  annotate("text", x = c(2), y=-22.5, label = c("136"), fontface =2, size=8) + 
  annotate("text", x = c(2), y=30, label = c("100"), fontface =2, size=8) +
  annotate("text", x = c(3), y=-24, label = c("137"), fontface =2, size=8) + 
  annotate("text", x = c(3), y=30, label = c("51"), fontface =2, size=8)+
  annotate("text", x = c(4), y=-24.5, label = c("131"), fontface =2, size=8) +
  annotate("text", x = c(4), y=30, label = c("142"), fontface =2, size=8) +
  annotate("text", x = c(5), y=-22.5, label = c("31"), fontface =2, size=8) + 
  annotate("text", x = c(5), y=30, label = c("151"), fontface =2, size=8) +
  annotate("text", x = c(6), y=-22.5, label = c("8"), fontface =2, size=8) + 
  annotate("text", x = c(6), y=30, label = c("11"), fontface =2, size=8) +
  annotate("text", x = c(7), y=-22.5, label = c("190"), fontface =2, size=8) + 
  annotate("text", x = c(7), y=30, label = c("159"), fontface =2, size=8) +
  annotate("text", x = c(8), y=-24, label = c("140"), fontface =2, size=8) + 
  annotate("text", x = c(8), y=30, label = c("93"), fontface =2, size=8)

########################################################################################################


plot_info4 = data.frame("Dataset" = rep("Dataset 1: Controlled Mixture - DDA - Progenesis", 8),
                        "Version" = rep(c("MSstats","DeqMS","MSqRob2","pmartR","DEP","ProDA","Prolfqua","limma"),2),
                        "Stat" = c("TP", "TP", "TP", "TP","TP","TP","TP","TP","FP", "FP", "FP", "FP","FP","FP","FP","FP"),
                        "Value" = c(-193,-194,-181,-201,-203,-194,-211,-193,50,460,134,550,470,417,398,394))


## Method Comp
plot_info4<-plot_info4 %>% mutate(Version = factor(Version, levels=c("MSstats","DeqMS","MSqRob2","pmartR","DEP","ProDA","Prolfqua","limma"))) %>% 
  ggplot(aes(x=Version, y=Value, fill=Version)) + 
  geom_bar(stat="identity", position="identity") +
  scale_fill_manual(values=c("violet","#33FFFF","yellow","green","purple","magenta","orange","blue"))+
  geom_hline(yintercept = 0) +ylim(-300,550)+
  labs(y = " TP                                           FP",
       x="", title="Dataset 1: Controlled Mixture - DDA - Progenesis") + theme(axis.text.y=element_blank(),
                                                                               axis.title.y=element_text(size=10, face='bold'),
                                                                               axis.text.x=element_text(size=10, face='bold',angle = 90, vjust = 0.5, hjust=1))+
  annotate("text", x = c(1), y=-50, label = c("193"), fontface =2, size=8) +
  annotate("text", x = c(1), y=55, label = c("50"), fontface =2, size=8) +
  annotate("text", x = c(2), y=-50, label = c("194"), fontface =2, size=8) + 
  annotate("text", x = c(2), y=80, label = c("460"), fontface =2, size=8) +
  annotate("text", x = c(3), y=-50, label = c("181"), fontface =2, size=8) + 
  annotate("text", x = c(3), y=55, label = c("134"), fontface =2, size=8)+
  annotate("text", x = c(4), y=-50, label = c("201"), fontface =2, size=8) +
  annotate("text", x = c(4), y=80, label = c("938"), fontface =2, size=8) +
  annotate("text", x = c(5), y=-50, label = c("203"), fontface =2, size=8) + 
  annotate("text", x = c(5), y=55, label = c("470"), fontface =2, size=8) +
  annotate("text", x = c(6), y=-50, label = c("194"), fontface =2, size=8) + 
  annotate("text", x = c(6), y=80, label = c("417"), fontface =2, size=8) +
  annotate("text", x = c(7), y=-50, label = c("211"), fontface =2, size=8) + 
  annotate("text", x = c(7), y=55, label = c("398"), fontface =2, size=8) +
  annotate("text", x = c(8), y=-50, label = c("193"), fontface =2, size=8) + 
  annotate("text", x = c(8), y=80, label = c("394"), fontface =2, size=8) 

#################################################################################################


plot_info5 = data.frame("Dataset" = rep("Dataset 1: Controlled Mixture - DDA - PD", 8),
                        "Version" = rep(c("MSstats","DeqMS","MSqRob2","pmartR","DEP","ProDA","Prolfqua","limma"),2),
                        "Stat" = c("TP", "TP", "TP", "TP","TP","TP","TP","TP","FP", "FP", "FP", "FP","FP","FP","FP","FP"),
                        "Value" = c(-199,-72,-56,-41,-100,-48,-61,-68,205,24,22,133,77,25,29,17))


## Method Comp
plot_info5<-plot_info5 %>% mutate(Version = factor(Version, levels=c("MSstats","DeqMS","MSqRob2","pmartR","DEP","ProDA","Prolfqua","limma"))) %>% 
  ggplot(aes(x=Version, y=Value, fill=Version)) + 
  geom_bar(stat="identity", position="identity") +
  scale_fill_manual(values=c("violet","#33FFFF","yellow","green","purple","magenta","orange","blue"))+
  geom_hline(yintercept = 0) +ylim(-200,300)+
  labs(y = " TP                                         FP",
       x="", title="Dataset 1: Controlled Mixture - DDA - PD") + theme(axis.text.y=element_blank(),
                                                                       axis.title.y=element_text(size=10, face='bold'),
                                                                       axis.text.x=element_text(size=10, face='bold',angle = 90, vjust = 0.5, hjust=1))+
  annotate("text", x = c(1), y=-24.5, label = c("199"), fontface =2, size=8) +
  annotate("text", x = c(1), y=30, label = c("205"), fontface =2, size=8) +
  annotate("text", x = c(2), y=-24, label = c("72"), fontface =2, size=8) + 
  annotate("text", x = c(2), y=30, label = c("24"), fontface =2, size=8) +
  annotate("text", x = c(3), y=-24, label = c("56"), fontface =2, size=8) + 
  annotate("text", x = c(3), y=30, label = c("22"), fontface =2, size=8)+
  annotate("text", x = c(4), y=-24.5, label = c("41"), fontface =2, size=8) +
  annotate("text", x = c(4), y=30, label = c("133"), fontface =2, size=8) +
  annotate("text", x = c(5), y=-24, label = c("100"), fontface =2, size=8) + 
  annotate("text", x = c(5), y=30, label = c("77"), fontface =2, size=8) +
  annotate("text", x = c(6), y=-24, label = c("48"), fontface =2, size=8) + 
  annotate("text", x = c(6), y=30, label = c("25"), fontface =2, size=8) +
  annotate("text", x = c(7), y=-24, label = c("61"), fontface =2, size=8) + 
  annotate("text", x = c(7), y=30, label = c("29"), fontface =2, size=8) +
  annotate("text", x = c(8), y=-24, label = c("68"), fontface =2, size=8) + 
  annotate("text", x = c(8), y=30, label = c("17"), fontface =2, size=8) 
  
  
  
#################################################################################################


plot_info6 = data.frame("Dataset" = rep("Dataset 2: Choi2017 - DDA - MAxQuant", 6),
                        "Version" = rep(c("MSstats","DeqMS","MSqRob2","pmartR","DEP","ProDA"),2),
                        "Stat" = c("TP", "TP", "TP", "TP","TP","TP", "FP", "FP", "FP", "FP","FP","FP"),
                        "Value" = c(-21,-8,-19,-20,-21,-20,27,3,9,227,216,30))


## Method Comp
plot_info6 %>% mutate(Version = factor(Version, levels=c("MSstats","DeqMS","MSqRob2","pmartR","DEP","ProDA"))) %>% 
  ggplot(aes(x=Version, y=Value, fill=Version)) + 
  geom_bar(stat="identity", position="identity") +
  scale_fill_manual(values=c("violet","#33FFFF","yellow","green","purple","magenta"))+
  geom_hline(yintercept = 0) + ylim(-25,230)+
  labs(y = " TP                                          FP",
       x="", title="Dataset 2: Choi2017 - DDA - MAxQuant") + theme(axis.text.y=element_blank(),
                                                                   axis.title.y=element_text(size=10, face='bold'),
                                                                   axis.text.x=element_text(size=10, face='bold',angle = 90, vjust = 0.5, hjust=1))+
  annotate("text", x = c(1), y=-24.5, label = c("21"), fontface =2, size=8) +
  annotate("text", x = c(1), y=30, label = c("27"), fontface =2, size=8) +
  annotate("text", x = c(2), y=-22.5, label = c("8"), fontface =2, size=8) + 
  annotate("text", x = c(2), y=30, label = c("3"), fontface =2, size=8) +
  annotate("text", x = c(3), y=-15, label = c("19"), fontface =2, size=8) + 
  annotate("text", x = c(3), y=30, label = c("9"), fontface =2, size=8)+
  annotate("text", x = c(4), y=-24.5, label = c("20"), fontface =2, size=8) +
  annotate("text", x = c(4), y=30, label = c("227"), fontface =2, size=8) +
  annotate("text", x = c(5), y=-22.5, label = c("21"), fontface =2, size=8) + 
  annotate("text", x = c(5), y=30, label = c("216"), fontface =2, size=8) +
  annotate("text", x = c(6), y=-22.5, label = c("20"), fontface =2, size=8) + 
  annotate("text", x = c(6), y=30, label = c("30"), fontface =2, size=8) 

#################################################################################################


plot_info7 = data.frame("Dataset" = rep("Dataset 2: Choi2017 - DDA - Skyline", 6),
                        "Version" = rep(c("MSstats","DeqMS","MSqRob2","pmartR","DEP","ProDA"),2),
                        "Stat" = c("TP", "TP", "TP", "TP","TP","TP", "FP", "FP", "FP", "FP","FP","FP"),
                        "Value" = c(-23,-10,-22,-24,-5,-15,6,35,12,245,98,100))


## Method Comp
plot_info7 %>% mutate(Version = factor(Version, levels=c("MSstats","DeqMS","MSqRob2","pmartR","DEP","ProDA"))) %>% 
  ggplot(aes(x=Version, y=Value, fill=Version)) + 
  geom_bar(stat="identity", position="identity") +
  scale_fill_manual(values=c("violet","#33FFFF","yellow","green","purple","magenta"))+
  geom_hline(yintercept = 0) + ylim(-25,360)+
  labs(y ="TP                                                                             FP",
       x="", title="Dataset 2: Choi2017 - DDA - Skyline") + theme(axis.text.y=element_blank(),
                                                                  axis.title.y=element_text(size=10, face='bold'),
                                                                  axis.text.x=element_text(size=10, face='bold',angle = 90, vjust = 0.5, hjust=1))+
  annotate("text", x = c(1), y=-20, label = c("23"), fontface =2, size=6) +
  annotate("text", x = c(1), y=30, label = c("6"), fontface =2, size=6) +
  annotate("text", x = c(2), y=-20, label = c("10"), fontface =2, size=6) + 
  annotate("text", x = c(2), y=30, label = c("35"), fontface =2, size=6) +
  annotate("text", x = c(3), y=-20, label = c("22"), fontface =2, size=6) + 
  annotate("text", x = c(3), y=30, label = c("12"), fontface =2, size=6)+
  annotate("text", x = c(4), y=-20, label = c("24"), fontface =2, size=6) +
  annotate("text", x = c(4), y=30, label = c("245"), fontface =2, size=6) +
  annotate("text", x = c(5), y=-20, label = c("5"), fontface =2, size=6) + 
  annotate("text", x = c(5), y=30, label = c("98"), fontface =2, size=6) +
  annotate("text", x = c(6), y=-20, label = c("15"), fontface =2, size=6) + 
  annotate("text", x = c(6), y=30, label = c("100"), fontface =2, size=6)  

#################################################################################################
plot_info8 = data.frame("Dataset" = rep("Dataset 2: Choi2017 - DDA - Progenesis", 6),
                        "Version" = rep(c("MSstats","DeqMS","MSqRob2","pmartR","DEP","ProDA"),2),
                        "Stat" = c("TP", "TP", "TP", "TP","TP","TP", "FP", "FP", "FP", "FP","FP","FP"),
                        "Value" = c(-22,-19,-24,-22,-19,-17,110,350,179,431,676,335))

## Method Comp
plot_info8 %>% mutate(Version = factor(Version, levels=c("MSstats","DeqMS","MSqRob2","pmartR","DEP","ProDA"))) %>% 
  ggplot(aes(x=Version, y=Value, fill=Version)) + 
  geom_bar(stat="identity", position="identity") +
  scale_fill_manual(values=c("violet","#33FFFF","yellow","green","purple","magenta"))+
  geom_hline(yintercept = 0) + ylim(-30,600)+
  labs(y = " TP                                             FP",
       x="", title="Dataset 2: Choi2017 - DDA - Progenesis") + theme(axis.text.y=element_blank(),
                                                                     axis.title.y=element_text(size=10, face='bold'),
                                                                     axis.text.x=element_text(size=10, face='bold',angle = 90, vjust = 0.5, hjust=1))+
  annotate("text", x = c(1), y=-20, label = c("22"), fontface =2, size=6) +
  annotate("text", x = c(1), y=50, label = c("110"), fontface =2, size=6) +
  annotate("text", x = c(2), y=-20, label = c("19"), fontface =2, size=6) + 
  annotate("text", x = c(2), y=50, label = c("350"), fontface =2, size=6) +
  annotate("text", x = c(3), y=-20, label = c("24"), fontface =2, size=6) +
  annotate("text", x = c(3), y=20, label = c("179"), fontface =2, size=6) +
  annotate("text", x = c(4), y=-20, label = c("22"), fontface =2, size=6) + 
  annotate("text", x = c(4), y=20, label = c("431"), fontface =2, size=6) +
  annotate("text", x = c(5), y=-20, label = c("19"), fontface =2, size=6) + 
  annotate("text", x = c(5), y=20, label = c("676"), fontface =2, size=6) +
  annotate("text", x = c(6), y=-20, label = c("17"), fontface =2, size=6) + 
  annotate("text", x = c(6), y=50, label = c("335"), fontface =2, size=6)

