library(dplyr)
library(ggplot2)

full_dat = read.csv("AOF.csv",header=T)
auc_fun <- function(risk,status){
  mean(sapply(risk[status==1],function(r){
    mean(risk[status==0]<=r)
  }))
}
ap_fun <- function(risk,status){
  ppv = sapply(risk[status==1],function(r){
    mean(status[risk>=r])
  })
  ppv[is.na(ppv)] = 0
  mean(ppv)
}

acc_dat = full_dat %>% group_by(model) %>% summarise(
  auc=auc_fun(risk,status),
  ap = ap_fun(risk,status))
acc_dat

uu = seq(0,1,by=0.01)
curve_dat = full_dat %>% split(.,.$model) %>% 
  lapply(.,function(dat){
    data.frame(
      "u" = uu,
      "TPR" = sapply(uu,function(u){mean(dat$risk[dat$status==1]>=u)}),
      "FPR" = sapply(uu,function(u){mean(dat$risk[dat$status==0]>=u)}),
      "PPV" = sapply(uu,function(u){mean(dat$status[dat$risk>=u]==1)}),
      "model" = unique(dat$model)
    )
  }) %>% do.call(rbind,.)


p1 = ggplot(curve_dat,aes(x=FPR,y=TPR,group=model)) + geom_step(aes(linetype=model)) + 
  scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,1)) + 
  ggtitle("ROC") + 
  geom_text(mapping=aes(x=0.5,y=0.5,label="Prescribed: AUC=0.96"))+
  geom_text(mapping=aes(x=0.5,y=0.42,label="Ovarian: AUC=0.94")) + 
  geom_text(mapping=aes(x=0.5,y=0.34,label="Delta~AUC==-0.02"),parse=TRUE) + 
  theme(legend.position = "none",plot.margin = unit(c(0.1,1,0.1,1), "cm")) 

p2 = ggplot(curve_dat,aes(x=TPR,y=PPV,group=model)) + geom_line(aes(linetype=model)) + 
  scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,1)) + 
  ggtitle("PR") + 
  geom_text(mapping=aes(x=0.5,y=0.25,label="Prescribed: AP=0.46"))+
  geom_text(mapping=aes(x=0.5,y=0.17,label="Ovarian: AP=0.68")) + 
  geom_text(mapping=aes(x=0.5,y=0.09,label="Delta~AP==0.22"),parse=T)
pdf("AOF_ROC_PR.pdf",height=6,width=12)
gridExtra::grid.arrange(p1,p2,nrow=1)
dev.off()

qq = seq(0.01,0.99,by=0.01)
delta_dat = full_dat %>% split(.,.$model) %>% 
  lapply(.,function(dat){
    cc = quantile(dat$risk[dat$status==1],qq)
    data.frame(
      "quantile" = qq,
      "F0" =  sapply(cc,function(c){mean(dat$risk[dat$status==0]<=c)}),
      model=unique(dat$model)
    )
  }) %>% do.call(rbind,.) %>% 
  spread(key="model",value="F0")

pi1 = mean(ovarian.dat$status)
tmpdat = delta_dat %>% 
  mutate(
    Delta = Ovarian - Prescribed,
    wgt = (1/pi1-1)/(1-quantile),
    wgt0 = wgt/((1+wgt*(1-`Prescribed`))*(1+wgt*(1-`Ovarian`))),
    wDelta = Delta*wgt0,
    wgt = max(wDelta)*wgt0/max(wgt0)
  ) %>% select(-Ovarian,-Prescribed) %>% 
  gather(key="type",value="values",-quantile) 

max.wDelta = max(tmpdat$values[tmpdat$type %in% c("wDelta","Delta")])
max.wgt = max(tmpdat$values[tmpdat$type %in% c("wgt0")])
pdf("AOF_Delta.pdf",height=6,width=10)
tmpdat %>% filter(type!="wgt0") %>% 
  mutate(type=factor(type,levels=c("Delta","wDelta","wgt"))) %>% 
  ggplot(.,aes(x=quantile,y=values)) + geom_line(aes(linetype=type)) +
  scale_y_continuous(bquote(Delta(alpha)~"or"~w[AP](alpha)*Delta(alpha)),
                     sec.axis = sec_axis(~ .*max.wgt/max.wDelta, name = expression(w[AP](alpha)))) + 
  scale_linetype_manual(values=c("Delta"="solid","wDelta"="dashed","wgt"="dotted"), labels=expression(Delta(alpha),w[AP](alpha)*Delta(alpha),w[AP](alpha)))   +
  labs(x=bquote(alpha))
dev.off()


