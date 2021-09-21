library(reshape2)
`%notin%` <- Negate(`%in%`)
source("PRI_function.R")
library(readr)
library(ggplot2)
library(lubridate)
library(tidyverse)

# METABOLOMICS DATA
df.metabo <- read.csv("./data/df.metabo.15.04.2021.csv", header = T, sep = ";")
df.metabo <- df.metabo[, -1]

notiam.id <- c('Control1', 'Control2', 'Control3', 'Control4', 'Control5', 'Control6', 'Control7', 'Control8', 'Pool1', 'Pool2')
df.metabo <- df.metabo[df.metabo$person_id %notin% notiam.id, ]
repeats.id <- c('IAM02_3_batch 4', 'IAM02_5_batch 4', 'IAM02_11_batch 4', 'IAM03_7_batch 4', 'IAM04_3_batch 4', 'IAM05_3_batch 4', 'IAM05_11_batch 4', 'IAM06_1_batch 4', 'IAM06_7_batch 4', 'IAM06_9_batch 4', 'IAM07_3_batch 4', 'IAM07_5_batch 4', 'IAM07_11_batch 4', 'IAM08_3_batch 4', 'IAM09_5_batch 4', 'IAM09_9_batch 4', 'IAM11_3_batch 4', 'IAM11_5_batch 4', 'IAM11_9_batch 4', 'IAM12_3_batch 4', 'IAM12_5_batch 4', 'IAM12_11_batch 4', 'IAM15_3_batch 4', 'IAM16_1_batch 4', 'IAM16_9_batch 4', 'IAM21_1_batch 4', 'IAM21_5_batch 4', 'IAM21_9_batch 4', 'IAM23_3_batch 4', 'IAM23_5_batch 4', 'IAM25_3_batch 4', 'IAM25_7_batch 4', 'IAM25_9_batch 4', 'IAM26_3_batch 4', 'IAM26_5_batch 4', 'IAM26_9_batch 4', 'IAM27_3_batch 4', 'IAM28_1_batch 4', 'IAM28_5_batch 4', 'IAM28_9_batch 4')
df.metabo <- df.metabo[df.metabo$batch_id %notin% repeats.id, ]

m.cvd<-c('Sphingomyelins','Tyr','Phe','Ile','Leu')
df.metabo <- df.metabo[df.metabo$metabolite_label %in% m.cvd & df.metabo$month<11,]
#df.metabo <- df.metabo[df.metabo$metabolite_label %in% m.cvd,]
table(df.metabo$month, df.metabo$metabolite_label)

df.metabo$collecting_date <- substr(df.metabo$collecting_date, 1, 10)
df.metabo$age.c <- floor(interval(dmy(df.metabo$date_of_birth), dmy(df.metabo$collecting_date)) / years(1))

colnames(df.metabo)

metabo.wide<-dcast(df.metabo,person_id+age.c+gender+month ~ metabolite_label, 
                   value.var = "value")
metabo.wide <- metabo.wide %>% mutate(sex=ifelse(gender=="M",1,0)) %>%
  rename(time=month, age=age.c)
metabo.wide$subject <- parse_number(metabo.wide$person_id)
colnames(metabo.wide)
metabo.wide <- metabo.wide[,c(11,4,2,10,5:9)]

jqm.res.metabo<-data.frame(NA)
for (i in 1:5) {
  db <- metabo.wide[,c(1:4,i+4)]
  db<-db%>%rename(y=colnames(db)[5])
  db<-subset(db, is.na(y)==FALSE)
  subjects<-unique(db$subject)
  for (j in subjects) {
    db[db$subject==j,]$time <- seq(1:length(db[db$subject==j,]$time))
  }
  source("JQM_Function_08072021.R")
  alpha=0.05
  res<-jqm(db=db,
           alpha=alpha,
           lambda.u.seq = seq(0.02,0.1,0.02),
           lambda.z.seq = seq(0.5,5,0.5))
  
  uz <- cbind.data.frame(res$u, res$z, db$age[db$time==1], db$sex[db$time==1])
  colnames(uz)<-c("u","z","age","sex")
  # uz$low <- res$beta0 + res$u + res$z*res$beta1 + res$beta3*uz$age + res$beta4*uz$sex
  # uz$up  <- res$beta0 + res$u + res$z*res$beta2 + res$beta3*uz$age + res$beta4*uz$sex
  uz$low <- res$beta0 + res$u + res$z*res$beta1
  uz$up  <- res$beta0 + res$u + res$z*res$beta2
  
  
  jqm.res.metabo<-cbind(jqm.res.metabo, uz[,5:6])
  
}
jqm.res.metabo<-jqm.res.metabo[,-1]
colnames(jqm.res.metabo) <- rep(colnames(metabo.wide[,5:9]),each=2)


# CLINICAL DATA
df.clinical <- read.csv("./data/clinical_27.05.21.csv", header = T, sep = ",")
df.clinical <- df.clinical[, -1]
iamt<-c( "iamtest01", "iamtest02", "iamtest03")
c.cvd<-c('Vitamine B12','Creatinine','Vitamine D (25-OH-vit D)','Triglyceriden','Cholesterol','Foliumzuur in serum','HDL-cholesterol','IJzer','Cortisol','Aldosteron')
df.clinical <- df.clinical[df.clinical$label %in% c.cvd & df.clinical$person_id %notin% iamt &
                             df.clinical$sample_type!="fasting urine" & df.clinical$month<11,]
# df.clinical <- df.clinical[df.clinical$label %in% c.cvd & df.clinical$person_id %notin% iamt & 
#                              df.clinical$sample_type!="fasting urine",]
table(df.clinical$month, df.clinical$label)

colnames(df.clinical)

df.clinical <- df.clinical[,c("person_id", "month","age","gender","label","value")]
clin.wide <- dcast(df.clinical, person_id+age+gender+month ~ label, 
                   value.var = "value")
clin.wide <- clin.wide %>% mutate(sex=ifelse(gender=="M",1,0)) %>%
  select(-age) %>%
  rename(time=month)
clin.wide$subject <- parse_number(clin.wide$person_id)

#clin.wide$month_id <- paste0(clin.wide$person_id, "_", clin.wide$month)

clin.wide <- clin.wide %>% left_join(metabo.wide %>% select(c(subject,time,age)),
                                     by=c("subject","time"))
for (i in 1:nrow(clin.wide)) {
  clin.wide$age[i] <- ifelse(is.na(clin.wide$age[i]),clin.wide$age[i-1],clin.wide$age[i])  
}
colnames(clin.wide)
clin.wide <- clin.wide[,c(15,3,16,14,4:13)]

jqm.res.clin<-data.frame(NA)
for (i in 1:10) {
  db <- clin.wide[,c(1:4,i+4)]
  db<-db%>%rename(y=colnames(db)[5])
  db<-subset(db, is.na(y)==FALSE)
  subjects<-unique(db$subject)
  for (j in subjects) {
    db[db$subject==j,]$time <- seq(1:length(db[db$subject==j,]$time))
  }
  source("JQM_Function_08072021.R")
  alpha=0.05
  res<-jqm(db=db,
           alpha=alpha,
           lambda.u.seq = seq(0.02,0.1,0.02),
           lambda.z.seq = seq(0.5,5,0.5))
  
  uz <- cbind.data.frame(res$u, res$z, db$age[db$time==1], db$sex[db$time==1])
  colnames(uz)<-c("u","z","age","sex")
  # uz$low <- res$beta0 + res$u + res$z*res$beta1 + res$beta3*uz$age + res$beta4*uz$sex
  # uz$up  <- res$beta0 + res$u + res$z*res$beta2 + res$beta3*uz$age + res$beta4*uz$sex
  uz$low <- res$beta0 + res$u + res$z*res$beta1
  uz$up  <- res$beta0 + res$u + res$z*res$beta2
  
  jqm.res.clin<-cbind(jqm.res.clin, uz[,5:6])
  
}
jqm.res.clin<-jqm.res.clin[,-1]
colnames(jqm.res.clin) <- rep(colnames(clin.wide[,5:14]),each=2)

jqm.res<-cbind(jqm.res.clin, jqm.res.metabo)
jqm.res.noconf<-jqm.res
save(jqm.res.noconf, file="./data/jqm.res.noconf.RData")

cvd.iri.low <- cbind(jqm.res.clin[,c(seq(1,20,by=2))], jqm.res.metabo[,c(1,3,5,7,9)])
cvd.iri.up <- cbind(jqm.res.clin[,c(seq(2,20,by=2))], jqm.res.metabo[,c(2,4,6,8,10)])
save(cvd.iri.low, file="./data/cvd.iri.low.RData")
save(cvd.iri.up, file="./data/cvd.iri.up.RData")

cvd.data <- clin.wide %>% left_join(select(metabo.wide, -c(age, sex)),
                                    by=c("subject", "time"))
save(cvd.data, file="./data/cvd.data.RData")

my_theme <-   theme(legend.position = "none",
                    legend.text = element_text(size=12),
                    legend.title = element_text(size=12),
                    axis.text.x = element_text(size = 14),
                    axis.text.y = element_text(size = 14, face = "bold"),
                    axis.title.x = element_text(size=14),
                    axis.title.y = element_text(size=14),
                    axis.ticks = element_blank(),
                    panel.border = element_rect(NA),
                    strip.text = element_text(size = 11, face = "bold"))  
j<-1
for (i in seq(1,30,by=2)) {
  uz<-jqm.res[,c(i,i+1)]
  colnames(uz)<-c("low","up")
  db<-cvd.data[,c(1:4,j+4)]
  colnames(db)<-c(colnames(cvd.data)[1:4],"y")
  
  plt<-ggplot(uz, aes(x=as.factor(1:nrow(uz)))) +
    geom_errorbar(aes(ymin = low, ymax = up), width=0.75, size=1,
                  position=position_dodge(width=0.7), color="darkblue")+
    geom_point(data=db, aes(x =as.factor(subject), y = y, group=sex, color=sex),
               size=2.5, position=position_dodge(width=0.7)) +
    geom_vline(xintercept=seq(1.5, length(unique(db$subject))-0.5, 1),
               lwd=0.5, colour="grey") +
    labs(x="Individual", title = paste0(colnames(jqm.res)[i])) +
    theme_classic()+
    my_theme
  
  ggsave(paste0("iri_", colnames(jqm.res)[i], ".png"))

  j<-j+1
}
plot.this <- ggplot(uz.this, aes(x=as.factor(id2))) +
  geom_errorbar(aes(ymin = low, ymax = up), width=0.75, size=1.5,
                position=position_dodge(width=0.7), color="darkblue") +
  geom_point(data = df.this, aes(x =as.factor(id2), y = y, group=Gender, color=Gender),
             size=3.5, position=position_dodge(width=0.7)) +
  geom_vline(xintercept=seq(1.5, length(unique(df.this$id2))-0.5, 1), 
             lwd=0.5, colour="grey") +
  scale_color_manual(values=c("seagreen", "red3")) +
  annotate("rect", xmin = 0, xmax =  31, ymin = 60, ymax = 80,
           alpha = .1, fill="#92D050") +
  annotate("segment", x = 0, xend =  31, y = 80, yend = 80, color = "#746f6e",
           alpha=1, linetype="dashed") +
  annotate("segment", x = 0, xend =  31, y = 60, yend = 60, color = "#746f6e",
           alpha=1, linetype="dashed") +
  theme_classic()+
  my_theme
plot.this 

load("./data/cvd.iri.low.RData")
load("./data/cvd.iri.up.RData")
cvd.iri<-cbind(cvd.iri.low,cvd.iri.up)
cvd.iri<-cvd.iri[,order(colnames(cvd.iri))]

load("./data/jqm.res.noconf.RData")
cvd.iri<-jqm.res.noconf
cvd.iri<-cvd.iri[,order(colnames(cvd.iri))]

load("./data/cvd.data.RData")
data2<-cvd.data[,5:19]
data2<-data2[,order(colnames(data2))]
cvd.data<-cbind(cvd.data[,1:4],data2)

colnames(cvd.data)
t.11<-cvd.data[cvd.data$time==11,]
rownames(t.11)<-NULL
t.11.iri2<-cbind.data.frame(t.11$subject, matrix(NA, nrow=30,ncol=15))
colnames(t.11.iri2)<-c("subject",colnames(t.11)[5:19])
j<-1 #for columns sequence in t.11
for (k in seq(1,30,by=2)) {
  for (i in 1:nrow(t.11)) {
    t.11.iri2[i,j+1]<-ifelse(t.11[i,j+4]>cvd.iri[i,k] & t.11[i,j+4]<cvd.iri[i,k+1], 1, 0)
  }
  j<-j+1
}


load("./data/cvd.data.RData")
ri.cvd<-read.csv("./data/CVD_pri.csv", sep=";", header = T)

seq1<-seq(1,11,by=2)
t.prot<-cvd.data[cvd.data$time %in% seq1,]

rownames(t.prot)<-NULL
t.stat<-cbind.data.frame(t.prot[,c(1,2)], matrix(NA, nrow=180,ncol=15))

colnames(t.stat)<-c("subject","time",colnames(t.prot)[5:19])

#0 for inside, 1 for outside
for (j in 1:nrow(ri.cvd)) {
  for (i in 1:nrow(t.prot)) {
    t.stat[i,j+2]<-ifelse(t.prot[i,j+4]>ri.cvd[j,2] & t.prot[i,j+4]<ri.cvd[j,3], 0, 1)
  }
}

#1 for below 'low', 2 for within, 3 for above 'up'
for (j in 1:nrow(ri.cvd)) {
  for (i in 1:nrow(t.prot)) {
    t.stat[i,j+2]<-ifelse(t.prot[i,j+4]<ri.cvd[j,2], 1, 
                          ifelse(t.prot[i,j+4]>ri.cvd[j,3], 3, 2))
  }
}

#wide<-spread(t.stat[,c(1:3)], subject, Aldosteron)
# from long to wide
res.wide<-NA
for (i in 3:ncol(t.stat)) {
  db<-t.stat[,c(1:2,i)]
  colnames(db)<-c("subject","time","y")
  
  wide<-spread(db, subject, y)
  wide$prot<-colnames(t.stat)[i]
  
  res.wide<-rbind(res.wide, wide)
}
res.wide<-res.wide[-1,]
res.wide<-res.wide[,c(32, 1:31)]
write.csv(res.wide, file="./data/int_clin.csv", row.names = F)


my_theme <-   theme(legend.position = "bottom",
                    legend.text = element_text(size=16),
                    legend.title = element_text(size=16),
                    axis.text.x = element_text(size=13),
                    axis.text.y = element_text(size = 14, face = "bold"),
                    axis.title.x = element_text(size=16),
                    axis.title.y = element_text(size=20),
                    axis.ticks = element_blank(),
                    panel.border = element_rect(NA),
                    strip.text = element_text(size = 11, face = "bold"))  
color.month = c('black','forestgreen', 'red2',
                'orange', 'cornflowerblue', 'magenta',
                'darkolivegreen4','indianred1','tan4',
                'darkblue', 'mediumorchid1', 'firebrick4',
                'yellowgreen')

color.month = c("grey","grey","grey","grey","grey",
                "grey","grey","grey","grey","grey",
                "forestgreen", 'magenta', 'orange')


cvd.data$time<-as.factor(cvd.data$time)
for (i in 1:15) {
  db<-cvd.data[,c(1:4,i+4)]
  colnames(db)<-c(colnames(cvd.data)[1:4],"y")
  
  plt<-ggplot(db) +
     geom_point(aes(x=factor(subject), y=y, group=time, color=time), size=3.5) +
     geom_vline(xintercept=seq(1.5, length(unique(db$subject))-0.5, 1), 
                lwd=0.5, colour="grey") +
     labs(x="Individual", title = paste0(colnames(db)[i])) +
     scale_color_manual(name="Months", values= color.month)+
     theme_classic()+ my_theme
  
  ggsave(paste0("ind_", colnames(cvd.data)[i+4], ".png"))
   
}


load("./data/cvd.data.RData")
load("./data/jqm.res.noconf.RData")
jqm.res<-jqm.res.noconf
data<-cvd.data
data$time<-as.factor(data$time)
ri.cvd<-read.csv("./data/CVD_pri.csv", sep=";", header = T)

# iri with confounding
load("./data/cvd.iri.low.RData")
load("./data/cvd.iri.up.RData")
cvd.iri<-cbind(cvd.iri.low,cvd.iri.up)
cvd.iri<-cvd.iri[,order(colnames(cvd.iri))]
jqm.res<-cvd.iri

load("./data/cvd.data.RData")
dat<-cvd.data[cvd.data$time<12,5:19]
data<-cbind(cvd.data[cvd.data$time<12,1:4], dat[,order(colnames(dat))])
data$time<-as.factor(data$time)

ri.cvd<-ri.cvd[order(ri.cvd$CVD.compound),]

my_theme <-   theme(legend.position = "bottom",
                    legend.text = element_text(size=12),
                    legend.title = element_text(size=12),
                    axis.text.x = element_text(size = 12),
                    axis.text.y = element_text(size = 14, face = "bold"),
                    axis.title.x = element_text(size=14),
                    axis.title.y = element_text(size=14),
                    axis.ticks = element_blank(),
                    panel.border = element_rect(NA),
                    strip.text = element_text(size = 11, face = "bold")) 
color.month = c("grey","grey","grey","grey","grey",
                "forestgreen")
#, 'magenta', 'orange'
color.month = c('black','forestgreen', 'red2',
                'orange', 'cornflowerblue', 'magenta',
                'darkolivegreen4','indianred1','tan4',
                'darkblue', 'mediumorchid1', 'firebrick4',
                'yellowgreen')

#jqm.res<-jqm.res.prot.noconf
j<-1
for (i in seq(1,30,by=2)) {
  uz<-jqm.res[,c(i,i+1)]
  colnames(uz)<-c("low","up")
  db<-data[,c(1:4,j+4)]
  colnames(db)<-c(colnames(data)[1:4],"y")
  db<-db[is.na(db$y)==F,]
  
  pri<-ri.cvd[j,]
  
  plt<-ggplot(uz, aes(x=as.factor(1:nrow(uz)))) +
    geom_errorbar(aes(ymin = low, ymax = up), width=0.75, size=1,
                  position=position_dodge(width=0.7), color="darkblue")+
    geom_point(data=db, aes(x =as.factor(subject), y = y, group=time, color=time),
               size=2.5, position=position_dodge(width=0.7)) +
    geom_vline(xintercept=seq(1.5, length(unique(db$subject))-0.5, 1),
               lwd=0.5, colour="grey") +
    annotate("rect", xmin = 0, xmax =  30, ymin = pri[,2], ymax = pri[,3],
             alpha = .1, fill="#92D050") +
    annotate("segment", x = 0, xend =  30, y = pri[,2], yend = pri[,2], color = "#746f6e",
             alpha=1, linetype="dashed") +
    annotate("segment", x = 0, xend =  30, y = pri[,3], yend = pri[,3], color = "#746f6e",
             alpha=1, linetype="dashed") +
    
    scale_color_manual(name="Month", values= color.month)+
    labs(x="Individual", title = paste0(colnames(jqm.res)[i])) +
    theme_classic()+
    my_theme
  
  ggsave(paste0("./IRI_CVD/iri_conf/iri_conf_", colnames(jqm.res)[i], ".png"))
  
  j<-j+1
}

