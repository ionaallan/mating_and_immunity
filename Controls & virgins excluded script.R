library(survival)
library(survminer)
library(survivalAnalysis)
library(tidyverse)
library(tidytidbits)
library(DT)
library(coxme)
library(gtsummary)
library(multcomp)
library(scales)
library(gridExtra)
library(grid)
library(ggplot2)
library(kableExtra)
library(ggpubr)



SuperSmallfont= 6
Smallfont= 16
Mediumfont= 18
Largefont= 20
verylargefont = 22
pointsize= 0.7
linesize=1.1
meansize = 1.5
Margin=c(0,0,0,0)
fontsizeaxes = 12
fontsizeaxes2 = 10




#Function to extract coxme table
extract_coxme_table <- function (mod){
  beta <- fixef(mod)
  nvar <- length(beta)
  nfrail <- nrow(mod$var) - nvar
  se <- sqrt(diag(mod$var)[nfrail + 1:nvar])
  z<- round(beta/se, 2)
  p<- format(as.numeric(1 - pchisq((beta/se)^2, 1)), 4)
  table=data.frame(cbind(beta,se,z,p))
  return(table)
}



#Function to extract survival data to plot
ggplotprep2 <- function(x, times){
  d <- data.frame(condition=rep(names(x$strata), x$strata), time=x$time, survival=x$surv, upper=x$upper, lower=x$lower)
  fillup0 <- function(s) rbind(c(condition=s, time=0, survival=1, upper=1, lower=1), d[d$condition==s, ], deparse.level = 0)
  
  indexes <- function(x, time) {
    if(x%in%time) return(x)
    return(floor(time[which.min(abs(time[time<x]-x))]))
  }
  
  fillup <- function(s) {
    d.temp <- d[d$condition==s, ]
    time <- as.numeric(d.temp$time)
    id <- sapply(times, indexes, time=time)
    d.temp <- d.temp[match(id, time), ]
    d.temp$time <- times
    return(d.temp)
  }
  
  if(times[1]==0) d <- do.call("rbind", sapply(names(x$strata), fillup0, simplify=F))
  d <- do.call("rbind", sapply(names(x$strata), fillup, simplify=F))
  clean.name <- function(name) unlist(lapply(strsplit(as.character(name), split="="), function(x) x[2]))
  d <- data.frame(Condition=clean.name(d$condition), Time=as.numeric(d$time), Survival=as.numeric(d$survival), upper=as.numeric(d$upper), lower=as.numeric(d$lower))
  return(d)
}




# Data import (MK's method - doesn't work for me)
# path.to.data <- "/Users/Documents/honours project regan lab/final version CoxPH analysis (minus controls and virgins)/CoxPH input data.csv"

# d <- list()
# path <- list()
# model <- list()
# for(f in list.files(path=path.to.data,pattern="*.csv$",recursive=T,full.names=T)) {
  # nom <- gsub(".*/(.*).csv","\\1",f)    
  # cat(nom,"\n")
  # path[[nom]] <- gsub("(.*)/.*csv","\\1/",f)
  # d[[nom]] <- read.table(f,header=T,sep=",",dec=",")
  
# }

# Mating_Immunity_Survival= d[["CoxPH_input_data"]]



# Data import (my method)
fly_survival_data <- read_csv("/Users/ionaa/Documents/honours project regan lab/final version CoxPH analysis (minus controls and virgins)/CoxPH input data.csv")




# Making Age another factor
Mating_Immunity_Survival_Young <-
  fly_survival_data %>% 
  filter(Group %in%  c("Group5_Trt_Once_Proximal","Group5_Control_Once_Proximal")) %>% 
  mutate(Age = paste("Young"))

Mating_Immunity_Survival_Old <-
  fly_survival_data %>% 
  filter(!Group %in%  c("Group5_Trt_Once_Proximal","Group5_Control_Once_Proximal")) %>% 
  mutate(Age = paste("Aged"))

## Combine the two datasets to make up original again 
fly_survival_data <- rbind(Mating_Immunity_Survival_Old, Mating_Immunity_Survival_Young)




## Making Controlless and Virginless dataframes:

# Remove controls and make a new dataframe
fly_survival_data_Controlless <-
  fly_survival_data %>% 
  filter(!Infection_Status == "Sham") %>% 
  droplevels()
# view(fly_survival_data_Controlless) # commented out so that table doesn't pop up every time code is run


# Remove virgins and make a new dataframe
fly_survival_data_Virginless <-
  fly_survival_data_Controlless %>% 
  filter(!Mating_Status == "Virgin") %>% 
  droplevels()
# view(fly_survival_data_Virginless) # commented out so that table doesn't pop up every time code is run





## CoxPH analyis on Controlless first:

# Making a survival object (controlless)
attach(fly_survival_data_Controlless)
Surv_Object <- Surv(time = Time_to_death_, event = Censor)

# Unordering factors (controlless)
fly_survival_data_Controlless$Mating_Status <- factor(fly_survival_data_Controlless$Mating_Status, ordered = FALSE)
fly_survival_data_Controlless$Group <- factor(fly_survival_data_Controlless$Group, ordered = FALSE)
fly_survival_data_Controlless$Age <- factor(fly_survival_data_Controlless$Age, ordered = FALSE)

# Set the reference level for each factor (controlless)
fly_survival_data_Controlless$Mating_Status <- relevel(fly_survival_data_Controlless$Mating_Status, ref = "Once")
#fly_survival_data_Controlless$Group <- relevel(fly_survival_data_Controlless$Group, ref = "Group1_Trt_Virgin") # Not needed because we're not including group as a covariate
fly_survival_data_Controlless$Age <- relevel(fly_survival_data_Controlless$Age, ref = "Aged")

# Perform CoxPH analysis (controlless)
CoxPH_analysis <- coxph(Surv_Object ~ Mating_Status + 
                          Age + 
                          # Mating_Status*Age +     # Can't be included because not all combos present (get NAs)
                          cluster(Vial), 
                        data = fly_survival_data_Controlless)
CoxPH_analysis
summary(CoxPH_analysis)

# Perform CoxME analysis (Controlless)
coxME <-  coxme(Surv_Object ~ Mating_Status + Age +  (1|Vial) , data = fly_survival_data_Controlless)
summary(coxME)




## CoxPH analysis on Virginless second:

# Making a survival object (virginless)
attach(fly_survival_data_Virginless)
Surv_Object2 <- Surv(time = Time_to_death_, event = Censor)

# Unordering factors (virginless)
fly_survival_data_Virginless$Mating_Status <- factor(fly_survival_data_Virginless$Mating_Status, ordered = FALSE)
fly_survival_data_Virginless$Group <- factor(fly_survival_data_Virginless$Group, ordered = FALSE)
fly_survival_data_Virginless$Age <- factor(fly_survival_data_Virginless$Age, ordered = FALSE)
fly_survival_data_Virginless$Mating_Infection_Distance <- factor(fly_survival_data_Virginless$Mating_Infection_Distance, ordered = FALSE)

# Set the reference level for each factor (virginless)
fly_survival_data_Virginless$Mating_Infection_Distance <- relevel(fly_survival_data_Virginless$Mating_Infection_Distance, ref = "Proximal")
fly_survival_data_Virginless$Mating_Status <- relevel(fly_survival_data_Virginless$Mating_Status, ref = "Once")
fly_survival_data_Virginless$Age <- relevel(fly_survival_data_Virginless$Age, ref = "Aged")

# Perform CoxPH analysis (virginless)
CoxPH_analysis_2 <- coxph(Surv_Object2 ~ Mating_Status + 
                            Age + 
                            Mating_Infection_Distance + 
                            # Mating_Status*Mating_Infection_Distance +   ## these covariates can't be included because not all combinations are present
                            # Mating_Status*Age +                         ## e.g. no distal twice-mated flies
                            # Mating_Infection_Distance*Age +             ## hence get NAs in the analysis
                            cluster(Vial), 
                          data = fly_survival_data_Virginless)
CoxPH_analysis_2
summary(CoxPH_analysis_2)

# Perform CoxME analysis (virginless)
coxME_2 <-  coxme(Surv_Object2 ~ Mating_Status + Age + Mating_Infection_Distance +  (1|Vial) , data = fly_survival_data_Virginless)
summary(coxME)





# Prepare Tukey pairwise comparison of groups: 
CoxmeGroup <- coxme(Surv(Time_to_death_,Censor) ~ Group + (1|Vial), data= fly_survival_data_Controlless)

multcomp = glht(CoxmeGroup, linfct=mcp(Group="Tukey"))
Comp = cld(multcomp)

unlist(Comp$mcletters$Letters)%>%
  kable(col.names = "Sign.group") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = F)


# Anova comparison of groups followed by Tukey comparison:
# anova_of_groups <- aov(fly_survival_data_Controlless$Group)
# TukeyHSD(CoxmeGroup, "Group", ordered = FALSE, conf.level = 0.95)





## Plot Survival curve using CoxME analysis of controlless data: 
survdata <- survfit(Surv_Object~Group, data=fly_survival_data_Controlless)

# pairwise comparison of groups: 
Pairwise <- pairwise_survdiff(Surv(Time_to_death_, Censor) ~ Group, 
                              data = fly_survival_data_Controlless)
Pairwise

#time checked
time_check = sort(unique(fly_survival_data_Controlless$Time_to_death_))

toplot <- ggplotprep2(survdata, times=c(0,time_check))

Name= ""
Limits=c("Group1_Trt_Virgin",
         "Group2_Trt_Once_Distal",
         "Group3_Trt_Twice_Proximal",
         "Group4_Trt_Once_Proximal",
         "Group5_Trt_Once_Proximal")

n1 = survdata$n[1]
n2 = survdata$n[2]
n3 = survdata$n[3]
n4 = survdata$n[4]
n5 = survdata$n[5]



Labels = c(paste("Virgin (n=",n1,") ab",sep=""),
           paste("Once-Mated Distal (n=",n2,") ab",sep=""),
           paste("Twice-Mated (n=",n3,") ab",sep=""),
           paste("Once-Mated Proximal (n=",n4, ") b", sep = ""),
           paste("Young-Mated Proximal (n=",n5,") c",sep="")
)

p1=
  ggplot(subset(toplot,Survival!=0), aes(x=Time,y=Survival,colour=as.factor(Condition)))+
  #  ggtitle("Canton S mated & P. rettgeri 0.02")+
  geom_line(aes(linetype=as.factor(Condition)),size=linesize)+
  geom_point(size = pointsize) +    
  # geom_errorbar(data=subset(toplot,Time==4),aes(ymin=lower, ymax=upper), width=.1, alpha=0.4, size=1, show.legend=FALSE)+
  ggtitle("Survival Curves")+
  scale_colour_manual(Name,
                      limits=Limits,
                      values=c("#00B050","#FF0000", "#1D78FF","#FF9900", "#9900CC"),
                      labels=Labels)+
  scale_linetype_manual(Name,
                        limits=Limits,
                        values=c("solid","solid","solid","solid","solid"),
                        labels=Labels)+
  scale_x_continuous("Days post-infection",
                     limits=c(0, 6),
                     breaks=c(seq(0,6,by=1)))+
  scale_y_continuous("Proportion of survivors",
                     limits=c(0, 1),breaks=c(0,0.2,0.4,0.6,0.8,1))+
  theme(title =element_text(size=Mediumfont, face='bold'),
        plot.title = element_text(hjust = 0.5,size = Largefont),
        panel.background = element_blank(),
        strip.text.x = element_text(size =Mediumfont, colour = "black",face="italic"),
        strip.text.y = element_text(size =Mediumfont, colour = "black",face="italic"),
        axis.text.x = element_text(size=Mediumfont,colour="black"),
        axis.text.y = element_text(size=Mediumfont,colour="black"),
        axis.title.x = element_text(size=Mediumfont,colour="black"),
        axis.title.y = element_text(size=Mediumfont,colour="black"), 
        axis.line.x = element_line(colour="black",size=0.75),
        axis.line.y = element_line(colour="black",size=0.75),
        axis.ticks.x = element_line(size = 0.75),
        axis.ticks.y = element_line(size = 0.75),
        legend.direction = "vertical", 
        legend.box = "horizontal",
        legend.position = c(0.3,0.35),
        legend.key.height = unit(1.4, "cm"),
        legend.key.width= unit(2.6, "cm"),
        legend.title = element_text(face="italic",size=Smallfont), 
        legend.key = element_rect(colour = 'white', fill = "white", linetype='dashed'),
        legend.text = element_text(size=Smallfont),
        legend.background = element_rect(fill=NA),
        plot.margin = unit(c(0,0,1.2,0), "cm"))+
  guides(shape=guide_legend(ncol=1),
         fill=guide_legend(ncol=1),
         col=guide_legend(ncol=1))

p1






# Plot hazard ratio graph using CoxME analysis of Controlless data:
# log(hazard ratio) (i.e. beta) is plotted rather than just hazard ratio because you can get a SE of it
# plotting hazard ratio and calculating 95%CI is also an option

CoxmeGroup <- coxme(Surv(Time_to_death_,Censor) ~ Group + (1|Vial), data= fly_survival_data_Controlless)

tab_res_surv_int = extract_coxme_table(CoxmeGroup)

rownames(tab_res_surv_int) <- sub("GroupG","G",rownames(tab_res_surv_int))

tab_res_surv_int$Treatment = as.factor(rownames(tab_res_surv_int))
tab_res_surv_int$beta = as.numeric(as.character(tab_res_surv_int$beta))
tab_res_surv_int$se = as.numeric(as.character(tab_res_surv_int$se))
tab_res_surv_int$z = as.numeric(as.character(tab_res_surv_int$z))
tab_res_surv_int$p = as.numeric(as.character(tab_res_surv_int$p))
tab_res_surv_int$HazardRatio = exp(tab_res_surv_int$beta)
tab_res_surv_int[nrow(tab_res_surv_int)+1,] = NA

tab_res_surv_int[5,1] = 0
tab_res_surv_int[5,2] = 0
tab_res_surv_int$Treatment = as.character(tab_res_surv_int$Treatment)
tab_res_surv_int[5,5] = "Group_Trt_Virgin"
tab_res_surv_int$Treatment = as.factor(tab_res_surv_int$Treatment)
tab_res_surv_int[5,6] = 1


Limits=c("Group1_Trt_Virgin",
         "Group2_Trt_Once_Distal",
         "Group3_Trt_Twice_Proximal",
         "Group4_Trt_Once_Proximal",
         "Group5_Trt_Once_Proximal")

p2=
  ggplot(tab_res_surv_int, aes(x=Treatment, y=beta)) + 
  geom_errorbar(aes(ymin=beta-se, ymax=beta+se, col=Treatment),
                width=.4,  size=1, show.legend=FALSE)+
  scale_color_manual("Treatment", breaks=c(1,2,3,4),values=c("#00B050","#FF0000", "#1D78FF","#FF9900", "#9900CC", "Grey0"))+
  geom_point(shape="circle", 
             aes(colour=Treatment,
                 fill=Treatment),
             alpha=0.4, stat="identity", show.legend=FALSE,  size=5) +
  scale_shape_manual(values=seq(0,9))+
  scale_y_continuous("log(Hazard ratio) \u00B1se",
                     limits = NULL)+
  geom_hline(yintercept = 0,colour="black",linetype=4)+
  ggtitle("log(Hazard Ratio) relative to Virgins")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust=0.5, size = Largefont),
        axis.title.x = element_text(size=Mediumfont,colour="black"),
        axis.title.y = element_text(size=Mediumfont,colour="black"),
        axis.line.x = element_line(colour="black",size=0.75),
        axis.line.y = element_line(colour="black",size=0.75),
        axis.ticks.x = element_line(size = 0.75),
        axis.ticks.y = element_line(size = 0.75),
        axis.text.x = element_text(size=Smallfont,colour="black",angle=30,hjust=1),
        axis.text.y = element_text(size=Smallfont,colour="black"),
        plot.margin = unit(Margin, "cm"),
        legend.direction = "horizontal", 
        legend.box = "vertical",
        legend.position = "bottom",
        legend.key.height = unit(0.4, "cm"),
        legend.key.width= unit(0.3, "cm"),
        legend.title = element_text(face="italic",size=Smallfont), 
        legend.key = element_rect(colour = 'white', fill = "white", linetype='dashed'),
        legend.text = element_text(size=Smallfont),
        legend.background = element_rect(fill=NA))+
  guides(shape=guide_legend(ncol=5),
         fill=guide_legend(ncol=5),
         col=guide_legend(ncol=5))

p2





# Export the two figures:
ggsave("/Users/ionaa/Documents/honours project regan lab/final version CoxPH analysis (minus controls and virgins)/survival curves.png", 
       p1, width = 15, height = 10, dpi = 500)

ggsave("/Users/ionaa/Documents/honours project regan lab/final version CoxPH analysis (minus controls and virgins)/hazard ratio plot.png", 
       p2, width = 15, height = 10, dpi = 500)




## check assumptions of CoxPH model: 

# Check proportional hazards assumption first
# if all covariates as well as global test are not statistically significant,
# proportional hazards can be assumed
test.ph <- cox.zph(CoxPH_analysis)
test.ph

test.ph_2 <- cox.zph(CoxPH_analysis_2)
test.ph_2

# plot graphs of the scaled Schoenfeld residuals against the transformed time for each covariate: 
# systematic departures from a horizontal line are indicative of non-proportional hazards, 
# since proportional hazards assumes that estimates β1,β2,β3 do not vary much over time.
ggcoxzph(test.ph)
ggcoxzph(test.ph_2)

# Next, check influential observations or outliers by visualising dfbeta values and/or deviance residuals: 
ggcoxdiagnostics(CoxPH_analysis, type = "dfbeta", linear.predictions = FALSE)
ggcoxdiagnostics(CoxPH_analysis_2, type = "dfbeta", linear.predictions = FALSE)

# deviance residuals should be roughly symmetrically distributed about zero, 
# with a standard deviation of one
ggcoxdiagnostics(CoxPH_analysis, type = "deviance", linear.predictions = FALSE)
ggcoxdiagnostics(CoxPH_analysis_2, type = "deviance", linear.predictions = FALSE)

# No need to check non-linearity, as all covariates are categorical
# (None are continuous)


