library(tidyverse)
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

#  (!requireNamespace("BiocManager", quietly = TRUE))
  # install.packages("BiocManager")

# BiocManager::install("survcomp")

# library(survcomp)



SuperSmallfont= 6
Smallfont= 14
Mediumfont= 16
Largefont= 18
verylargefont = 20
pointsize= 0.7
linesize=1.1
meansize = 1.5
Margin=c(0,0,0,0)
fontsizeaxes = 12
fontsizeaxes2 = 10



# Function to extract coxme table
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



# Functions to extract survival data to plot
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



# Data import
survival_data <- read_csv("/Users/ionaa/Documents/honours project regan lab/R for survival curve newest/CoxPH & ME input.csv")
survival_data
# view(survival_data)



# Distinguishing treatment groups 4 and 5. 
# Then combining the two datasets to make up original again
Mating_Immunity_Survival_Young <-
  survival_data %>% 
  filter(Group %in%  c("Group5_Trt_Once_Proximal","Young_Control")) %>% 
  mutate(Age = paste("Young"))

Mating_Immunity_Survival_Old <-
  survival_data %>% 
  filter(!Group %in%  c("Group5_Trt_Once_Proximal","Young_Control")) %>% 
  mutate(Age = paste("Aged"))

survival_data <- rbind(Mating_Immunity_Survival_Old, Mating_Immunity_Survival_Young)



# removing group 3 controls
survival_data <-
  survival_data %>% 
  filter(!Group == "Group3_Control_Twice_Proximal") %>% 
  droplevels()


# Selecting column headings
Mating.stat <- c("Virgin", "Once", "Twice")
Never.mate<- c("Works", "FALSE", "FALSE")
Proximal.mate <- c("FALSE", "Works", "Works")
Distal.mate <- c("FALSE", "Works", "Works")
Illustrative.grid <- data.frame(Mating.stat, Never.mate,Proximal.mate,Distal.mate)


# Making a survival object
attach(survival_data)
Surv_Object <- Surv(time = Time_to_death_, event = Death)



# Unordering factors
survival_data$Mating_Status <- factor(survival_data$Mating_Status, ordered = FALSE)
survival_data$Group <- factor(survival_data$Group, ordered = FALSE)
survival_data$Infection_Status <- factor(survival_data$Infection_Status, ordered = FALSE)
survival_data$Age <- factor(survival_data$Age, ordered = FALSE)



# CoxPH analysis:
# First, set the reference leve for each factor, then run. Analysis output printed in console.
survival_data$Mating_Status <- relevel(survival_data$Mating_Status, ref = "Once")
survival_data$Infection_Status <- relevel(survival_data$Infection_Status, ref = "Sham")
survival_data$Age <- relevel(survival_data$Age, ref = "Aged")

CoxPH_analysis_All <- coxph(Surv_Object ~ Infection_Status + 
                              Mating_Status + 
                              Age + 
                              Infection_Status * Age + 
                              # Infection_Status * Mating_Status + 
                              # Mating_Status * Age + 
                              cluster(Vial) , 
                            data = survival_data)
CoxPH_analysis_All



# CoxME analysis: 
CoxME <- coxme(Surv_Object ~ Mating_Status + Infection_Status + Age +  (1|Vial),
               data = survival_data)
print(CoxME)



# ANOVA to compare models: 
Fit_Compare <- anova(CoxPH_analysis_All, CoxME)
print(Fit_Compare)



# Making survival curves from CoxMe data: 
# First, make a CoxME object, 
# Then use Tukey comparison to compare differences between Sham and Infected
# letter combinations represent significantly different groups, 
# These letters will be visible on downstream plots
CoxmeGroup <- coxme(Surv(Time_to_death_,Death) ~ Group + (1|Vial), data= survival_data)

multcomp = glht(CoxmeGroup, linfct=mcp(Group="Tukey"))
Comp = cld(multcomp)
Comp

unlist(Comp$mcletters$Letters)%>%
  kable(col.names = "Sign.group") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = F)


# The plot the data: 
survdata <- survfit(Surv_Object~Group, data=survival_data)

# pairwise comparison of groups: 
Pairwise <- pairwise_survdiff(Surv(Time_to_death_, Death) ~ Group, 
                              data = survival_data)
Pairwise


time_check = sort(unique(survival_data$Time_to_death_))

toplot <- ggplotprep2(survdata, times=c(0,time_check))

Name= ""
Limits=c("Young_Control",
         "Aged_Control",
         "Group1_Trt_Virgin",
         "Group2_Trt_Once_Distal",
         "Group3_Trt_Twice_Proximal",
         "Group4_Trt_Once_Proximal",
         "Group5_Trt_Once_Proximal")

n1 = survdata$n[1]
n2 = survdata$n[2]
n3 = survdata$n[3]
n4 = survdata$n[4]
n5 = survdata$n[5]
n6 = survdata$n[6]
n7 = survdata$n[7]


Labels = c(paste("Young Control (n=",n7,") d",sep=""),
           paste("Aged Control (n=",n1,") cd",sep=""), 
           paste("Virgin (n=",n2,") ab",sep=""),
           paste("Once-Mated Distal (n=",n3,") ab",sep=""),
           paste("Twice-Mated (n=",n4,") ab",sep=""),
           paste("Once-Mated Proximal (n=",n5,") bc",sep=""),
           paste("Young-Mated Proximal (n=",n6,") a",sep="")
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
                      values=c("Grey0","Grey60", "#00B050","#FF0000", "#1D78FF","#FF9900", "#9900CC", "#9900CC", "#9900CC"),
                      labels=Labels)+
  scale_linetype_manual(Name,
                        limits=Limits,
                        values=c("dotted","dotted", "solid","solid","solid","solid","solid","dotted","solid"),
                        labels=Labels)+
  scale_x_continuous("Days post-infection",
                     limits=c(0, 6),
                     breaks=c(seq(0,6,by=1)))+
  scale_y_continuous("Proportion of survivors",
                     limits=c(0, 1),breaks=c(0,0.2,0.4,0.6,0.8,1))+
  theme(title =element_text(size=Mediumfont, face='bold'),
        plot.title = element_text(hjust = 0.5,size = verylargefont),
        panel.background = element_blank(),
        strip.text.x = element_text(size =Mediumfont, colour = "black",face="italic"),
        strip.text.y = element_text(size =Mediumfont, colour = "black",face="italic"),
        axis.text.x = element_text(size=Mediumfont,colour="black"),
        axis.text.y = element_text(size=Mediumfont,colour="black"),
        axis.title.x = element_text(size=Largefont,colour="black"),
        axis.title.y = element_text(size=Largefont,colour="black"), 
        axis.line.x = element_line(colour="black",size=0.75),
        axis.line.y = element_line(colour="black",size=0.75),
        axis.ticks.x = element_line(size = 0.75),
        axis.ticks.y = element_line(size = 0.75),
        legend.direction = "vertical", 
        legend.box = "horizontal",
        legend.position = c(0.25,0.35),
        legend.key.height = unit(1.4, "cm"),
        legend.key.width= unit(2.6, "cm"),
        legend.title = element_text(face="italic",size=Mediumfont), 
        legend.key = element_rect(colour = 'white', fill = "white", linetype='dashed'),
        legend.text = element_text(size=Mediumfont),
        legend.background = element_rect(fill=NA),
        plot.margin = unit(c(0,0,1.2,0), "cm"))+
  guides(shape=guide_legend(ncol=1),
         fill=guide_legend(ncol=1),
         colour=guide_legend(ncol=1))

p1



# Extract info from CoxME model to plot hazard ratio
CoxmeGroup <- coxme(Surv(Time_to_death_,Death) ~ Group + (1|Vial), data= survival_data)

tab_res_surv_int = extract_coxme_table(CoxmeGroup)

rownames(tab_res_surv_int) <- sub("GroupG","G",rownames(tab_res_surv_int))
rownames(tab_res_surv_int) <- sub("GroupY","Y",rownames(tab_res_surv_int))

tab_res_surv_int$Treatment = as.factor(rownames(tab_res_surv_int))
tab_res_surv_int$beta = as.numeric(as.character(tab_res_surv_int$beta))
tab_res_surv_int$se = as.numeric(as.character(tab_res_surv_int$se))
tab_res_surv_int$z = as.numeric(as.character(tab_res_surv_int$z))
tab_res_surv_int$p = as.numeric(as.character(tab_res_surv_int$p))
tab_res_surv_int$HazardRatio = exp(tab_res_surv_int$beta)
tab_res_surv_int[nrow(tab_res_surv_int)+1,] = NA

tab_res_surv_int[7,1] = 0
tab_res_surv_int[7,2] = 0
tab_res_surv_int$Treatment = as.character(tab_res_surv_int$Treatment)
tab_res_surv_int[7,5] = "Aged_Control"
tab_res_surv_int$Treatment = as.factor(tab_res_surv_int$Treatment)
tab_res_surv_int[7,6] = 1




# Plot log(hazard ratio)
Limits=c("Aged_Control",
         "Group1_Trt_Virgin",
         "Group2_Trt_Once_Distal",
         "Group3_Trt_Twice_Proximal",
         "Group4_Trt_Once_Proximal",
         "Young_Control",
         "Group5_Trt_Once_Proximal")

p2=
  ggplot(data=tab_res_surv_int, aes(x=Treatment, y=beta)) + 
  geom_errorbar(aes(ymin=beta-se, ymax=beta+se, colour=Treatment),
                width=.4,  size=1, show.legend=FALSE,)+
  scale_color_manual("Treatment", breaks=c(1,2,3,4),values=c("Grey60", "#00B050","#FF0000", "#1D78FF","#FF9900", "#9900CC", "Grey0"))+
  geom_point(shape="circle", 
             aes(colour=Treatment,
                 fill=Treatment),
             alpha=0.4, stat="identity", show.legend=FALSE,  size=5) +
  scale_shape_manual(values=seq(0,9))+
  scale_y_continuous("log(Hazard ratio) \u00B1se",
                     limits = NULL)+
  geom_hline(yintercept = 0,colour="black",linetype=4)+
  ggtitle("log(Hazard Ratio) Relative to Aged Controls")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust=0.5, size = Largefont),
        axis.title.x = element_text(size=Largefont,colour="black"),
        axis.title.y = element_text(size=Largefont,colour="black"),
        axis.line.x = element_line(colour="black",size=0.75),
        axis.line.y = element_line(colour="black",size=0.75),
        axis.ticks.x = element_line(size = 0.75),
        axis.ticks.y = element_line(size = 0.75),
        axis.text.x = element_text(size=Mediumfont,colour="black",angle=35,hjust=1),
        axis.text.y = element_text(size=Smallfont,colour="black"),
        plot.margin = unit(Margin, "cm"),
        legend.direction = "horizontal", 
        legend.box = "vertical",
        legend.position = "bottom",
        legend.key.size = unit(4, "line"),
        #legend.key.height = unit(0.4, "cm"),
        #legend.key.width= unit(0.8, "cm"),
        legend.title = element_text(face="italic",size=Smallfont), 
        legend.key = element_rect(colour = 'white', fill = "white", linetype='dashed'),
        legend.text = element_text(size=Smallfont),
        legend.background = element_rect(fill=NA))+
  guides(shape=guide_legend(ncol=5),
         fill=guide_legend(ncol=5),
         colour=guide_legend(ncol=5))

p2



# Export Figures
ggsave("/Users/ionaa/Documents/honours project regan lab/R for survival curve newest/survival plot.png", p1, width = 15, height = 10, dpi = 500)

ggsave("/Users/ionaa/Documents/honours project regan lab/R for survival curve newest/hazard ratio plot.png", p2, width = 15, height =10, dpi = 500)


# print a more extensive summary of the coxPH data:
summary(CoxPH_analysis_All)



## check assumptions of CoxPH model: 

# Check proportional hazards assumption first
# if all covariates as well as global test are not statistically significant,
# proportional hazards can be assumed
test.ph <- cox.zph(CoxPH_analysis_All)
test.ph

# plot graphs of the scaled Schoenfeld residuals against the transformed time for each covariate: 
# systematic departures from a horizontal line are indicative of non-proportional hazards, 
# since proportional hazards assumes that estimates β1,β2,β3 do not vary much over time.
ggcoxzph(test.ph)

# Next, check influential observations or outliers by visualising dfbeta values and/or deviance residuals: 
ggcoxdiagnostics(CoxPH_analysis_All, type = "dfbeta", linear.predictions = FALSE)

# deviance residuals should be roughly symmetrically distributed about zero, 
# with a standard deviation of one
ggcoxdiagnostics(CoxPH_analysis_All, type = "deviance", linear.predictions = FALSE)

# No need to check non-linearity, as all covariates are categorical
# (None are continuous)


