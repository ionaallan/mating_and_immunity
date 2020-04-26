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



# Data Import 


# path.to.data <- "/Users/ionaa/Documents/honours project regan lab/MK Cox analysis download 2"

# d <- list()
# path <- list()
# model <- list()
# for(f in list.files(path=path.to.data,pattern="*.csv$",recursive=T,full.names=T)) {
  # nom <- gsub(".*/(.*).csv","\\1",f)    
  # cat(nom,"\n")
  # path[[nom]] <- gsub("(.*)/.*csv","\\1/",f)
  # d[[nom]] <- read.table(f,header=T,sep=",",dec=",")
  
#}

Mating_Immunity_Survival= d[["CoxPH input"]]

Mating_Immunity_Survival <- read_csv("/Users/ionaa/Documents/honours project regan lab/MK Cox analysis download 2/CoxPH input (1).csv")
Mating_Immunity_Survival
# view(Mating_Immunity_Survival)

## Datatable Visualisation

#datatable(Mating_Immunity_Survival, rownames = FALSE, filter="top", options = list(pageLength = 5, scrollX=T) )



## Data organisation 

  
#  Here I am modifying the dataset, as opposed to the spreadsheet, so that we have a column for age.
# As you can see, without such a ditinction between groups 4 and 5, R will consider them the same
# bar their group numbers.




Mating_Immunity_Survival_Young <-
  Mating_Immunity_Survival %>% 
  filter(Group %in%  c("Group5_Trt_Once_Proximal","Group5_Control_Once_Proximal")) %>% 
  mutate(Age = paste("Young"))

Mating_Immunity_Survival_Old <-
  Mating_Immunity_Survival %>% 
  filter(!Group %in%  c("Group5_Trt_Once_Proximal","Group5_Control_Once_Proximal")) %>% 
  mutate(Age = paste("Aged"))

## Combine the two datasets to make up original again 
Mating_Immunity <- rbind(Mating_Immunity_Survival_Old, Mating_Immunity_Survival_Young)


  
 # Now we want to select Group, Mating_status, Infection_status, age, time to death and Censor. We also want to drop the group3 controls because of contamination.
# We don't select mating to infection distance because it interupts the model as we can see, the empty grids in this table represent spaces where the Cox Model can't compute term interactions, eg the effect of being a virgin and having mated distally to infection (its not possible).



Mating.stat <- c("Virgin", "Once", "Twice")
Never.mate<- c("Works", "FALSE", "FALSE")
Proximal.mate <- c("FALSE", "Works", "Works")
Distal.mate <- c("FALSE", "Works", "Works")
Illustrative.grid <- data.frame(Mating.stat, Never.mate,Proximal.mate,Distal.mate)



## Data select, remove group 3 control

Mating_Immunity <-
  Mating_Immunity %>% 
  filter(!Group == "Group3_Control_Twice_Proximal") %>% 
  droplevels()



  
  # Survival Analysis 
  ## Making a Survival Object 


attach(Mating_Immunity)
Surv_Object <- Surv(time = Time_to_death_, event = Censor)



## Feed the Survival Object into a CoxPH formula 


  
#  We can compare the fits, or outputs, of two models, a multivariate CoxPH analysis 
# that accounts for random effects and a Cox mixed effects (coxme) model. First we
# use the function to assign the reference level.


# Unordering factors: 
Mating_Immunity$Mating_Status <- factor(Mating_Immunity$Mating_Status, ordered = FALSE)
Mating_Immunity$Group <- factor(Mating_Immunity$Group, ordered = FALSE)
Mating_Immunity$Age <- factor(Mating_Immunity$Age, ordered = FALSE)
Mating_Immunity$Infection_Status <- factor(Mating_Immunity$Infection_Status, ordered = FALSE)

## Set the reference level for each factor

Mating_Immunity$Mating_Status <- relevel(Mating_Immunity$Mating_Status, ref = "Once")
Mating_Immunity$Infection_Status <- relevel(Mating_Immunity$Infection_Status, ref = "Sham")
Mating_Immunity$Age <- relevel(Mating_Immunity$Age, ref = "Aged")


CoxPH_analysis_All <- coxph(Surv_Object ~ Infection_Status + 
                              Mating_Status + 
                              Age + 
                              Infection_Status * Age + 
                              Infection_Status * Mating_Status + 
                              Mating_Status * Age + 
                              cluster(Vial) , data = Mating_Immunity)
CoxPH_analysis_All





## Coxme Analysis 



CoxME <- coxme(Surv_Object ~ Mating_Status + Infection_Status + Age +  (1|Vial),
               data = Mating_Immunity)

summary(CoxME)

print(CoxME)




## Lets use ANOVA to compare the model fits.

Fit_Compare <- anova(CoxPH_analysis_All, CoxME)

print(Fit_Compare)



  
 # You should include the Coxme output in your report.


print(CoxME)




# Survival

  
  # **Survival Curves Prep** 
  
# Here we create a coxme object on the first line od code. 
# The subsequent code draws on "Tukey" comparison to compare differences between sham and infected, 
# while created letter combinations to represent significantly different groups. 
# These letters will be visible on downstream plots. 



  
  ## Survival to *P. rettgeri*
  
  


CoxmeGroup <- coxme(Surv(Time_to_death_,Censor) ~ Group + (1|Vial), data= Mating_Immunity)

multcomp = glht(CoxmeGroup, linfct=mcp(Group="Tukey"))
Comp = cld(multcomp)

unlist(Comp$mcletters$Letters)%>%
  kable(col.names = "Sign.group") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = F)




#  **Survival Curves** 
  
 # Here we use `survfit` to create a survival plot and later use the function at the beginning of this code "ggplotprep2" to extract our survival matrix from that, which is called topplot.

# We provide the limits, or levels of our group (ie. group names). And we can then call those when it comes to plotting. 

# ggplot` is a widely used ploting package in r and it has many customisation options. I have reandomly chosen colours for the groups but you can change this, colour palette can be found [here](http://sape.inf.usi.ch/quick-reference/ggplot2/colour).

# p1 represents the survival curves while p2 is the comparison of the log of the harzard ratios between the groups. These were extracted from the survival curve. 



## plot
survdata <- survfit(Surv_Object~Group, data=Mating_Immunity)

# pairwise comparison of groups: 
Pairwise <- pairwise_survdiff(Surv(Time_to_death_, Censor) ~ Group, 
                              data = Mating_Immunity)
Pairwise

###time checked
time_check = sort(unique(Mating_Immunity$Time_to_death_))

toplot <- ggplotprep2(survdata, times=c(0,time_check))

Name= ""
Limits=c("Group1_Control_Virgin",
         "Group1_Trt_Virgin",
         "Group2_Control_Once_Distal",
         "Group2_Trt_Once_Distal",
         "Group3_Trt_Twice_Proximal",
         "Group4_Control_Once_Proximal",
         "Group4_Trt_Once_Proximal",
         "Group5_Control_Once_Proximal",
         "Group5_Trt_Once_Proximal")

n1 = survdata$n[1]
n2 = survdata$n[2]
n3 = survdata$n[3]
n4 = survdata$n[4]
n5 = survdata$n[5]
n6 = survdata$n[6]
n7 = survdata$n[7]
n8 = survdata$n[8]
n9 = survdata$n[9]


Labels = c(paste("Virgin_Cntrl (n=",n1,") bcd",sep=""), 
           paste("Virgin_Trt (n=",n2,") ab",sep=""),
           paste("Once_Mated_Distal_Cntrl (n=",n3,") bcd",sep=""), 
           paste("Once_Mated_Distal_Trt (n=",n4,") a",sep=""),
           paste("Twice_Mated_Trt (n=",n5,") ac",sep=""),
           paste("Once_Mated_Proximal_Cntrl (n=",n6,") cd",sep=""),
           paste("Once_Mated_Proximal_Trt (n=",n7,") ac",sep=""),
           paste("Young_Mated_Proximal_Cntrl (n=",n8,") d",sep=""),
           paste("Young_Mated_Proximal_Trt (n=",n9,") a",sep="")
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
                      values=c("#00B050","#00B050",
                               "#FF0000", "#FF0000",
                               "#1D78FF",
                               "#FF9900", "#FF9900",
                               "#9900CC", "#9900CC"),
                      labels=Labels)+
  scale_linetype_manual(Name,
                        limits=Limits,
                        values=c("dotted","solid",
                                 "dotted","solid",
                                 "solid",
                                 "dotted","solid",
                                 "dotted","solid"),
                        labels=Labels)+
  scale_x_continuous("Days post-injection",
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
         col=guide_legend(ncol=1))

p1

#plot hazard ratio
## Extract info from coxme model

CoxmeGroup <- coxme(Surv(Time_to_death_,Censor) ~ Group + (1|Vial), data= Mating_Immunity)

tab_res_surv_int = extract_coxme_table(CoxmeGroup)

rownames(tab_res_surv_int) <- sub("GroupG","G",rownames(tab_res_surv_int))

tab_res_surv_int$Treatment = as.factor(rownames(tab_res_surv_int))
tab_res_surv_int$beta = as.numeric(as.character(tab_res_surv_int$beta))
tab_res_surv_int$se = as.numeric(as.character(tab_res_surv_int$se))
tab_res_surv_int$z = as.numeric(as.character(tab_res_surv_int$z))
tab_res_surv_int$p = as.numeric(as.character(tab_res_surv_int$p))
tab_res_surv_int$HazardRatio = exp(tab_res_surv_int$beta)
tab_res_surv_int[nrow(tab_res_surv_int)+1,] = NA

tab_res_surv_int[9,1] = 0
tab_res_surv_int[9,2] = 0
tab_res_surv_int$Treatment = as.character(tab_res_surv_int$Treatment)
tab_res_surv_int[9,5] = "Group1_Control_Virgin"
tab_res_surv_int$Treatment = as.factor(tab_res_surv_int$Treatment)
tab_res_surv_int[9,6] = 1

# Plot log(hazard ratio) - better than hazard ratio has you can get a standard deviation of it. (note that one could show the HR and calculate the 95%CI)

Limits=c("Group1_Control_Virgin",
         "Group1_Trt_Virgin",
         "Group2_Control_Once_Distal",
         "Group2_Trt_Once_Distal",
         "Group3_Trt_Twice_Proximal",
         "Group4_Control_Once_Proximal",
         "Group4_Trt_Once_Proximal",
         "Group5_Control_Once_Proximal",
         "Group5_Trt_Once_Proximal")

p2=
  ggplot(tab_res_surv_int, aes(x=Treatment, y=beta)) + 
  geom_errorbar(aes(ymin=beta-se, ymax=beta+se, col=Treatment),
                width=.4,  size=1, show.legend=FALSE)+
  scale_color_manual("Treatment", breaks=c(1,2,3,4),values=c("#7ADEA7", "#00B050",
                                                             "#FC9292", "#FF0000", 
                                                             "#1D78FF",
                                                             "#FCD79F", "#FF9900",
                                                             "#B978CF", "#9900CC")) +
  geom_point(shape="circle",aes(col=Treatment, fill=c("empty","solid",
                                                 "empty","solid",
                                                 "solid",
                                                 "empty","solid",
                                                 "empty","solid")),
             alpha=0.4, stat="identity", show.legend=FALSE,  size=3) +
  scale_shape_manual(values=seq(0,9))+
  scale_y_continuous("log(Hazard ratio) \u00B1se",
                     limits = NULL)+
  geom_hline(yintercept = 0,colour="black",linetype=4)+
  ggtitle("log(Hazard Ratio) relative to control Virgins")+
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust=0.5, size = Largefont),
        axis.title.x = element_text(size=Largefont,colour="black"),
        axis.title.y = element_text(size=Largefont,colour="black"),
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




# Export figures


  
#  Use `ggsave` to save and export the figures to a specified folder. This is done by giving the file name as a character string while specifying the file type - "users/R/Figure_output/Plot1.png". You can play around with dimensions and resolution as well.




ggsave("/Users/ionaa/Documents/honours project regan lab/MK Cox analysis download 2/survival curves.png", p1, width = 15, height = 10, dpi = 300)

ggsave("/Users/ionaa/Documents/honours project regan lab/MK Cox analysis download 2/hazard ratio plot.png", p2, width = 19, height = 10, dpi = 300)



summary(CoxPH_analysis_All)


