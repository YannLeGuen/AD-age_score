################################################################################
## Simulate power for AD phenotype score versus case-control association test ##
## *Simulation for determining Type I Error rate ###############################
################################################################################
##################################################################### 02/17/2021
################################################################################
## Libraries and Input

## Libraries
library(simpleboot)
library(ggplot2)
library(survival)

# # This script is meant to be run on a cluster
# For example the submission script below can be used on a cluster with SLURM
# ##SBATCH --job-name=SimuTypeIError
# #SBATCH --array=1-60
# #SBATCH --time=2:00:00
# #SBATCH --ntasks=1
# #SBATCH --cpus-per-task=2
# #SBATCH --mem-per-cpu=6G
# #SBATCH --partition=XXX
# 
# # ml load R
# # Rscript statistical_power_simulation_Type_I_Error.R ${SLURM_ARRAY_TASK_ID}

## Enable command line arguments
args = commandArgs(TRUE)
job_index = as.numeric(args[1])

## Output folder
outdir<- # set the output directory

## Parameter settings: Odds ratios (OR), minor allele frequencies (MAF), sample
## sizes for cases and controls (n), median age differences cases-controls (ageDiff),
## Prevalence estimates at age 100 (PREV100)
df <- data.frame(Doubles=double(),Doubles=double(),Doubles=double(),Doubles=double(),
                 Doubles=double())
colnames(df) <- c("OR","MAF","N", 'ageDiff',"PREV100")
ORs <- c(1)
MAFs <- c(0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45)
ageDiffs <- c(0, 5, 10)
Ns <- c(1000, 5000)
PREV100s <- c(0.80)

## Parameter settings: create combinations
c <- 1
for (id1 in 1:length(ORs)) {
  for (id2 in 1:length(MAFs)) {
    for (id3 in 1:length(Ns)) {
      for (id4 in 1:length(ageDiffs)) {
        for (id5 in 1:length(PREV100s)) {
          df[c,] <- c(ORs[id1], MAFs[id2], Ns[id3], ageDiffs[id4], PREV100s[id5])
          c <- c + 1
        }
      }
    }
  }
}

## Parameter settings: set alpha values per sample size
MAF <- df[job_index, 'MAF']; OR <- df[job_index, 'OR']; n <- df[job_index, 'N']; 
ageDiff <- df[job_index, 'ageDiff']; PREV100 <- df[job_index, 'PREV100'] 
if (n == 1000) {
  alpha1 = 0.10; alpha2 = 0.05; alpha3 = 1;# trend, nominal
} else if (n == 5000) {
  alpha1 = 1e-5; alpha2 = 5e-7; alpha3 = 5e-8; # suggestive, exome-wide
} else {
  alpha1 = 1e-5; alpha2 = 5e-8; alpha3 = 1;# suggestive, genome-wide
}

################################################################################
## Functions 

# Cauchy function for aggregating p-values
Get.cauchy<-function(p){
  p[p>0.99]<-0.99
  is.small<-(p<1e-16) & !is.na(p)
  is.regular<-(p>=1e-16) & !is.na(p)
  temp<-rep(NA,length(p))
  temp[is.small]<-1/p[is.small]/pi
  temp[is.regular]<-as.numeric(tan((0.5-p[is.regular])*pi))
  
  cct.stat<-mean(temp,na.rm=T)
  if(is.na(cct.stat)){return(NA)}
  if(cct.stat>1e+15){return((1/cct.stat)/pi)}else{
    return(1-pcauchy(cct.stat))
  }
}

################################################################################
## A population draw and related power simulations will be performed 10 times

p<-c() # variable to store p-values over different simulations 
d<-c() # variable to store case-control age differences over different simulations 
CASE.AGE<-c(); CNTRL.AGE<-c() # variables to store age distributions over different simulations 
for (i in 1:10){
  
  ##############################################################################
  ## Generate population data for Alzheimer's disease case-control status
  
  # Supporting references:
  # Van Der Lee et al 2018; APOE genotype prevalence estimates across age 60-100y
  
  # Number of individuals in the population considered
  N<-10000000
  
  # Get P.age
  P.age <- rnorm(N,75,8)
  P.age <- P.age[which(P.age>60 & P.age<100)]
  
  # Number of individual in the population considered
  N<-length(P.age) 
  
  # Simulate genotype (frequencies)
  beta.G <- log(OR)
  P.G<-rbinom(N,2,MAF)
  
  # Set prevalence at baseline age (60y) and intercept
  Prev60<-0.01; 
  b0_pheno<-log(Prev60/(1-Prev60))
  
  # Set prevalence estimates at age 100 
  Prev100 <- PREV100
  
  # Set beta.age assuming a linear increase of prevalence with age
  beta.age<-(log(Prev100/(1-Prev100))-log(Prev60/(1-Prev60)))/(100-60)
  
  ##############################################################################
  ## Initial case-control draw from population
  beta.noise <- log(1.05)
  P.noise <- rnorm(N, 0, 1)
  # Define the model
  mu<-b0_pheno+(P.age-60)*beta.age+P.G*beta.G+beta.noise*P.noise
  # Draw the AD status for each individual
  P.y<-as.matrix(rbinom(N, 1, exp(mu)/(1+exp(mu))))
  # Case-control indexing
  index.cases = which(P.y == 1)
  index.controls = which(P.y == 0)
  
  # df <- as.data.frame(cbind(P.age[index.controls],P.y[index.controls]))
  # colnames(df) <- c("age","y")
  # p<-ggplot(data=df, aes(x=age)) +
  #   geom_histogram(binwidth=1)
  # 
  # df <- as.data.frame(cbind(P.age[index.cases],P.y[index.cases]))
  # colnames(df) <- c("age","y")
  # p<-ggplot(data=df, aes(x=age)) +
  #   geom_histogram(binwidth=1)
  
  ##############################################################################
  # Determine age modifier term to enforce age differences between cases and controls
  
  # basic settings 
  # Note that ageMod is set to start of with a negative median age difference for controls minus cases
  ageMod <- -70; ageModInc <- 0.5; c <- 1
  # allow start of while loop by setting a negative age difference for controls minus cases
  aca <- 79; aco <- 78
  
  # while loop to identify optimal ageMod value for requested ageDiff value
  while ((aco-aca)<(ageDiff)) {
    
    # set starting condition
    if (c==1) {
      ageMod <- ageMod
      ageModInc <- 1;
    } else {
      ageMod <- ageMod + ageModInc
    }
    
    # if ageMod reaches extreme values (100 is max and produces NA) then adapt increment
    if (ageMod==95) {
      ageModInc <- 0.25
    }
    if (ageMod==97.5) {
      ageModInc <- 0.125
    }
    if (ageMod==99) {
      ageModInc <- 0.01
    }
    
    # Increase probability to be selected among controls if older
    Pselected0 <- 0.1; b0_sel <- log(Pselected0/(1-Pselected0))
    Pselected100 <- 0.90; offset <- 85 
    beta_sel <-(log(Pselected100/(1-Pselected100))-log(Pselected0/
                                                         (1-Pselected0)))/(100-ageMod) 
    mu_sampling_controls <- (P.age[index.controls]-offset)*beta_sel
    
    # Select the controls
    N.controls = length(index.controls)
    # Draw the controls
    P.yc<-as.matrix(rbinom(N.controls, 1, exp(mu_sampling_controls+b0_sel)/
                             (1+exp(mu_sampling_controls+b0_sel))))
    index.selected.controls = which(P.yc == 1)
    P.age.controls = P.age[index.controls][index.selected.controls]
    
    # extract median age for controls
    aco <- summary(P.age.controls)[3] 
    
    # # visualize distribution
    # df <- as.data.frame(cbind(P.age[index.controls][index.selected.controls],
    #                           P.y[index.controls][index.selected.controls]))
    # colnames(df) <- c("age","y")
    # p<-ggplot(data=df, aes(x=age)) +
    #   geom_histogram(binwidth=1)
    
    # Increase probability to be selected among cases if younger
    Pselected0 <- 0.9; b0_sel <- log(Pselected0/(1-Pselected0))
    Pselected100 <- 0.10; offset <- 75
    beta_sel <-(log(Pselected100/(1-Pselected100))-log(Pselected0/
                                                         (1-Pselected0)))/(100-ageMod) 
    mu_sampling_cases <- (P.age[index.cases]-offset)*beta_sel
    
    # Select the cases
    N.cases = length(index.cases)
    # Draw the cases
    P.yc<-as.matrix(rbinom(N.cases, 1, exp(mu_sampling_cases+b0_sel)/
                             (1+exp(mu_sampling_cases+b0_sel))))
    index.selected.cases = which(P.yc == 1)
    P.age.cases = P.age[index.cases][index.selected.cases]
    
    # extract median age for cases
    aca <- summary(P.age.cases)[3] 
    
    # # visualize distribution
    # df <- as.data.frame(cbind(P.age[index.cases][index.selected.cases],
    #                           P.y[index.cases][index.selected.cases]))
    # colnames(df) <- c("age","y")
    # p<-ggplot(data=df, aes(x=age)) +
    #   geom_histogram(binwidth=1)
    
    # iteration counter
    c <- c + 1
  }
  
  # collect output
  AGE_MOD <- ageMod; AGE_DIFF <- aco-aca
  
  ##############################################################################
  # Second control draw using automatically generated age_mod term
  
  # Increase probability to be selected among controls if older
  Pselected0 <- 0.1; b0_sel <- log(Pselected0/(1-Pselected0))
  Pselected100 <- 0.90; offset <- 85 
  beta_sel <- (log(Pselected100/(1-Pselected100))-log(Pselected0/
                                                        (1-Pselected0)))/(100-AGE_MOD) 
  mu_sampling_controls <- (P.age[index.controls]-offset)*beta_sel
  
  # Select the controls
  N.controls = length(index.controls)
  # Draw the controls
  P.yc<-as.matrix(rbinom(N.controls,1, exp(mu_sampling_controls+b0_sel)/
                           (1+exp(mu_sampling_controls+b0_sel))))
  index.selected.controls = which(P.yc == 1)
  
  ##############################################################################
  # Second case draw using automatically generated age_mod term
  
  # Increase probability to be selected among cases if younger
  Pselected0 <- 0.9; b0_sel <- log(Pselected0/(1-Pselected0))
  Pselected100 <- 0.10; offset <- 75
  beta_sel <-(log(Pselected100/(1-Pselected100))-log(Pselected0/
                                                       (1-Pselected0)))/(100-AGE_MOD) 
  mu_sampling_cases <- (P.age[index.cases]-offset)*beta_sel
  
  # Select the cases
  N.cases = length(index.cases)
  # Draw the cases
  P.yc<-as.matrix(rbinom(N.cases,1, exp(mu_sampling_cases+b0_sel)/
                           (1+exp(mu_sampling_cases+b0_sel))))
  index.selected.cases = which(P.yc == 1)
  
  ##############################################################################
  ## Aggegrate data for power simulations
  
  P.age.cases = P.age[index.cases][index.selected.cases]
  P.age.controls = P.age[index.controls][index.selected.controls]
  
  P.G.cases <- P.G[index.cases][index.selected.cases]
  P.G.controls <- P.G[index.controls][index.selected.controls]
  
  P.y.cases <- P.y[index.cases][index.selected.cases]
  P.y.controls <- P.y[index.controls][index.selected.controls]
  
  # Merge cases and controls
  P.age <- c(P.age.cases, P.age.controls)
  P.G <- c(P.G.cases, P.G.controls)
  P.y <- c(P.y.cases, P.y.controls)
 
  ##############################################################################
  ## Simulate power for 100-fold case-control random sampling
  
  for(k in 1:100){
    
    # random case sampling
    index.case <- sample(which(P.y==1),n); Case.y <- P.y[index.case];
    Case.age <- P.age[index.case]; Case.G <- P.G[index.case];
    
    # random control sampling
    index.control <- sample(which(P.y==0),n); Control.y <- P.y[index.control];
    Control.age <- P.age[index.control]; Control.G <- P.G[index.control];
    
    # aggegate data
    y<-c(Case.y,Control.y); age<-c(Case.age,Control.age); G<-c(Case.G,Control.G); 
    
    # m1. logistic regression without adjusting for age
    fit.unadjust<-summary(glm(y~G,family='binomial'))
    
    # m2. logistic regression adjusting by age
    fit.age_adjust<-summary(glm(y~G+age,family='binomial'))
    
    # create 1st phenotype score that weighs case status by age 
    Case.score1<--log((Case.age-59.5)/(100.5-59.5))+0.5
    Control.score1<-log(1-(Control.age-59.5)/(100.5-59.5))-0.5
    y.score1<-c(Case.score1,Control.score1)
    
    # m3. linear regression on the first phenotype score
    fit.score1<-summary(lm(y.score1~G))
    # m4. bootstrapped linear regression on the first phenotype score
    fit.boot1<-summary(lm.boot(lm(y.score1~G),R=100))
    
    # m5. Cox regression
    tdf <- list(time=age, status=y, geno=G)
    fit.cox1 <- summary(coxph(Surv(time, status) ~ geno, tdf))
    
    # Gather the p-values from each model
    temp.p<-c(fit.unadjust$coefficients[2,4],
              fit.age_adjust$coefficients[2,4],
              fit.score1$coefficients[2,4],
              as.double(2*pnorm(abs(fit.boot1$orig.lm$coefficients[2]/
                                    fit.boot1$stdev.params[2]),lower.tail=F)),
              fit.cox1$coefficients[1,5])
    
    # Add Cauchy p-value aggregate 
    temp.p <- append(temp.p,Get.cauchy(temp.p[c(1:5)])) # m6. evaluate with regard to first score
    
    # Append p-values over iterations (10-fold population by 100-fold random sampling)
    p<-rbind(p,temp.p)
    
    # Append case-control age difference over iterations (10-fold population by 100-fold random sampling)
    d<-rbind(d,AGE_DIFF)
    
    ############################################################################
    
    # visualize distribution
    df <- data.frame("age"=Case.age,"y"=1)
    pp<-ggplot(data=df, aes(x=age)) +
      geom_histogram(binwidth=1,breaks=c(60:100))
    pg <- ggplot_build(pp)
    
    # collect data
    CASE.AGE<-rbind(CASE.AGE,pg$data[[1]]$count)
    
    # visualize distribution
    df <- data.frame("age"=Control.age,"y"=1)
    pp<-ggplot(data=df, aes(x=age)) +
      geom_histogram(binwidth=1,breaks=c(60:100))
    pg <- ggplot_build(pp)
    
    # collect data
    CNTRL.AGE<-rbind(CNTRL.AGE,pg$data[[1]]$count)
    
  }
}

################################################################################
## Collect current results

results1 <- cbind(rbind(apply(p<alpha1, 2, mean, na.rm=TRUE)), mean(d))
results2 <- cbind(rbind(apply(p<alpha2, 2, mean, na.rm=TRUE)), mean(d))

colnames(results1) <- c('Log_reg', 'Log_reg_age_adj', 'Lin_reg_age_sc', 
                        'Lin_reg_boot_age_sc', "cox_reg", "Cauchy_age_sc", "age_diff")
colnames(results2) <- colnames(results1); rownames(results1) <- ""; rownames(results2) <- ""

output <- paste(outdir, "PREV", PREV100, '_MAF', MAF, '_OR', OR, '_n', n,
                '_ageDiff_', ageDiff, '_thr', alpha1, '.csv' ,sep='')
write.csv(results1, output, row.names=FALSE)
output <- paste(outdir, "PREV", PREV100, '_MAF', MAF, '_OR', OR, '_n', n,
                '_ageDiff_', ageDiff, '_thr', alpha2, '.csv' ,sep='')
write.csv(results2, output, row.names=FALSE)

if (alpha3 != 1){
  results3 <- cbind(rbind(apply(p<alpha3, 2, mean, na.rm=TRUE)), mean(d))
  colnames(results3) <- colnames(results1); rownames(results3) <- ""
  output <- paste(outdir, "PREV", PREV100, '_MAF', MAF, '_OR', OR, '_n', n,
                  '_ageDiff_', ageDiff, '_thr', alpha3, '.csv' ,sep='')
  write.csv(results3, output, row.names=FALSE)
}


# figure of age distributions
df <- as.data.frame(cbind(colMeans(CNTRL.AGE),colMeans(CASE.AGE)))/n
colnames(df) <- c("CNTRL","CASE")
ma <- max(c(max(df$CNTRL),max(df$CNTRL)))
# save figure of age distributions
output <- paste(outdir, "PREV", PREV100, '_MAF', MAF, '_OR', OR, '_n', n,
                '_ageDiff_', ageDiff, '_age_sampling_offset.png' ,sep='')
png(output,width=6,height=6,units="in",res=1200)

plot(df$CNTRL,type="l",lty=1,lwd=1,col="blue",ylim=c(0,ma+ma/10),ylab="",xlab="")
par(new=T)
plot(df$CASE,type="l",lty=1,lwd=1,col="red",ylim=c(0,ma+ma/10),ylab="density",xlab="Years above 60")
dev.off()