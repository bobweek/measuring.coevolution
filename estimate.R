rm(list=ls())
require(mvtnorm)
require(matrixcalc)
require(ggplot2)
require(rgl)
source("~/Research/ecology letters/rscripts/functions.R")

#      #############      ########################################################
#      #           #
#      # Toju data #      ########################################################
#      #           #
#      #############      ########################################################

# data for abiotic optimum
male_snout <- c(6.19,5.57,5.5,5.45,5.52,5.71,6.01,5.95,6.33,6.42,6.38)
male_snout_location <- c("Nagaoka","Kyoto","Jurinji","Taiji","Muroto","Ashizuri","Reihoku","Takahama",
                       "Yahazu","Shitoko","Hanyama")
male_snout_island <- c(rep("HNS",4),rep("SKK",2),rep("AMK",2),rep("YKI",3))
male_snout_df <- data.frame(snout=male_snout,location=male_snout_location,island=male_snout_island)

weevil_absence <- c(5.22,5.27,5.48,7.54,4.68,5.19,4.91,4.99,5.55,7.65,6.82,6.42,7.36,7.69,7.28,7.56,6.61)
weevil_absence_location <- c("Miyatsuka","Hatake","Mukai","Izuhara","Kamitsuki","Izu","Igatani","AkoN2",
                             "Tsubota","Jinoshima","Kishiku","Tomie","Kuwanoura","Ichinoura","Nakakoshiki",
                             "Kashima","Teuchi")
weevil_absence_island <- c(rep("NIJ",3),"TSM",rep("MYK",5),"JNS","FKE","FKE","KKS","KKS","NKS","SKS","SKS")
weevil_absence_df <- data.frame(pericarp=weevil_absence,location=weevil_absence_location,island=weevil_absence_island)

# core data
female_snout <- c(9.63,7.98,9.13,10.42,10.79,9.89,10.31,10.05,11.91,10.66,9.12,9.21,9.61,13.63,10.06,12.98,11.68,
                  11.48,14.54,16.92,12.99,19.48,21.11,18.28,20.59)
snout_var <- mean(c(0.72,0.61,1.02,0.83,0.87,1.17,0.74,0.69,0.59,0.91,0.66,0.93,0.94,0.70,1.87,1.114,1.87,1.85,1.55,1.61,0.46)^2)
female_body <- c(8.4,6.95,7.82,8.19,7.87,7.86,7.78,7.53,8.05,7.75,7.35,7.42,7.73,8.21,7.5,8.05,7.87,7.77,8.3,8.83,
                 8.16,9.31,9.6,8.87,9.34)
female_body_var <- mean(c(.42,.18,.55,.27,.45,.56,.51,.26,.37,.38,.32,.43,.42,.34,.5,.48,.57,.68,.73,.52,.16)^2)
pericarp_thickness <- c(5.42,4.32,3.88,6.07,4.90,6.30,6.13,6.66,7.95,7.73,6.76,6.42,7.52,11.65,7.77,12.80,
                        11.89,11.13,12.49,17.83,12.97,20.41,19.69,21.21,19.35)
pericarp_var <- mean(c(0.85,1.19,0.59,0.87,0.71,1.65,1.20,1.38,2.09,2.47,1.14,1.08,7.52,2.82,1.69,2.45,2.61,2.68,3.39,1.8,3.6,3.99,3.39,2.38,2.68)^2)
loc <- c("Nodzumi","Kamo","Nagaoka","Kutsuki","Hiratsuka","Kyoto","Jurinji","Nara","Shimonoseki",
              "Kiikatsuura","Taiji","Arafune","Kiioshima","Usa","Muroto","Ashizuri","Reihoku","Takahama",
              "Yahazu","Shitoko","Fukagawa","Hanyama","Shiratani","Kawahara","Ohko-rindoh")
isl <- c(rep("HNS",13),rep("SKK",3),rep("AMK",2),rep("YKI",7))
toju_df <- data.frame(rostrum=female_snout,body=female_body,pericarp=pericarp_thickness,location=loc,island=isl)
save(toju_df,file="~/Research/ecology letters/toju_pairs.Rda")
write.csv(toju_df,file="~/Research/ecology letters/toju_pairs.csv")

th1 <- mean(male_snout)
th2 <- mean(weevil_absence)
n_weevil   <- c(26389,41944,19167,66389,44167,23333,30278)
n_camellia <- c(1178,1334,2155,1841,3112,2389,1770)
n1 <- 1/mean(1/n_weevil)
n2 <- 1/mean(1/n_camellia)
G1 <- snout_var*mean(c(0.496,0.353,0.206,0.595))
G2 <- pericarp_var*mean(c(0.8,0.82,0.7,0.63))

par_df <- data.frame(thetas=c(th1,th2),eff_ns=c(n1,n2),species=c("C. camellia","C. japonica"))
save(par_df,file="~/Research/ecology letters/toju_par.Rda")
write.csv(par_df,file="~/Research/ecology letters/toju_par.csv")

data <- toju_df[c("rostrum","pericarp")]
N <- dim(data)[1]
m <- colMeans(data)
S <- (N-1)*var(data)/N

######
# ML #
######

ML <- ml_sol(m[1],m[2],S[1,1],S[2,2],S[1,2],n1,n2,th1,th2,G1,G2)
lnL <- likelihood(data,m,S)

##########
# null 1 #
##########

mu_n1 <- c(th1,m[2])
lnL1  <- likelihood(data,mu_n1,S)

##########
# null 2 #
##########

mu_n2 <- c(m[1],th2)
lnL2  <- likelihood(data,mu_n2,S)

# chisqrs
p1_toju <- 1-pchisq(2*(lnL-lnL1),1)
p2_toju <- 1-pchisq(2*(lnL-lnL2),1)

# storing ML sols
A1_toju <- ML$A1
A2_toju <- ML$A2
B1_toju <- ML$B1
B2_toju <- ML$B2
k_toju <- ML$k

#############
# conf ints #
#############

PTS <- 1000
cnt <- 1
solutions <- matrix(numeric(PTS*4),ncol=4,nrow=PTS)
for(cnt in 1:PTS){
  fdata <- rmvnorm(N,m,S)
  fm <- colMeans(fdata)
  fS <- (N-1)*var(fdata)/N
  ML <- ml_sol(fm[1],fm[2],fS[1,1],fS[2,2],fS[1,2],n1,n2,th1,th2,G1,G2)
  
  solutions[cnt,1] <- ML$A1
  solutions[cnt,2] <- ML$A2
  solutions[cnt,3] <- ML$B1
  solutions[cnt,4] <- ML$B2
  
}

A1_95upper_toju <- q_fct(0.975,solutions[,1])
A1_95lower_toju <- q_fct(0.025,solutions[,1])

A2_95upper_toju <- q_fct(0.975,solutions[,2])
A2_95lower_toju <- q_fct(0.025,solutions[,2])

B1_95upper_toju <- q_fct(0.975,solutions[,3])
B1_95lower_toju <- q_fct(0.025,solutions[,3])

B2_95upper_toju <- q_fct(0.975,solutions[,4])
B2_95lower_toju <- q_fct(0.025,solutions[,4])

###############
# effect size #
###############

A1_toju_noco <- 1/(2*S[1,1]*n1)
A2_toju_noco <- 1/(2*S[2,2]*n2)
mu_noco_toju <- mu(A1_toju_noco,A2_toju_noco,0,0,k_toju,n1,n2,th1,th2)
S_noco_toju  <- SIGMA(A1_toju_noco,A2_toju_noco,0,0,k_toju,n1,n2,th1,th2,G1,G2)

print("effect size for toju data")
print(mahalanobis(m,mu_noco_toju,(S+S_noco_toju)/2))
print("percent change for toju data")
print(abs(m-mu_noco_toju)/mu_noco_toju)

l <- mu_noco_toju[1]-sqrt(S_noco_toju[1,1])-3
b <- mu_noco_toju[2]-sqrt(S_noco_toju[2,2])
r <- colMeans(data)[1]+sqrt(var(data)[1,1])
t <- colMeans(data)[2]+sqrt(var(data)[2,2])+3
pts <- NULL
for(i in seq(b,t,0.1)){
  p <- NULL
  for(j in seq(l,r,0.1)){
    p <- rbind(p,c(i,j))
  }
  pts <- rbind(pts,p)
}

zs_co <- dmvnorm(pts,m,S)
cont_dat_co <- data.frame(X1=pts[,1],X2=pts[,2],X3=zs_co)
zs_noco <- dmvnorm(pts,mu_noco_toju,S_noco_toju)
cont_dat_noco <- data.frame(X1=pts[,1],X2=pts[,2],X3=zs_noco)

cont_dat <- data.frame(cbind(rbind(cont_dat_co,cont_dat_noco),c(rep("Coevolution",length(zs_co)),rep("No coevolution",length(zs_noco)))))
colnames(cont_dat) <- c("X1","X2","X3","type")
coev_effect_toju <- ggplot(cont_dat,aes(x=X1,y=X2,z=X3,color=type))+geom_contour()+theme_bw()+theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values=c("#f8766d","#1ec56f"))+
  ylab(expression(paste("Pericarp thickness ", (mm))))+xlab(expression(paste("Rostrum length ", (mm))))+ggtitle("Camellia-Weevil (CW) System")

ggsave(coev_effect_toju,filename = "~/Research/ecology letters/effect_toju.eps",width=8,height=5,device="eps")


#      #############      ########################################################
#      #           #
#      # Pauw data #      ########################################################
#      #           #
#      #############      ########################################################

snout <- c(79.3,54.8,57.0,72.8,42.9,58.4,85.8,75.6)
tubes <- c(77.0,mean(c(50.0,27.5)),49.1,52.6,41.1,58.6,60.8,68.0)
data <- cbind(snout,tubes)

G1 <- mean(c(4.6,4.4,5.9,5.8,1.6,6.8,6.5,7.4)^2)
G2 <- mean(c(5.9,4.9,4.3,4.7,2.5,7.5,4.6,5.5)^2)

siblings <- c(11.5,32,(38+44)/2)

ancepsAsolo <- c(31.2,51.7)
Moegestorynchus <- c(siblings,mean(siblings)) # brevirostris, perplexus, braunsi

th1s <- Moegestorynchus
th2 <- mean(ancepsAsolo)
n1 <- 100
n2 <- 1000

N <- dim(data)[1]
m <- colMeans(data)
S <- (N-1)*var(data)/N

######
# ML #
######

B1s <- NULL
B2s <- NULL
A1s <- NULL
A2s <- NULL
ks <- NULL
p1s <- NULL
p2s <- NULL
for(th1 in th1s){
  
  ML  <- ml_sol(m[1],m[2],S[1,1],S[2,2],S[1,2],n1,n2,th1,th2,G1,G2)
  lnL <- likelihood(data,m,S)
  
  A1s <- c(A1s,ML$A1)
  A2s <- c(A2s,ML$A2)
  B1s <- c(B1s,ML$B1)
  B2s <- c(B2s,ML$B2)
  ks <- c(ks,ML$k)
  
  ##########
  # null 1 #
  ##########
  
  mu_n1 <- c(th1,m[2])
  S_n1 <- S
  lnL1  <- likelihood(data,mu_n1,S_n1)
  
  ##########
  # null 2 #
  ##########
  
  mu_n2 <- c(m[1],th2)
  S_n2 <- S
  lnL2  <- likelihood(data,mu_n2,S_n2)
  
  # chisqrs
  p1s <- c(p1s,1-pchisq(2*(lnL-lnL1),1))
  p2s <- c(p2s,1-pchisq(2*(lnL-lnL2),1))
  
}

# reporting results for abiotic optimum averaged over sister species
th1 <- th1s[4]
p1_pauw <- p1s[4]
p2_pauw <- p2s[4]
B1_pauw <- B1s[4]
B2_pauw <- B2s[4]
A1_pauw <- A1s[4]
A2_pauw <- A2s[4]
k_pauw <- ks[4]

mu(ML$A1,ML$A2,ML$B1,ML$B2,ML$k,n1,n2,th1,th2)
SIGMA(ML$A1,ML$A2,ML$B1,ML$B2,ML$k,n1,n2,th1,th2,G1,G2)

# conf ints
PTS <- 1000
cnt <- 1
solutions <- matrix(numeric(PTS*4),ncol=4,nrow=PTS)
for(cnt in 1:PTS){
  fdata <- rmvnorm(N,m,S)
  fm <- colMeans(fdata)
  fS <- (N-1)*var(fdata)/N
  ML <- ml_sol(fm[1],fm[2],fS[1,1],fS[2,2],fS[1,2],n1,n2,th1,th2,G1,G2)
  
  solutions[cnt,1] <- ML$A1
  solutions[cnt,2] <- ML$A2
  solutions[cnt,3] <- ML$B1
  solutions[cnt,4] <- ML$B2
  
}

A1_95upper_pauw <- q_fct(0.975,solutions[,1])
A1_95lower_pauw <- q_fct(0.025,solutions[,1])

A2_95upper_pauw <- q_fct(0.975,solutions[,2])
A2_95lower_pauw <- q_fct(0.025,solutions[,2])

B1_95upper_pauw <- q_fct(0.975,solutions[,3])
B1_95lower_pauw <- q_fct(0.025,solutions[,3])

B2_95upper_pauw <- q_fct(0.975,solutions[,4])
B2_95lower_pauw <- q_fct(0.025,solutions[,4])

###############
# effect size #
###############

A1_pauw_noco <- 1/(2*S[1,1]*n1)
A2_pauw_noco <- 1/(2*S[2,2]*n2)
th1 <- mean(th1s)
mu_noco_pauw <- mu(A1_pauw_noco,A2_pauw_noco,0,0,k_pauw,n1,n2,th1,th2)
S_noco_pauw  <- SIGMA(A1_pauw_noco,A2_pauw_noco,0,0,k_pauw,n1,n2,th1,th2,G1,G2)

print("effect size for pauw data")
print(mahalanobis(m,mu_noco_pauw,(S+S_noco_pauw)/2))
print("percent change for pauw data")
print(abs(m-mu_noco_pauw)/mu_noco_pauw)

l <- mu_noco_pauw[1]-sqrt(S_noco_pauw[1,1])
b <- mu_noco_pauw[2]-sqrt(S_noco_pauw[2,2])-30
r <- colMeans(data)[1]+sqrt(var(data)[1,1])
t <- colMeans(data)[2]+sqrt(var(data)[2,2])+10
pts <- NULL
for(i in seq(b,t,1)){
  p <- NULL
  for(j in seq(l,r,1)){
    p <- rbind(p,c(i,j))
  }
  pts <- rbind(pts,p)
}

zs_co <- dmvnorm(pts,m,S)
cont_dat_co <- data.frame(X1=pts[,1],X2=pts[,2],X3=zs_co)
zs_noco <- dmvnorm(pts,mu_noco_pauw,S_noco_pauw)
cont_dat_noco <- data.frame(X1=pts[,1],X2=pts[,2],X3=zs_noco)

cont_dat <- data.frame(cbind(rbind(cont_dat_co,cont_dat_noco),c(rep("Coevolution",length(zs_co)),rep("No coevolution",length(zs_noco)))))
colnames(cont_dat) <- c("X1","X2","X3","type")
coev_effect_pauw <- ggplot(cont_dat,aes(x=X1,y=X2,z=X3,color=type))+geom_contour()+theme_bw()+theme(legend.position = "none",plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values=c("#f8766d","#1ec56f"))+
  ylab(expression(paste("Nectar tube depth ", (mm))))+xlab(expression(paste("Proboscis length ", (mm))))+ggtitle("Fly-Flower (FF) System")

ggsave(coev_effect_pauw,filename = "~/Research/ecology letters/effect_pauw.eps",width=6,height=5,device="eps")

upper <- c(B1_95upper_toju,B2_95upper_toju,B1_95upper_pauw,B2_95upper_pauw)
lower <- c(B1_95lower_toju,B2_95lower_toju,B1_95lower_pauw,B2_95lower_pauw)
yupp <- max(upper)
ylow <- min(lower)

strengths <- data.frame(strength=c(B1_toju,B2_toju,B1_pauw,B2_pauw),
                        species=c("C. camellia","C. japonica","M. longirostris","L. anceps"),
                        system=c("Toju & Sota 2006","Toju & Sota 2006","Pauw et al 2008","Pauw et al 2008"),
                        up=upper,low=lower)

bigPic <- ggplot(strengths,aes(x=reorder(species,as.numeric(system)),y=strength,colour=system))+
  geom_point(size=5)+
  geom_errorbar(data=strengths,mapping=aes(x=species,ymin=low,ymax=up))+
  scale_color_manual(values=c("#f8766d","#1ec56f"))+
  scale_y_continuous(trans='log10',limits=c(ylow,1e-2))+
  xlab("Species")+
  ylab(expression(paste("Strength of biotic selection ",(mm^-2))))+
  labs(colour="System")+
  theme_bw()

bigPic

ggsave(filename = "~/Research/ecology letters/strengths.png",bigPic,width=10,height=5,device="png")

