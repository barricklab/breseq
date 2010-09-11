X<-read.table("cnv_cov.tab", header=T)
X<-subset(X, cov!=0)
X<-subset(SX, red==0)
Y<-SSX$cov;
GC<-SSX$gc;

m<-mean(Y)
v<-var(Y)
D<-v/m

library(vcd)
f_nb<-goodfit(Y, type="nbinomial", method="ML")
print(f_nb$par);

counts<-c();
for(i in 1:max(Y)) {
  counts[i] = 0;
}

for(i in 1:length(Y)) {
  counts[Y[i]] = counts[Y[i]]+1;
}

pdf("t.pdf")
plot(1:length(counts), counts, bg="white", cex=1.2)
lines(f_nb$count, f_nb$fitted, lwd=3, col="black");
dev.off()

pdf("gc.pdf")
plot(GC, Y, bg="white", cex=1.2)
dev.off()

fit<-lm(cov ~ gc, X)

cor<-c();
fitted<-fitted.values(fit)

cor<- X$cov/fitted*m;
cor<-trunc(cor);

cor_counts<-c();
for(i in 1:max(cor)) {
  cor_counts[i] = 0;
}

for(i in 1:length(cor)) {
  cor_counts[cor[i]] = cor_counts[cor[i]]+1;
}

f_nb_cor<-goodfit(cor, type="nbinomial", method="ML")
print(f_nb_cor$par);

pdf("cor.pdf")
plot(1:length(cor_counts), cor_counts, bg="white", cex=1.2)
lines(f_nb_cor$count, f_nb_cor$fitted, lwd=3, col="black");
dev.off()

#barely perceptible improvement by doing correction

X<-read.table("ZDB172.2/cnv_cov.tab", header=T)
Y<-read.table("JEB559_cnv_2/cnv_cov.tab", header=T)
Z<-c();
for(i in 1:max(X$pos,Y$pos)) {
  Z<-rbind(Z, c(i,0,0));
}

X<-subset(X, cov!=0)
X<-subset(X, red==0)
for(i in 1:length(X$pos)) {
  Z[X$pos[i],2] = X$cov[i]; 
}

Y<-subset(Y, cov!=0)
Y<-subset(Y, red==0)
for(i in 1:length(Y$pos)) {
  Z[Y$pos[i],3] = Y$cov[i]; 
}

ZZ<-c();
for(i in 1:length(Z[,1])) {
  if ( (Z[i,2] != 0) && (Z[i,3] != 0) ) {
    ZZ<-rbind(ZZ, Z[i,]); 
  }
}

pdf("cor.pdf")
plot(ZZ[,2], ZZ[,3]);
dev.off()
