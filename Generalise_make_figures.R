#genotype = "WT"
#age = "7_dpf"
#Rearing_conditions = "NORM"



atab<-vector("list",length(prefix_list))
seltab<-vector("list",length(prefix_list))
ttab<-vector("list",length(prefix_list))
enseltab<-vector("list",length(prefix_list))
prefix_list=list.files(pattern = pattern)
dflist<-vector("list",length(prefix_list))
cellXYZ_list<-vector("list",length(prefix_list))

cpalette=colorRampPalette(colors=c("blue","red"))(5)
#pchseq=c(0,1,2,15,16,17,18)
pchseq=rep(16,length(prefix_list))

exclusions<-vector("list",length(prefix_list))
side<-vector("list",length(prefix_list))
sample_size<-matrix(NA,length(prefix_list),2)
orig_labels<-vector("list",length(prefix_list))
explained_cells<-matrix(0,length(prefix_list),2)

exclusions<-vector("list",length(prefix_list))

for (i in 1:length(exclusions)) {
  exclusions [[i]] = c(-1)
}



pmu_thresh=0.005
lambda1_thresh=0.05
lambda0_thresh=0.05 
gs_thresh=5

counter=0; 
for(pref in prefix_list){
  folder=paste(pref,sep="") 
  cat(folder,"\n")
  counter=counter+1; 
  centers<-read.table(paste(data_location,pref,centers_suffix,sep=""));
  #centers[,1:2]<-read.table(paste("~/SeG/results/Baysian_network_inference/data_for_Gibbs/7_dpf/",pref,"_registered_cell_centers.dat",sep=""));
  selection<-read.table(paste(folder,"/selection.dat",sep=""))$V1
  centers<-centers[selection+1,]
  cellXYZ_list[[counter]]=centers;
  
  load(paste(folder,"/pars.RData",sep="")); 
  atab[[counter]]=A;
  
  seltab[[counter]] = as.vector(A$gs.mean>gs_thresh & 
                                  A$lambda1.mean>lambda1_thresh & 
                                  A$pmu.mean>pmu_thresh &
                                  A$lambda0.mean<lambda0_thresh &
                                  !as.numeric(gsub("V","",names(A$gs.mean)))%in% exclusions[[counter]])
  
  enseltab[[counter]] = as.numeric(gsub("V","",names(A$gs.mean)))[seltab[[counter]]]
  counter2=0
  side[[counter]]<-rep(NA,sum(seltab[[counter]]))
  sidefile<-data.frame((1:length(seltab[[counter]]))[seltab[[counter]]],1,1)
  names(sidefile)<-c("enID","RL","AP")
  for(en in enseltab[[counter]]){
    counter2=counter2+1
    side[[counter]][counter2]=1
  }
  
  dflist[[counter]] = data.frame(ID=counter,lambda1=A$lambda1.mean,lambda0=A$lambda0.mean,gs=A$gs.mean,pmu=A$pmu.mean,centers=A$centers)[seltab[[counter]],]
  dflist[[counter]]$side=side[[counter]]
  dflist[[counter]]$LI <- unlist(read.table(file = paste(folder,"/LI.dat",sep="")))
  dflist[[counter]]$wcc <- unlist(read.table(file = paste(folder,"/wcc.dat",sep="")))
  dflist[[counter]]$Assembly_freq <- unlist(read.table(file = paste(folder,"/Assembly_Freq.dat",sep="")))
  dflist[[counter]]$abs_LI <- abs(unlist(read.table(file = paste(folder,"/LI.dat",sep=""))))
  dflist[[counter]]$areas <- unlist(read.table(file = paste(folder,"/ensem_area.dat",sep=""))[1])
  dflist[[counter]]$density <- unlist(read.table(file = paste(folder,"/ensem_area.dat",sep=""))[2])
  dflist[[counter]]$tectum_total_area <- unlist(read.table(file = paste(folder, "/eigen_decomp_of_cov.dat", sep =""))$X)
  dflist[[counter]]$norm_eigen_decomp <- unlist(read.table(file = paste(folder, "/eigen_decomp_of_cov.dat", sep =""))$norm_eigen_decomp)
  
  sample_size[counter,] = as.numeric(read.table(paste(folder,"/dims.dat",sep="")))[1:2]
  orig_labels[[counter]]=as.numeric(gsub("V","",names(A$gs.mean)[seltab[[counter]]]))
  cat(orig_labels[[counter]],"\n")
  explained_cells[counter,1] = sum(table(A$membership[A$membership %in% orig_labels[[counter]]]));
  explained_cells[counter,2] = sample_size[counter,1] - explained_cells[counter,1]
}


df<-do.call(rbind,dflist)

plot_cen1<- function(lis=1:7){
  
  i=lis[1];
  
  plot(atab[[i]]$centers[seltab[[i]],1:2],
       col=cpalette[as.numeric(cut(atab[[i]]$centers[seltab[[i]],3],breaks=seq(1,5,l=6)))],
       pch=pchseq[i],
       cex=atab[[i]]$gs.mean[seltab[[i]]]/20,
       xlim=c(-140,100),ylim=c(-120,150));
  points(cellXYZ_list[[i]][,1:2],cex=0.1)
  if(length(lis)>1){			
    for(i in lis[-1]){
      points(atab[[i]]$centers[seltab[[i]],1:2],col=cpalette[as.numeric(cut(atab[[i]]$centers[seltab[[i]],3],breaks=seq(1,5,l=6)))],pch=pchseq[i],cex=atab[[i]]$gs.mean[seltab[[i]]]/20);
    }
  }
}

plot_cen2<- function(lis){
  
  i=lis;
  
  plot(cellXYZ_list[[i]][,1:2],cex=0.1)
  counter=0;
  for(id in as.numeric(gsub("V","",names(atab[[i]]$gs.mean)))[seltab[[i]]]){
    counter=counter+1;
    points(cellXYZ_list[[i]][atab[[i]]$membership==id,1:2],col=rainbow(20)[counter],pch=19)
  }
}

plot_explained<-function(){
  p=barplot(t(explained_cells),ylim=c(0,2500))
  text(p,explained_cells[,1],label=paste(round(explained_cells[,1]/sample_size[,1]*100,digits=1),"%"),pos=3)
}

plot_params_across <- function(){
  layout(matrix(c(1,2),1,2)); 
  par(cex=1.2,cex.lab=1.5)
  plot(df$lambda1~df$pmu,col=c("#ff9500","#0059ff")[df$side+1],pch=19,xlab="activity",ylab="coherence"); 
  plot(df$lambda1~df$lambda0,col=c("#ff9500","#0059ff")[df$side+1],pch=19,xlab="noise",ylab="coherence")
  points(df2$lambda1~df2$lambda0,col=c("#0059ff", "#ff9500")[df2$side+1],pch=19)
}


plot_params_across <- function(){
  layout(matrix(c(1,2,3),1,3)); 
  par(cex=1.2,cex.lab=1.5)
  plot(df$gs~df$lambda0,col=c("#ff9500","#0059ff")[df$side+1],pch=19,xlab="activity",ylab="coherence"); 
  points(df2$gs~df2$lambda0,col=c("#0059ff", "#ff9500")[df2$side+1],pch=19)
  plot(df$lambda1~df$pmu,col=c("#ff9500","#0059ff")[df$side+1],pch=19,xlab="activity",ylab="coherence", log = "xy"); 
  points(df2$lambda1~df2$pmu,col=c("#0059ff", "#ff9500")[df2$side+1],pch=19)
  plot(df$lambda1~df$lambda0,col=c("#ff9500","#0059ff")[df$side+1],pch=19,xlab="noise",ylab="coherence")
  points(df2$lambda1~df2$lambda0,col=c("#0059ff", "#ff9500")[df2$side+1],pch=19)
}

make.boxplots <- function(){
  layout(matrix(1:4,2,2))
  par(cex=1.2,cex.lab=1.5)
  boxplot(df$lambda1 ~ df$ID,ylim=c(0,.7),ylab="coherence",outline=F)
  boxplot(df$pmu ~ df$ID,ylim=c(0,.3),ylab="activity",outline=F)
  boxplot(df$lambda0 ~ df$ID,ylim=c(0,.05),ylab="noise",outline=F)
  boxplot(df$gs~df$ID,outline=F,ylim=c(0,100),ylab="ensemble size")
}

plot_n_gs<- function(){
  par(cex=1.2,cex.lab=1.5)
  barplot(table(df$ID), ylab="number of ensembles",ylim=c(0,55))
}

plot_AP_corr_all<- function(print=FALSE){
  AP=c();
  AA=c();
  PP=c();
  for(i in 1:7){
    pref=prefix_list[i]
    folder=paste(pref,"_S1",sep="")
    AP=c(AP,read.table(paste(folder,"/cormatAP.dat",sep=""))$V1)
    AA=c(AA,read.table(paste(folder,"/cormatAA.dat",sep=""))$V1)
    PP=c(PP,read.table(paste(folder,"/cormatPP.dat",sep=""))$V1)
  }
  
  ind=seq(0,max(c(AA,PP)),0.01)
  if(print) pdf(paste("figures2/AP_correlations_combined.pdf",sep=""),5,5)
  plot(ind,1-ecdf(AA)(ind),type='l',col="red",xlab="time correlation",ylab="frequency",cex.lab=1.5,cex.axis=1.5)
  lines(ind,1-ecdf(PP)(ind),col="blue")
  lines(ind,1-ecdf(AP)(ind),col="magenta")
  legend("topright",legend=c("anterior-anterior","posterior-posterior","anterior-posterior"),lty=1,col=c("red","blue","magenta"))
  if(print) dev.off();
}

plot_AP_corr_across_all<- function(print=FALSE){
  AP=c();
  AA=c();
  PP=c();
  for(i in 1:7){
    pref=prefix_list[i]
    folder=paste(pref,"_S1",sep="")
    AP=c(AP,read.table(paste(folder,"/cormatAP_across.dat",sep=""))$V1)
    AA=c(AA,read.table(paste(folder,"/cormatAA_across.dat",sep=""))$V1)
    PP=c(PP,read.table(paste(folder,"/cormatPP_across.dat",sep=""))$V1)
  }
  
  ind=seq(0,max(c(AA,PP)),0.01)
  if(print) pdf(paste("figures2/AP_correlations_across_combined.pdf",sep=""),5,5)
  plot(ind,1-ecdf(AA)(ind),type='l',col="red",xlab="time correlation",ylab="frequency",cex.lab=1.5,cex.axis=1.5)
  lines(ind,1-ecdf(PP)(ind),col="blue")
  lines(ind,1-ecdf(AP)(ind),col="magenta")
  legend("topright",legend=c("AA","PP","AP"),lty=1,col=c("red","blue","magenta"))
  if(print) dev.off();
}

plot_corrdist <- function(print=FALSE){
  df<-vector("list",7)
  for(i in 1:7){
    pref=prefix_list[i]
    df[[i]]<-read.table(paste("figures2/corr_vs_dist_",pref,".dat",sep=""))
  }
  df_glob<- do.call(rbind,df)
  colnames(df_glob)=c("dist","correlation")
  if(print) pdf("figures2/corr_vs_dist.pdf",5,5)
  plot(df_glob$dist*0.8,df_glob$correlation,pch=19,cex=0.5)
  if(print) dev.off()
}

df$norm_eigen <- df$norm_eigen_decomp/df$tectum_total_area
dat <- apply(df,2, function(x) by(x,df$ID,mean))



colnames(explained_cells) <- c("ensembled_neurons", "free_neurons")
numb_ensembles <- data.frame(table(df$ID)) [,2]
dat <- as.data.frame(cbind(dat,explained_cells, numb_ensembles)) 

dat$Genotype <-  rep(paste(genotype),1, length(prefix_list))
dat$age <- rep(paste(age), 1, length(prefix_list))
dat$Rearing_conditions <- rep(paste(Rearing_conditions), 1, length(prefix_list))
dat$mean_corr_total <-  unlist(as.vector(read.table(paste("mean_across_fish/mean_corr_", genotype,"_", age,"_", Rearing_conditions, ".dat", sep = ""))))
dat$sn_freq <- unlist(as.vector(read.table(paste("mean_across_fish/sn_freq_", genotype,"_", age,"_", Rearing_conditions, ".dat", sep = ""))))
dat$cell_num <- unlist(as.vector(read.table(paste("mean_across_fish/cell_num_", genotype,"_", age,"_", Rearing_conditions, ".dat", sep = ""))))


write.table(dat, paste("mean_across_fish/", genotype,"_", age,"_", Rearing_conditions, "_mean_datatable.dat", sep =""))

df$Genotype <-  rep(paste(genotype),1, nrow(df))
df$age <- rep(paste(age), 1, nrow(df))
df$Rearing_conditions <- rep(paste(Rearing_conditions), 1, nrow(df))
write.table(df, paste("mean_across_fish/", genotype,"_", age,"_", Rearing_conditions, "_all_assemblies_datatable.dat", sep =""))





calculate_difference_mat <- function(x = "gs",  y = "lambda1", cond_a, cond_b) {
  par(mfrow= c(2,2))
  max_a <- max(max(cond_a [x]), max(cond_b [x]))
  max_b <- max(max(cond_a [y]), max(cond_b[y]))
  grav <- kde2d(unlist(cond_a[x]), unlist(cond_a [y]), n= 100, lims = c(0,max_a, 0, max_b))
  norm <- kde2d(unlist(cond_b[x]), unlist(cond_b [y]), n= 100, lims = c(0, max_a, 0, max_b))
  grav$z <- grav$z/sum(grav$z)
  norm$z <- norm$z/sum(norm$z)
  diff <- norm
  diff$z <- norm$z - grav$z
  library(RColorBrewer)
  rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
  r <- rf(32)
  absval=max(abs(range(diff$z)))
  image(diff,col = r, xlab = x, ylab = y,breaks=seq(-absval,absval,l=33))
  image(norm, col = r, xlab = x, ylab = y,breaks=seq(0,max(norm$z),l=33))
  title("norm")
  image(grav, col = r, xlab = x, ylab = y,breaks=seq(0,max(grav$z),l=33))
  title("grav")
  return(diff)
}