atab<-vector("list",7)
seltab<-vector("list",7)
ttab<-vector("list",7)
enseltab<-vector("list",7)
prefix_list=list.files(pattern = "WT_grav_h2b_gc6s")
dflist<-vector("list",7)
cellXYZ_list<-vector("list",7)

cpalette=colorRampPalette(colors=c("blue","red"))(5)
#pchseq=c(0,1,2,15,16,17,18)
pchseq=rep(16,7)

exclusions<-vector("list",7)
side<-vector("list",7)
sample_size<-matrix(NA,7,2)
orig_labels<-vector("list",7)
explained_cells<-matrix(0,7,2)

exclusions[[1]]=c(-1);
exclusions[[2]]=c(-1);
exclusions[[3]]=c(-1)
exclusions[[4]]=c(-1)
exclusions[[5]]=c(-1)
exclusions[[6]]=c(-1)
exclusions[[7]]=c(-1)


pmu_thresh=0.005
lambda1_thresh=0.05
lambda0_thresh=0.05 
gs_thresh=5

counter=0; 
for(pref in prefix_list){
  folder=paste(pref,sep="") 
  cat(folder,"\n")
  counter=counter+1; 
  centers<-read.table(paste("~/SeG/results/Baysian_network_inference/data_for_Gibbs/Grav/",pref,"_scaled_aligned_all_cells_centers_cut_ordered.dat",sep=""));
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
  plot(df$lambda1~df$lambda0,col=c("#ff9500","#0059ff")[df$side+1],pch=19,xlab="noise",ylab="coherence", log ="xy")
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


grav_dpf7_df <- apply(df,2, function(x) by(x,df$ID,mean))

colnames(explained_cells) <- c("ensembled_neurons", "free_neurons")
numb_ensembles <- data.frame(table(df$ID)) [,2]
grav_dpf7_df <- as.data.frame(cbind(grav_dpf7_df,explained_cells, numb_ensembles)) 

grav_dpf7_df$group <-  rep("GRAV",1, 7)
grav_dpf7_df$age <- rep("7_dpf", 1, 7)
write.table(grav_dpf7_df, "/Users/MeyerLab/SeG/results/Baysian_network_inference/BPTS_Tom/grav_dpf7_df_mean_datatable.dat")
