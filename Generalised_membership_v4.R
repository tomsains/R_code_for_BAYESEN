# attempt to make the member ship script general


#data_location = "../data_for_Gibbs/7_dpf/"
#centers_suffix = "_sa_aligned_all_cells_centers_cut_ordered.dat"
#spikes_suffix = "_sa_aligned_all_cells_spikes.dat"

#midline <- data.frame(cbind(intercept = c(10, 10, 10, 10, 10, 10, 10, 10), slope = c(-0.97, -0.97, -0.97, -0.97, -0.97, -0.97, -0.97, -0.97)))
#pattern <- "*WT_h2b_gc6s_7dpf_*"


library(splancs)

prefix_list = list.files(pattern = pattern)

for (i in 1:length(prefix_list)){
  # folder containing the results of BPTS
  prefix = prefix_list [i]
  folder=paste(prefix,sep="")
  
  # graph plots
  library(igraph)
  # color-blind color scheme
  library(viridis)
  
  # load the prefix list and determine the index of the prefix under consideration in the list
  list.files(pattern = pattern)
  prefix_ID=which(prefix==prefix_list)
  
  # excluded ensembles because only active at the beginning of the recording
  exclusions<-vector("list",length(prefix_list))
  
 for (i in 1:length(exclusions)) {
   exclusions [[i]] = c(-1)
 }
  
  exc=exclusions[[which(prefix==prefix_list)]]
  
  # file containing infos about side and anterior-posterior location of the ensemble
  #side_AP<-read.table(paste(folder,"/side_AP.dat",sep=""))
  #names(side_AP)<-c("enID","RL","AP")
  
  # Thresholds used for the analysis
  pmu_thresh=0.005
  lambda1_thresh=0.05
  lambda0_thresh=0.05 
  gs_thresh=5
  
  # number of samples extracted from the MCMC
  nsamples=20
  
  # load input binary matrix S
  loadS=T
  if(loadS) s<-read.table(paste(data_location, prefix, spikes_suffix ,sep=""));
  
  # load results
  mem_traj<-read.table(paste(folder,"/membership_traj.dat",sep=""))
  centers<-read.table(paste(data_location ,prefix,centers_suffix,sep=""));
  
  if (reg == TRUE) centers[,1:2]<-read.table(paste(data_location, prefix,"_registered_cell_centers.dat",sep=""));
  
  
  selection<-read.table(paste(folder,"/selection.dat",sep=""))$V1
  centers<-centers[selection+1,]
  FreeEnergy<-read.table(paste(folder,"/F.dat",sep=""))
  tot_samples=nrow(mem_traj)
  mem<-apply(mem_traj[(tot_samples-nsamples):tot_samples,],2,function(x) which.max(tabulate(x+1)))
  tab=sort(table(mem),decreasing=T)
  enIDs=as.numeric(names(tab))
  #mem=as.numeric(mem_traj[3000,]+1)
  probs<-apply(mem_traj[(tot_samples-nsamples):tot_samples,],2,function(x) max(tabulate(x+1)/length(x)))
  
  # write centers and membership
  write.table(centers,file=paste(folder,"/centers.dat",sep=""),row.names=FALSE,col.names=FALSE)
  write.table(mem,file=paste(folder,"/mem.dat",sep=""),row.names=FALSE,col.names=FALSE)
  maxcl=max(mem);
  
  plot_par <- function(P,plot=FALSE){
    
    assembly_centers<-matrix(NA,length(enIDs),3)
    pmu=as.matrix(read.table(paste(folder,"/pmu.dat",sep="")))[(tot_samples-nsamples+1):tot_samples,enIDs]
    lambda0=as.matrix(read.table(paste(folder,"/lambda0.dat",sep="")))[(tot_samples-nsamples+1):tot_samples,enIDs]
    lambda1=as.matrix(read.table(paste(folder,"/lambda1.dat",sep="")))[(tot_samples-nsamples+1):tot_samples,enIDs]
    gs=as.matrix(read.table(paste(folder,"/n.dat",sep="")))[(tot_samples-nsamples+1):tot_samples,enIDs]
    
    counter=0;
    for(i in enIDs){
      counter=counter+1;
      assembly_centers[counter,] <- apply(centers[mem==i,],2,mean);
    }
    
    if(plot){
      X11();layout(matrix(1:4,1,4))
      boxplot(pmu[,sel],main="pmu")
      boxplot(lambda0[,sel],main="lambda0")
      boxplot(lambda1[,sel],main="lambda1")
      boxplot(gs[,sel],main="group sizes")
    }
    
    return(list(pmu.mean=apply(pmu,2,mean),
                lambda0.mean=apply(lambda0,2,mean),
                lambda1.mean=apply(lambda1,2,mean),
                centers=assembly_centers,
                gs.mean=apply(gs,2,mean),
                membership=mem))
    
  }
  
  # define the list of parameters nad save it
  A=plot_par();
  save(A, file=paste(folder,"/pars.RData",sep=""))
  # define the selection to be applied to tab
  sel=as.vector(A$gs.mean>gs_thresh & 
                  A$lambda1.mean>lambda1_thresh & 
                  A$pmu.mean>pmu_thresh &
                  A$lambda0.mean<lambda0_thresh &
                  !as.numeric(gsub("V","",names(A$gs.mean)))%in% exc)
  # ensemble selection
  ensel=enIDs[sel]
  
  # calculate side_AP_ordered
  #side_AP_ord=data.frame(enID=ensel,RL=NA,AP=NA)
  #counter=0;
  #for(i in ensel){
  #	counter=counter+1
  #	j=which(side_AP$enID==i)
  #	if(length(j)==0) cat("problems in finding ensemble ID!!!!\n")
  #	side_AP_ord[counter,]=side_AP[j,]
  #}
  
  plot_F <- function(){
    plot(apply(FreeEnergy,1,mean),type='l')
  }
  
  plot_P <- function(){
    ptrace<-read.table(paste(folder,"/P.dat",sep=""))$V1
    ptrace<-ptrace[(length(ptrace)-nsamples):length(ptrace)]
    plot(ptrace,type='l')
  }
  
  plot_assemblies_simple <- function(minprob=0){
    layout(matrix(1:35,nrow=5))
    par(mar=c(0,0,0,0));
    for(i in 1:35){
      plot(centers[,1:2],pch=19,cex=1,xaxt='n',yaxt='n')
      title(paste(i),line=-1)
      points(centers[mem==i & probs>minprob,1:2],col="red",pch=19)
    }
    plot(apply(FreeEnergy,1,mean),type='l')
  }
  
  
  plot_assemblies_custom <- function(print=F,snc=FALSE){
    
    subnet_mem<-wc$membership
    minprob=0.99
    
    if(print) png(paste("figures2/assemblies_",prefix,".png",sep=""),width=6*2,height=7*2,units='in',res=300)
    layout(matrix(1:42,ncol=6,byrow=T))
    par(mar=rep(0.3,4));
    counter=0;
    
    for(i in ensel[order(subnet_mem)]){
      counter=counter+1;			
      plot(centers[,1:2],pch=18,cex=.3,xaxt='n',yaxt='n',col=rgb(0,0,0,1))
      text(x=-100,y=100,paste(i),cex=2.5)
      #text(x=-100,y=100,i,cex=1.5)
      points(centers[mem==i & probs>minprob,1:2],col=viridis(max(subnet_mem))[sort(subnet_mem)[counter]],pch=19)
    }
    if(print) dev.off();
  }
  
  plot_cormat<- function(print=F){
    if(print) png(paste("figures2/cormat_",prefix,".png",sep=""),width=5*2,height=5*2,units='in',res=300)
    
    subnet_mem=wc$membership
    layout(matrix(1:4,2,2),widths=c(.5,5),heights=c(5,.5))
    par(mar=c(0,1,1,0))
    image(matrix(1:sum(sel),1,sum(sel)),col=viridis(max(subnet_mem))[subnet_mem[order(subnet_mem)]],xlab='',ylab='',xaxt='n',yaxt='n')
    par(mar=c(1,1,0,0))
    plot.new()
    par(mar=c(0,0,1,1))
    cmat=cor(t(omega[order(wc$membership),]))
    image(cmat,col=grey.colors(100),frame=F,axes=F)
    par(mar=c(1,0,0,1))
    image(matrix(1:sum(sel),sum(sel),1),col=viridis(max(subnet_mem))[subnet_mem[order(subnet_mem)]],xlab='',ylab='',xaxt='n',yaxt='n')
  }
  
  
  
  plot_assemblies <- function(x,minprob=0){
    if(x<0){
      layout(matrix(1:55,nrow=5))
      par(mar=c(0,0,0,0));
      for(i in ensel){
        
        plot(centers[,1:2],pch=18,cex=.3,xaxt='n',yaxt='n',col=rgb(0,0,0,1))
        title(paste(i),line=-1)
        points(centers[mem==i & probs>minprob,1:2],col="red",pch=19)
      }
      plot(apply(FreeEnergy,1,mean),type='l')
    }
    
    if(x>0) {
      plot(centers[,1:2],col="black",pch=19)
      points(centers[mem==x & probs>minprob,1:2],col="red",pch=19)
    }
    
    
  }
  
  get_omega <- function(print=FALSE,P,PMAX=200,nav=nsamples){
    system(paste("tail -n ",PMAX*nav," ",folder,"/omega_traj.dat > ",folder,"/omega_traj_tmp.dat",sep=""))
    omegatraj<-as.matrix(read.table(paste(folder,"/omega_traj_tmp.dat",sep="")))
    
    times=ncol(omegatraj)
    samples=nrow(omegatraj)/PMAX
    
    omegaprob=matrix(0,ncol=times,nrow=PMAX)
    for(i in (samples-nav):(samples-1)){
      omegaprob=omegaprob+omegatraj[ i*PMAX + 1:PMAX,]
    }
    
    omegaprob=omegaprob/nav
    omegaprob=omegaprob[ensel,]
    return(omegaprob)
  }
  
  get_omega0 <- function(print=FALSE,P,PMAX=200,nav=nsamples){
    system(paste("tail -n ",PMAX*nav," ",folder,"/omega_traj.dat > ",folder,"/omega_traj_tmp.dat",sep=""))
    omegatraj<-as.matrix(read.table(paste(folder,"/omega_traj_tmp.dat",sep="")))
    
    times=ncol(omegatraj)
    samples=nrow(omegatraj)/PMAX
    
    omegaprob=matrix(0,ncol=times,nrow=PMAX)
    for(i in (samples-nav):(samples-1)){
      omegaprob=omegaprob+omegatraj[ i*PMAX + 1:PMAX,]
    }
    
    omegaprob=omegaprob/nav
    return(omegaprob)
  }
  
  
  plot_omega <- function(print=FALSE,PMAX=200,subset=NA,nav=nsamples){
    system(paste("tail -n ",PMAX*nav," ",folder,"/omega_traj.dat > ",folder,"/omega_traj_tmp.dat",sep=""))
    omegatraj<-as.matrix(read.table(paste(folder,"/omega_traj_tmp.dat",sep="")))
    
    times=ncol(omegatraj)
    samples=nrow(omegatraj)/PMAX
    
    omegaprob=matrix(0,ncol=times,nrow=PMAX)
    for(i in (samples-nav):(samples-1)){
      omegaprob=omegaprob+omegatraj[ i*PMAX + 1:PMAX,]
    }
    
    par(mar=c(1,4,1,1))
    omegaprob=omegaprob/nav
    if(is.na(subset)){
      omegaprob=omegaprob[ensel,]
      image(x=1:times,y=1:sum(sel),z=t(omegaprob),col=grey.colors(100),yaxt='n',xlab="",ylab="",xaxt='n',cex.axis=2,cex.lab=2,breaks=seq(0,1,length=101))
      axis(2,at=1:sum(sel),labels=ensel)
    } else {
      omegaprob=omegaprob[subset,]
      image(x=1:times,y=1:length(subset),z=t(omegaprob),col=grey.colors(100),yaxt='n',xlab="synchronous frames",ylab="assemblies",xaxt='n',cex.axis=2,cex.lab=2,breaks=seq(0,1,length=101))
      axis(2,at=1:length(subset),labels=subset)
    }
    
  }
  
  plot_omega_custom <- function(print=FALSE,P,PMAX=100,subset=NA){
    
    if(print) png(paste("figures2/omega_",prefix,".png",sep=""),width=5*2,height=7*2,units='in',res=300)
    
    subnet_mem=wc$membership
    times=ncol(omega)
    
    par(mar=c(1,4,1,0))
    omegaprob=as.matrix(omega[order(subnet_mem),])
    layout(matrix(1:2,1,2),widths=c(5,.3))
    
    image(x=1:times,y=1:sum(sel),z=1-t(omegaprob),col=grey.colors(100),yaxt='n',xlab="",ylab="",xaxt='n',cex.axis=2,cex.lab=2,breaks=seq(0,1,length=101))
    
    #axis(2,at=1:sum(sel),labels=rep("",sum(sel)))
    axis(2,at=1:sum(sel),labels=ensel[order(subnet_mem)],las=2,cex.axis=2)
    
    par(mar=c(1,0,1,1))
    image(matrix(1:sum(sel),1,sum(sel)),col=viridis(max(subnet_mem))[subnet_mem[order(subnet_mem)]],xlab='',ylab='',xaxt='n',yaxt='n')
    if(print) dev.off()
    
  }
  
  autocor<- function(s){
    ls=length(s)
    tmat<-matrix(0,2,2)
    for(i in 1:(ls-1)){
      tmat[s[i+1]+1,s[i]+1]=tmat[s[i+1]+1,s[i]+1]+1;
    }
    tmat = t(t(tmat)/colSums(tmat))
    return(tmat[2,2])
  }
  
  cordist<- function(use="L"){
    
    if(use=="L"){
      sel_side=(side_AP$enID %in% ensel) & side_AP$RL==1 & !is.na(side_AP$RL);
    } else {
      sel_side=(side_AP$enID %in% ensel) & side_AP$RL==0 & !is.na(side_AP$RL);
    }
    
    ng=sum(sel_side)
    
    assembly_centers<-matrix(NA,ng,2);
    counter=0;
    for(i in side_AP$enID[sel_side]){
      counter=counter+1;
      assembly_centers[counter,]<-apply(centers[mem==i,1:2],2,mean)
    }
    omegamat=get_omega0()[side_AP$enID[sel_side],];
    cormat=cor(t(omegamat))[lower.tri(matrix(NA,ng,ng))];
    dist=as.matrix(dist(assembly_centers,upper=T))[lower.tri(matrix(NA,ng,ng))]
    
    return(list(cormat=cormat,dist=dist))
    
  }
  
  get_corr_AP <- function(print=F){
    
    cormat=cor(t(get_omega()));
    diag(cormat)=0;
    sel_A=(side_AP_ord$AP==0 & !is.na(side_AP_ord$AP))
    sel_P=(side_AP_ord$AP==1 & !is.na(side_AP_ord$AP))
    
    # use Left side 
    sel_side=(side_AP_ord$RL==1 & !is.na(side_AP_ord$RL));
    cormat_APL=as.vector(cormat[sel_side & sel_A, sel_side & sel_P]);
    cormat_AAL=as.vector(cormat[sel_side & sel_A, sel_side & sel_A]);
    cormat_PPL=as.vector(cormat[sel_side & sel_P, sel_side & sel_P]);
    
    # use right side 
    sel_side=(side_AP_ord$RL==0 & !is.na(side_AP_ord$RL));
    cormat_APR=as.vector(cormat[sel_side & sel_A, sel_side & sel_P]);
    cormat_AAR=as.vector(cormat[sel_side & sel_A, sel_side & sel_A]);
    cormat_PPR=as.vector(cormat[sel_side & sel_P, sel_side & sel_P]);
    if(print){
      write.table(c(cormat_APR,cormat_APL),file=paste(folder,"/cormatAP.dat",sep=""),col.names=F,row.names=F) 
      write.table(c(cormat_AAR,cormat_AAL),file=paste(folder,"/cormatAA.dat",sep=""),col.names=F,row.names=F) 
      write.table(c(cormat_PPR,cormat_PPL),file=paste(folder,"/cormatPP.dat",sep=""),col.names=F,row.names=F) 
    }
    
    return(list(AP=c(cormat_APR,cormat_APL),AA=c(cormat_AAR,cormat_AAL), PP=c(cormat_PPR,cormat_PPL)))
    
  }
  
  get_corr_AP_across <- function(print=F){
    
    cormat=cor(t(get_omega()));
    diag(cormat)=0;
    sel_A=(side_AP_ord$AP==0 & !is.na(side_AP_ord$AP))
    sel_P=(side_AP_ord$AP==1 & !is.na(side_AP_ord$AP))
    left_side=(side_AP_ord$RL==1 & !is.na(side_AP_ord$RL));
    right_side=(side_AP_ord$RL==0 & !is.na(side_AP_ord$RL));
    
    # use Left/right side 
    cormat_APL=as.vector(cormat[right_side & sel_A, left_side & sel_P]);
    cormat_AAL=as.vector(cormat[right_side & sel_A, left_side & sel_A]);
    cormat_PPL=as.vector(cormat[right_side & sel_P, left_side & sel_P]);
    
    # use right side 
    cormat_APR=as.vector(cormat[left_side & sel_A, right_side & sel_P]);
    cormat_AAR=as.vector(cormat[left_side & sel_A, right_side & sel_A]);
    cormat_PPR=as.vector(cormat[left_side & sel_P, right_side & sel_P]);
    if(print){
      write.table(c(cormat_APR,cormat_APL),file=paste(folder,"/cormatAP_across.dat",sep=""),col.names=F,row.names=F) 
      write.table(c(cormat_AAR,cormat_AAL),file=paste(folder,"/cormatAA_across.dat",sep=""),col.names=F,row.names=F) 
      write.table(c(cormat_PPR,cormat_PPL),file=paste(folder,"/cormatPP_across.dat",sep=""),col.names=F,row.names=F) 
    }
    
    return(list(AP=c(cormat_APR,cormat_APL),AA=c(cormat_AAR,cormat_AAL), PP=c(cormat_PPR,cormat_PPL)))
    
  }
  
  plot_corr_AP <- function(print=FALSE){
    u<-get_corr_AP(print);
    ind=seq(0,max(c(u$AA,u$PP)),0.01)
    if(print) pdf(paste("figures2/AP_correlations_",prefix,".pdf",sep=""),5,5)
    plot(ind,1-ecdf(u$AA)(ind),type='l',col="red",xlab="time correlation",ylab="frequency",cex.lab=1.5,cex.axis=1.5)
    lines(ind,1-ecdf(u$PP)(ind),col="blue")
    lines(ind,1-ecdf(u$AP)(ind),col="magenta")
    legend("topright",legend=c("AA","PP","AP"),lty=1,col=c("red","blue","magenta"))
    if(print) dev.off();
  }
  
  plot_corr_AP_across <- function(print=FALSE){
    u<-get_corr_AP_across(print);
    ind=seq(0,max(c(u$AA,u$PP)),0.01)
    if(print) pdf(paste("figures2/AP_correlations_across",prefix,".pdf",sep=""),5,5)
    plot(ind,1-ecdf(u$AA)(ind),type='l',col="red",xlab="time correlation",ylab="frequency",cex.lab=1.5,cex.axis=1.5)
    lines(ind,1-ecdf(u$PP)(ind),col="blue")
    lines(ind,1-ecdf(u$AP)(ind),col="magenta")
    legend("topright",legend=c("AA","PP","AP"),lty=1,col=c("red","blue","magenta"))
    if(print) dev.off();
  }
  
  plot_corr <- function(print=FALSE){
    uL=cordist("L")
    uR=cordist("R")
    if(print){
      write.table(cbind(c(uL$dist,uR$dist),c(uL$cormat,uR$cormat)),file=paste("figures2/corr_vs_dist_",prefix,".dat",sep=""),col.names=F,row.names=F)
    }
    if(print) pdf(paste("figures2/corr_vs_dist_",prefix,".pdf",sep=""),height=5,width=10)
    layout(matrix(1:2,1,2))
    plot(uL$dist,uL$cormat,main="OTL",xlab="distance",ylab="correlation",cex.lab=1.5)
    plot(uR$dist,uR$cormat,main="OTR",xlab="distance",ylab="correlation",cex.lab=1.5)
    if(print) dev.off();
  }
  
  plot_par_means<- function(pmu_thresh=0,gs_thresh=1,xlim=c(0,1),ylim=c(0,1)) {
    plot(A$lambda1.mean[A$pmu.mean>pmu_thresh & A$gs.mean>gs_thresh],A$lambda0.mean[A$pmu.mean>pmu_thresh & A$gs.mean>gs_thresh],cex=2,cex.axis=3,cex.lab=2,xlim=xlim,ylim=ylim)
  }
  
  omega=get_omega();
  plot_cells <- function(i,thresh=0){
    layout(matrix(1:2,2,1),heights=c(7,1))
    id=which(ensel==i)
    sum_act=apply(s[selection+1,],2,sum)
    par(mar=c(0,2,2,1))
    image(t(as.matrix(s)[selection+1,sum_act>thresh][mem==i & probs>0.9,]),axes=F,frame=F)
    par(mar=c(2,2,0,1))
    image(as.matrix(omega[id,]),frame=F,axes=F,col=gray.colors(100),breaks=seq(0,1,length=101))
  }
  
  plot_s <- function(thresh=15){
    mem2=mem;
    cls=as.numeric(gsub("V","",names(A$gs.mean)))[sel]
    for(i in 1:length(cls)){
      mem2[mem==cls[i]]=i
    }
    
    mem2[!mem%in%cls]=-1
    sum_act=apply(s[selection+1,],2,sum)
    ord2=order(mem2[mem2!=-1])
    
    layout(matrix(1:2,2,1))
    par(mar=c(0.1,0.1,0.1,0.1));
    image(1-t(as.matrix(s)[selection[mem2!=-1][ord2]+1,sum_act>thresh]),frame=F,xaxt='n',yaxt='n',col=gray.colors(100))
    image(1-t(as.matrix(s)[selection[mem2==-1]+1,sum_act>thresh]),frame=F,xaxt='n',yaxt='n',col=gray.colors(100))
  }
  
  plot2_s <- function(print=F,width=200){
    if(print) png(paste("figures2/",prefix,"_smatrix.png",sep=""),15*width,5*200) 
    stmp1=as.matrix(s)
    stmp1_vec = as.numeric(stmp1)
    u1=rep(1:ncol(stmp1),1,each=nrow(stmp1))
    v1=rep(1:nrow(stmp1),ncol(stmp1))
    plot(u1[stmp1_vec>0],v1[stmp1_vec>0],pch=19,cex=0.1,col=viridis(10)[4],axes=F,xlab="",ylab="")
    if(print) dev.off();
  }
  
  plot2_cells_beh <- function(i,b="rs",print=F,width=200,len=NA,
                              col=viridis(10)[4]){
    if(print) png(paste("figures2/",prefix,"_enIDs_",i,".png",sep=""),15*width,5*200*11/19.)
    layout(matrix(1:3,3,1),heights=c(5,1,3))
    
    if(is.na(len)) len=ncol(s)
    btraj=colMeans(s[selection+1,][mem==i & probs>0.99,]);
    
    tabsel=as.numeric(names(sort(table(mem),decreasing=T))[sel])
    enID=which(tabsel==i)
    stmp1=as.matrix(s)[selection+1,][mem==i & probs>0.99,]
    stmp1_vec=as.numeric(stmp1)
    u1=rep(1:ncol(stmp1),1,each=nrow(stmp1))
    v1=rep(1:nrow(stmp1),ncol(stmp1))
    omega_row=omega_extended[enID,1:len]
    
    par(mar=c(0,5,2,1))
    image(xlim=c(0,len),ylim=c(1,nrow(stmp1)),z=matrix(NA,2,2),xlab="",ylab="",xaxt='n',yaxt='n')
    points(u1[stmp1_vec>0],v1[stmp1_vec>0],pch=19,cex=0.1,col=col)
    par(mar=c(0,5,0,1))
    image(xlim=c(0,len),ylim=c(0,1.2),z=matrix(NA,2,2),xlab="",ylab="",xaxt='n',yaxt='n')
    abline(v=which(omega_row>0.8)[1:len],col=gray.colors(100))
    par(mar=c(4,5,0,1))
    image(xlim=c(0,len),ylim=range(btraj[1:len]),z=matrix(NA,1,2),cex.axis=2,xlab="",ylab="",xaxt='n',las=2)
    lines(btraj[1:len])
    if(print) dev.off()
  }
  
  group_plots <- function(start){
    if(any(is.na(start))){
      X11(); plot_assemblies(-1,.99);
      X11(); plot_omega();
      X11(); P<-read.table(paste(folder,"/P.dat",sep=""))$V1; barplot(table(P)); 
    } else {
      dev.set(start[1]); plot_assemblies(-1,.99);
      dev.set(start[2]); plot_omega();
      dev.set(start[3]); P<-read.table(paste(folder,"/P.dat",sep=""))$V1; barplot(table(P));
    }
  }
  
  get_props <- function(x) return(c(A$pmu.mean[x],A$lambda0.mean[x],A$lambda1.mean[x]))
  
  plot_subfig <- function(){
    plot(A$pmu.mean[A$gs.mean>gs_thresh],A$lambda1.mean[A$gs.mean>gs_thresh],
         pch=19,xlab="activity", ylab="coherence",cex.lab=1.2,cex.axis=1.2,cex.lab=1.5); 
    lines(x=rep(pmu_thresh,2), y=c(0,1),lty=2,col="green",lwd=2)
    lines(x=c(0,1), y=rep(lambda1_thresh,2),lty=2,col="green",lwd=2)
    
  }
  
  plot_graph <- function(print=F,th=0,snc=FALSE,nav=nsamples){
    corrmat=cor(t(get_omega(nav=nav)))
    corrmat[corrmat<th]=0
    diag(corrmat)=0
    if(snc){
      graph<-graph.adjacency(corrmat[order(subnet_col),order(subnet_col)],weighted=TRUE,mode="lower")
    } else {
      graph<-graph.adjacency(corrmat,weighted=TRUE,mode="lower")
    }
    
    # 1 is left, 2 is right
    shapes=c("circle","square")
    if(snc){
      V(graph)$color=viridis(3)[subnet_col[order(subnet_col)]]
      V(graph)$shape=shapes[side_AP_ord$RL[order(subnet_col)]+1]
    } else {
      V(graph)$color[!is.na(side_AP_ord$RL)]=c("#ff9500","#0059ff")[side_AP_ord$RL[!is.na(side_AP_ord$RL)]+1]
    }
    V(graph)$label=ensel
    
    if(print) pdf(paste(folder,"/graph.pdf",sep=""))
    plot(graph)
    if(print) dev.off();
  }
  
  #cordata=get_cormat_bayes(nav=200,th=0.05)
  #graph <- cordata$graph 
  #wc <- edge.betweenness.community(graph, weight=get.data.frame(graph)$weight, directed = FALSE, bridges=TRUE)
  omega_extended = matrix(0,nrow=nrow(omega),ncol=ncol(s))
  for(i in 1:nrow(omega)) omega_extended[i,colSums(s[selection+1,])>15]=omega[i,]
  
  #mem2=mem
  #mem2[ !mem2 %in% ensel] = -1
  #mem3=mem2; for(i in ensel) mem3[mem2==i]=(which(ensel[order(subnet_mem)]==i))
  #write.table(mem3,file=paste(folder,"/mem3.dat",sep=""),col.names=F,row.names=F)
  
  plot_subnet <- function(print=F,nav=200,th=0.05){
    subnet_mem<-wc$membership
    
    # 1 is left, 2 is right
    shapes=c("circle","square")
    V(graph)$shape[!is.na(side_AP_ord$RL)]=shapes[side_AP_ord$RL[!is.na(side_AP_ord$RL)]+1]
    V(graph)$shape[is.na(side_AP_ord$RL)]="none"
    V(graph)$label.color="black"
    
    if(print) pdf(paste("figures2/subnets_",prefix,".pdf",sep=""))
    
    plot(wc,
         graph,
         mark.col=viridis(max(wc$membership),alpha=0.15), #shaded community colors
         mark.border=viridis(max(wc$membership)), # border colors of the communities
         col=(viridis(max(wc$membership),alpha=0.4)[wc$membership])) # vertex colors
    
    if(print) dev.off();
  }
  
  get_feature_df<- function(){
    return(data.frame(gs=A$gs.mean[sel],pmu=A$pmu.mean[sel],lambda0=A$lambda0.mean[sel],lambda1=A$lambda1.mean[sel]))
  }
  
  get_cormat_bayes<- function(nav=10,PMAX=200,th,disp=FALSE){
    system(paste("tail -n ",PMAX*nav," ",folder,"/omega_traj.dat > ",folder,"/omega_traj_tmp.dat",sep=""))
    
    corsel=matrix(0,length(ensel),length(ensel))
    cortab=array(0, dim=c(length(ensel),length(ensel),nav))
    omegatraj<-as.matrix(read.table(paste(folder,"/omega_traj_tmp.dat",sep="")))
    
    times=ncol(omegatraj)
    samples=nrow(omegatraj)/PMAX
    
    omegaprob=matrix(0,ncol=times,nrow=PMAX)
    counter=0;
    for(i in (samples-nav):(samples-1)){
      counter=counter+1
      omegasample=omegatraj[ i*PMAX + 1:PMAX,]
      omegasample=omegasample[ensel,]
      corsample=cor(t(omegasample))
      corsample[is.na(corsample)]=0
      cortab[,,counter]=corsample
      corsel[corsample>th]=corsel[corsample>th]+1
      if(disp) image(corsel,col=gray.colors(100))
    }
    
    cormean=apply(cortab,c(1,2),mean)
    ad=cormean;
    ad[corsel/nav<0.9]=0
    diag(corsel)=0
    diag(ad)=0
    
    graph <- graph.adjacency(ad,weighted=TRUE,mode="lower")
    V(graph)$color[!is.na(side_AP_ord$RL)]=c("#ff9500","#0059ff")[side_AP_ord$RL[!is.na(side_AP_ord$RL)]+1]
    V(graph)$label=ensel
    
    return(list(cortab=cortab,corsel=corsel/nav,ad=ad,graph=graph))
  }
  
  
  
  cha<-function(x,y){
    chull(x,y)->i
    return(areapl(cbind(x[i],y[i])))
  }
  
  total_area <- function(minprob = 0.99){
    areas = vector("numeric", length = length(as.numeric(names(tab)[sel])))
    densities = vector("numeric", length = length(as.numeric(names(tab)[sel])))
    counter = 1
    for(i in as.numeric(ensel)){
      areas [counter] = cha(centers[mem==i & probs>minprob,1], centers[mem==i & probs>minprob,2])
      densities [counter] = length(centers[mem==i & probs>minprob,1])/(areas [counter])
      counter = counter +1
    }
    return(data.frame(areas, densities))
  }
  
  eigen_decomp_of_cov <- function() {
    whole_tectal_area <- cha(centers [,1], centers [,2])
    counter = 1
    norm_eigen_decomp <- vector("numeric", length = length(as.numeric(names(tab)[sel])))
    for (i in as.numeric(ensel)) {
      assembly_mem <- centers[mem==i & probs>.99, 1:2]
      norm_eigen_decomp [counter] <- (eigen(cov(assembly_mem))$values [1]*eigen(cov(assembly_mem))$values [2])/whole_tectal_area
      counter = counter +1
    }
    return(cbind(rep(whole_tectal_area, length(as.numeric(ensel))), norm_eigen_decomp))
  }
  
  Lateral_index <- function(){
    LI <- vector("numeric", length = length(ensel))
    counter = 1
    for (i in as.numeric(ensel)){
      assembly_mem <- centers[mem==i & probs>.99, 1:2]
      LH = sum(assembly_mem [,2] > (midline$slope [1]*assembly_mem [,1] +midline$intercept [1]))
      RH <- sum(assembly_mem [,2] < (midline$slope [1]*assembly_mem [,1] +midline$intercept [1]))
      ALL <- nrow(assembly_mem)
      LI [counter] <- ((LH-RH)/ALL)
      counter = counter +1
      
    }
    return(LI)
  }
  
  within_cluster_corr <- function(){
    ensemble_mean <- vector("numeric", length = length(ensel))
    counter = 1
    spikes <- s [selection +1,]
    for (i in as.numeric(ensel)){
      assembly_mem <- centers[mem==i & probs>.99, 1:2]
      cor_mat <- cor(t(spikes[mem == i,]))
      lower <- cor_mat[lower.tri(cor_mat)]
      ensemble_mean [counter] <- mean(lower)
      counter = counter+1
    }
    return(ensemble_mean)
  }
  
  
  
  find_midline <- function(){
    plot_points <- function(data, a, b){
      smoothScatter(data[,1],data[,2], transformation = function(x) x^.25)
      points(data[,1],data[,2], cex= 0.2)
      abline(a = a, b = b, col ="red", lwd =5)
    }
    manipulate(plot_points(all_centers [,1:2], a, b), a =slider(-100,300), b = slider(-3,-0.5))
  }
  
  
  plot_midline <- function() {
    smoothScatter(centers[,1],centers[,2], transformation = function(x) x^.25)
    points(centers[,1],centers[,2])
    abline(a = midline [prefix_ID,1], b= midline [prefix_ID,2], col = "red", lwd = 2 )
  }
  
  raster_plot <- function(stmp1, m = "WT_NORM") {
    # stmp1 is the binary activity matrix
    circular_permutation <- function(t) {
      permute_by <- round(runif(1,1, length(t)))
      pt <- c(tail(t, -permute_by), head(t, permute_by))
      return(pt)
    }
    stmp1<- as.matrix(stmp1)
    stmp1_vec = as.numeric(unlist(stmp1))
    u1=rep(1:ncol(stmp1),1,each=nrow(stmp1))
    v1=rep(1:nrow(stmp1),ncol(stmp1))
    layout(matrix(c(1,1,1,1,1,1,2,2), nrow = 4, ncol = 2, byrow = TRUE))
    par(mar = c(0,6,5,2))
    plot(u1[stmp1_vec>0],v1[stmp1_vec>0],pch=19,cex=0.0001,col="black",cex.main=3, xaxt='n', xaxs = "i", yaxs ="i", ylim = c(0,nrow(stmp1)), ylab =  "Cell ID", cex.axis = 2, cex.lab = 3)
    title(main = list( m, cex = 2.5), adj = 0.5, line = 2)
    par(mar = c(6,6,0,2))
    plot((seq(1,ncol(s))/4.85/60), colSums(stmp1), type ="h", axes=T, xaxs = "i", yaxs ="i", xlab = "Time (mins)", ylim = c(0,100), cex.lab=2.5, cex.axis=2, ylab = "Cell #" )
    abline(h = quantile(colSums(t(apply(stmp1, 1, circular_permutation))), .95), col ="red", lwd = 3)
  }
  
  png(paste(folder, "/maps.png", sep =""), width = 1500, height = 1000)
  plot_assemblies(-1, .99)
  dev.off()
  
  
  png(paste(folder, "/omega.png", sep =""), width = 1500, height = 1000)
  plot_omega()
  dev.off()
  
  
  png(paste(folder, "/raster.png", sep =""), width = 1500, height = 1000)
  raster_plot(s, paste(folder))
  dev.off()
  
  
  
  write.table(omega, file = paste(folder, "/omega.dat", sep =""))
  write.table(omega_extended, file = paste(folder, "/omega_extended.dat", sep =""))
  write.table(eigen_decomp_of_cov(), file = paste(folder, "/eigen_decomp_of_cov.dat", sep ="") )
  write.table(total_area(), file = paste(folder, "/ensem_area.dat", sep =""))
  write.table(within_cluster_corr(), file = paste(folder, "/wcc.dat", sep ="")) 
  write.table(apply(omega_extended, 1, function(x) {sum(x)/ncol(omega_extended)}), file = paste(folder, "/Assembly_Freq.dat", sep ="")) 
  write.table(x = Lateral_index(), file = paste(folder, "/LI.dat", sep =""))
}



find_ant_post <- function(){
  plot_points <- function(data, a, b){
    smoothScatter(data[,1],data[,2], transformation = function(x) x^.25)
    points(data[,1],data[,2])
    abline(a = a, b = b, col ="red", lwd =5)
  }
  manipulate(plot_points(centers [,1:2], a, b), a =slider(-100,300), b = slider(0.5,2))
}

rotate_tect_by_line <- function(centers = all_centers) {
  angle <- atan(-1) +pi
  rot_angle <- angle - (pi/2)
  trans_mat <- matrix(c(cos(rot_angle), sin(rot_angle), -sin(rot_angle), cos(rot_angle)), 2, 2)
  rot_cent <- as.matrix(centers [,1:2]) %*% trans_mat
  shift_cent <- 10* cos(angle)
  rot_cent [,1] <- rot_cent [,1] - shift_cent
  return(rot_cent)
}


mean_by_bin <- function(centers) {
  new_cent <- rotate_tect_by_line(centers = centers )
  RL <- list(new_cent [new_cent [,1] > 0,], new_cent [new_cent [,1] < 0,])
  theta = c(pi/-3, pi/3)
  rotRL <- RL
  means <- list("numeric", length(theta))
  line <- list("numeric", length(theta))
  for (i in 1:length(theta)) {
    rotRL [[i]] <- cbind(Rotate(RL [[i]], theta = theta [i])$x, Rotate(RL [[i]], theta = theta [i])$y)
    cuts <- cut(x = rotRL [[i]] [,1], breaks = 10)
    means [[i]] <- cbind(by(data = rotRL [[i]] [,1], INDICES = cuts, FUN =function(x) {mean(x)}), by(data = rotRL [[i]] [,2], INDICES = cuts, FUN =function(x) {mean(x)}))
    line [[i]] <- cbind(Rotate(approx(means [[i]], n = 100), theta = theta [-i])$x,Rotate(approx(means [[i]], n = 100), theta = theta [-i])$y)
  }
  
  return(list(new_cent, line, RL))
} 

