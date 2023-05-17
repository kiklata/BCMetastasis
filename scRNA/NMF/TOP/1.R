
step1=function(dir_input="count_data",dir_output="res1",k=3:10,iteration=200,feature=2000){
  
  library(tidyverse)

  files=dir(dir_input)
  for (i in files) {
    
    one.file.path=paste(dir_input,"/",i,sep = "")
    prefix=strsplit(i,"\\.")[[1]][1]
    k.choices=paste(k,collapse = " ")
    junk.file=paste(dir_output,"/",prefix,"/","cnmf_tmp/",prefix,".spectra.k_*.iter_*.df.npz",sep = "")
    
    system(paste("python -W ignore cnmf.py prepare --output-dir ",dir_output," --name ",prefix," -c ",one.file.path," -k ",k.choices," --n-iter ",iteration," --total-workers 1 --numgenes ",feature,sep=""))
    system(paste("python -W ignore cnmf.py factorize --output-dir ",dir_output," --name ",prefix," --worker-index 0",sep = ""))
    system(paste("python -W ignore cnmf.py combine --output-dir ",dir_output," --name ",prefix,sep = ""))
    system(paste("rm -f ",junk.file,sep = ""))
    system(paste("MPLBACKEND='Agg' python -W ignore cnmf.py k_selection_plot --output-dir ",dir_output," --name ",prefix,sep = ""))
    
  }
  
  dir.create( paste(dir_output,"/k_selection",sep = ""))
  file.create(paste(dir_output,"/k_selection.txt",sep = ""))
  path3=paste(dir_output,"/k_selection.txt",sep = "")
  
  for (i in str_replace(files,"\\..*$","")) {
    path1=paste(dir_output,"/",i,"/",i,".k_selection.png",sep = "")
    path2=paste(dir_output,"/k_selection/",i,".k_selection.png",sep = "")
    system(paste("cp ",path1," ",path2,sep = ""))
    
    cat(i,",\n",sep = "",file=path3,append = T)
  }
  
}