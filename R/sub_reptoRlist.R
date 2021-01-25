#function to read cout output from hydra_sim into an R list


sub_reptoRlist = function(fn) {
  # reads in filename, fn, and interprets evrything as a character
  ifile=scan(fn,what="character",flush=T,blank.lines.skip=F,quiet=T)
  
  # tries to convert to double."Real" characters wont. 
  # We want these since they are variable names
  idx=sapply(as.double(ifile),is.na)
  vnam=ifile[idx] #list the variable names
  nv=length(vnam) #number of names
  A=list()
  ir=0
  # loop through the list of varoable names and get the data associated with each
  for(i in 1:nv)
  {
    ir=match(vnam[i],ifile)
    if(i!=nv){ # not on last element.count lines of data between variable names
      irr=match(vnam[i+1],ifile) 
    } else {
      irr=length(ifile)+1
    }
    dum<-NA
    # if a single line, read into a vector
    if(irr-ir==2) dum<-as.double(scan(fn,skip=ir,nlines=1,quiet=T,what=""))
    # otherwise read into a matrix
    if(irr-ir>2) dum<-as.matrix(read.table(fn,skip=ir,nrow=irr-ir-1,fill=T))
    #Logical test to ensure dealing with numbers. if so add to list
    if(is.numeric(dum)){
      A[[ vnam[i ]]]<-dum
    }
  }
  # a list of objects, each accessible by its variable name
  return(A)
}