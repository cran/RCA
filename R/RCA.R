RCA <-
function(inputMatrix, bootstrap=1000, p_value=0.05) { 
    # If the library hasn't been loaded yet, load it 
    PKG_PATH=system.file("libs", .Platform$r_arch, package="RCA", mustWork=TRUE)
    PKG_OBJ=paste(PKG_PATH,"RCA.so",sep="")
    if (!is.loaded('RCA_R_wrapper')) { 
        dyn.load(PKG_OBJ) 
    } 
    # Call the C function. A list of parameters values after the function is called 
    # is returned, assigned to the same names as are given before the 
    # = signs in the arguments. 
    # Change 2D matrix data to 1D data
	Dim=dim(inputMatrix) # Dim=[row,col]
	inputData=double(Dim[1]*Dim[2])
	for(i in 1:Dim[1]){
        for(j in 1:Dim[2]){
            inputData[j+(i-1)*Dim[2]]=inputMatrix[i,j]
        }
    }
    # Calling RCA-R function
    returned_data = .C('RCA_R_wrapper', inputData=inputData, numObs=as.integer(Dim[1]), numVars=as.integer(Dim[2]), bootstrap=as.integer(bootstrap), z_score=as.double(qnorm(1-p_value/2)), member=integer(Dim[1]), mergeDim=integer(2), mod=as.double(0), stats=double(2), merge=matrix(0,1,1), PACKAGE="RCA")
    # get membership tag one by one because it is too big an array for R to handle at the same time. 
    for(i in 0:(Dim[1]-1)){
        each=.C('get', i=as.integer(i), element=as.integer(1), PACKAGE="RCA")
        returned_data$member[i+1]=each$element
    }
    
    # C only accepts 1D input. We have to transform merge result here.
    row=returned_data$mergeDim[1]
    col=returned_data$mergeDim[2]
    merges2D=matrix(0,row,col)
    mergeResult=.C('getMerges', merges=integer(row*col), PACKAGE="RCA")
    for(i in 1:row){
        for(j in 1:col){
            merges2D[i,j]=mergeResult$merges[(j-1)*row+i]
        }
    }
    result=list(member=returned_data$member, mod=returned_data$mod, merge=merges2D)
    # Return the value of the result parameter 
    return(result) 
}

