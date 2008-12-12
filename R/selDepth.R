`selDepth`<- function(x= NULL, notAt= NULL, at= NULL, baseline= NULL, out= "df", ties=TRUE){

    if(is.null(x)){
        stop("Vector of allele ratios is null!")
    }
    if(min(x)<0 | max(x)>1){
        warning("\nFound allele ratios below 0 or above 1!")
    }
    if(!is.null(notAt)){
		if(length(notAt)>length(x)){
	        stop("More markers to test in 'notAt' than markers present in 'x'.")
	    }
	    if(max(notAt)>length(x)){
	        stop("The index of a marker to test in 'notAt' exceeds the length of 'x'.")
	    }
    }
    if(!is.null(baseline)){
		if(length(baseline)>length(x)){
	        stop("More markers to baseline than markers present in 'x'.")
	    }
	    if(max(baseline)>length(x)){
	        stop("The index of a marker to baseline exceeds the length of 'x'.")
	    }
    }
    if(out %in% c("df", "pos", "depth")== FALSE){
        stop("Argument to 'out' should be one of 'df', 'pos', 'depth'.")
    }
        
    if(is.null(baseline)){
        meanAR<- mean(x)
    }   else {
            xx<- x[baseline]
            meanAR<- mean(xx)
        }
    
    if(!is.null(notAt)){
        x[c(notAt)]<- NA
    }
    if(!is.null(at)){
        x[-c(at)]<- NA
    }
    axv<- abs(x - meanAR)
    maxDepth<- max(axv, na.rm=TRUE)
    if(ties==TRUE){
        atPos<- which(axv==maxDepth)
        }else{
            atPos<- which(axv==maxDepth)[1]
        }
    depth<- x[atPos] - meanAR
    if(out=="df"){
        return(data.frame(maxDepth= depth, atPos= atPos))
    }else if(out=="pos"){            
        return(atPos)  
    }
    else if(out=="depth"){
        return(depth)
    }
}

        
