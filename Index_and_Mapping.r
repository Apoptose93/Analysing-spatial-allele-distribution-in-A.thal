
# Generalized All Marker Script

require ("ks")
require ("Peacock.test")
require ("geosphere")
require ("astrolibR")
require ("doMC")
require ("ff")

# setting the number of cores used with %dopar%
registerDoMC(2)

# X: SNP data matrix of Ecotypes
# W: Location data matrix of Ecotypes
# max_mark is optionally cutting the matrix smaller for testing purpose
# Ecoschnitt: vector of ecotypes both present in imputed and unimputed datasets for comparison

# The generate_raw function first cuts down the SNP matrix to the Ecotypes present in Ecoschnitt vector and 
# up to the number of max_mark. Then binding together matrices of location and SNP for easier calculation.
# Output is the foundation for Marker_test

generate_mat <- function(x=X, w =W, max_mark = ncol(X)){
  matrix <- X[,1:max_mark]
  #matrix <- matrix[rownames(matrix)%in%Ecoschnitt,]
  m <- paste(W[,3],"-", W[,4])
  Wx <- W[which(duplicated(m) != T),]
  matrix <- matrix[rownames(matrix)%in%Wx[,1],]
  names <- rownames(matrix)
  matrix <- cbind(Wx[Wx[,1]%in%rownames(matrix),3:7], matrix)
  rownames(matrix) <- names
  remove(names)
  
  return(matrix)
}

# matrix : output of generate_mat function
# calc.pval: turns off/on calculation of p-values, this step is by far most time consuming
# snp: foundation index including general information about location in gene, syn or non-syn mutation

Marker_test <- function(matrix = Matrix, calc.pval=T, snp = SNPs){
  
  ### makes sure, the snp index and matrix have the same intersection and order to combine later
  
  rownames(snp) <- paste(snp[,2],"- ",snp[,3], sep = "")
  snp <- snp[rownames(snp)%in%colnames(matrix),]
  
  matrix <- cbind(matrix[,1:5], matrix[,colnames(matrix)%in%rownames(snp)])

  ### First calculates the minor allele count (MAC) of the alleles by summing mutation and comparing to total 	   ### count
  MAC <- foreach (x = 6:ncol(matrix), .combine = "c") %dopar% {
    
    # MAC  / AC0, AC1
    Count1 <- sum(matrix[matrix[,x]%in%1,x])
    Count0 <- length(matrix[matrix[,x]%in%c(1,0),x]) - Count1
    
    if (Count1 > Count0) {
      Count0
    } else {Count1}
    
  }
   gc()
   
   ### Calculates the mean geographical location for all marker having an MAC greater 9, otherwise NA
   
     mittelpunkte <- foreach (x = 6:ncol(matrix), .combine = "rbind") %dopar% {

    if (MAC[x-5] > 9) {
      c(geomean(matrix[matrix[,x]%in%1,1:2])[1:2],
      geomean(matrix[matrix[,x]%in%0,1:2])[1:2])
    } else {
      rep(NA,4)
    }
     }  
   gc()
   
   ### Calculates the distance between the mean locations for all MAC > 9
   
    distanzen01 <- foreach (x = 6:ncol(matrix), .combine = "c") %dopar% {
      if (MAC[x-5] > 9) {
        distm(x=as.numeric( mittelpunkte[x-5,3:4]), y = as.numeric(mittelpunkte[x-5,1:2]))
        
      } else{ 
        NA
      }
    }
  gc()
  
  ### Calculates the p-value(kde.test) / value (peacock) for all marker with MAC > 9
  if (calc.pval == T) {
  values <- foreach (x = 6:ncol(matrix), .combine = "c") %dopar% {
    if(MAC[x-5] > 9 ) {
      #kde.test(x1 = matrix[matrix[,x]%in%1, 3:4],
      #         x2 = matrix[matrix[,x]%in%0, 3:4])$pvalue
               
      ## at this point, either kde.test, peacock2 or peacock3 can be used for different tests
      peacock2(x = matrix[matrix[,x]%in%1, 1:2],
                y = matrix[matrix[,x]%in%0, 1:2])
      #peacock3(x = matrix[matrix[,x]%in%1, 3:5],
      #          y = matrix[matrix[,x]%in%0, 3:5])
      
    } else{
      NA
      }
  }}
  gc()
  
  ### merging the results with the snp matrix
  snpindex <- cbind(snp[,1:8], values,distanzen01, MAC, mittelpunkte)
  remove(mittelpunkte, MAC, distanzen01,values)
  
  return (snpindex)
}

### Loading ressources
load("European_lat_long_eci.rda")
load("snps_2029.rda")


# generates the Matrix with generate_mat and creates an index using Marker_test()
Matrix <- generate_mat()
Index <- Marker_test()



locate_snp <- function(matrix = Matrix, snpindex = Index){
  
  G0lat <- mean(matrix[matrix[,snp]%in%0, 1])
  G0long <- mean(matrix[matrix[,snp]%in%0, 2])
  G0lat <- c(G0lat,mean(matrix[matrix[,snp]%in%1, 1]))
  G0long <- c(G0long, mean(matrix[matrix[,snp]%in%1, 2]))
  geomatrix <- cbind(G0lat, G0long)
  geomatrix <- data.frame(geomatrix)
  require(ggmap)
  europe_map<-get_map(location=c(left = min(matrix[matrix[,snp]%in%c(1,0), 2])-5, bottom = min(matrix[matrix[,snp]%in%c(1,0), 1])-5, 
                                 right = max(matrix[matrix[,snp]%in%c(1,0), 2])+5, top= max(matrix[matrix[,snp]%in%c(1,0), 1])+5),maptype = "satellite")
  ggmap(europe_map)+
    labs(x = "Longitude", y = "Latitude") +
    geom_point(aes(x = longitude, y = latitude), data = matrix[matrix[,snp]%in%1,], colour = "red", size = 1) +
    geom_point(aes(x = longitude, y = latitude), data = matrix[matrix[,snp]%in%0,], colour = "blue", size = 1) +
    geom_point(aes(x = G0long, y = G0lat), data = geomatrix, colour = "white", size = 3) +
    geom_path(aes(x = G0long, y = G0lat), data = geomatrix, colour = "white", size = .5) +
    annotate("text", x = min(matrix[matrix[,snp]%in%c(1,0), 2])-2, y = max(matrix[matrix[,snp]%in%c(1,0), 1])+1, size = 4, label = snp, colour = "white") +
    annotate("text", x = min(matrix[matrix[,snp]%in%c(1,0), 2])-1, y = max(matrix[matrix[,snp]%in%c(1,0), 1]), size = 4, label = paste("MAC:", snpindex[snp,"MAC"], "/", length(matrix[,snp]), sep = ""), colour = "white") +
    annotate("text", x = min(matrix[matrix[,snp]%in%c(1,0), 2])-1, y = max(matrix[matrix[,snp]%in%c(1,0), 1])-1, size = 4, label = paste("pvalue:", format(snpindex[snp,"pvalue"],digits = 4), "/", length(matrix[,snp]), sep = ""), colour = "white") 
  
  
}

# select snp and plot it using the previously generated Matrix of generate_mat
snp <- "1- 73"
locate_snp()