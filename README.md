# Analysing-spatial-allele-distribution-in-A.thal

"SampleX" contains a short sample SNP-matrix for tests of functionality

"SampleIndex" contains the corresponding index to SampleX when functioning properly. 

"SNPs1 - SNPs6" contains information about positon, genes, non-syn/syn mutation. 
for uploading, SNPs data has been splitted into  parts, which needs to be merged to SNPs again.
To load it into R, you can use following lines:


SNP <- vector()
SNPs <- vector()
for ( i in 1:6) {
  load(paste("SNPs", i , ".rda", sep = ""))
  SNP <- rbind(SNP, SNPs)
}
SNPs <- SNP
remove(SNP)





"European_lat_long_eci.rda" contains location data of most sequences ecotypes in latitude/longitude for 2-dimensional analysis, but also ECI (earth centered inertial) for 3-dimensional analysis. 


"Index_and_Mapping.r" is the R-script containing the following 3 functions:

generate_mat() filters the SNP-matrix by overlap with location data and merged them together

Marker_test() calculates MAC(minor allele count), mean locations of 0(reference allele) and 1(mutated allele) group, geographical distancs between these groups aswell as the values of peacock-test and the corresponding p-value as indicator for geographical clustering.

locate_snp() uses the ggmap package and google maps service to plot all ecotypes of an selected snp on the map colorized in blue for 0(reference allele) and red for 1(mutated allele)
