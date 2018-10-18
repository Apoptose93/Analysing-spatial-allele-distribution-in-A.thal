# Analysing-spatial-allele-distribution-in-A.thal

"SampleX" contains a short sample SNP-matrix for tests of functionality

"SampleIndex" contains the corresponding index to SampleX when functioning properly. 

"snps_2029.rda" contains information about positon, genes, non-syn/syn mutation. 

"European_lat_long_eci.rda" contains location data of most sequences ecotypes in latitude/longitude for 2-dimensional analysis, but also ECI (earth centered inertial) for 3-dimensional analysis. 




"Index_and_Mapping.r" is the R-script containing the following 3 functions:

generate_mat() filters the SNP-matrix by overlap with location data and merged them together

Marker_test() calculates MAC(minor allele count), mean locations of 0(reference allele) and 1(mutated allele) group, geographical distancs between these groups aswell as the values of peacock-test and the corresponding p-value as indicator for geographical clustering.

locate_snp() uses the ggmap package and google maps service to plot all ecotypes of an selected snp on the map colorized in blue for 0(reference allele) and red for 1(mutated allele)
