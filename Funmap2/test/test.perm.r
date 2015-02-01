library(Funmap2);

pheno_csv <- "LC.populus.pheno.csv";
geno_csv  <- "LC.populus.geno.csv";
marker_csv <- "LC.populus.marker.csv";
curve_type <- CURVE_LC;
cross_type <- CROSS_BC;

dat <- FM2.load_data(pheno_csv, geno_csv, marker_csv, curve_type, cross_type);
r <- FM2.permutation(dat, options=list(cluster_count=7, permu_loop=200))

