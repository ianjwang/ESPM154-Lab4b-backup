#' Calculate genetic distances
#'
#' @param gen path to vcf file, a `vcfR` type object, or a dosage matrix
#' @param dist_type the type of genetic distance to calculate (options: `"euclidean"` (default), `"bray_curtis"`, `"dps"` for proportion of shared alleles (requires vcf), `"plink"`, or `"pc"` for PC-based)
#' @param plink_file if `"plink"` dist_type is used, path to plink distance file (typically ".dist"; required only for calculating plink distance). File must be a **square** distance matrix.
#' @param plink_id_file if `"plink"` dist_type is used, path to plink id file (typically ".dist.id"; required only for calculating plink distance)
#' @param npc_selection if `dist_type = "pc"`, how to perform K selection (options: `"auto"` for automatic selection based on significant eigenvalues from Tracy-Widom test (default), or `"manual"` to examine PC screeplot and enter no. PCs into console)
#' @param criticalpoint if `dist_type = "pc"` used with `npc_selection = "auto"`, the critical point for the significance threshold for the Tracy-Widom test within the PCA (defaults to 2.0234 which corresponds to an alpha of 0.01)
#'
#' @details
#' Euclidean and Bray-Curtis distances calculated using the ecodist package: Goslee, S.C. and Urban, D.L. 2007. The ecodist package for dissimilarity-based analysis of ecological data. Journal of Statistical Software 22(7):1-19. DOI:10.18637/jss.v022.i07.
#' Proportions of shared alleles calculated using the adegenet package: Jombart T. and Ahmed I. (2011) adegenet 1.3-1: new tools for the analysis of genome-wide SNP data. Bioinformatics. doi:10.1093/bioinformatics/btr521.
#' For calculating proportions of shared alleles, missing values are ignored (i.e., prop shared alleles calculated from present values; no scaling performed)
#'
#' @return pairwise distance matrix for given distance metric
#' @export
gen_dist <- function(gen = NULL, dist_type = "euclidean", plink_file = NULL, plink_id_file = NULL,
                     npc_selection = "auto", criticalpoint = 2.0234) {

  # Process input data ------------------------------------------------------
  # Read in vcf if path provided
  if (is.character(gen)) gen <- vcfR::read.vcfR(gen)
  # Convert vcf to dosage matrix
  if (inherits(gen, "vcfR") & (dist_type == "euclidean" | dist_type == "bray_curtis" | dist_type == "pc")) {
    gen <- vcf_to_dosage(gen)
    # Perform imputation with warning
    if (any(is.na(gen))) {
      gen <- simple_impute(gen, median)
      warning("NAs found in genetic data, imputing to the median (NOTE: this simplified imputation approach is strongly discouraged. Consider using another method of removing missing data)")
    }
  }

  # Calculate Euclidean distances -------------------------------------------
  if (dist_type == "euclidean") {
    # Check for NAs
    if (any(is.na(gen))) {
      stop("NA values found in genetic data")
    }
    dists <- ecodist::distance(gen, method = "euclidean")
    dists <- as.matrix(dists)
    return(as.data.frame(dists))
  }

  # Calculate Bray-Curtis distances -----------------------------------------
  if (dist_type == "bray_curtis") {
    # Check for NAs
    if (any(is.na(gen))) {
      stop("NA values found in genetic data")
    }

    dists <- ecodist::distance(gen, method = "bray-curtis")
    dists <- as.matrix(dists)
    return(as.data.frame(dists))
  }

  # Calculate proportion of shared alleles ----------------------------------
  if (dist_type == "dps") {
    if (!inherits(gen, "vcfR")) stop("VCF file required for calculating DPS distances")
    # Convert to genind
    genind <- vcfR::vcfR2genind(gen)
    # Show DPS warning about previously incorrect calculation
    dps_warning()
    dists <- 1 - adegenet::propShared(genind)
    return(as.data.frame(dists))
  }

  # Process Plink distance output files -------------------------------------
  if (dist_type == "plink") {
    if (is.null(plink_file)) stop("No plink distance file provided")
    dists <- as.data.frame(readr::read_tsv(plink_file, col_names = FALSE))
    plink_names <- readr::read_tsv(plink_id_file, col_names = FALSE) %>%
      dplyr::select(-`X1`) %>%
      as.matrix()
    # Assign row and col names according to sampleID
    rownames(dists) <- plink_names
    colnames(dists) <- plink_names
    return(dists)
  }

  # PC-based dist -----------------------------------------------------------
  if (dist_type == "pc") {
    # Check for NAs
    if (any(is.na(gen))) {
      stop("NA values found in genetic data")
    }

    gl <- adegenet::as.genlight(gen)

    # Perform PCA
    pc <- stats::prcomp(gl)

    # Get eig
    eig <- pc$sdev^2

    # Automatic npc selection based on number of significant eigenvalues
    if (npc_selection == "auto") {
      # Run Tracy-Widom test
      # NOTE: critical point corresponds to significance level.
      # If the significance level is 0.05, 0.01, 0.005, or 0.001,
      # the criticalpoint should be set to be 0.9793, 2.0234, 2.4224, or 3.2724, respectively.
      # The default is 2.0234.
      tw_result <- tw(eig, eigenL = length(eig), criticalpoint = criticalpoint)
      npc <- tw_result$SigntEigenL
    }

    # Manual npc selection: screeplot printout and selecting no. PCs to retain
    if (npc_selection == "manual") {
      stats::screeplot(pc, type = "barplot", npcs = 10, main = "PCA Eigenvalues")
      npc <- as.numeric(readline(prompt = "Number of PC axes to retain:"))
    }

    # Calculate PC-based distance
    dists <- as.matrix(dist(pc$x[, 1:npc], diag = TRUE, upper = TRUE))

    return(as.data.frame(dists))
  }
}

#' Plot the relationship between two distance metrics
#'
#' @param dist_x df containing square distance matrix for x axis
#' @param dist_y df containing square distance matrix for y axis
#' @param metric_name_x name of distance metric for x axis; if DPS used, must be `"dps"`
#' @param metric_name_y name of distance metric for y axis; if DPS used, must be `"dps"`
#'
#' @return scatterplot comparing two user-defined genetic distance metrics
#' @export
gen_dist_corr <- function(dist_x, dist_y, metric_name_x, metric_name_y) {
  if (!is.null(dist_x)) if (!inherits(dist_x, "data.frame")) dist_x <- as.data.frame(dist_x)
  if (!is.null(dist_y)) if (!inherits(dist_y, "data.frame")) dist_y <- as.data.frame(dist_y)

  # Check to ensure sample IDs match ----------------------------------------
  if (all(rownames(dist_x) == rownames(dist_y)) == FALSE) {
    stop("Sample IDs do not match")
  }

  # Melt data from square to long -------------------------------------------
  # Assign NAs to upper triangle of square matrix
  dist_x[upper.tri(dist_x, diag = FALSE)] <- NA

  melt_x <- dist_x %>%
    tibble::rownames_to_column(var = "comparison") %>%
    tidyr::pivot_longer(cols = -(comparison)) %>%
    na.omit() %>%
    dplyr::filter(comparison != name) %>%
    dplyr::rename(!!metric_name_x := value)

  melt_y <- dist_y %>%
    tibble::rownames_to_column(var = "comparison") %>%
    tidyr::pivot_longer(cols = -(comparison)) %>%
    na.omit() %>%
    dplyr::filter(comparison != name) %>%
    dplyr::rename(!!metric_name_y := value)

  # Build plots -------------------------------------------------------------
  if (metric_name_x == "dps" || metric_name_y == "dps") {
    joined <- dplyr::full_join(melt_x, melt_y) %>%
      dplyr::mutate(rev_dps = (1 - dps))
    joined %>%
      ggplot2::ggplot(ggplot2::aes_string(x = metric_name_x, y = metric_name_y)) +
      ggplot2::geom_abline(ggplot2::aes(intercept = 0.0, slope = 1), color = "gray") +
      ggplot2::geom_point(color = "black", size = .2, alpha = .5)
  } else {
    joined <- dplyr::full_join(melt_x, melt_y)
    joined %>%
      ggplot2::ggplot(ggplot2::aes_string(x = metric_name_x, y = metric_name_y)) +
      ggplot2::geom_abline(ggplot2::aes(intercept = 0.0, slope = 1), color = "gray") +
      ggplot2::geom_point(color = "black", size = .2, alpha = .5)
  }
}

#' Make heatmap of genetic distances
#'
#' @param dist Matrix of genetic distances
#'
#' @return heatmap of genetic distances
#' @export
gen_dist_hm <- function(dist) {
  if (!is.null(dist)) if (!inherits(dist, "data.frame")) dist <- as.data.frame(dist)

  dist %>%
    tibble::rownames_to_column("sample") %>%
    tidyr::gather("sample_comp", "dist", -"sample") %>%
    ggplot2::ggplot(ggplot2::aes(x = sample, y = sample_comp, fill = dist)) +
    ggplot2::geom_tile() +
    ggplot2::coord_equal() +
    viridis::scale_fill_viridis(option = "inferno") +
    ggplot2::xlab("Sample") +
    ggplot2::ylab("Sample") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))
}


#' Convert a vcf to a dosage matrix
#'
#' @param x can either be an object of class `vcfR` or a path to a .vcf file
#'
#' @return dosage matrix
#' @export
vcf_to_dosage <- function(x) {
  if (!inherits(x, "vcfR")) {
    x <- vcfR::read.vcfR(x)
  }
  genlight <- vcfR::vcfR2genlight(x)
  gen <- as.matrix(genlight)
  return(gen)
}

#' Impute NA values
#' NOTE: use extreme caution when using this form of simplistic imputation. We mainly provide this code for creating test datasets and highly discourage its use in analyses.
#' @param x matrix
#' @param f function to use for imputation (defaults to median)
#'
#' @return matrix of values with missing values imputed
#' @export
#' @family Imputation functions
simple_impute <- function(x, FUN = median) {
  x_noNA <- apply(x, 2, impute_helper, FUN)
  return(x_noNA)
}

#' Helper function for imputation
#' @export
#' @noRd
#' @family Imputation functions
impute_helper <- function(i, FUN = median) {
  i[which(is.na(i))] <- FUN(i[-which(is.na(i))], na.rm = TRUE)
  return(i)
}

#' Imputation of missing values using population structure inferred with `LEA::snmf`
#'
#' @param gen a dosage matrix, an object of class 'vcfR', or an object of type snmfProject
#' @param quiet whether to operate quietly and suppress the results of cross-entropy scores (defaults to TRUE; only does so if more than one K-value); only displays run with minimum cross-entropy
#' @param save_output if TRUE, saves SNP GDS and ped (plink) files with retained SNPs in new directory; if FALSE returns object (defaults to FALSE)
#' @param output_filename if `save_output = TRUE`, name prefix for saved .geno file, sNMF project file, and sNMF output file results (defaults to FALSE, in which no files are saved)
#'
#' @inheritParams LEA::snmf
#'
#' @return dosage matrix with imputed missing values
#' @export
#' @family Imputation functions
str_impute <- function(gen, K, entropy = TRUE, repetitions = 10, project = "new", quiet = TRUE, save_output = FALSE, output_filename = NULL) {
  if (is.null(output_filename)) filename <- "tmp" else filename <- output_filename
  
  # Convert vcf to dosage
  if (inherits(gen, "vcfR")) gen <- vcf_to_dosage(gen)
  
  # Load sNMF project
  if (inherits(gen, "snmfProject")) snmf_proj <- gen
  
  # Convert gen to .geno type file unless snmfProject provided
  if (!inherits(gen, "snmfProject")) {
    geno <- gen_to_geno(gen)
    # sNMF requires an input file saved to file (cannot accept an R object)
    LEA::write.geno(geno, paste0(filename, ".geno"))
    
    # Run sNMF
    snmf_quiet <- purrr::quietly(LEA::snmf)
    quiet_result <- snmf_quiet(input.file = paste0(filename, ".geno"), K = K, entropy = entropy, repetitions = repetitions, project = project)
    snmf_proj <- quiet_result[[1]]
  }
  
  # Look through directories
  bestK <- snmf_bestK(snmf_proj, K = K, quiet = quiet)
  
  if (!quiet & length(K) > 1) {
    print(plot_crossent(bestK$ce_values))
  }
  
  # Impute missing values based on sNMF groupings
  impute_quiet <- purrr::quietly(LEA::impute)
  result <- impute_quiet(object = snmf_proj, input.file = paste0(filename, ".geno"), method = "mode", K = as.integer(bestK$K_value), run = as.integer(bestK$bestrun))
  
  # Read .lfmm file in
  imputed <- LEA::read.lfmm(paste0(filename, ".lfmm_imputed.lfmm"))
  
  # Add individual and variant names back in
  rownames(imputed) <- rownames(gen)
  colnames(imputed) <- colnames(gen)
  
  imputed <- geno_to_dosage(imputed)
  
  # Remove created files
  if (!save_output) {
    # Removes snmf project and associated snmf files
    LEA::remove.snmfProject(paste0(filename, ".snmfProject"))
    
    # Delete other associated files
    unlink(paste0(filename, ".geno"))
    unlink(paste0(filename, ".lfmm"))
    unlink(paste0(filename, ".lfmm_imputed.lfmm"))
  }
  return(imputed)
}

#' Helper function to select "best" K based on minimizing cross-entropy criteria from sNMF results
#'
#' @param snmf_proj object of type snmfProject
#' @param K integer corresponding to K-value
#' @param quiet whether to operate quietly and suppress the results of cross-entropy scores (defaults to TRUE; only does so if more than one K-value); only displays run with minimum cross-entropy
#'
#' @return list with best K-value and run number and all cross-entropy scores
#' @export
#' @family Imputation functions
#' @keywords internal
snmf_bestK <- function(snmf_proj, K, quiet) {
  if (length(K) == 1) {
    bestrun <- which.min(LEA::cross.entropy(snmf_proj, K = K))
    results <- list(K_value = K, bestrun = bestrun)
  }
  
  if (length(K) > 1) {
    ce_values <- as.data.frame(purrr::map(K, snmf_crossent_helper, snmf_proj = snmf_proj, select_min = FALSE))
    ce_values <-
      ce_values %>%
      tibble::rownames_to_column(var = "run") %>%
      tidyr::pivot_longer(names_to = "K_value", values_to = "cross_entropy", -run)
    ce_values$run <- stringr::str_replace_all(ce_values$run, "run ", "")
    ce_values$K_value <- stringr::str_replace_all(ce_values$K_value, "K...", "")
    
    best <- ce_values %>% dplyr::slice(which.min(cross_entropy))
    results <- list(K_value = best$K_value, bestrun = best$run, ce_values = ce_values)
    
    if (!quiet) {
      print(plot_crossent(ce_values))
    }
  }
  return(results)
}

#' Helper function to retrieve cross entropy scores from SNMF project
#'
#' @param snmf_proj object of type snmfProject
#' @param K K-value(s)
#' @param select_min whether to return minimum
#'
#' @return cross entropy scores for given K
#' @export
#' @family Imputation functions
#' @keywords internal
snmf_crossent_helper <- function(snmf_proj, K, select_min = TRUE) {
  if (select_min) results <- which.min(LEA::cross.entropy(snmf_proj, K = K))
  if (!select_min) results <- LEA::cross.entropy(snmf_proj, K = K)
  return(results)
}

#' Helper function to plot cross entropy scores from SNMF
#'
#' @param ce_values df with run, K-value, and cross entropy created in \link[algatr]{snmf_bestK}
#'
#' @return ggplots of cross entropy values compared to K-values (and across runs)
#' @export
#' @keywords internal
plot_crossent <- function(ce_values) {
  if (length(unique(ce_values$run)) == 1) {
    plt <-
      ce_values %>%
      dplyr::group_by(K_value) %>%
      dplyr::slice(which.min(cross_entropy)) %>%
      ggplot2::ggplot(ggplot2::aes(x = K_value, y = cross_entropy)) +
      ggplot2::geom_point(size = 3, color = "red") +
      ggplot2::theme_bw() +
      ggplot2::ylab("Cross entropy value") +
      ggplot2::xlab("K value") +
      ggplot2::geom_vline(xintercept = best$K_value, linetype = "dashed", color = "blue")
    print(plt)
  }
  
  if (length(unique(ce_values$run)) > 1) {
    plt <-
      ce_values %>%
      ggplot2::ggplot(ggplot2::aes(x = run, y = cross_entropy)) +
      ggplot2::geom_point(size = 3, color = "red") +
      ggplot2::theme_bw() +
      ggplot2::ylab("Cross entropy value") +
      ggplot2::xlab("Run") +
      ggplot2::facet_grid(~K_value)
    print(plt)
  }
}

#' Convert dosage matrix or vcf to geno type object (N.B.: this only works for diploids!)
#'
#' @inheritParams str_impute
#'
#' @return matrix encoded as geno type object
#' @export
#' @family Imputation functions
gen_to_geno <- function(gen) {
  if (!is.matrix(gen)) if (!inherits(gen, "vcfR")) {
    gen <- vcfR::read.vcfR(gen)
    gen <- vcf_to_dosage(gen)
  }
  if (inherits(gen, "vcfR")) gen <- vcf_to_dosage(gen)
  
  # Recode data for geno type object
  gen[gen == 2] <- 555 # tmp placeholder so 0s aren't overwritten
  gen[gen == 0] <- 2 # two ref alleles (vcf 0/0 or dosage 0)
  gen[gen == 555] <- 0 # no ref alleles (vcf 1/1 or dosage 2)
  gen[is.na(gen)] <- 9 # missing data encoded as 9
  
  return(gen)
}

#' Convert lfmm/geno matrix to dosage matrix (N.B.: this only works for diploids!)
#'
#' @param geno matrix of LEA geno or lfmm format (i.e., 0 corresponds to zero reference alleles)
#'
#' @return matrix encoded as dosage type object (0 corresponds to two reference alleles)
#' @export
#' @family Imputation functions
geno_to_dosage <- function(geno) {
  # Recode data for geno type object
  geno[geno == 2] <- 555 # tmp placeholder so 2s aren't overwritten
  geno[geno == 0] <- 2 # two ref alleles (vcf 0/0 or dosage 0)
  geno[geno == 555] <- 0 # no ref alleles (vcf 1/1 or dosage 2)
  geno[geno == 9] <- NA # missing data encoded as NA
  
  return(geno)
}

#' Calculate distance between environmental vars
#'
#' @param env dataframe or vector of environmental variables for locations
#' @param stdz if TRUE then environmental values will be standardized (default = TRUE)
#'
#' @return list of environmental distances between samples (for each environmental variable)
#' @export
env_dist <- function(env, stdz = TRUE) {
  if (!is.null(dim)) distmat <- dplyr::as_tibble(env) %>% purrr::map(env_dist_helper, stdz = stdz)
  if (is.null(dim)) distmat <- env_dist_helper(env, stdz)
  return(distmat)
}

#' Helper function to convert an environmental vector to a distance matrix
#'
#' @inheritParams env_dist
#'
#' @export
#' @noRd
env_dist_helper <- function(env, stdz = TRUE) {
  # Standardize environmental variables
  if (stdz) env <- scale(env, center = TRUE, scale = TRUE)
  
  distmat <- as.matrix(dist(env, diag = TRUE, upper = TRUE))
  
  return(distmat)
}

#' Calculate geographic distance between coordinates
#'
#' @param coords dataframe with x and y coordinates
#' @param type the type of geographic distance to be calculated; options are "Euclidean" for direct distance, "topographic" for topographic distances, and "resistance" for resistance distances.
#' @param lyr SpatRaster or Raster* DEM for calculating topographic distances or resistance raster for calculating resistance distances (RasterLayer or SpatRaster object)
#' @details
#' Euclidean, or linear, distances are calculated using the geodist package: Padgham M, Sumner M (2021). geodist: Fast, Dependency-Free Geodesic Distance Calculations. R package version 0.0.7, Available: https://CRAN.R-project.org/package=geodist.
#' Topographic distances are calculated using the topoDistance package: Wang I.J. (2020) Topographic path analysis for modeling dispersal and functional connectivity: calculating topographic distances using the TOPODISTANCE R package. Methods in Ecology and Evolution, 11: 265-272.
#' Resistance distances are calculated using the gdistance package: van Etten, J. (2017). R package gdistance: Distances and routes on geographical grids. Journal of Statistical Software, 76(1), 1–21.
#'
#' @return geographic distance matrix
#' @export
geo_dist <- function(coords, type = "Euclidean", lyr = NULL) {
  if (type == "Euclidean" | type == "euclidean" | type == "linear") {
    # Format coordinates
    coords <- coords_to_sf(coords)
    # Calculate geodesic distance between points
    distmat <- sf::st_distance(coords)
  } else if (type == "topo" | type == "topographic") {
    # Format coordinates
    coords <- coords_to_df(coords)
    
    if (is.null(lyr)) stop("Calculating topographic distances requires a DEM layer for argument lyr.")
    message("Calculating topo distances... This can be time consuming with many points and large rasters.")
    
    # Convert to RasterLayer if SpatRaster object
    if (inherits(lyr, "SpatRaster")) lyr <- raster::raster(lyr)
    
    distmat <- topoDistance::topoDist(lyr, coords, paths = FALSE)
  } else if (type == "resistance" | type == "cost" | type == "res") {
    if (is.null(lyr)) stop("Calculating resistance distances requires a resistance surface for argument lyr.")
    message("Calculating resistance distances... This can be time consuming with many points and large rasters.")
    
    # Format coordinates
    coords <- coords_to_df(coords)
    
    # Convert to RasterLayer if SpatRaster object
    if (inherits(lyr, "SpatRaster")) lyr <- raster::raster(lyr)
    
    # Convert resistance surface to conductance surface
    cond.r <- 1 / lyr
    trSurface <- gdistance::transition(cond.r, transitionFunction = mean, directions = 8) # Create transition surface
    trSurface <- gdistance::geoCorrection(trSurface, type = "c", scl = FALSE)
    sp <- sp::SpatialPoints(coords = coords)
    distmat <- as.matrix(gdistance::commuteDistance(trSurface, sp)) # Calculate circuit distances
  }
  
  return(distmat)
}

#' convert coordinates to sf
#' @noRd
coords_to_sf <- function(coords) {
  if (inherits(coords, "sf")) {
    return(coords)
  }
  if (inherits(coords, "SpatVector")) {
    return(sf::st_as_sf(coords))
  }
  if (is.matrix(coords)) coords <- data.frame(coords)
  if (is.data.frame(coords)) colnames(coords) <- c("x", "y")
  return(sf::st_as_sf(coords, coords = c("x", "y")))
}

#' Helper function to create statistic vector for gt() tables
#'
#' @param stat_name name of statistic
#' @param stat value of statistic
#' @param df dataframe you will eventually add statistic to
#'
#' @return vector with statistic value and name and NA fillers for remaining columns
#' @export
#' @noRd
#'
make_stat_vec <- function(stat_name, stat, df) {
  stat_vec <- rep(NA, ncol(df))
  stat_vec[1] <- stat_name
  stat_vec[2] <- stat
  names(stat_vec) <- colnames(df)
  return(stat_vec)
}

# convert from matrix/data.frame/sf to formatted df
coords_to_df <- function(coords) {
  if (inherits(coords, "sf")) coords <- sf::st_coordinates(coords)
  if (is.matrix(coords)) coords <- data.frame(coords)
  colnames(coords) <- c("x", "y")
  return(coords)
}