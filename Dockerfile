FROM rocker/shiny:4.3.3

# System dependencies for R package compilation
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    libglpk-dev \
    libhdf5-dev \
    libmagick++-dev \
    libgsl-dev \
    cmake \
    libnlopt-dev \
    libudunits2-dev \
    && rm -rf /var/lib/apt/lists/*

# Install renv and BiocManager
RUN R -e "install.packages(c('renv', 'BiocManager'), repos = 'https://cran.r-project.org')"

# Install Bioconductor packages (renv can't resolve archived Bioc URLs)
RUN R -e "BiocManager::install(version = '3.18', ask = FALSE)" && \
    R -e "BiocManager::install(c( \
      'Biobase', 'BiocGenerics', 'BiocParallel', 'BiocVersion', \
      'DelayedArray', 'DelayedMatrixStats', 'GenomeInfoDb', 'GenomeInfoDbData', \
      'GenomicRanges', 'IRanges', 'MatrixGenerics', \
      'rhdf5', 'rhdf5filters', 'Rhdf5lib', \
      'S4Arrays', 'S4Vectors', 'SingleCellExperiment', \
      'SparseArray', 'SpatialExperiment', 'SummarizedExperiment', \
      'UCell', 'XVector', 'beachmat', 'BiocFileCache', \
      'BiocNeighbors', 'glmGamPoi', 'HDF5Array', \
      'sparseMatrixStats', 'zlibbioc' \
    ), ask = FALSE, update = FALSE)"

# Copy lockfile for renv restore (CRAN + GitHub packages)
WORKDIR /pkg
COPY renv.lock renv.lock

# Restore exact CRAN + GitHub package versions from lockfile
RUN R -e " \
  bioc_pkgs <- c( \
    'Biobase', 'BiocGenerics', 'BiocParallel', 'BiocVersion', \
    'DelayedArray', 'DelayedMatrixStats', 'GenomeInfoDb', 'GenomeInfoDbData', \
    'GenomicRanges', 'IRanges', 'MatrixGenerics', \
    'rhdf5', 'rhdf5filters', 'Rhdf5lib', \
    'S4Arrays', 'S4Vectors', 'SingleCellExperiment', \
    'SparseArray', 'SpatialExperiment', 'SummarizedExperiment', \
    'UCell', 'XVector', 'beachmat', 'BiocFileCache', \
    'BiocNeighbors', 'glmGamPoi', 'HDF5Array', \
    'sparseMatrixStats', 'zlibbioc' \
  ); \
  renv::restore(lockfile = 'renv.lock', library = .libPaths()[1], prompt = FALSE, exclude = bioc_pkgs) \
"

# Copy package source and install
COPY . /pkg
RUN R -e "install.packages('remotes', repos = 'https://cran.r-project.org')" && \
    R -e "remotes::install_local('.', dependencies = FALSE)"

# Expose port
EXPOSE 3838

CMD ["R", "-e", "matrispace.app::run_app(host='0.0.0.0', port=3838)"]
