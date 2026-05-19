FROM rocker/r-ver:4.5.3 AS builder

ARG TARGETARCH

LABEL maintainer="Felipe Marques de Almeida <almeidafmarques@outlook.com>"
LABEL description="Asc-Seurat v3: Interactive scRNA-seq analysis"
LABEL version="3.0.6"

ENV DEBIAN_FRONTEND=noninteractive
ENV RETICULATE_PYTHON=/opt/venv/bin/python

RUN set -eux; \
    if [ "$(dpkg --print-architecture)" = "arm64" ]; then \
        find /etc/apt -type f \( -name '*.list' -o -name '*.sources' \) -exec sed -i 's|http://ports.ubuntu.com/ubuntu-ports|https://mirrors.mit.edu/ubuntu-ports/ubuntu-ports|g; s|https://ports.ubuntu.com/ubuntu-ports|https://mirrors.mit.edu/ubuntu-ports/ubuntu-ports|g' {} +; \
    fi; \
    for attempt in 1 2 3 4 5; do \
        apt-get -o Acquire::Retries=5 -o Acquire::http::Timeout=30 -o Acquire::https::Timeout=30 update \
        && apt-get -o Acquire::Retries=5 -o Acquire::http::Timeout=30 -o Acquire::https::Timeout=30 install -y --no-install-recommends \
            build-essential \
            ca-certificates \
            cmake \
            curl \
            gfortran \
            git \
            libabsl-dev \
            libcairo2-dev \
            libcurl4-openssl-dev \
            libfftw3-dev \
            libfontconfig1-dev \
            libfreetype6-dev \
            libfribidi-dev \
            libgit2-dev \
            libgdal-dev \
            libgeos-dev \
            libglpk-dev \
            libharfbuzz-dev \
            libhdf5-dev \
            libnode-dev \
            libpng-dev \
            libproj-dev \
            libpython3.12-dev \
            libssl-dev \
            libtiff5-dev \
            libudunits2-dev \
            libx11-dev \
            libxml2-dev \
            pandoc \
            python3 \
            python3-pip \
            python3-venv \
        && break; \
        if [ "$attempt" = "5" ]; then exit 1; fi; \
        rm -rf /var/lib/apt/lists/*; \
        sleep $((attempt * 15)); \
    done; \
    rm -rf /var/lib/apt/lists/*

RUN echo 'options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/noble/latest"), timeout = 600)' >> /usr/local/lib/R/etc/Rprofile.site

RUN python3 -m venv /opt/venv \
    && /opt/venv/bin/pip install --no-cache-dir --upgrade pip uv \
    && /opt/venv/bin/uv pip install --python /opt/venv/bin/python \
        scanpy anndata numpy scipy pandas leidenalg igraph \
    && /opt/venv/bin/python -c "import scanpy, anndata, numpy, scipy, pandas, leidenalg, igraph"

RUN R -q -e 'pkgs <- c("pak", "cli", "fs", "httpuv"); for (attempt in seq_len(5)) { try(install.packages(pkgs, repos = getOption("repos"), type = "source"), silent = TRUE); missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]; if (!length(missing)) break; if (attempt == 5) stop("Failed to install packages: ", paste(missing, collapse = ", ")); Sys.sleep(15 * attempt) }'

RUN R -q -e 'options(Ncpus = 1, pkg.sysreqs = FALSE); pak::pkg_install(c( \
    "shiny", "bslib", "Seurat", "SeuratObject", "harmony", \
    "ggplot2", "dplyr", "patchwork", "reactable", "DT", \
    "shinyWidgets", "shinyFeedback", "shinycssloaders", "reticulate", "sass", "metap", \
    "bioc::SingleCellExperiment", "bioc::SummarizedExperiment", "bioc::DelayedMatrixStats", "bioc::slingshot", "bioc::hdf5r", \
    "bioc::scDblFinder", "bioc::SingleR", \
    "bioc::celldex", "bioc::tradeSeq", \
    "bioc::anndataR"), ask = FALSE)'

RUN R -q -e 'options(Ncpus = 1, pkg.sysreqs = FALSE); pak::pkg_install("samuel-marsh/scCustomize", ask = FALSE)'
RUN R -q -e 'options(Ncpus = 1, pkg.sysreqs = FALSE); pak::pkg_install("BioBam/scMaSigPro", ask = FALSE)'
ARG BPCELLS_REF=adc4a3c30f60a03522f58947d733d7d77a6eb2cf
RUN R -q -e 'options(Ncpus = 1, pkg.sysreqs = FALSE); pak::pkg_install(c("hexbin", "RcppEigen", "readr", "scattermore", "ggrepel", "RColorBrewer", "tidyr", "tibble", "vctrs", "lifecycle", "stringr", "magrittr", "scales"), ask = FALSE)'
RUN git init /tmp/BPCells \
    && cd /tmp/BPCells \
    && git remote add origin https://github.com/bnprks/BPCells.git \
    && git fetch --depth 1 origin "${BPCELLS_REF}" \
    && git checkout --detach FETCH_HEAD \
    && CI=true R CMD INSTALL r \
    && cd / \
    && rm -rf /tmp/BPCells
ARG MONOCLE3_REF=536f1033d6de7c957f26a1f403f81efbd825e0db
RUN R -q -e 'options(Ncpus = 1, pkg.sysreqs = FALSE); pak::pkg_install(c("bioc::batchelor", "bioc::ResidualMatrix", "classInt", "distributional", "furrr", "ggdist", "ggforce", "leidenbase", "lme4", "minqa", "nloptr", "pbmcapply", "pscl", "quadprog", "reformulas", "rsample", "s2", "sf", "slam", "slider", "spData", "spdep", "cole-trapnell-lab/speedglm", "tweenr", "units", "warp", "wk"), ask = FALSE)'
RUN git init /tmp/monocle3 \
    && cd /tmp/monocle3 \
    && git remote add origin https://github.com/cole-trapnell-lab/monocle3.git \
    && git fetch --depth 1 origin "${MONOCLE3_REF}" \
    && git checkout --detach FETCH_HEAD \
    && CI=true R CMD INSTALL . \
    && cd / \
    && rm -rf /tmp/monocle3
RUN R -q -e 'options(Ncpus = 1, pkg.sysreqs = FALSE); pak::pkg_install("dsong-lab/PseudotimeDE", ask = FALSE, upgrade = FALSE)'

RUN R -q -e 'library(reticulate); use_python(Sys.getenv("RETICULATE_PYTHON"), required = TRUE); py_run_string("import scanpy, anndata"); print(py_config())'
RUN R -q -e 'stopifnot(requireNamespace("PseudotimeDE", quietly = TRUE)); stopifnot(requireNamespace("monocle3", quietly = TRUE)); stopifnot(requireNamespace("tradeSeq", quietly = TRUE))'
RUN R -q -e 'pkgs <- c("shiny", "bslib", "Seurat", "SeuratObject", "scCustomize", "slingshot", "harmony", "SingleCellExperiment", "SummarizedExperiment", "DelayedMatrixStats", "ggplot2", "dplyr", "patchwork", "reactable", "DT", "shinyWidgets", "shinyFeedback", "shinycssloaders", "reticulate", "sass", "hdf5r", "metap", "monocle3", "PseudotimeDE", "scMaSigPro", "tradeSeq", "anndataR", "scDblFinder", "SingleR", "celldex", "BPCells"); missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]; if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))'

COPY . /tmp/ascseurat
RUN R -q -e 'pak::local_install("/tmp/ascseurat", dependencies = FALSE)' \
    && rm -rf /tmp/ascseurat


FROM rocker/r-ver:4.5.3 AS runtime

ARG TARGETARCH

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="/opt/venv/bin:${PATH}"
ENV RETICULATE_PYTHON=/opt/venv/bin/python

# Runtime-only shared libraries needed by the installed R packages.
# Many runtime libs (libcairo2, libcurl4t64, libfontconfig1, libfreetype6,
# libfribidi0, libgfortran5, libharfbuzz0b, libpng16-16t64, libssl3t64,
# libtiff6, libx11-6, libxml2) ship in rocker/r-ver:4.5.3 already and are
# not listed here. Everything below is added on top of that base.
RUN set -eux; \
    if [ "$(dpkg --print-architecture)" = "arm64" ]; then \
        find /etc/apt -type f \( -name '*.list' -o -name '*.sources' \) -exec sed -i 's|http://ports.ubuntu.com/ubuntu-ports|https://mirrors.mit.edu/ubuntu-ports/ubuntu-ports|g; s|https://ports.ubuntu.com/ubuntu-ports|https://mirrors.mit.edu/ubuntu-ports/ubuntu-ports|g' {} +; \
    fi; \
    for attempt in 1 2 3 4 5; do \
        apt-get -o Acquire::Retries=5 -o Acquire::http::Timeout=30 -o Acquire::https::Timeout=30 update \
        && apt-get -o Acquire::Retries=5 -o Acquire::http::Timeout=30 -o Acquire::https::Timeout=30 install -y --no-install-recommends \
            ca-certificates \
            curl \
            libfftw3-double3 \
            libgdal34t64 \
            libgeos-c1t64 \
            libgit2-1.7 \
            libglpk40 \
            libhdf5-103-1t64 \
            libnode109 \
            libproj25 \
            libpython3.12t64 \
            libudunits2-0 \
            libuv1t64 \
            libwebpmux3 \
            pandoc \
            python3 \
        && break; \
        if [ "$attempt" = "5" ]; then exit 1; fi; \
        rm -rf /var/lib/apt/lists/*; \
        sleep $((attempt * 15)); \
    done; \
    rm -rf /var/lib/apt/lists/*

COPY --from=builder /usr/local/lib/R/site-library /usr/local/lib/R/site-library
COPY --from=builder /opt/venv /opt/venv

RUN useradd -m -u 1001 ascseurat \
    && mkdir -p /home/ascseurat/data /home/ascseurat/RDS_files \
    && chown -R ascseurat:ascseurat /home/ascseurat

USER ascseurat
WORKDIR /home/ascseurat

EXPOSE 3838

HEALTHCHECK --interval=30s --timeout=10s --start-period=60s \
    CMD curl -f http://localhost:3838/ || exit 1

CMD ["R", "-q", "-e", "ascseurat::run_app(host = '0.0.0.0', port = 3838, launch.browser = FALSE)"]
