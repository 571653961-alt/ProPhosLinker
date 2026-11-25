# =============================================================================
# ProPhosLinker Docker Image
# Multi-omics integration pipeline: protein–phosphosite correlation,
# pattern discovery, differential network & functional enrichment.
#
# Build:
#   docker build -t prophoslinker .
#
# Run (mount data + output, expose Neo4j if needed):
#   docker run -it --rm \
#     -v ./casedata:/data \
#     -v ./results:/results \
#     prophoslinker \
#     --pro_file /data/protein_abundance.tsv \
#     --phos_file /data/phosphoprotein_abundance.tsv \
#     --sample_group /data/compare_groups.tsv \
#     --mapping_file /data/protein_phosphoproSite.tsv \
#     --outdir /results \
#     --group_comparing "T:N"
# =============================================================================

FROM ubuntu:22.04

LABEL maintainer="ppy1222"
LABEL description="ProPhosLinker: Multi-omics proteomics & phosphoproteomics analysis pipeline"
LABEL version="1.0"

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Asia/Shanghai

# =============================================================================
# 1. System dependencies
# =============================================================================
RUN apt-get update && apt-get install -y --no-install-recommends \
    # Build tools & utilities
    build-essential \
    wget \
    curl \
    ca-certificates \
    gnupg \
    git \
    unzip \
    # R system dependencies
    r-base=4.1.2-1ubuntu2 \
    r-base-dev \
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
    libudunits2-dev \
    libv8-dev \
    libnode-dev \
    libgmp-dev \
    libmpfr-dev \
    libcairo2-dev \
    libxt-dev \
    libgit2-dev \
    librsvg2-dev \
    # Python
    python3.10 \
    python3.10-dev \
    python3-pip \
    python3-venv \
    # Neo4j requires Java 17+
    openjdk-17-jre-headless \
    # Misc
    procps \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# =============================================================================
# 2. Python setup (pip + virtualenv)
# =============================================================================
RUN python3 -m pip install --no-cache-dir --upgrade pip setuptools wheel

# Install Python dependencies — layer cached separately from project code
RUN python3 -m pip install --no-cache-dir \
    pandas \
    numpy \
    scipy \
    scikit-learn \
    matplotlib \
    statsmodels \
    networkx \
    python-igraph \
    leidenalg \
    neo4j \
    pyyaml \
    openpyxl

# =============================================================================
# 3. R setup — install CRAN + Bioconductor packages
# =============================================================================
# Use a standalone R script to avoid shell escaping issues in RUN
COPY install_r_packages.R /tmp/install_r_packages.R
RUN Rscript /tmp/install_r_packages.R && rm /tmp/install_r_packages.R

# =============================================================================
# 4. Neo4j Community Edition
# =============================================================================
ENV NEO4J_VERSION=5.26.4
ENV NEO4J_HOME=/opt/neo4j
ENV NEO4J_DATA=/var/lib/neo4j/data
ENV NEO4J_LOGS=/var/log/neo4j

RUN wget -q "https://dist.neo4j.org/neo4j-community-${NEO4J_VERSION}-unix.tar.gz" -O /tmp/neo4j.tar.gz \
    && tar xzf /tmp/neo4j.tar.gz -C /opt \
    && mv /opt/neo4j-community-${NEO4J_VERSION} ${NEO4J_HOME} \
    && rm /tmp/neo4j.tar.gz \
    && mkdir -p ${NEO4J_DATA} ${NEO4J_LOGS}

# Configure Neo4j for headless/container use
RUN echo "dbms.security.auth_enabled=false" >> ${NEO4J_HOME}/conf/neo4j.conf \
    && echo "dbms.connector.bolt.listen_address=0.0.0.0:7687" >> ${NEO4J_HOME}/conf/neo4j.conf \
    && echo "dbms.connector.http.listen_address=0.0.0.0:7474" >> ${NEO4J_HOME}/conf/neo4j.conf \
    && echo "server.memory.heap.initial_size=512m" >> ${NEO4J_HOME}/conf/neo4j.conf \
    && echo "server.memory.heap.max_size=1g" >> ${NEO4J_HOME}/conf/neo4j.conf \
    && echo "server.memory.pagecache.size=512m" >> ${NEO4J_HOME}/conf/neo4j.conf

# Add Neo4j to PATH
ENV PATH="${NEO4J_HOME}/bin:${PATH}"

# =============================================================================
# 5. Copy project code
# =============================================================================
WORKDIR /app
COPY . /app/

# Add project to Python path
ENV PYTHONPATH="/app"

# =============================================================================
# 6. Entrypoint — start Neo4j, then run the pipeline
# =============================================================================
COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]
CMD ["--help"]