FROM rocker/r-ver:4.6.1

RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    libpng-dev libjpeg-dev cmake && \
    apt-get clean

RUN R -e "install.packages('BiocManager', repos='https://cran.r-project.org')" && \
    R -e "BiocManager::install(c('TCGAbiolinks','maftools','org.Hs.eg.db'), ask=FALSE, update=FALSE)" && \
    R -e "install.packages(c('survival','survminer','ggplot2','data.table','shiny','plotly','DT'), repos='https://cran.r-project.org')"

WORKDIR /pipeline
COPY . /pipeline

CMD ["Rscript", "run_analysis.R"]
