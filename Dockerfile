FROM nfcore/base
LABEL authors="Bilge Sürün" \
      description="Docker image containing all requirements for nf-core/clinvap pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-clinvap-1.0dev/bin:$PATH
