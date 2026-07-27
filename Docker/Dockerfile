FROM python:3.11-slim

ARG DSF_REPO_URL=https://github.com/MethodsDev/Diff-Splice-Finder.git
ARG DSF_REF=master

# Install system dependencies. samtools is required by util/depth_windows.py,
# which shells out to `samtools depth` for read-depth denominators.
RUN apt-get update && apt-get install -y \
    git \
    samtools \
    && rm -rf /var/lib/apt/lists/*

# Install the repository from the requested ref. The default tracks the latest
# code on the main development branch used by this repository.
WORKDIR /opt
RUN git clone --depth 1 --branch "${DSF_REF}" "${DSF_REPO_URL}" Diff-Splice-Finder \
    && cd Diff-Splice-Finder \
    && pip install --no-cache-dir -r requirements.txt

# Set up PATH to include base directory and util/ subdirectory
ENV PATH="/opt/Diff-Splice-Finder:/opt/Diff-Splice-Finder/util:${PATH}"

# Set working directory
WORKDIR /data

CMD ["/bin/bash"]
