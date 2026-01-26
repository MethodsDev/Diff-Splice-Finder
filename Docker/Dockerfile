FROM python:3.11-slim

# Install system dependencies
RUN apt-get update && apt-get install -y \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
RUN pip install --no-cache-dir pysam

# Clone the repository
WORKDIR /opt
RUN git clone https://github.com/MethodsDev/Diff-Splice-Finder.git

# Set up PATH to include base directory and util/ subdirectory
ENV PATH="/opt/Diff-Splice-Finder:/opt/Diff-Splice-Finder/util:${PATH}"

# Set working directory
WORKDIR /data

CMD ["/bin/bash"]
