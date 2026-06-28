FROM python:3.11-slim

ARG DSF_REPO_URL=https://github.com/MethodsDev/Diff-Splice-Finder.git
ARG DSF_COMMIT=c88c80ccad01aede97974027489527042fcbf7f8

# Install system dependencies
RUN apt-get update && apt-get install -y \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install the repository at a pinned commit.
WORKDIR /opt
RUN git clone "${DSF_REPO_URL}" Diff-Splice-Finder \
    && cd Diff-Splice-Finder \
    && git checkout "${DSF_COMMIT}" \
    && pip install --no-cache-dir -r requirements.txt

# Set up PATH to include base directory and util/ subdirectory
ENV PATH="/opt/Diff-Splice-Finder:/opt/Diff-Splice-Finder/util:${PATH}"

# Set working directory
WORKDIR /data

CMD ["/bin/bash"]
