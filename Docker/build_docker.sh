#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
version="$(tr -d '[:space:]' < "${repo_root}/Docker/VERSION.txt")"
image_base="${1:-us-central1-docker.pkg.dev/methods-dev-lab/diff-splice-finder/diff-splice-finder}"
image="${image_base}:${version}"
latest="${image_base}:latest"
dsf_ref="${DSF_REF:-master}"

docker build --pull --no-cache --build-arg DSF_REF="${dsf_ref}" -f "${repo_root}/Docker/Dockerfile" -t "${image}" -t "${latest}" "${repo_root}/Docker"
echo "Built ${image} and ${latest} from ${dsf_ref}"
