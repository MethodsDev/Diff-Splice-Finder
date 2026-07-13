#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
version="$(tr -d '[:space:]' < "${repo_root}/Docker/VERSION.txt")"
image_base="${1:-us-central1-docker.pkg.dev/methods-dev-lab/diff-splice-finder/diff-splice-finder}"
image="${image_base}:${version}"
dsf_ref="${DSF_REF:-master}"

docker build --pull --no-cache --build-arg DSF_REF="${dsf_ref}" -f "${repo_root}/Docker/Dockerfile" -t "${image}" "${repo_root}/Docker"
echo "Built ${image} from ${dsf_ref}"
