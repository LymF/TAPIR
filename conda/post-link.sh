#!/bin/bash
# Install pip-only dependencies not available in bioconda/conda-forge
"${PREFIX}/bin/pip" install --no-deps cobra-meta viralquest 2>&1 || true
