FROM antsx/ants:latest

ENV DEBIAN_FRONTEND=noninteractive \
    PIP_NO_CACHE_DIR=1 \
    PYTHONUNBUFFERED=1 \
    VIRTUAL_ENV=/opt/nifreeze-venv

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      python3 \
      python3-pip \
      python3-venv \
      ca-certificates && \
    rm -rf /var/lib/apt/lists/*

COPY requirements.txt /tmp/nifreeze-benchmark-requirements.txt

RUN python3 -m venv "$VIRTUAL_ENV" && \
    python3 -m pip install --upgrade pip setuptools wheel --break-system-packages && \
    "$VIRTUAL_ENV/bin/pip" install --upgrade pip setuptools wheel && \
    "$VIRTUAL_ENV/bin/pip" install -r /tmp/nifreeze-benchmark-requirements.txt && \
    rm -f /tmp/nifreeze-benchmark-requirements.txt

ENV PATH="/opt/nifreeze-venv/bin:/opt/ants/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"

WORKDIR /work

ENTRYPOINT ["python", "inst/benchmarks/run_nifreeze_realdata.py"]
