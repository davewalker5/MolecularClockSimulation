FROM python:3.13-slim-bullseye

ARG APPLICATION_VERSION

COPY strictclockexplorer-${APPLICATION_VERSION}.0 /opt/strictclockexplorer

WORKDIR /opt/strictclockexplorer

RUN apt-get update -y \
    && apt-get install -y --no-install-recommends graphviz \
    && rm -rf /var/lib/apt/lists/*
RUN pip install --no-cache-dir /opt/strictclockexplorer/molecular_clock_simulation-${APPLICATION_VERSION}-py3-none-any.whl

ENTRYPOINT [ "strict-clock-explorer" ]
CMD ["--server.address=0.0.0.0", "--server.port=8501", "--runner.fastReruns=false"]
