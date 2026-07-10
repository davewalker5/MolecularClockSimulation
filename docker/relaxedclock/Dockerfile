FROM python:3.13-slim-bullseye

ARG APPLICATION_VERSION

COPY relaxedclockexplorer-${APPLICATION_VERSION}.0 /opt/relaxedclockexplorer

WORKDIR /opt/relaxedclockexplorer

RUN apt-get update -y
RUN pip install --no-cache-dir /opt/relaxedclockexplorer/molecular_clock_simulation-${APPLICATION_VERSION}-py3-none-any.whl

ENTRYPOINT [ "relaxed-clock-explorer" ]
CMD ["--server.address=0.0.0.0", "--server.port=8501", "--runner.fastReruns=false"]
