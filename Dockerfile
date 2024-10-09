FROM python:3.10 AS base

# Update pip to latest version
RUN python -m pip install --upgrade pip

# Install dependencies
RUN apt-get update \
    && apt-get upgrade -y \
    && apt-get install -y curl \
    && apt-get install -y rsync

# Needed to install requirements,
# in devcontainer a local mounted version of flowcell_parser is used
RUN git clone https://github.com/SciLifeLab/flowcell_parser.git

RUN cd flowcell_parser \
    && python -m pip install -r requirements.txt \
    && pip3 install -e .

# Install dependencies
COPY requirements.txt requirements.txt
COPY requirements-dev.txt requirements-dev.txt
RUN python -m pip install -r requirements-dev.txt

RUN mkdir /root/.taca/

FROM base AS testing
COPY . /taca
RUN python -m pip install -e /taca
WORKDIR /taca/tests
CMD ["python", "-m", "unittest"]