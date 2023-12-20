FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y && \
    apt-get install -y gcc libblas-dev liblapacke-dev make vim gnuplot

WORKDIR /app
COPY . .

ENV HOSTNAME=docker
RUN mkdir bin && \
    make all && \
    make run
