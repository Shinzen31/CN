# Use a base image with the necessary dependencies
FROM ubuntu:latest

# Install build dependencies
RUN apt-get update -y && \
    apt-get install -y build-essential liblapacke-dev

# Set the working directory
WORKDIR /src

# Copy the source code and header files into the container
COPY src/tp_env.c .
COPY src/tp_poisson1D_direct.c .
COPY src/lib_poisson1D.c .
COPY include ./include

# Copy the Makefile into the container
COPY Makefile .

# Compile the program
RUN make all

# Set the default command to run the program
CMD ["./bin/tp_poisson1d_direct"]
