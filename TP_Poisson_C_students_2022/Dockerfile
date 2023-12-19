# Use a base image with the necessary dependencies
FROM ubuntu:latest

# Install build dependencies
RUN apt-get update && \
    apt-get install -y build-essential liblapacke-dev

# Set the working directory
WORKDIR /src

# Copy the source code and header files into the container
#COPY src/tp_env.c .
COPY src/tp_poisson1D_direct.c .
COPY src/lib_poisson1D.c .
COPY include ./include

# Compile the program
#RUN gcc -o tp_env tp_env.c -I./include -llapacke -llapack -lblas -lm
RUN gcc -o tp_poisson1d_direct tp_poisson1D_direct.c lib_poisson1D.c -I./include -llapacke -llapack -lblas -lm

# Set the default command to run the program
#CMD ["./tp_env"]
CMD ["./tp_poisson1d_direct"]