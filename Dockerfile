# Build arguments
ARG BASE_DISTRO="debian"
ARG CMAKE_CONFIG="Release"

## THALASSA runtime images
# Alpine
FROM alpine:latest as thalassa_runtime_alpine
# Update, upgrade, and install packages
RUN apk update && \
    apk upgrade && \
    apk add libgfortran libgcc libquadmath

# Debian
FROM debian:bullseye-slim AS thalassa_runtime_debian
# Update, upgrade, and install packages
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y libgfortran5 libgcc-s1 libquadmath0


## THALASSA build images
# Alpine
FROM thalassa_runtime_alpine AS thalassa_builder_alpine
# Update, upgrade, and install packages
RUN apk update && \
    apk upgrade && \
    apk add gcc g++ gfortran libc-dev make cmake git

# Debian
FROM thalassa_runtime_debian AS thalassa_builder_debian
# Update, upgrade, and install packages
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y gcc g++ gfortran make cmake git

# Generic
FROM thalassa_builder_${BASE_DISTRO} as thalassa_builder
# Import arguments
ARG CMAKE_CONFIG
# Change workdirectory
WORKDIR /thalassa
# Copy THALASSA source
COPY ./app/ /thalassa/app/
COPY ./external/ /thalassa/external/
COPY ./interface/ /thalassa/interface/
COPY ./src/ /thalassa/src/
COPY ./CMakeLists.txt /thalassa/CMakeLists.txt
# Configure CMake files for THALASSA
RUN cmake -B /thalassa/build
# Build THALASSA with CMake
RUN cmake --build /thalassa/build --config ${CMAKE_CONFIG} --target thalassa_main


## Release image
# TODO: addMetadata
FROM thalassa_runtime_${BASE_DISTRO} AS thalassa
# Change workdirectory
WORKDIR /thalassa
# Copy THALASSA files from BUILDER_THALASSA
COPY --from=thalassa_builder /thalassa/thalassa_main /thalassa/thalassa_main
COPY ./data /thalassa/data
# Create input directory
RUN mkdir /thalassa/in
# Create output directory
RUN mkdir /thalassa/out
# Run THALASSA
ENTRYPOINT ./thalassa_main