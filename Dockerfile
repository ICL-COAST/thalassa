# Build arguments
ARG BASE_DISTRO="debian"
ARG CMAKE_CONFIG="Release"

## THALASSA Builder images
# Alpine
FROM alpine:latest AS thalassa_builder_alpine
# Import arguments
ARG CMAKE_CONFIG
# Update, upgrade, and install packages
RUN apk update && \
    apk upgrade && \
    apk add gcc g++ gfortran libc-dev make cmake
# Change workdirectory
WORKDIR /thalassa
# Copy THALASSA source
COPY . .
# Configure CMake files for THALASSA
RUN cmake -B /thalassa/build
# Build THALASSA with CMake
RUN cmake --build /thalassa/build --config ${CMAKE_CONFIG} --target thalassa_main

# Debian
FROM debian:bullseye-slim AS thalassa_builder_debian
# Import arguments
ARG CMAKE_CONFIG
# Update, upgrade, and install packages
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y gcc g++ gfortran make cmake 
# Change workdirectory
WORKDIR /thalassa
# Copy THALASSA source
COPY . .
# Configure CMake files for THALASSA
RUN cmake -B /thalassa/build
# Build THALASSA with CMake
RUN cmake --build /thalassa/build --config ${CMAKE_CONFIG} --target thalassa_main


## Release images
# Alpine
FROM alpine:latest AS thalassa_release_alpine
# Update, upgrade, and install packages
RUN apk update && \
    apk upgrade && \
    apk add libgfortran libgcc libquadmath
# Change workdirectory
WORKDIR /thalassa
# Copy THALASSA files from BUILDER_THALASSA
COPY --from=thalassa_builder_alpine /thalassa/thalassa_main /thalassa/thalassa_main
COPY --from=thalassa_builder_alpine /thalassa/data /thalassa/data
COPY --from=thalassa_builder_alpine /thalassa/in /thalassa/in
# Create output directory
RUN mkdir /thalassa/out

# Debian
FROM debian:bullseye-slim AS thalassa_release_debian
# Update, upgrade, and install packages
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y libgfortran5 libgcc-s1 libquadmath0
# Change workdirectory
WORKDIR /thalassa
# Copy THALASSA files from BUILDER_THALASSA
COPY --from=thalassa_builder_debian /thalassa/thalassa_main /thalassa/thalassa_main
COPY --from=thalassa_builder_debian /thalassa/data /thalassa/data
COPY --from=thalassa_builder_debian /thalassa/in /thalassa/in
# Create output directory
RUN mkdir /thalassa/out


## Final image
# TODO: addMetadata
FROM thalassa_release_${BASE_DISTRO} AS thalassa
# Run THALASSA
ENTRYPOINT ./thalassa_main