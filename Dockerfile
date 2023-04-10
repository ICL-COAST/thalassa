# Argument for base image
ARG BASEIMAGE="alpine"

## THALASSA Builder images
# Alpine
FROM alpine:latest AS thalassa_builder_alpine

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
RUN cmake --build /thalassa/build --config Release --target thalassa_main


## Release images
# Alpine
FROM alpine:latest AS thalassa_release_alpine

# Update, upgrade, and install packages
RUN apk update && \
    apk upgrade && \
    apk add libgfortran libgcc libquadmath

# Change workdirectory
WORKDIR /thalassa

# Arguments for user and group IDs
ARG UID=1000
ARG GID=1000

# Create THALASSA group
RUN addgroup -g "${GID}" thalassagroup

# Create THALASSA user
RUN adduser -u "${UID}" -G thalassagroup -D thalassa

# Copy THALASSA files from BUILDER_THALASSA
COPY --from=thalassa_builder_alpine /thalassa/thalassa_main /thalassa/thalassa_main
COPY --from=thalassa_builder_alpine /thalassa/data /thalassa/data
COPY --from=thalassa_builder_alpine /thalassa/in /thalassa/in

# Change ownership of the THALASSA directory
RUN chown -R thalassa:thalassagroup /thalassa

# Switch to THALASSA user
USER thalassa

# Create output directory
RUN mkdir /thalassa/out


## Final image
# TODO: addMetadata
FROM thalassa_release_${BASEIMAGE} AS thalassa

# Run THALASSA
ENTRYPOINT ./thalassa_main