## THALASSA Builder
# Base on Alpine Linux
FROM alpine:latest AS BUILDER_THALASSA

# Update and upgrade packages
RUN apk update && apk upgrade

# Install packages
RUN apk add gcc g++ gfortran libc-dev make cmake

# Change workdirectory
WORKDIR /thalassa

# Copy THALASSA source
COPY . .

# Configure CMake files for THALASSA
RUN cmake -B /thalassa/build

# Build THALASSA with CMake
RUN cmake --build /thalassa/build --config Release --target thalassa_main


## Release
# Base on Alpine Linux
FROM alpine:latest AS RELEASE

# TODO: addMetadata

# Update and upgrade packages
RUN apk update && apk upgrade

# Install packages
RUN apk add libgfortran libgcc libquadmath

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
COPY --from=BUILDER_THALASSA /thalassa/thalassa_main /thalassa/thalassa_main
COPY --from=BUILDER_THALASSA /thalassa/data /thalassa/data
COPY --from=BUILDER_THALASSA /thalassa/in /thalassa/in

# Change ownership of the THALASSA directory
RUN chown -R thalassa:thalassagroup /thalassa

# Switch to THALASSA user
USER thalassa

# Create output directory
RUN mkdir /thalassa/out

# Run THALASSA
ENTRYPOINT ./thalassa_main