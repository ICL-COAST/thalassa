## SPICE builder
# Base on Alpine Linux
FROM alpine:latest AS BUILDER_SPICE

# Update and upgrade packages
RUN apk update && apk upgrade

# Install packages
RUN apk add gfortran curl tar tcsh

# Change work directory
WORKDIR /spice

# Copy SPICE script
COPY external/spice.sh spice.sh

# Download, unpack, and compile SPICE
RUN sh spice.sh


## SOFA builder
# Base on Alpine Linux
FROM alpine:latest AS BUILDER_SOFA

# Update and upgrade packages
RUN apk update && apk upgrade

# Install packages
RUN apk add gfortran curl tar make

# Change work directory
WORKDIR /sofa

# Copy SOFA script
COPY external/sofa.sh sofa.sh

# Download, unpack, and compile SOFA
RUN sh sofa.sh


## THALASSA Builder
# Base on Alpine Linux
FROM alpine:latest AS BUILDER_THALASSA

# Update and upgrade packages
RUN apk update && apk upgrade

# Install packages
RUN apk add gfortran libc-dev make cmake

# Change workdirectory
WORKDIR /thalassa

# Copy THALASSA source
COPY . .

# Copy SPICE from BUILDER_SPICE
COPY --from=BUILDER_SPICE /spice/spicelib.a /thalassa/external/spicelib.a

# Copy SOFA from BUILDER_SOFA
COPY --from=BUILDER_SOFA /sofa/libsofa.a /thalassa/external/libsofa.a

# Configure CMake files for THALASSA
RUN cmake -B/thalassa/build

# Build THALASSA with CMake
RUN cmake --build /thalassa/build --config Release


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

# Copy THALASSA files from BUILDER_THALASSA
COPY --from=BUILDER_THALASSA /thalassa/thalassa_main /thalassa/thalassa_main
COPY --from=BUILDER_THALASSA /thalassa/data /thalassa/data
COPY --from=BUILDER_THALASSA /thalassa/in /thalassa/in

# Set folder permissions
RUN chmod -R 777 /thalassa

# Create THALASSA user
RUN adduser -D thalassa

# Switch to THALASSA user
USER thalassa

# Run THALASSA
ENTRYPOINT ./thalassa_main