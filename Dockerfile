## SPICE builder
# Base on Alpine Linux
FROM alpine:latest AS BUILDER_SPICE

# Update and upgrade packages
RUN apk update && apk upgrade

# Change workdirectory
WORKDIR /spice

# Install packages
RUN apk add gfortran curl tar tcsh

# Download and unpack SPICE
RUN curl -O https://naif.jpl.nasa.gov/pub/naif/toolkit//FORTRAN/PC_Linux_gfortran_64bit/packages/toolkit.tar.Z
RUN tar -zxvf toolkit.tar.Z

# Compile SPICE
WORKDIR /spice/toolkit
RUN tcsh -f makeall.csh

## SOFA builder
# Base on Alpine Linux
FROM alpine:latest AS BUILDER_SOFA

# Update and upgrade packages
RUN apk update && apk upgrade

# Change workdirectory
WORKDIR /sofa

# Install packages
RUN apk add gfortran curl unzip make

# Download and unpack SOFA
# WARNING: downloaded using http due to SSL certificate issues
RUN curl -O http://www.iausofa.org/2021_0512_F/sofa_f.zip
RUN unzip sofa_f.zip

# Compile SOFA
WORKDIR /sofa/sofa/20210512/f77/src
RUN make


## THALASSA Builder
# Base on Alpine Linux
FROM alpine:latest AS BUILDER_THALASSA

# Update and upgrade packages
RUN apk update && apk upgrade

# Change workdirectory
WORKDIR /thalassa

# Install packages
RUN apk add gfortran libc-dev make

# Copy THALASSA source
COPY . .

# Copy SPICE from BUILDER_SPICE
# TODO: handle SPICE kernels
COPY --from=BUILDER_SPICE /spice/toolkit/lib/spicelib.a /thalassa/lib/spicelib.a

# Copy SOFA from BUILDER_SOFA
COPY --from=BUILDER_SOFA /sofa/sofa/20210512/f77/src/libsofa.a /thalassa/lib/libsofa.a

# Build THALASSA
RUN make


## Release
# Base on Alpine Linux
FROM alpine:latest AS RELEASE

# TODO: addMetadata

# Update and upgrade packages
RUN apk update && apk upgrade

# Change workdirectory
WORKDIR /thalassa

# Install packages
RUN apk add libgfortran libgcc libquadmath

# Change workdirectory
WORKDIR /thalassa

# Copy THALASSA files from BUILDER_THALASSA
COPY --from=BUILDER_THALASSA /thalassa/thalassa.x /thalassa/thalassa.x
COPY --from=BUILDER_THALASSA /thalassa/data /thalassa/data
COPY --from=BUILDER_THALASSA /thalassa/in /thalassa/in

# Set folder permissions
RUN chmod -R 777 /thalassa

# Create THALASSA user
RUN adduser -D thalassa

# Switch to THALASSA user
USER thalassa

# Run THALASSA
ENTRYPOINT ./thalassa.x