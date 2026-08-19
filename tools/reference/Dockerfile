# syntax=docker/dockerfile:1
FROM debian:bookworm-slim@sha256:abd67ffcfa541b485a3dff59865ab629aa048a6c613e639d36e7456b0b229241

ARG AT_ARCHIVE_SHA256=f8a7a2c1e80a73431cd230a10bef5fcfc996c88889a0e1540771c3922ee2a21f

RUN apt-get update \
    && apt-get install --yes --no-install-recommends ca-certificates gfortran make \
    && rm -rf /var/lib/apt/lists/* \
    && test "$(gfortran -dumpfullversion)" = "12.2.0" \
    && test "$(dpkg-query -W -f='${Version}' gfortran-12)" = "12.2.0-14+deb12u1" \
    && test "$(dpkg-query -W -f='${Version}' libgfortran5)" = "12.2.0-14+deb12u1" \
    && test "$(dpkg-query -W -f='${Version}' libc6)" = "2.36-9+deb12u14"

COPY Acoustics-Toolbox-475108519289c6fb488b58980c644ea14eccc604.tar.gz /tmp/acoustics-toolbox.tar.gz
RUN echo "${AT_ARCHIVE_SHA256}  /tmp/acoustics-toolbox.tar.gz" | sha256sum --check \
    && mkdir /opt/acoustics-toolbox \
    && tar --extract --gzip --strip-components=1 \
        --file /tmp/acoustics-toolbox.tar.gz --directory /opt/acoustics-toolbox \
    && rm /tmp/acoustics-toolbox.tar.gz

WORKDIR /opt/acoustics-toolbox
RUN make -C misc FC=gfortran \
        FFLAGS='-O1 -ffast-math -funroll-all-loops -fomit-frame-pointer -std=gnu' \
    && make -C Bellhop bellhop.exe FC=gfortran \
        FFLAGS='-O1 -ffast-math -funroll-all-loops -fomit-frame-pointer -std=gnu -I../misc' \
    && install -m 0755 Bellhop/bellhop.exe /usr/local/bin/bellhop-fortran \
    && printf '%s\n' \
        'Acoustics Toolbox v2023.5' \
        'commit 475108519289c6fb488b58980c644ea14eccc604' \
        'GNU Fortran 12.2.0 (Debian 12.2.0-14+deb12u1)' \
        'GNU Fortran runtime 12.2.0-14+deb12u1' \
        'glibc 2.36-9+deb12u14' \
        'flags -O1 -ffast-math -funroll-all-loops -fomit-frame-pointer -std=gnu' \
        > /opt/reference-build.txt

LABEL org.opencontainers.image.title="bellhop-rs numerical reference" \
      org.opencontainers.image.source="https://github.com/leihuang1/bellhop-rs" \
      org.opencontainers.image.version="Acoustics Toolbox v2023.5"

WORKDIR /work
ENTRYPOINT ["/usr/local/bin/bellhop-fortran"]
