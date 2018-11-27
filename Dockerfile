# See ../triqs/packaging for other options
FROM flatironinstitute/triqs:master-ubuntu-clang

RUN apt-get install -y libgsl-dev || yum install -y gsl-devel

ARG APPNAME
COPY . $SRC/$APPNAME
WORKDIR $BUILD/$APPNAME
RUN chown build .
USER build
RUN cmake $SRC/$APPNAME -DTRIQS_ROOT=${INSTALL} && make -j2 && CTEST_OUTPUT_ON_FAILURE=1 make test
USER root
RUN make install
