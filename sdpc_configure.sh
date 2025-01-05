#! /bin/sh

compilers=gnu

# -- HEAD network
OS_LABEL=rh8
DEVROOT=/tempo/nas0

# -- ITS network
# OS_LABEL=rocky8
# DEVROOT=/proj

build_arch="${compilers}-${OS_LABEL}"

soft_home="$DEVROOT/sdpc_soft"
archive_home="$DEVROOT/sdpc_archive"
pipe_home="$DEVROOT/sdpc/liveroot"

refdata_dir="$soft_home/refdata"
ancdata_dir="$archive_home/ancillary"

install_root="$soft_home/install/$build_arch"
inrroot="$install_root/inr_r3.0.2"
otsroot="$install_root/ots"
s6root="$install_root/skarnet"

prefix="$install_root/sdpc/devel"

# --with-iocsdpc=$DEVROOT/sdpc_soft/src/ots_sdpc.git/iocsdpc_test \

./configure --prefix=$prefix \
            --with-compilers=$compilers \
            --with-otsroot=$otsroot \
            --with-slang=$otsroot \
            --with-inrroot=$inrroot \
            --with-s6root=$s6root \
            --with-refdata=$refdata_dir \
            --with-ancdata=$ancdata_dir \
            --with-atlasblas \
            --with-archive-home=$archive_home \
            --with-pipe-home=$pipe_home
