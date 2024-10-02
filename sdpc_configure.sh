#! /bin/sh

compilers=gnu

build_arch="${compilers}-rh8"

soft_home="/tempo/nas0/sdpc_soft"
archive_home="/tempo/nas0/sdpc_archive"
pipe_home="/tempo/nas0/sdpc/liveroot"

refdata_dir="$soft_home/refdata"
ancdata_dir="$archive_home/ancillary"

install_root="$soft_home/install/$build_arch"
inrroot="$install_root/inr_r3.0.2"
otsroot="$install_root/ots"
s6root="$install_root/skarnet"

prefix="$install_root/sdpc/v4p4p3"

# --with-iocsdpc=/tempo/nas0/sdpc_soft/src/ots_sdpc.git/iocsdpc_test \

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
