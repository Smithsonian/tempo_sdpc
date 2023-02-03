#! /bin/sh

compilers=gnu

build_arch="${compilers}"
#build_arch="${compilers}-rh8"

soft_home="/tempo/nas0/sdpc_soft"
archive_home="/tempo/nas0/sdpc_archive"
pipe_home="/tempo/nas0/sdpc/liveroot"

refdata_dir="$soft_home/refdata"
ancdata_dir="$archive_home/ancillary"

install_root="$soft_home/install/$build_arch"
prefix="$install_root/sdpc/v4"
inrroot="$install_root/inr_r2.3.6"
otsroot="$install_root/ots"
s6root="$install_root/skarnet"

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
