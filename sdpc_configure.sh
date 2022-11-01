#! /bin/sh

compilers=gnu

install_root="/tempo/nas0/sdpc_soft/install/$compilers"
refdata_dir="/tempo/nas0/sdpc_soft/refdata"
ancdata_dir="/tempo/nas0/sdpc_archive/ancillary"

archive_home="/tempo/nas0/sdpc_archive"
pipe_home="/tempo/nas0/sdpc/liveroot"

prefix="$install_root/sdpc/v3"
inrroot="$install_root/inr_r2.3.5"
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
