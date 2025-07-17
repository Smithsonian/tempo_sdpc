#! /bin/sh

compilers=gnu

site=HEAD
#site=ITS

case "$site" in
  ITS )
     OS_LABEL=rocky8
     DEVROOT=/proj
     localroot=/h
     ;;

  HEAD )
     OS_LABEL=rh8
     DEVROOT=/tempo/nas0
     localroot=/scratch
     ;;

  * )
     echo "Unsupported site: $site"
     exit 1
     ;;
esac

build_arch="${compilers}-${OS_LABEL}"

soft_home="$DEVROOT/sdpc_soft"
pipe_home="$DEVROOT/sdpc/liveroot"
archive_home="$DEVROOT/sdpc_archive"

refdata_dir="$soft_home/refdata"
ancdata_dir="$archive_home/ancillary"

install_root="$soft_home/install/$build_arch"
inrroot="$install_root/inr_r3.0.4"
otsroot="$install_root/ots"
s6root="$install_root/skarnet"

pythonbindir="$install_root/python-3.12.8/bin"

prefix="$install_root/sdpc/jch_devel"

# --with-iocsdpc=$DEVROOT/sdpc_soft/src/ots_sdpc.git/iocsdpc_test \

./configure --prefix=$prefix \
            --with-compilers=$compilers \
            --with-otsroot=$otsroot \
            --with-slang=$otsroot \
            --with-inrroot=$inrroot \
            --with-s6root=$s6root \
            --with-localrootdir=$localroot \
            --with-pythonbindir=$pythonbindir \
            --with-refdata=$refdata_dir \
            --with-ancdata=$ancdata_dir \
            --with-atlasblas \
            --with-archive-home=$archive_home \
            --with-pipe-home=$pipe_home
