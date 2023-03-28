#! /bin/sh

if [ $# -ne 1 ]; then
  echo "Usage: pack.sh TARGET_DIR"
  exit 0
fi

SHA=HEAD
version=$(git describe --tag $SHA | sed -e s,TEMPO_SDPC_v,,)
if test -z "$version" ; then
   echo "HEAD is untagged -- using sha1 label"
   version=$(git rev-parse HEAD | cut -c 1-8)
fi

target="sdpc-$version"

target_dir=$1

if ! test -d $target_dir ; then
   mkdir -p $target_dir
fi

if test -d ${target_dir}/$target ; then
   echo "WARNING: ${target_dir}/$target exists"
   exit 1
fi

git archive --format=tar --prefix="$target/" $SHA | (cd ${target_dir} && tar xf - )
git --no-pager log -1 --pretty=format:%H%n $SHA > ${target_dir}/$target/git-commit-hash
./utils/get_sdpc_version.sh | sed -e s,TEMPO_SDPC_v,, > ${target_dir}/$target/sdpc-release

cd ${target_dir} || exit 1

tar czf $target.tar.gz --remove-files $target

echo "Created ${target_dir}/$target.tar.gz"
