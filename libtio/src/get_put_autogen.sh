#! /bin/sh

template="get_put_tmpl.in"
out_code="get_put_code.inc"
out_decl="get_put_decl.inc"

expand_template(){
  act=$1
  typ=$2
  dim=$3
  colons=$4
  case $act in
    get )
      ioaction=read
      iodir=from
      arraydir=out
      ;;
    put )
      arraydir=in
      ioaction=write
      iodir=to
      ;;
  esac
  case $typ in
    i1 ) ntype=nf90_byte
         ftype=integer
         fkind=i1
         ;;
    i2 ) ntype=nf90_short
         ftype=integer
         fkind=i2
         ;;
    i4 ) ntype=nf90_int
         ftype=integer
         fkind=i4
         ;;
    r4 ) ntype=nf90_float
         ftype=real
         fkind=4
         ;;
    r8 ) ntype=nf90_double
         ftype=real
         fkind=8
         ;;
  esac
  sed -e s/@type@/$typ/g \
      -e s/@ftype@/$ftype/g \
      -e s/@fkind@/$fkind/g \
      -e s/@ntype@/$ntype/g \
      -e s/@action@/$act/g \
      -e s/@dim@/$dim/g \
      -e s/@colons@/$colons/g \
      -e s/@arraydir@/$arraydir/g \
      -e s/@ioaction@/$ioaction/g \
      -e s/@iodir@/$iodir/g \
      $template >> $out_code

  echo "public tiof_${act}${dim}d_${typ}" >> $out_decl
}

echo '! Auto-generated file -- do not edit.' > $out_decl
echo '! Auto-generated file -- do not edit.' > $out_code

dim_list="1 2 3"
type_list="i1 i2 i4 r4 r8"

colons=":"
for dim in $dim_list ; do
  for typ in $type_list; do
    expand_template "get" $typ $dim $colons
    expand_template "put" $typ $dim $colons
  done
  colons="${colons},:"
done

