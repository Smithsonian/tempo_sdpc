#! /bin/sh

template_scalar="get_put_scalar.in"
template_array="get_put_array.in"
out_code="get_put_code.inc"
out_decl="get_put_decl.inc"

expand_action(){
  act=$1
  case $act in
    get )
      ioaction=read
      iodir=from
      intent=out
      ;;
    put )
      intent=in
      ioaction=write
      iodir=to
      ;;
  esac
}

expand_type(){
  typ=$1
  case $typ in
    i1 ) ntype=nf90_byte
         ftype=integer
         fkind=i1
         ;;
   ui1 ) ntype=nf90_ubyte
         ftype=integer
         fkind=i1
         ;;
    i2 ) ntype=nf90_short
         ftype=integer
         fkind=i2
         ;;
   ui2 ) ntype=nf90_ushort
         ftype=integer
         fkind=i2
         ;;
    i4 ) ntype=nf90_int
         ftype=integer
         fkind=i4
         ;;
   ui4 ) ntype=nf90_uint
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
}

expand_template_scalar(){
  act=$1
  typ=$2
  expand_action $act
  expand_type $typ
  sed -e s/@type@/$typ/g \
      -e s/@ftype@/$ftype/g \
      -e s/@fkind@/$fkind/g \
      -e s/@ntype@/$ntype/g \
      -e s/@action@/$act/g \
      -e s/@intent@/$intent/g \
      -e s/@ioaction@/$ioaction/g \
      -e s/@iodir@/$iodir/g \
      $template_scalar >> $out_code

  echo "public tiof_${act}_${typ}" >> $out_decl
}

expand_template_array(){
  act=$1
  typ=$2
  dim=$3
  colons=$4
  expand_action $act
  expand_type $typ
  sed -e s/@type@/$typ/g \
      -e s/@ftype@/$ftype/g \
      -e s/@fkind@/$fkind/g \
      -e s/@ntype@/$ntype/g \
      -e s/@action@/$act/g \
      -e s/@dim@/$dim/g \
      -e s/@colons@/$colons/g \
      -e s/@intent@/$intent/g \
      -e s/@ioaction@/$ioaction/g \
      -e s/@iodir@/$iodir/g \
      $template_array >> $out_code

  echo "public tiof_${act}${dim}d_${typ}" >> $out_decl
}

echo '! Auto-generated file -- do not edit.' > $out_decl
echo '! Auto-generated file -- do not edit.' > $out_code

dim_list="1 2 3"
type_list="i1 ui1 i2 ui2 i4 ui4 r4 r8"

for typ in $type_list; do
  expand_template_scalar "get" $typ
  expand_template_scalar "put" $typ
done

colons=":"
for dim in $dim_list ; do
  for typ in $type_list; do
    expand_template_array "get" $typ $dim $colons
    expand_template_array "put" $typ $dim $colons
  done
  colons="${colons},:"
done

