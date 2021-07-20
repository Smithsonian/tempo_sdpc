#! /bin/bash

out_decl="ndinterp_decl.inc"
out_code="ndinterp_code.inc"

expand_type(){
  typ=$1
  case $typ in
    r4 ) fkind=r4
         typchar=f
         ;;
    r8 ) fkind=r8
         typchar=d
         ;;
  esac
}

expand_index_string(){
  aname=$1
  dim=$2
  s=""
  for d in $(seq $dim) ; do
     s+="$aname($d)"
     if test $d -lt $dim ; then
        s+=","
     fi
  done
  index_string=$s
}

expand_colons(){
  num=$1
  colons=":"
  for d in $(seq $(($num-1)) ) ; do
     colons+=",:"
  done
}

expand_code_template(){
  template=$1
  typ=$2
  dim=$3

  expand_type $typ
  expand_colons $dim
  expand_index_string "dimlens" $dim
  idxdimlens="$index_string"

  expand_index_string "id_k" $dim
  idxid="$index_string"

  sed -e s/@type@/$typ/g \
      -e s/@fkind@/$fkind/g \
      -e s/@typchar@/$typchar/g \
      -e s/@dim@/$dim/g \
      -e s/@colons@/$colons/g \
      -e s/@idxdimlens@/$idxdimlens/g \
      -e s/@idxid@/$idxid/g \
      $template >> $out_code
}

expand_modproc_decl(){
  type_list="$1"
  dim_list="$2"
  indent="    "
  printf "%s\n" "${indent}interface ndi_table_interp"
  for typ in $type_list; do
     expand_type $typ
     for dim in $dim_list ; do
        printf "%s\n" "${indent}   module procedure ndi_table_${typchar}${dim}d_interp"
     done
  done
  printf "%s\n" "${indent}end interface $prefix"
}

dim_list="1 2 3 4 5 6 7"
type_list="r4 r8"

echo "!> Auto-generated file -- do not edit." > $out_decl
echo "!> Auto-generated file -- do not edit." > $out_code

for dim in $dim_list ; do
  for typ in $type_list; do
    expand_code_template "ndinterp_module_code.in" $typ $dim
  done
done

(
expand_modproc_decl "$type_list" "$dim_list"
) >>  $out_decl
