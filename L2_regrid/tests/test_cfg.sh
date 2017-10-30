#! /bin/bash

FILENAME_LIST=".test_input_files.lst"

test_syntax(){
outfile=$1
cat >$outfile <<EOF
target_grid: {
   longitude = ;
   latitude  = {delta = 0.05; min =   17.0; num =  940;};
};
data_products: {
   "nonexistent_file.nc"
};
EOF
}

test_missing_target_grid(){
outfile=$1
cat >$outfile <<EOF
data_products: {
};
EOF
}

test_defining_target_grid(){
outfile=$1
cat >$outfile <<EOF
target_grid: {
   longitude = {};
   latitude  = {};
};
data_products: {
};
EOF
}

test_missing_latitude(){
outfile=$1
cat >$outfile <<EOF
target_grid: {
   longitude = {delta = 0.05; min = -155.0; num = 2611;};
   latitude  = {};
};
data_products: {
};
EOF
}

test_missing_data_products(){
outfile=$1
cat >$outfile <<EOF
target_grid: {
   longitude = {delta = 0.05; min = -155.0; num = 2611;};
   latitude  = {delta = 0.05; min =   17.0; num =  940;};
};
EOF
}

test_missing_file_list(){
outfile=$1
cat >$outfile <<EOF
target_grid: {
   longitude = {delta = 0.05; min = -155.0; num = 2611;};
   latitude  = {delta = 0.05; min =   17.0; num =  940;};
};
data_products: {
   prod:{
      name = "prod";
   };
};
EOF
}

default_input_files_list(){
echo "foo.nc\nfoo.nc" > .test_input_files.lst
}

test_missing_input_files(){
outfile=$1
echo "foo.nc" > .test_input_files.lst
cat >$outfile <<EOF
target_grid: {
   longitude = {delta = 0.05; min = -155.0; num = 2611;};
   latitude  = {delta = 0.05; min =   17.0; num =  940;};
};
data_products: {
   prod:{
     name="prod";
     filename_list = "$FILENAME_LIST";
   };
};
EOF
}

test_missing_vars(){
outfile=$1
default_input_files_list
cat >$outfile <<EOF
target_grid: {
   longitude = {delta = 0.05; min = -155.0; num = 2611;};
   latitude  = {delta = 0.05; min =   17.0; num =  940;};
};
data_products: {
   prod:{
     name="prod";
     filename_list = "$FILENAME_LIST";
   };
};
EOF
}

test_zero_vars(){
outfile=$1
default_input_files_list
cat >$outfile <<EOF
target_grid: {
   longitude = {delta = 0.05; min = -155.0; num = 2611;};
   latitude  = {delta = 0.05; min =   17.0; num =  940;};
};
data_products: {
   prod:{
     name="prod";
     filename_list = "$FILENAME_LIST";
     vars = ();
   };
};
EOF
}

test_var_missing_in(){
outfile=$1
default_input_files_list
cat >$outfile <<EOF
target_grid: {
   longitude = {delta = 0.05; min = -155.0; num = 2611;};
   latitude  = {delta = 0.05; min =   17.0; num =  940;};
};
data_products: {
   prod:{
     name="prod";
     filename_list = "$FILENAME_LIST";
     vars = (
       {out="foo";}
     );
   };
};
EOF
}

test_bad_input_files(){
outfile=$1
echo "foo.nc\n7" > .test_input_files.lst
cat >$outfile <<EOF
target_grid: {
   longitude = {delta = 0.05; min = -155.0; num = 2611;};
   latitude  = {delta = 0.05; min =   17.0; num =  940;};
};
data_products: {
   prod:{
     name="prod";
     filename_list = "$FILENAME_LIST";
     vars = (
       {in="foo";}
     );
   };
};
EOF
}

test_var_bad_bitfield_type(){
outfile=$1
default_input_files_list
cat >$outfile <<EOF
target_grid: {
   longitude = {delta = 0.05; min = -155.0; num = 2611;};
   latitude  = {delta = 0.05; min =   17.0; num =  940;};
};
data_products: {
   prod:{
     name="prod";
     filename_list = "$FILENAME_LIST";
     vars = (
       {in="bar"; bitfield_type=0;},
       {in="foo"; bitfield_type=7;}
     );
   };
};
EOF
}

test_missing_longlat_group(){
outfile=$1
default_input_files_list
cat >$outfile <<EOF
target_grid: {
   longitude = {delta = 0.05; min = -155.0; num = 2611;};
   latitude  = {delta = 0.05; min =   17.0; num =  940;};
};
data_products: {
   prod:{
     name="prod";
     filename_list = "$FILENAME_LIST";
     vars = (
       {in="foo";}
     );
   };
};
EOF
}

test_have_longlat_group(){
outfile=$1
default_input_files_list
cat >$outfile <<EOF
target_grid: {
   longitude = {delta = 0.05; min = -155.0; num = 2611;};
   latitude  = {delta = 0.05; min =   17.0; num =  940;};
};
data_products: {
   prod:{
     name="prod";
     filename_list = "$FILENAME_LIST";
     longlat_group = {in="geolocation"};
     vars = (
       {in="foo";}
     );
   };
};
EOF
}

ARCHOBJS=$1
TMPDIR=tmp_cfgtest
/bin/mkdir -p $TMPDIR || exit 1

perform_test(){
 func=$1
 cfg_file=$TMPDIR/${1}.cfg
 log_file=$TMPDIR/log_${1}
 $func $cfg_file
 ../src/${ARCHOBJS}/L2_regrid $cfg_file > $log_file 2>&1
 exit_status=$?
 if test $exit_status -eq 0; then
    echo "FAIL[$1]:  expected non-zero exit status from ${cfg_file}"
    exit 1
 fi
}

perform_test test_syntax
perform_test test_missing_target_grid
perform_test test_missing_latitude
perform_test test_defining_target_grid
perform_test test_missing_data_products
perform_test test_missing_file_list
perform_test test_missing_input_files
perform_test test_missing_vars
perform_test test_zero_vars
perform_test test_var_missing_in
perform_test test_bad_input_files
perform_test test_var_bad_bitfield_type
perform_test test_missing_longlat_group
perform_test test_have_longlat_group
