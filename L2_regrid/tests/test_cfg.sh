#! /bin/bash

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

test_missing_output_file(){
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

test_missing_input_files(){
outfile=$1
cat >$outfile <<EOF
target_grid: {
   longitude = {delta = 0.05; min = -155.0; num = 2611;};
   latitude  = {delta = 0.05; min =   17.0; num =  940;};
};
data_products: {
   prod:{
     name="prod";
     output_file = "foo.nc";
   };
};
EOF
}

test_zero_input_files(){
outfile=$1
cat >$outfile <<EOF
target_grid: {
   longitude = {delta = 0.05; min = -155.0; num = 2611;};
   latitude  = {delta = 0.05; min =   17.0; num =  940;};
};
data_products: {
   prod:{
     name="prod";
     output_file = "foo.nc";
     input_files = [
     ];
   };
};
EOF
}

test_missing_vars(){
outfile=$1
cat >$outfile <<EOF
target_grid: {
   longitude = {delta = 0.05; min = -155.0; num = 2611;};
   latitude  = {delta = 0.05; min =   17.0; num =  940;};
};
data_products: {
   prod:{
     name="prod";
     output_file = "foo.nc";
     input_files = [
        "foo.nc"
     ];
   };
};
EOF
}

test_zero_vars(){
outfile=$1
cat >$outfile <<EOF
target_grid: {
   longitude = {delta = 0.05; min = -155.0; num = 2611;};
   latitude  = {delta = 0.05; min =   17.0; num =  940;};
};
data_products: {
   prod:{
     name="prod";
     output_file = "foo.nc";
     vars = ();
     input_files = [
        "foo.nc"
     ];
   };
};
EOF
}

test_var_missing_in(){
outfile=$1
cat >$outfile <<EOF
target_grid: {
   longitude = {delta = 0.05; min = -155.0; num = 2611;};
   latitude  = {delta = 0.05; min =   17.0; num =  940;};
};
data_products: {
   prod:{
     name="prod";
     output_file = "foo.nc";
     vars = (
       {out="foo";}
     );
     input_files = [
        "foo.nc"
     ];
   };
};
EOF
}

test_bad_input_files(){
outfile=$1
cat >$outfile <<EOF
target_grid: {
   longitude = {delta = 0.05; min = -155.0; num = 2611;};
   latitude  = {delta = 0.05; min =   17.0; num =  940;};
};
data_products: {
   prod:{
     name="prod";
     output_file = "foo.nc";
     vars = (
       {in="foo";}
     );
     input_files = [
       7
     ];
   };
};
EOF
}

test_var_bad_bitfield_type(){
outfile=$1
cat >$outfile <<EOF
target_grid: {
   longitude = {delta = 0.05; min = -155.0; num = 2611;};
   latitude  = {delta = 0.05; min =   17.0; num =  940;};
};
data_products: {
   prod:{
     name="prod";
     output_file = "foo.nc";
     vars = (
       {in="bar"; bitfield_type=0;},
       {in="foo"; bitfield_type=7;}
     );
     input_files = [
        "foo.nc"
     ];
   };
};
EOF
}

test_missing_longlat_group(){
outfile=$1
cat >$outfile <<EOF
target_grid: {
   longitude = {delta = 0.05; min = -155.0; num = 2611;};
   latitude  = {delta = 0.05; min =   17.0; num =  940;};
};
data_products: {
   prod:{
     name="prod";
     output_file = "foo.nc";
     vars = (
       {in="foo";}
     );
     input_files = [
       "bar.nc"
     ];
   };
};
EOF
}

test_have_longlat_group(){
outfile=$1
cat >$outfile <<EOF
target_grid: {
   longitude = {delta = 0.05; min = -155.0; num = 2611;};
   latitude  = {delta = 0.05; min =   17.0; num =  940;};
};
data_products: {
   prod:{
     name="prod";
     output_file = "foo.nc";
     longlat_group = {in="geolocation"};
     vars = (
       {in="foo";}
     );
     input_files = [
       "bar.nc"
     ];
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
perform_test test_missing_input_files
perform_test test_missing_output_file
perform_test test_zero_input_files
perform_test test_missing_vars
perform_test test_zero_vars
perform_test test_var_missing_in
perform_test test_bad_input_files
perform_test test_var_bad_bitfield_type
perform_test test_missing_longlat_group
perform_test test_have_longlat_group
