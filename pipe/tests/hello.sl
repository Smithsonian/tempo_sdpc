#! /usr/bin/env slsh

prepend_to_slang_load_path ("../src");
require ("daemon");

daemonize ("hello", "$PWD/hello.log"$);
loop (3)
{
   daemon_log (LOG_INFO, "Hello.");
}

