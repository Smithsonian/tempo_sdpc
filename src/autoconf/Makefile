
../L1b_L2/configure: aclocal.m4 L1b_L2_configure.ac
	/bin/rm -rf autom4te.cache
	autoconf -o ../L1b_L2/configure L1b_L2_configure.ac 
update: config.sub config.guess
config.guess: /usr/share/misc/config.guess
	/bin/cp -f /usr/share/misc/config.guess config.guess
config.sub: /usr/share/misc/config.sub
	/bin/cp -f /usr/share/misc/config.sub config.sub
