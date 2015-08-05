#
config_file = Doxyfile

docs:
	doxygen -s $(config_file)

clean:
	-/bin/rm -rf html/* latex/*
	-/bin/rmdir html latex

distclean:
	$(MAKE) clean
