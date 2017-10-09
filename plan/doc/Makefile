#
doxy_file = Doxyfile

docs:
	doxygen -s $(doxy_file)

pngs:
	-cd figures ; make pngs

clean:
	-/bin/rm -rf html/* latex/*
	-/bin/rmdir html latex
#	-cd figures ; $(MAKE) clean

distclean:
	$(MAKE) clean
