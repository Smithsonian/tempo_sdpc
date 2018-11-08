include ../../Makefile.defines

objdir=../src/$(archobjs)

CACHEMON_TESTS=$(shell grep -v '\#' tests.lis)

cachemon-tests:
	@./run_tests.sh $(CACHEMON_TESTS)

export FILEDB_TEST_ROOT_DIR=/tmp/filedb_test.$(USER)

filedb-test:
	./create_filedb_test.sl 7 $(FILEDB_TEST_ROOT_DIR)
	./test_filedb.sh $(objdir)/filedb

check:	filedb-test cachemon-tests

GCDA_FILES_DIR = "../src/gnuobjs"

coverage: 
	# zero all execution counters
	lcov -z -d $(GCDA_FILES_DIR) -d .
	# create baseline coverage data file
	lcov -c -i -d $(GCDA_FILES_DIR) -d . -o base.info
	# perform test
	$(MAKE) check
	# create test coverage data file
	lcov -c -d $(GCDA_FILES_DIR) -d . -o test_cov.info
	# combine baseline and test coverage data
	lcov -a base.info -a test_cov.info -o total.info
	# remove selected paths from coverage statistics
	lcov -r total.info $(OTS_ROOT)/include/cfortran.h -o filtered.info
	# generate summary web pages
	genhtml filtered.info -o lcov.dir
