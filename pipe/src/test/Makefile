include ../../../Makefile.defines

objdir=../$(archobjs)

CACHEMON_TESTS=$(shell grep -v '\#' tests.lis)

cachemon-tests:
	@./run_tests.sh $(CACHEMON_TESTS)

export FILEDB_TEST_ROOT_DIR=/tmp/filedb_test.$(USER)

filedb-test:
	./create_filedb_test.sl 7 $(FILEDB_TEST_ROOT_DIR)
	./test_filedb.sh $(objdir)/filedb

check:	filedb-test cachemon-tests
