# -*- sh -*-

export PATH := ../src/gnuobjs:$(PATH)

EPOCH="1980-01-06T00:00:00Z"
MANEUVERFILE="tempo_d15479m53112000u000_r0_maneuver.csv"

PLAN_CONFIG="plan_opt1_panhandle_4.cfg"

#DEBUG = gdb --args
#DEBUG = valgrind --tool=memcheck --leak-check=yes --error-limit=no --num-callers=25 --log-file=valgrind.log

check-plan-std:
	@printf "Checking plan: scan=std\n"
	@source ../../pipe/etc/sdpc_env.sh && \
	$(DEBUG) plan --config $(PLAN_CONFIG) \
	     --epoch $(EPOCH) \
	     --maneuver $(MANEUVERFILE) \
	     --date 2022-05-27 --ndays 14 \
	     --scan std \
	     --output test_earthscan_std.csv
	@diff --ignore-matching-lines "^# Created" \
             earthscan_std.csv test_earthscan_std.csv
	@printf "OK\n"

check-plan-opt1:
	@printf "Checking plan: scan=opt1\n"
	@source ../../pipe/etc/sdpc_env.sh && \
	$(DEBUG) plan --config $(PLAN_CONFIG) \
	     --epoch $(EPOCH) \
	     --maneuver $(MANEUVERFILE) \
	     --date 2022-05-27 --ndays 14 \
	     --scan opt1 \
	     --output test_earthscan_opt1.csv
	@diff --ignore-matching-lines "^# Created" \
             earthscan_opt1.csv test_earthscan_opt1.csv
	@printf "OK\n"

check-plan-split-opt1-cbm:
	@printf "Checking plan: scan=split-opt1 with custom CBM\n"
	@source ../../pipe/etc/sdpc_env.sh && \
	$(DEBUG) plan --config $(PLAN_CONFIG) \
	     --epoch $(EPOCH) \
	     --maneuver $(MANEUVERFILE) \
	     --date 2022-05-27 --ndays 14 \
	     --scan split-opt1-Panhandle-4 \
	     --output test_earthscan_opt1_panhandle_4.csv
	@diff --ignore-matching-lines "^# Created" \
             earthscan_opt1_panhandle_4.csv test_earthscan_opt1_panhandle_4.csv
	@printf "OK\n"

check-sza:
	(. ../../pipe/etc/sdpc_env.sh ; ./check_sza.sl)

check:	check-sza check-plan-std check-plan-opt1 check-plan-split-opt1-cbm

TEST_OUTPUT = test_earthscan_std.csv test_earthscan_opt1.csv test_earthscan_opt1_panhandle_4.csv

clean:
	/bin/rm -f $(TEST_OUTPUT)
