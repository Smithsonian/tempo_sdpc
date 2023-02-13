# -*- sh -*-

export PATH := ../src/gnuobjs:$(PATH)
export LD_LIBRARY_PATH := ../../libtio/src/gnuobjs:$(LD_LIBRARY_PATH)

PLAN_EXEC = ../src/gnuobjs/plan

EPOCH="1980-01-06T00:00:00Z"

#DEBUG = gdb --args
#DEBUG = valgrind --tool=memcheck --leak-check=yes --error-limit=no --num-callers=25 --log-file=valgrind.log

check-plan-std:
	@printf "Checking plan: scan=std\n"
	@source ../../pipe/etc/sdpc_env.sh && \
	$(DEBUG) $(PLAN_EXEC) --config "plan_opt1_panhandle_4.cfg" \
	     --epoch $(EPOCH) \
	     --maneuver "tempo_d15479m53112000u000_r0_maneuver.csv" \
	     --date 2022-05-27 --ndays 14 \
	     --scan std \
	     --output test_earthscan_std.csv
	@diff --ignore-matching-lines "^#" \
             earthscan_std.csv test_earthscan_std.csv
	@printf "OK\n"

check-plan-opt1:
	@printf "Checking plan: scan=opt1\n"
	@source ../../pipe/etc/sdpc_env.sh && \
	$(DEBUG) $(PLAN_EXEC) --config "plan_opt1_panhandle_4.cfg" \
	     --epoch $(EPOCH) \
	     --maneuver "tempo_d15479m53112000u000_r0_maneuver.csv" \
	     --date 2022-05-27 --ndays 14 \
	     --scan opt1 \
	     --output test_earthscan_opt1.csv
	@diff --ignore-matching-lines "^#" \
             earthscan_opt1.csv test_earthscan_opt1.csv
	@printf "OK\n"

check-plan-split-opt1-cbm:
	@printf "Checking plan: scan=split-opt1 with custom CBM\n"
	@source ../../pipe/etc/sdpc_env.sh && \
	$(DEBUG) $(PLAN_EXEC) --config "plan_opt1_panhandle_4.cfg" \
	     --epoch $(EPOCH) \
	     --maneuver "tempo_d15479m53112000u000_r0_maneuver.csv" \
	     --date 2022-05-27 --ndays 14 \
	     --scan split-opt1-Panhandle-4 \
	     --output test_earthscan_opt1_panhandle_4.csv
	@diff --ignore-matching-lines "^#" \
             earthscan_opt1_panhandle_4.csv test_earthscan_opt1_panhandle_4.csv
	@printf "OK\n"

check-plan-twilight:
	@printf "Checking plan: scan=opt1 + twilight\n"
	@source ../../pipe/etc/sdpc_env.sh && \
	$(DEBUG) $(PLAN_EXEC) --config "plan_citylights.cfg" \
	     --epoch $(EPOCH) \
	     --maneuver "tempo_d15682m00018000u000_r0_maneuver.csv" \
	     --date 2022-12-15 --ndays 14 \
	     --scan opt1 --nightlights \
	     --output test_earthscan_opt1_twilight.csv
	@diff --ignore-matching-lines "^#" \
             earthscan_opt1_twilight.csv test_earthscan_opt1_twilight.csv
	@printf "OK\n"

check-sza:
	(. ../../pipe/etc/sdpc_env.sh ; ./check_sza.sl)

check:	check-sza check-plan-std check-plan-opt1 check-plan-split-opt1-cbm check-plan-twilight

TEST_OUTPUT = test_earthscan_std.csv test_earthscan_opt1.csv test_earthscan_opt1_panhandle_4.csv \
              test_earthscan_opt1_twilight.csv

clean:
	/bin/rm -f $(TEST_OUTPUT)
