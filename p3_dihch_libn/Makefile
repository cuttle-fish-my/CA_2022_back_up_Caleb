CC = gcc
CFLAGS = -O2 -fopenmp -pthread -mavx -mavx2
LIBS = -lm
all : base_test fast_test check_test

base_fast_test: base_test fast_test check_test
base_tmp_test: base_test tmp_test check_base_tmp
fast_tmp_test: fast_test tmp_test check_fast_tmp

clean:
	@-rm -f *.o gbfloat_base
	@-rm -f *.o gbfloat_fast
	@-rm -f *.o test_accuracy
	
base:
	@-rm -f *.o gbfloat_base
	$(CC) $(CFLAGS) -o gbfloat_base gbfloat_base.c $(LIBS)

fast:
	@-rm -f *.o gbfloat_fast
	$(CC) $(CFLAGS) -o gbfloat_fast gbfloat_fast.c $(LIBS)

tmp:
	@@-rm -f *.o gbfloat_tmp
	$(CC) $(CFLAGS) -o gbfloat_tmp gbfloat_tmp.c $(LIBS)

check:
	@-rm -f *.o test_accuracy
	$(CC) $(CFLAGS) -o test_accuracy test_accuracy.c $(LIBS)

base_test: base
	./gbfloat_base test.jpg test_base.jpg 0.6 -2.0 2.0 1001 201

fast_test: fast
	./gbfloat_fast test.jpg test_fast.jpg 0.6 -2.0 2.0 1001 201

tmp_test: tmp
	./gbfloat_tmp test.jpg test_m256d.jpg 0.6 -2.0 2.0 1001 201

check_test: check
	./test_accuracy test_base.jpg test_fast.jpg

check_fast_tmp: check
	./test_accuracy test_fast.jpg test_tmp.jpg

check_base_tmp: check
	./test_accuracy test_base.jpg test_tmp.jpg