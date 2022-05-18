CC = gcc-9
CFLAGS = -O2 -fopenmp -pthread -mavx -mavx2 -mfma
LIBS = -lm
all : base_test fast_test check_test

test_fast_tmp: fast_test tmp_test check_fast_tmp

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
	@-rm -f *.o gbfloat_m128
	$(CC) $(CFLAGS) -o gbfloat_tmp gbfloat_m128.c $(LIBS)

check:
	@-rm -f *.o test_accuracy
	$(CC) $(CFLAGS) -o test_accuracy test_accuracy.c $(LIBS)

base_test: base
	./gbfloat_base test1.jpg test_base.jpg 0.6 -2.0 2.0 1001 501

fast_test: fast
	./gbfloat_fast test1.jpg test_fast.jpg 0.6 -2.0 2.0 1001 501

tmp_test: tmp
	./gbfloat_tmp test1.jpg test_tmp.jpg 0.6 -2.0 2.0 1001 501

check_test: check
	./test_accuracy test_base.jpg test_fast.jpg

check_fast_tmp: check
	./test_accuracy test_fast.jpg test_tmp.jpg