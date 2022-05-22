#gcc -shared -c -Wpedantic -Wall -fPIC -Wextra -Werror -std=c89 ringbuffer.c -o libringbuffer.so
#gcc -c -Wpedantic -Wall -fPIC -Wextra -Werror -std=c89 test.c
#gcc -g -shared -Wl,-soname,libringbuffer.so test.o -o dynamicringbuffer
gcc -c -fPIC ringbuffer.c
gcc -shared -o libringbuffer.so ringbuffer.o
gcc -c -fPIC test.c
gcc -L. -Wl,-rpath=. -o dynamicringbuffer test.o libringbuffer.so -lringbuffer
