gcc -c -Wpedantic -Wall -Wextra -Werror -std=c89 ringbuffer.c test.c
ar rcs libringbuffer.a ringbuffer.o
ld -o staticringbuffer test.o libringbuffer.a -lc -e main