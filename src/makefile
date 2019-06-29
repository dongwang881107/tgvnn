ifeq ($(PREFIX),)
	PREFIX := /usr/local
endif

CC=nvcc -Wno-deprecated-gpu-targets
CFLAGS=-O3 -D_FORCE_INLINES -DCUDA_HOST_MALLOC -g
LFLAGS=-lcufft -lcublas -lcusolver

objects = runtgv.o tgv.o ra.o

all: main

main: $(objects)
	$(CC) $(LFLAGS) $(objects) -o runtgv
%.o: %.cu
	$(CC) $(CFLAGS) -dc $< -o $@
clean:
	rm -f *.o runtgv
install:
	install -m 0755 runtgv $(PREFIX)/bin
