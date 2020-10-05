

bin/dap: bin/cmain.o bin/probdata_dap.o bin/reader_dap.o
	gcc -o bin/dap bin/cmain.o bin/probdata_dap.o bin/reader_dap.o -lscip

bin/cmain.o: src/cmain.c
	gcc -c -o bin/cmain.o src/cmain.c

bin/probdata_dap.o: src/probdata_dap.c src/probdata_dap.h
	gcc -c -o bin/probdata_dap.o src/probdata_dap.c

bin/reader_dap.o: src/reader_dap.c src/reader_dap.h
	gcc -c -o bin/reader_dap.o src/reader_dap.c

.PHONY: clean

clean:
	rm -f bin/*.o bin/dap

