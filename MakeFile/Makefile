BAMTOOLS=../Bamtools/include
BAMTOOLS_LD=../Bamtools/lib
#RootPath = ../ScaffoldGapFillerTest/
#NeedleManPath = ../Depends/needleman_wunsch/
#NormDependsPath = ../Depends/
#INCLUDES:= ../Depends ../Depends/needleman_wunsch ../Bamtools/include/  # 头文件目录

CC = g++
CFLAGS = -O3 -Wall -static -I$(BAMTOOLS) -L$(BAMTOOLS_LD) -Wl,-rpath,$(BAMTOOLS_LD)
LIBS += -lpthread
INCLUDEPATH += tr1/memory

OBJ = ./obj/

all: HuffmanCoding

recipe-577480-1.o: recipe-577480-1.cpp
	$(CC) $(CFLAGS) -c recipe-577480-1.cpp

clsbasealgorithm.o : clsbasealgorithm.cpp clsbasealgorithm.h
	$(CC) $(CFLAGS) -c clsbasealgorithm.cpp

clsbwa.o: clsbwa.cpp clsbwa.h
	$(CC) $(CFLAGS) -c clsbwa.cpp

clsfastareader.o: clsfastareader.cpp clsfastareader.h
	$(CC) $(CFLAGS) -c clsfastareader.cpp

clsvelvet.o: clsvelvet.cpp clsvelvet.h
	$(CC) $(CFLAGS) -c clsvelvet.cpp

KmerUtils.o :KmerUtils.cpp KmerUtils.h
	$(CC) $(CFLAGS) -c KmerUtils.cpp

clsmultithread.o :clsmultithread.cpp clsmultithread.h
	$(CC) $(CFLAGS) -c clsmultithread.cpp

output.o: ./DNAZip/output.h ./DNAZip/output.cc
	$(CC) $(CFLAGS) -c ./DNAZip/output.cc

input.o: ./DNAZip/input.h ./DNAZip/input.cc
	$(CC) $(CFLAGS) -c ./DNAZip/input.cc

bitfile.o: ./DNAZip/bitfile.h ./DNAZip/bitfile.cc
	$(CC) $(CFLAGS) -c ./DNAZip/bitfile.cc

huffman.o: ./DNAZip/huffman.h ./DNAZip/huffman.cc
	$(CC) $(CFLAGS) -c ./DNAZip/huffman.cc

clsmethod.o: clsmethod.cpp clsmethod.h
	$(CC) $(CFLAGS) -c clsmethod.cpp

main.o: main.cpp 
	$(CC) $(CFLAGS) -c main.cpp

HuffmanCoding: recipe-577480-1.o clsbasealgorithm.o clsbwa.o clsfastareader.o clsvelvet.o KmerUtils.o clsmultithread.o output.o input.o bitfile.o huffman.o clsmethod.o main.o 
	$(CC) $(CFLAGS) -o HuffmanCoding recipe-577480-1.o clsbasealgorithm.o clsbwa.o clsfastareader.o clsvelvet.o KmerUtils.o clsmultithread.o output.o input.o bitfile.o huffman.o clsmethod.o main.o \
	-lbamtools -lz

RM-F := rm -f
.PHONY : all deps objs clean veryclean rebuild info
clean :	
	@$(RM-F) *.o
	@$(RM-F) *.d

veryclean: clean
	@$(RM-F) CheckScafInfo


