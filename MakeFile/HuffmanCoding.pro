TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
LIBS += -lz
CONFIG += thread #This is to support Multi_Thread (pthread!!!)

SOURCES += main.cpp \
    recipe-577480-1.cpp \
    FastqFileParse/kseq_test.c \
    clsbasealgorithm.cpp \
    clsmethod.cpp \
    clsfastareader.cpp \
    clsvelvet.cpp \
    clsbwa.cpp \
    DNAZip/watsonFileGenerator.cc \
    DNAZip/perftest.cc \
    DNAZip/output.cc \
    DNAZip/kMerFreqGenerator.cc \
    DNAZip/input.cc \
    DNAZip/huffman.cc \
    DNAZip/dbSNPDeCompression.cc \
    DNAZip/dbSNPCompression.cc \
    DNAZip/bitfile.cc \
    KmerUtils.cpp \
    clsmultithread.cpp

OTHER_FILES += \
    ../../../ThirdPartyTools/DNAZip/makefile \
    DNAZip/makefile

HEADERS += \
    FastqFileParse/kseq.h \
    clsbasealgorithm.h \
    clsfastareader.h \
    clsvelvet.h \
    clsbwa.h \
    DNAZip/watsonFileGenerator.h \
    DNAZip/output.h \
    DNAZip/kMerFreqGenerator.h \
    DNAZip/input.h \
    DNAZip/huffman.h \
    DNAZip/fileGenerator.h \
    DNAZip/dbSNPDeCompression.h \
    DNAZip/dbSNPCompression.h \
    DNAZip/dbSNPCompression.cc.autosave \
    DNAZip/bitfile.h \
    clsmethod.h \
    KmerUtils.h \
    clsmultithread.h


unix:!macx: LIBS += -L$$PWD/../Bamtools/lib/ -lbamtools

INCLUDEPATH += $$PWD/../Bamtools/include
DEPENDPATH += $$PWD/../Bamtools/include

unix:!macx: PRE_TARGETDEPS += $$PWD/../Bamtools/lib/libbamtools.a
