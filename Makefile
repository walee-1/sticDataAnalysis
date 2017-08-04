CC=g++

CCFLAGS=`root-config --cflags`
CFLAGS=`root-config --cflags`

INCLUDE= -Iinclude/

LIBS=-lSpectrum -lrt
IMPORT_LIBS=
LDFLAGS=-Wl,--no-as-needed `root-config --glibs`

#binaries=CoinS
binaries=Coin16
obj=src/tdc_ch_values.o

all: $(binaries)

reader: src/EventDict.o $(obj)
	$(CC) $^ $(LDFLAGS) $(LIBS) -o reader

Coin16: src/EventDict.o $(obj) src/Coin16.o
	$(CC) $^ $(CFLAGS) $(LDFLAGS) $(LIBS) -o Coin16

CoinS: src/EventDict.o $(obj) src/CoinS.o
	$(CC) $^ $(CFLAGS) $(LDFLAGS) $(LIBS) -o CoinS

src/EventDict.cpp:
	@echo "Generating Dictionary $@..."
	cd include ; rootcint -f EventDict.cpp -c $(CFLAGS) -p EventType.h TSTiC2_Ana.h LinkDef.h ; mv EventDict.cpp ../src/
	@echo "Done generating, exit code: $!"

clear:
	rm -f src/*.o
	rm -f include/EventDict.h
	rm -f src/EventDict.cpp

clean:	clear
	rm -f $(binaries)

%.o: %.cpp
	$(CC) $(CFLAGS) $(INCLUDE) $(LIBS) -include sstream -o $@ -g -c $<

