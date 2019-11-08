pcmDir=/home/linux/pcm/
mysqlDir=/usr/include/mysql-cppconn/jdbc/

CC=g++
CFLAGS = -std=c++0x -I. -I/usr/include/cppconn
DEPS = $(pcmDir)client_bw.o $(pcmDir)cpucounters.o $(pcmDir)msr.o $(pcmDir)pci.o

all: quikoloSlave quikoloMaster
.PHONY: all

quikoloSlave: $(DEPS)
	$(CC) $(CFLAGS) -I$(pcmDir) -I$(mysqlDir) -o quikoloSlave $(DEPS) quikoloSlaveOneTable.cpp -lpthread -lmysqlcppconn

quikoloMaster:
	$(CC) $(CFLAGS) -I$(mysqlDir) -o quikoloMaster quikoloMasterPredictOneTable.cpp -lmysqlcppconn -fopenmp

clean:
	rm -f quikoloMaster
	rm -f quikoloSlave
