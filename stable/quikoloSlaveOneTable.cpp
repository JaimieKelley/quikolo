#include <cpucounters.h>
#include <mysql_connection.h>
#include <cppconn/driver.h>
#include <cppconn/exception.h>
#include <cppconn/resultset.h>
#include <cppconn/statement.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <fcntl.h>
#include <sys/types.h>
#include <unistd.h>
#include <errno.h>
#include <netdb.h>
#include <arpa/inet.h>
#include <sys/wait.h>
#include <signal.h>
#include <assert.h>
#include <unistd.h>
#include <errno.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
using namespace std;

#define PORT "1072"
#define MAX_LINE 16384
#define BACKLOG 10

void *get_in_addr(struct sockaddr *sa)
{
  if(sa->sa_family == AF_INET) {
    return &(((struct sockaddr_in*)sa)->sin_addr);
  }

  return &(((struct sockaddr_in6*)sa)->sin6_addr);
}

void writeToFile(string fName, int sec, int milli, string values)
{
  stringstream fNameSS;
  fNameSS << "traces/" << fName << "hardware.txt";
  string fileName = fNameSS.str();
  //cout << fileName << endl;
  fstream outStream(fileName, fstream::out | fstream::app);

  if(outStream){
    outStream << sec << " " << milli << " " << values;
  }
  else
  {
    cout << "write failed" << endl;
  }

  outStream.close();
}

int main(int argc, char * argv[]) {
  string line;
  string diskstr = "sda ";
  string diskstr2 = "sdb ";
  string cpustr = "cpu";
  string memFstr = "MemFree:";
  string memTstr = "MemTotal:";
  string iostr = ":";
  string table = "hardware";
  int prevReadms = 0;
  int prevWritems = 0;
  int totalreadms = 0;
  int totalwritems = 0;
  int netReadms = 0;
  int netWritems = 0;
  double diskUtilRead = 0.0;
  double diskUtilWrite = 0.0;
  int numDlines = 10;
  long int tprofile = 0;
  long int iseconds = 0;
  int imilliseconds = 0;
  int rseconds =0;
  long int seconds = 0;
  int numCPU = 16;
  int milliseconds = 0;
  int prevSeconds = 0;
  int prevMilliseconds = 0;
  int netMS;
  int cpu = 0;
  unsigned long long int prevCPU[numCPU+1];
  unsigned long long int prevIdle[numCPU+1];
  unsigned long long int totalCPU[numCPU+1];
  unsigned long long int idle[numCPU+1];
  int memTotal = 0;
  int memFree = 0;
  int prevReceived = 0;
  int prevTransmitted = 0;
  int received = 0;
  int transmitted = 0;
  int netReceived = 0;
  int netTransmitted = 0;
  double cpuUtil[numCPU+1];
  double memUtil = 0.0;
  unsigned long long int netTotalCPU[numCPU+1];
  unsigned long long int netCPUIdle[numCPU+1];
  struct timespec interval, ainterval;
  SystemCounterState after_system, before_system;
  vector<SocketCounterState> after_socket;
  vector<SocketCounterState> before_socket;
  vector<CoreCounterState> after_core;
  vector<CoreCounterState> before_core;
  double l3cacheHitRatio = 0.0; // miss ratio is calculated
  double l2cacheHitRatio = 0.0; // miss ratio is calculated
  bool subsequentRun = false;
  bool profile = false;
  int firsttime = 0;
  //SetTree times = new SetTree();
  long int lastprofile = 0;
  string fNameNew;
  string fileName = argv[1];

  fd_set master, read_fds;
  int fdmax;
  int listener;
  int new_fd;
  struct addrinfo hints, *ai, *p;
  struct sockaddr_storage remoteaddr; // connector's address info
  socklen_t addrlen;
  char buf[256];
  int nbytes;
  struct sigaction sa;
  int yes=1;
  char remoteIP[INET6_ADDRSTRLEN];
  int j, k, rv;

  FD_ZERO(&master);
  FD_ZERO(&read_fds);

  memset(&hints, 0, sizeof hints);
  hints.ai_family = AF_INET;
  hints.ai_socktype = SOCK_STREAM;
  hints.ai_flags = AI_PASSIVE; // use my IP

  if((rv = getaddrinfo(NULL, PORT, &hints, &ai)) != 0) {
	fprintf(stderr, "selectserver: %s\n", gai_strerror(rv));
	return 1;
  }

  // loop through all results and bind to first we can
  for(p = ai; p != NULL; p=p->ai_next) {
    listener = socket(p->ai_family, p->ai_socktype, p->ai_protocol);
    if(listener < 0) {
	continue;
    }

    setsockopt(listener, SOL_SOCKET, SO_REUSEADDR, &yes, sizeof(int));

    if(bind(listener, p->ai_addr, p->ai_addrlen) == -1) {
	close(listener);
	perror("selectserver: bind");
	continue;
    }

    break;
  }

  freeaddrinfo(ai); // all done with this structure

  if(p == NULL) {
	fprintf(stderr, "selectserver: failed to bind\n");
	exit(1);
  }

  if(listen(listener, 10) == -1) {
	perror("listen");
	exit(3);
  }

  FD_SET(listener, &master);

  fdmax = listener;

  printf("Server: waiting for connections...\n");

  for(int i = 0; i <= numCPU; i++)
  {
	prevCPU[i] = 0ULL;
	prevIdle[i] = 0ULL;
	totalCPU[i] = 0ULL;
	idle[i] = 0ULL;
	netTotalCPU[i] = 0ULL;
	netCPUIdle[i] = 0ULL;
  }

  PCM * m = PCM::getInstance();

  // program counters, and on a failure just exit

  if (m->program() != PCM::Success) return 1; 

  while(1) {
    firsttime++;

    //cout << "entered while loop\n";

    timeval tv_seconds;
    gettimeofday(&tv_seconds, NULL);
    iseconds = tv_seconds.tv_sec;
    imilliseconds = tv_seconds.tv_usec / 1000;
    tprofile = iseconds * 1000 + imilliseconds;

    timeval wait_seconds;
    wait_seconds.tv_sec = 0;
    wait_seconds.tv_usec = 000000100;

    read_fds = master;
    if(select(fdmax+1, &read_fds, NULL, NULL, &wait_seconds) == -1) {
	perror("select");
	exit(4);
    }

    for(k = 0; k <= fdmax; k++) {
	if(FD_ISSET(k, &read_fds)) {
		if(k == listener) {
			addrlen = sizeof remoteaddr;
			new_fd = accept(listener, (struct sockaddr *)&remoteaddr, &addrlen);

			if(new_fd == -1) {
				perror("accept");
			} else {
				FD_SET(new_fd, &master);
				if(new_fd > fdmax) {
					fdmax = new_fd;
				}
				printf("selectserver: new connection from %s on socket %d\n", inet_ntop(remoteaddr.ss_family, get_in_addr((struct sockaddr*)&remoteaddr), remoteIP, INET6_ADDRSTRLEN), new_fd);
			}
		} else {
			if ((nbytes = recv(k, buf, sizeof buf, 0)) <= 0) {
				if(nbytes == 0) {
					printf("selectserver: socket %d hung up\n", k);					} else {
					perror("server recv");
				}
				close(k);
				FD_CLR(k, &master);
			} else {
				if(strstr(buf, "start\n"))
				{
					profile = true;
					cout << "start tstamp " << iseconds << " " << iseconds * 1000 << " " << milliseconds << " " << tprofile << "\n";
					sprintf(buf, "%Ld\n", tprofile);
					if(send(k, buf, sizeof(buf), 0) == -1) {
						perror("server send timestamp failed");
					}
					stringstream ioFile;
					/*ioFile << buf;
					ioFile >> fNameNew;
					ioFile >> fNameNew;
					cout << "File Name: " << fNameNew;*/
				} else if(strstr(buf, "stop\n"))
				{
					profile = false;
					if(send(k, "-1\n", 3, 0) == -1) {
						perror("server send ok failed");
					}
					fNameNew = "";
				}
			}
		}
	}
    }



    if(profile)
    {
	if(tprofile < lastprofile+500)
	{
		struct timespec wait_seconds2;
        	wait_seconds2.tv_sec = 0;
		int millsec = 500 - (tprofile - lastprofile);
		if(millsec > 0)
		{
			cout << "wait " << millsec << endl;
        		wait_seconds2.tv_nsec = millsec * 1000000;
			nanosleep(&wait_seconds2, NULL);
		}
	}
    }

    if(profile && tprofile > lastprofile + 500)
    {
      // Here we are going to read the timestamp.
      gettimeofday(&tv_seconds, NULL);
      prevSeconds = seconds;
      prevMilliseconds = milliseconds;
      milliseconds = tv_seconds.tv_usec / 1000;
      seconds = tv_seconds.tv_sec;
      lastprofile = seconds*1000 + milliseconds;
      netMS = ((seconds*1000) + milliseconds) - ((prevSeconds * 1000) + prevMilliseconds);

      before_system = after_system;
      before_socket = after_socket;
      before_core = after_core;
      m->getAllCounterStates(after_system, after_socket, after_core);

      // Here we are going to read cpu utilization.
      ifstream cpufile ("/proc/stat");
      if(cpufile.is_open())
      {
        cpu = 0;
        while(getline(cpufile, line))
        {
	  if(line.find(cpustr) != string::npos)
	  {
	    string arr[10];
	    int i = 0;
	    prevCPU[cpu] = totalCPU[cpu];
	    prevIdle[cpu] = idle[cpu];
	    totalCPU[cpu] = 0ULL;
	    idle[cpu] = 0ULL;
	    stringstream ssin(line);
	    while(ssin.good() && i < 10)
	    {
	      ssin >> arr[i];
	      if(i > 0)
	      {
		char * pEnd;
		unsigned long long int temp = 0ULL;
		temp = strtoull((arr[i]).c_str(), &pEnd, 10);
		totalCPU[cpu] += temp;
	      }
	      i++;
	    }
	    char * iEnd;
	    idle[cpu] = strtoull((arr[4]).c_str(), &iEnd, 10);

	    cpu++;
	  }
        }

        cpufile.close();
      }

      for(int i = 0; i <= numCPU; i++)
      {
        netTotalCPU[i] = totalCPU[i] - prevCPU[i];
        netCPUIdle[i] = idle[i] - prevIdle[i];
        cpuUtil[i] = (double) (netTotalCPU[i] - netCPUIdle[i]) / (double) netTotalCPU[i];
      }

      // Here we are going to read memory.
      ifstream memfile ("/proc/meminfo");
      if(memfile.is_open())
      {
        memFree = 0;
        while(getline(memfile, line))
        {
	  if(memTotal == 0 && line.find(memTstr) != string::npos){
	    string arr[10];
	    int i = 0;
	    stringstream ssin(line);
	    while(ssin.good() && i < 3)
	    {
	      ssin >> arr[i];
	      i++;
	    }
	    memTotal += atoi((arr[1]).c_str());
	  }

	  if(line.find(memFstr) != string::npos){
	    string arr[10];
	    int i = 0;
	    stringstream ssin(line);
	    while(ssin.good() && i < 3)
	    {
	      ssin >> arr[i];
	      i++;
	    }
	    memFree += atoi((arr[1]).c_str());
	  }
        }

        memfile.close();
      }

      memUtil = (double) (memTotal - memFree) / (double) memTotal;

      // Here we are going to read network I/O.
      ifstream netfile("/proc/net/dev");
      if(netfile.is_open())
      {
        prevReceived = received;
        prevTransmitted = transmitted;
        received = 0;
        transmitted = 0;
        while(getline(netfile, line))
        {
	  if(line.find(iostr) != string::npos){
	    string arr[10];
	    int i = 0;
	    stringstream ssin(line);
	    while(ssin.good() && i < 10)
	    {
	      ssin >> arr[i];
	      i++;
	    }
	    received += atoi((arr[1]).c_str());
	    transmitted += atoi((arr[9]).c_str());
	  }
        }

        netfile.close();
      }

      netReceived = received - prevReceived;
      netTransmitted = transmitted - prevTransmitted;

      // Here we are going to read disk I/O.
      ifstream diskfile ("/proc/diskstats");
      if(diskfile.is_open())
      {
        prevReadms = totalreadms;
        prevWritems = totalwritems;
        totalreadms = 0;
        totalwritems = 0;
        numDlines = 0;
        while(getline(diskfile, line))
        {
	  if((line.find(diskstr) != string::npos) || (line.find(diskstr2) != string::npos)) {
	    string arr[15];
	    int i = 0;
	    int readms = 0;
	    int writems = 0;
	    stringstream ssin(line);
	    while (ssin.good() && i < 12)
	    {
	      ssin >> arr[i];
	      i++;
	    }
	    readms = atoi((arr[6]).c_str());
	    writems = atoi((arr[10]).c_str());
	    totalreadms += readms;
	    totalwritems += writems;
	    numDlines++;
	  }
        }	

        diskfile.close();
      }

      netReadms = totalreadms - prevReadms;
      netWritems = totalwritems - prevWritems;
      diskUtilRead = (double) netReadms / (double) (numDlines * netMS);
      diskUtilWrite = (double) netWritems / (double) (numDlines * netMS);

      if(subsequentRun)
      {
        try {
          sql::Driver *driver;
          sql::Connection *con;
          sql::Statement *stmt;

          /* Create a connection */
          driver = get_driver_instance();
          con = driver->connect("tcp://localhost:3306", "root", "root");

          /* Connect to the MySQL test database */
          con->setSchema("READINGS");
          stmt = con->createStatement();

          std::stringstream cpSS, intermedCPSS, toWrite;
          cpSS << "INSERT INTO " << table << " (seconds, milliseconds, ";
          intermedCPSS << " ) VALUES ( " << seconds << ", " << milliseconds << ", ";
          for(int i=1; i <= numCPU; i++)
          {
	    if(totalCPU[i] > 0) {
	      cpSS << "cpu" << (i-1) << ", ";
	      intermedCPSS << cpuUtil[i] << ", ";
	      toWrite << cpuUtil[i] << " ";
            }
          }

          for (int i = 0; i < after_core.size(); i++)
          {
            l2cacheHitRatio = 1 - getL2CacheHitRatio(before_core[i], after_core[i]);
            cpSS << "cpu" << i << "l2" << ", ";
            intermedCPSS << l2cacheHitRatio << ", ";
	    toWrite << l2cacheHitRatio << " ";
          }

          for (int i = 0; i < after_socket.size(); i++)
          {
            l3cacheHitRatio = 1 - getL3CacheHitRatio(before_socket[i], after_socket[i]);
            cpSS << "socket" << i << "l3, ";
            intermedCPSS << l3cacheHitRatio << ", ";
	    toWrite << l3cacheHitRatio << " ";
          }

          cpSS << "mem, received, transmitted, diskreads, diskwrites "; 
	  intermedCPSS << memUtil << ", " << netReceived << ", " << netTransmitted << ", " << diskUtilRead << ", " << diskUtilWrite << " )";
	  toWrite << memUtil << " " << netReceived << " " << netTransmitted << " " << diskUtilRead << " " << diskUtilWrite << endl;

	  writeToFile(fileName, seconds, milliseconds, toWrite.str());

          string cpS = cpSS.str() + intermedCPSS.str();
          stmt->execute(cpS);

          delete stmt;
          delete con;
        } catch (sql::SQLException &e) {
          cout << "#ERR: SQLException in " << __FILE__;
          cout << "(" << __FUNCTION__ << ") on line " << __LINE__ << endl;
          cout << "#ERR: " << e.what();
          cout << " (MySQL error code: " << e.getErrorCode();
          cout << ", SQLState: " << e.getSQLState() << " ) " << endl;
        }
      }
      else
      {
        subsequentRun = true;
      }
    }
  }

  return 0;
}

