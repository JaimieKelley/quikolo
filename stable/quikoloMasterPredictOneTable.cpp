#include <cmath>
#include <math.h>
#include <vector>
#include <algorithm>
#include <complex>
#include <mysql_connection.h>
#include <cppconn/driver.h>
#include <cppconn/exception.h>
#include <cppconn/resultset.h>
#include <cppconn/statement.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <unistd.h>
#include <errno.h>
#include <netdb.h>
#include <arpa/inet.h>
#include <fcntl.h>
#include <sys/select.h>
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

#define W 120
#define MAX_LINE 16384
#define MAXDATASIZE 100 // max number of bytes we can get at once

void *get_in_addr(struct sockaddr *sa)
{
  if(sa->sa_family == AF_INET) {
    return &(((struct sockaddr_in*)sa)->sin_addr);
  }

  return &(((struct sockaddr_in6*)sa)->sin6_addr);
}

long int initiateClient(char *argv[],int argc, int message)
{
        int sockfd, numbytes;
        char buf[MAXDATASIZE];
        char hostname[30];
        char portNum[6];
        struct addrinfo hints, *servinfo, *p;
        int rv;
        char s[INET6_ADDRSTRLEN];
        strcpy(hostname, "chris09");
        strcpy(portNum, "1072");

        memset(&hints, 0, sizeof hints);
        hints.ai_family = AF_INET;
        hints.ai_socktype = SOCK_STREAM;

        if((rv = getaddrinfo(hostname, portNum, &hints, &servinfo)) != 0) {
                fprintf(stderr, "getaddrinfo: %s\n", gai_strerror(rv));
                return 1;
        }


        // loop through all results and connect to the first
        for(p = servinfo; p != NULL; p = p->ai_next) {
                if((sockfd = socket(p->ai_family, p->ai_socktype, p->ai_protocol)) == -1) {
                        perror("client: socket");
                        continue;
                }

                if (connect(sockfd, p->ai_addr, p->ai_addrlen) == -1) {
                        close(sockfd);
                        perror("client: connect");
                        continue;
                }

                break;
        }

        if(p == NULL) {
                cout << "client: failed to connect\n";
                return 2;
        }

        inet_ntop(p->ai_family, get_in_addr((struct sockaddr *)p->ai_addr), s, sizeof s);
        printf("cli: connecting to %s\n", s);

        freeaddrinfo(servinfo); // all done with this structure

        if(message==1)
        {
		strcpy(buf, "start\n");
                /*sprintf(buf, "start %s\n", argv[1]);
		char * pos = strchr(buf, '\n');
		int number = pos-buf+1;*/
                if((numbytes = send(sockfd, buf, 6, 0)) == -1) {
                        perror("client send failed");
                        exit(1);
                }
                cout << "Message sent\n";
        }
        else
        {
                strcpy(buf, "stop\n");

                if((numbytes = send(sockfd, buf, 5, 0)) == -1) {
                        perror("client send failed");
                        exit(1);
                }
                cout << "Message sent\n";
        }

        if((numbytes = recv(sockfd, buf, MAXDATASIZE-1, 0)) == -1) {
                perror("client recv failed");
                exit(1);
        }
        cout << "Message received\n";

        buf[numbytes] = '\0';

        printf("Client: received '%s'\n", buf);

        close(sockfd);

        return atol(buf);
}

void sendToMYSQL(long int sTime, int mTime, double sloTime)
{
  try{
    sql::Driver *driver;
    sql::Connection *con;
    sql::Statement *stmt;

    /* Create a connection */
    driver = get_driver_instance();
    con = driver->connect("tcp://localhost:3306", "root", "root");

    /* Connect to the MySQL test database */
    con->setSchema("READINGS");
    stmt = con->createStatement();

    string cpS;
    stringstream memSS;
    memSS << "INSERT INTO slo(seconds, milliseconds, sloviolations) VALUES (" << sTime << ", "
	<< mTime << ", " << sloTime << ")";
    cpS = memSS.str();
    stmt->executeUpdate(cpS);

    delete stmt;
    delete con;
  } catch (sql::SQLException &e) {
    cout << "#ERR: SQLException in " << __FILE__;
    cout << "(" << __FUNCTION__ << ") on line " << __LINE__ << endl;
    cout << "#ERR: " << e.what();
    cout << " (MySQL error code: " << e.getErrorCode();
    cout << ", SQLState: " << e.getSQLState() << " ) " << endl;
  }

  cout << "end of sendToMYSQL" << endl;
}

int loadSLOs(int windowSTime, long int seconds[], int milliseconds[], double trace[])
{
  int count = 0;
  int total = 0;
  int secount = 0;

  try{
    sql::Driver *driver;
    sql::Connection *con;
    sql::Statement *stmt;
    sql::ResultSet *res;

    /* Create a connection */
    driver = get_driver_instance();
    con = driver->connect("tcp://localhost:3306", "root", "root");

    /* Connect to the MySQL test database */
    con->setSchema("READINGS");
    stmt = con->createStatement();

    string cpS;
    stringstream memSS;
    memSS << "SELECT seconds, milliseconds, sloviolations FROM slo WHERE seconds>=" << windowSTime << " ORDER BY seconds, milliseconds ASC";
    cpS = memSS.str();
    res = stmt->executeQuery(cpS);

    int win = W;
    count = 0;

    vector<long int> se;
    vector<int> me;
    vector<double> te;
    while (res->next()) {
        se.push_back(atol((res->getString(1)).c_str()));
        me.push_back(res->getInt(2));
        te.push_back(atof((res->getString(3)).c_str()));
	count++;
    }

    if(count < win) win = count;

    for(int i = 0 + (count - win); i < count; i++)
    {
        seconds[secount] = se.at(i);
        milliseconds[secount]=me.at(i);
        trace[secount] = te.at(i);
	secount++;
    }

    delete res;
    delete stmt;
    delete con;
  } catch (sql::SQLException &e) {
    cout << "#ERR: SQLException in " << __FILE__;
    cout << "(" << __FUNCTION__ << ") on line " << __LINE__ << endl;
    cout << "#ERR: " << e.what();
    cout << " (MySQL error code: " << e.getErrorCode();
    cout << ", SQLState: " << e.getSQLState() << " ) " << endl;
  }

  cout << "end of loadSLO" << endl;

  return secount;
}

void writeToFile(string fNameStr, long int seconds, int milliseconds, double traceVal)
{
  fstream myOutStream(fNameStr, fstream::out | fstream::app);

  if(myOutStream)
  {
    myOutStream << seconds << " " << milliseconds << " " << traceVal << endl;
  }

  myOutStream.close();
}

int loadArray(int windowSTime, long int seconds[], int milliseconds[], double (&trace)[39][W])
{
  int count = 0;
  int total = 0;
  int secount = 0;

  try{
    sql::Driver *driver;
    sql::Connection *con;
    sql::Statement *stmt;
    sql::ResultSet *res;

    /* Create a connection */
    driver = get_driver_instance();
    con = driver->connect("tcp://localhost:3306", "root", "root");

    /* Connect to the MySQL test database */
    con->setSchema("READINGS");
    stmt = con->createStatement();

    string cpS;
    stringstream memSS;
    memSS << "SELECT seconds, milliseconds, cpu0, cpu1, cpu2, cpu3, cpu4, cpu5, cpu6, cpu7, cpu8, cpu9, cpu10, cpu11, cpu12, cpu13, cpu14, cpu15, cpu0l2, cpu1l2, cpu2l2, cpu3l2, cpu4l2, cpu5l2, cpu6l2, cpu7l2, cpu8l2, cpu9l2, cpu10l2, cpu11l2, cpu12l2, cpu13l2, cpu14l2, cpu15l2, socket0l3, socket1l3, mem, received, transmitted, diskreads, diskwrites FROM hardware WHERE seconds>=" << windowSTime << " ORDER BY seconds, milliseconds ASC";
    cpS = memSS.str();
    res = stmt->executeQuery(cpS);

    int win = W;
    count = 0;

    vector<long int> se;
    vector<int> me;
    vector<double> values;
    while (res->next()) {
	se.push_back(atol((res->getString(1)).c_str()));
	me.push_back(res->getInt(2));
	for(int i = 3; i <= 41; i++)
	{
		values.push_back(atof((res->getString(i)).c_str()));
	}
	//if(count < 10) cout << count << " " << se.at(count) << " " << me.at(count) << " " << te.at(count) << endl;
	count++;
    }
    cout << count << " " << values.size() << endl;

    if(count < win) win = count;

    for(int i = 0 + (count - win); i < count; i++)
    {
	seconds[secount] = se.at(i);
        milliseconds[secount]=me.at(i);
	for(int j = 0; j < 39; j++)
	{
		int base = 39 * i;
		//cout << j << " " << secount << " " << i << " " << base << endl;
		trace[j][secount] = values.at(base + j);
        }
	//if(secount < 10) cout << secount << " " << seconds[secount] << " " << milliseconds[secount] << " " <<trace[secount] << endl;
        secount++;
    }

    delete res;
    delete stmt;
    delete con;
  } catch (sql::SQLException &e) {
    cout << "#ERR: SQLException in " << __FILE__;
    cout << "(" << __FUNCTION__ << ") on line " << __LINE__ << endl;
    cout << "#ERR: " << e.what();
    cout << " (MySQL error code: " << e.getErrorCode();
    cout << ", SQLState: " << e.getSQLState() << " ) " << endl;
  }

  cout << "end of loadArray" << endl;

  return secount;
}

template<class T>
class matrix
{
public:
	matrix(unsigned int nRows, unsigned int nCols) : 
		m_nRows( nRows ), 
		m_nCols( nCols ), 
		m_oData( nRows*nCols, 0 )
	{
		if ( !nRows || !nCols )
		{
			cout << "invalid matrix size\n";
		}
	}
 
	static matrix identity( unsigned int nSize )
	{
		matrix oResult( nSize, nSize );
 
		int nCount = 0;
		for(std::vector<double>::iterator it=oResult.m_oData.begin(); it!=oResult.m_oData.end(); ++it)
			*it = !(nCount++%(nSize + 1));
 
		return oResult;
	}
 
	inline T& operator()(unsigned int nRow, unsigned int nCol)
	{
		if ( nRow >= m_nRows || nCol >= m_nCols )
		{
			cout <<"position out of range\n";
		}
 
		return m_oData[nCol+m_nCols*nRow];
	}
 
	inline matrix operator*(matrix& other)
	{
		if ( m_nCols != other.m_nRows )
		{
			cout << "matrix dimensions are not multiplicable\n";
		}
 
		matrix oResult( m_nRows, other.m_nCols );
		for ( unsigned int r = 0; r < m_nRows; ++r )
		{
			for ( unsigned int ocol = 0; ocol < other.m_nCols; ++ocol )
			{
				for ( unsigned int c = 0; c < m_nCols; ++c )
				{
					oResult(r,ocol) += (*this)(r,c) * other(c,ocol);
				}
			}
		}
 
		return oResult;
	}
 
	inline matrix transpose()
	{
		matrix oResult( m_nCols, m_nRows );
		for ( unsigned int r = 0; r < m_nRows; ++r )
		{
			for ( unsigned int c = 0; c < m_nCols; ++c )
			{
				oResult(c,r) += (*this)(r,c);
			}
		}
		return oResult;
	}
 
	inline unsigned int rows() 
	{
		return m_nRows;
	}
 
	inline unsigned int cols() 
	{
		return m_nCols;
	}
 
	inline vector<T> data()
	{
		return m_oData;
	}
 
	void print()
	{
		for ( unsigned int r = 0; r < m_nRows; r++ )
		{
			for ( unsigned int c = 0; c < m_nCols; c++ )
			{
				std::cout << (*this)(r,c) << "\t";
			}
			std::cout << std::endl;
		}
	}
 
private:
	std::vector<T> m_oData;
 
	unsigned int m_nRows;
	unsigned int m_nCols;
};

template<typename T>
	class Givens
	{
	public:
		Givens() : m_oJ(2,2), m_oQ(1,1), m_oR(1,1)
		{
		}

		/*
			Calculate the inverse of a matrix using the QR decomposition.

			param:
				A	matrix to inverse
		*/
		const matrix<T> Inverse( matrix<T>& oMatrix )
		{
			if ( oMatrix.cols() != oMatrix.rows() )
			{
				cout << "matrix has to be square\n";
			}
			matrix<T> oIdentity = matrix<T>::identity( oMatrix.rows() );
			Decompose( oMatrix );
			return Solve( oIdentity );
		}

		/*
			Performs QR factorization using Givens rotations.
		*/
		void Decompose( matrix<T>& oMatrix )
		{
			int nRows = oMatrix.rows();
			int nCols = oMatrix.cols();


			if ( nRows == nCols )
			{
				nCols--;
			}
			else if ( nRows < nCols )
			{
				nCols = nRows - 1;
			}

			m_oQ = matrix<T>::identity(nRows);
			m_oR = oMatrix;

			for ( int j = 0; j < nCols; j++ )
			{
				for ( int i = j + 1; i < nRows; i++ )
				{
					GivensRotation( m_oR(j,j), m_oR(i,j) );
					PreMultiplyGivens( m_oR, j, i );
					PreMultiplyGivens( m_oQ, j, i );
				}
			}

			m_oQ = m_oQ.transpose();
		}
		
		/*
			Find the solution for a matrix.
			http://en.wikipedia.org/wiki/QR_decomposition#Using_for_solution_to_linear_inverse_problems
		*/
		matrix<T> Solve( matrix<T>& oMatrix )
		{
			matrix<T> oQtM( m_oQ.transpose() * oMatrix );
			int nCols = m_oR.cols();
			matrix<T> oS( 1, nCols );
			for (int i = nCols-1; i >= 0; i-- )
			{
				oS(0,i) = oQtM(i, 0);
				for ( int j = i + 1; j < nCols; j++ )
				{
					oS(0,i) -= oS(0,j) * m_oR(i, j);
				}
				oS(0,i) /= m_oR(i, i);
			}

			return oS;
		}

		const matrix<T>& GetQ()
		{
			return m_oQ;
		}

		const matrix<T>& GetR()
		{
			return m_oR;
		}

	private:
		/*
			Givens rotation is a rotation in the plane spanned by two coordinates axes.
			http://en.wikipedia.org/wiki/Givens_rotation
		*/
		void GivensRotation( T a, T b )
		{
			T t,s,c;
			if (b == 0)
			{
				c = (a >=0)?1:-1;
				s = 0; 
			}
			else if (a == 0)
			{
				c = 0;
				s = (b >=0)?-1:1;
			}
			else if (abs(b) > abs(a))
			{
				t = a/b;
				s = -1/sqrt(1+t*t);
				c = -s*t;
			}
			else
			{
				t = b/a;
				c = 1/sqrt(1+t*t);
				s = -c*t;
			}
			m_oJ(0,0) = c; m_oJ(0,1) = -s;
			m_oJ(1,0) = s; m_oJ(1,1) = c;
		}

		/*
			Get the premultiplication of a given matrix 
			by the Givens rotation.
		*/
		void PreMultiplyGivens( matrix<T>& oMatrix, int i, int j )
		{
			int nRowSize = oMatrix.cols();

			for ( int nRow = 0; nRow < nRowSize; nRow++ )
			{
				double nTemp = oMatrix(i,nRow) * m_oJ(0,0) + oMatrix(j,nRow) * m_oJ(0,1);
				oMatrix(j,nRow) = oMatrix(i,nRow) * m_oJ(1,0) + oMatrix(j,nRow) * m_oJ(1,1);
				oMatrix(i,nRow) = nTemp;
			}
		}

	private:
		matrix<T> m_oQ, m_oR, m_oJ;
	};

template<typename T>
std::vector<T> polyfit( const std::vector<T>& oX, const std::vector<T>& oY, int nDegree )
{
	if ( oX.size() != oY.size() )
		cout <<"X and Y vector sizes do not match\n";
 
	// more intuative this way
	nDegree++;
	
	size_t nCount =  oX.size();
	matrix<T> oXMatrix( nCount, nDegree );
	matrix<T> oYMatrix( nCount, 1 );
	
	// copy y matrix
	for ( size_t i = 0; i < nCount; i++ )
	{
		oYMatrix(i, 0) = oY[i];
	}
 
	// create the X matrix
	for ( size_t nRow = 0; nRow < nCount; nRow++ )
	{
		T nVal = 1.0f;
		for ( int nCol = 0; nCol < nDegree; nCol++ )
		{
			oXMatrix(nRow, nCol) = nVal;
			nVal *= oX[nRow];
		}
	}
 
	// transpose X matrix
	matrix<T> oXtMatrix( oXMatrix.transpose() );
	// multiply transposed X matrix with X matrix
	matrix<T> oXtXMatrix( oXtMatrix * oXMatrix );
	// multiply transposed X matrix with Y matrix
	matrix<T> oXtYMatrix( oXtMatrix * oYMatrix );
 
	Givens<T> oGivens;
	oGivens.Decompose( oXtXMatrix );
	matrix<T> oCoeff = oGivens.Solve( oXtYMatrix );
	// copy the result to coeff
	return oCoeff.data();
}

vector<complex<double> > fft(vector<double> inputData)
{
    int count = 0;
    unsigned long nn = inputData.size();
    double data[2*nn];
    int ni = 0;
    for(vector<double>::iterator di=inputData.begin(); di < inputData.end(); ++di)
    {
	data[ni] = *di;
	data[ni+1] = 0.0;
	ni+= 2;
    }

    unsigned long n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;
 
    // reverse-binary reindexing
    n = nn<<1;
    j=1;
    for (i=1; i<n; i+=2) {
        if (j>i){ // && j < n) {
            swap(data[j-1], data[i-1]);
            swap(data[j], data[i]);
        }
        m = nn;
	count++;
        while (m>=2 && j>m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    };

    // here begins the Danielson-Lanczos section
    mmax=2;
    while (n>mmax) {
        istep = mmax<<1;
        theta = -(2*M_PI/mmax);
        wtemp = sin(0.5*theta);
        wpr = -2.0*wtemp*wtemp;
        wpi = sin(theta);
        wr = 1.0;
        wi = 0.0;
        for (m=1; m < mmax; m += 2) {
            for (i=m; i <= n; i += istep) {
                j=i+mmax;
		if(j < n)
		{
               	  tempr = wr*data[j-1] - wi*data[j];
                  tempi = wr * data[j] + wi*data[j-1];
 
                  data[j-1] = data[i-1] - tempr;
                  data[j] = data[i] - tempi;
                  data[i-1] += tempr;
                  data[i] += tempi;
		}
            }
            wtemp=wr;
            wr += wr*wpr - wi*wpi;
            wi += wi*wpr + wtemp*wpi;
        }
        mmax=istep;
    }

    vector<complex<double> > outData;
    for(i = 1; i < n; i+=2)
    {
	outData.push_back(complex<double>(data[i-1],data[i]));
    }

    return outData;
}

vector<double> fftfreq(int n, double d)
{
  vector<double> ff;
  double f[n];
  f[0] = 0;
  for(int i = 1; i <= (n/2); i++)
  {
	f[i] = i / (d*n);
	f[n-i] = - i /(d*n);
  }

  for(int i =0; i< n; i++)
  {
    ff.push_back(f[i]);
  }

  return ff;
}

struct MyComparator
{
  const vector<double> & value_vector;

  MyComparator(const vector<double> & val_vec):
	value_vector(val_vec) {}

  bool operator()(int i1, int i2)
  {
	return (abs(value_vector[i1]) < abs(value_vector[i2]));
  }
};

vector<double> fftpredict(vector<double> x, int harmonies, int toPredict) {
  const unsigned long cn = x.size();
  int n = x.size();  // number in initial array
  int d = toPredict;	// number of predictions to make
  int num_Harmonies = harmonies;
  vector<double> indexes;
  vector<double> t_vec;
  vector<double> x_vec;
  double x_avg = 0.0;
  double pred_avg = 0.0;
  for(int i = 0; i < n; i++)
  {
	double t;
	t = (double)i;
	t_vec.push_back(t);
	x_vec.push_back(x[i]);
	indexes.push_back(t);
	x_avg += x[i];
  }
  x_avg = x_avg / n;
  vector<double> p = polyfit(t_vec, x_vec, 1);
  for(vector<double>::iterator it=p.begin(); it < p.end(); ++it)
  {
	cout << "P " << *it << "\n";
  }

  vector<double> x_notrend;
  for(int i = 0; i < n; i++)
  {
	double temp = x_vec[i] - p[1] * t_vec[i];
	x_notrend.push_back(temp);
  }

  vector<complex<double> > x_fft;
  x_fft = fft(x_notrend);

  vector<double> f = fftfreq(n, 1.0);
  std::sort(indexes.begin(), indexes.end(), MyComparator(f));

  int td[n+d];
  for(int i = 0; i < n + d; i++)
  {
	td[i] = i;
  }

  double restored_sig[n+d];
  for(int i = 0; i < n+d; i++)
  {
	restored_sig[i] = 0.0;
  }

  int max_Harm = 1+ 2*num_Harmonies;
  int ind = 0;
  for(vector<double>::iterator it=indexes.begin(); it < indexes.end(); ++it)
  {
	if(ind < max_Harm)
	{
		int ii = *it;
		double amplitude = abs(x_fft[ii]) / n;
		double phase = std::arg(x_fft[ii]);
		for(int j = 0; j < n+d; j++)
		{
			restored_sig[j] += amplitude * cos(2*M_PI*f[ii]*td[j] + phase);
		}
	}
	ind++;
  }

  double answer[n+d];
  vector<double> prediction;
  for(int i = 0; i < n+d; i++)
  {
	answer[i] = restored_sig[i] + p[1]*td[i];
	pred_avg += answer[i];
  }
  pred_avg = pred_avg / (n+d);

  for(int i = n; i < n+d; i++)
  {
	prediction.push_back(abs(answer[i]));
  }

  return prediction;
}


int main(int argc, char *argv[]) {
  string line;
  int numCPU = 16;
  int numSockets = 2;
  timeval start_time;
  long int start_seconds = 0;
  long int sSeconds[W];
  int sMilliseconds[W];
  double sloMetric[W];
  int start_millisec = 0;
  int halfW = W/2;  // window size 6000
  int ww = W;
  int D = 10; // number to predict
  double traceMetric[39][W];
  long int tseconds[W];
  int tmilliseconds[W];
  long int stime = 0L;
  long int rtime = 0L;
  vector<string> columns;
  size_t fSize;
  char command[100];
  bzero(command, 100);
  strcpy(command, "./universalSLO.sh ");
  strcat(command, argv[1]);
  strcat(command, " ");
  strcat(command, argv[2]);
  bool appdone = false;

  cout << "docker container: " << argv[1] << endl;
  cout << "args: " << argv[2] << endl;

  string filename = argv[1];
  char * founds = strstr(argv[1], "/");
  if(founds == NULL){
	founds = (char *)&argv[1][0];
  }
  else
  {
	founds += 1;
  }

  fSize = strlen(founds);

  char fileName[40];
  bzero(fileName, 40);
  strcpy(fileName, "traces/");
  strncat(fileName, founds, fSize);
  strcat(fileName, argv[2]);

  cout << "fileName " << fileName << endl;

  for(int i = 0; i < numCPU; i++)
  {
    stringstream cpSS;
    cpSS << "cpu" << i;
    string cpS = cpSS.str();
    columns.push_back(cpS);

    stringstream l2SS;
    l2SS << "cpu" << i << "l2";
    cpS = l2SS.str();
    columns.push_back(cpS);
  }

  for(int i = 0; i < numSockets; i++)
  {
    stringstream l3SS;
    l3SS << "socket" << i << "l3";
    string cpS = l3SS.str();
    columns.push_back(cpS);
  }

  columns.push_back("mem");
  columns.push_back("received");
  columns.push_back("transmitted");
  columns.push_back("diskreads");
  columns.push_back("diskwrites");

  stime = initiateClient(argv, argc, 1);
  start_seconds = stime / 1000;
  start_millisec = stime - (start_seconds*1000);
  long int slo_start = start_seconds;
  int slo_mill = start_millisec;

  cout << "start " << stime << " " << start_seconds << " " << start_millisec << "\n";

  stringstream fName;
  fName << fileName << ".txt";
  string fNameStr = fName.str();

  #pragma omp parallel sections
  {
    {
      FILE *pin;
      char buff[1024];

      if(!(pin = popen(command, "r"))){
	exit;
      }

      while(fgets(buff, sizeof(buff), pin) != NULL){
	char * rtimes;
	if((rtimes=strstr(buff,"runtime"))!= NULL)
	{
		cout << buff;
		long int bsec;
		int mbsec;
		double totaltime;

		bsec = atol(buff);
		totaltime = (double)bsec / 1000.0;
		rtimes = strstr(buff, " ");
		bsec = atol(rtimes);
		mbsec = bsec % 1000;
		bsec = bsec / 1000;
		cout << " " << bsec << " " << mbsec << " " << totaltime << endl;
		sendToMYSQL(bsec, mbsec, totaltime);
		writeToFile(fNameStr, bsec, mbsec, totaltime);
	}
      }
      pclose(pin);
      appdone = true;
    }
    #pragma omp section
    {
      sleep(halfW);

      for(int i = 0; i < W; i++)
      {
	for(int j = 0; j < 39; j++)
	{
          traceMetric[j][i] = 0.0;
	}
	sSeconds[i] = 0;
	sMilliseconds[i] = 0;
	sloMetric[i] = 0.0;
        tseconds[i] = 0;
        tmilliseconds[i] = 0;
      }

      while(!appdone)
      {
	struct timespec sltime;
	sltime.tv_sec = 0;
	sltime.tv_nsec = 500000000;
	nanosleep(&sltime, NULL);

	for(int i = 0; i < W; i++)
	{
	  tseconds[i] = 0;
	  tmilliseconds[i] = 0;
	  sSeconds[i] = 0;
	  sMilliseconds[i] = 0;
	  sloMetric[i] = 0.0;
	}

	for(int j = 0; j < 39; j++)
        {
	  for(int i = 0; i < W; i++)
          {
                traceMetric[j][i] = 0.0;
          }
	}

	cout << "Before load array" << endl;
        ww = loadArray(start_seconds, tseconds, tmilliseconds, traceMetric);
        cout << "After load array, before load slos" << endl;

	int vv = loadSLOs(slo_start, sSeconds, sMilliseconds, sloMetric);
	cout << "After load slos" << endl;

	for(int j = 0; j < 39; j++)
	{
          cout << columns.at(j) << endl;

          // do stuff with trace.
          vector<double> givenVals;
          for(int i = 0; i < ww; i++)
          {
                givenVals.push_back(traceMetric[j][i]);
                //cout << "sec " << tseconds[i] << ", msec " << tmilliseconds[i] << ", tval " << traceMetric[j][i] << endl;
          }

	  cout << "Predicted: " << endl;
          vector<double> predictedVals = fftpredict(givenVals, 10, D);
	  for(int i = 0; i < D; i++)
	  {
		cout << predictedVals[i] << ", ";
 	  }
	  cout << endl;
        }  // end for

	cout << "sloViolations" << endl;
	vector<double> givenVals;
        for(int i = 0; i < vv; i++)
        {
          givenVals.push_back(sloMetric[i]);
          cout << "sec " << sSeconds[i] << ", msec " << sMilliseconds[i] << ", tval " << sloMetric[i] << endl;
        }

        cout << "Predicted: " << endl;
        vector<double> predictedVals = fftpredict(givenVals, 10, D);
        for(int i = 0; i < D; i++)
        {
          cout << predictedVals[i] << ", ";
        }
        cout << endl;


	start_millisec = tmilliseconds[0];
	start_seconds = tseconds[0];
      } // end while
    } // end section
  }  // end parallel

  long throwaway = initiateClient(argv, argc, 0);

  return 0;
}

