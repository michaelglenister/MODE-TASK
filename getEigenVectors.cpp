//Run individually for each mode
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <math.h> 
#include <stdio.h>

#include <ap.h>
#include <time.h>

using namespace std;

// Get current date/time, format is YYYY-MM-DD HH:mm:ss
const std::string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);

    strftime(buf, sizeof(buf), "%Y-%m-%d %X", &tstruct);

    return buf;
}

int main(int argc, char *argv[])
{
	//Init vars
	string protomerVT, outdir = "output";
	bool hasVt = false;

	// Begin parameter handling
	// Add more else if statements for further parameters
	int i;
    for(i=0; i<argc; ++i)
    {
		if (strcmp(argv[i], "-h") == 0)
		{
			cout<<"Help"<<endl;
			return -1;
		}
		else if(strcmp(argv[i], "--vt_values") == 0)
		{
			protomerVT = argv[i+1];
			hasVt = true;
		}
		else if(strcmp(argv[i], "--outdir") == 0)
		{
			outdir = atof(argv[i+1]);
			//hasOutdir = true;
		}
    }

	if(!hasVt)
	{
		cout<<"A VT matrix file is required, use '-h' to view help"<<endl;
		return -1;
	}

	//string protomerVT = argv[1]; //Protomer3CG_VT.txt
	string protomerMode = outdir + "/ProtomerMode.txt"; //argv[2]; //ProtomerMode617.txt
	
	// change these varaibles as required
	int total = 627;
	//int mode1 = 12758;
	int mode1 = 617;
	//int mode1 = 12756;
	//int mode1 = 12755;
	//int mode1 = 12754;

	int direction = 1; //direction of overlap correction

	// End parameter handling

	// Start cronometer
	const int ONE_HOUR = 60 * 60;
	const int ONE_MINUE = 60;

	int hour;
	int min;
	int sec;
	std::cout << "Started at: " << currentDateTime() << std::endl;
	clock_t tStart = clock();

	
	double eigenVectors[total/3][3];
	int vectorIndex = 0;
	ifstream mfile (protomerVT.c_str());
	string line;
	string element;
	double ve;

	int countLines = 0;

	while (!mfile.eof())
	{
		getline (mfile,line);

		if(!line.empty())
		{
			countLines++;
		}
		if (countLines-1==mode1)
		{
			istringstream iss(line);
			for (int i = 0; i<total/3; i++)
			{
				for (int j= 0; j<3; j++)
				{
					iss>>element;
					ve = atof(element.c_str());
					eigenVectors[i][j]=ve;
				}//for
			}//for
			break;
		}//if
	}//while
	mfile.close();

	//convert all vectors to unit vectors
	double magnitude;

	for (int i = 0; i<total/3; i++)
	{
		magnitude=0;
		for (int j= 0; j<3; j++)
		{
			magnitude = magnitude+(eigenVectors[i][j]*eigenVectors[i][j]);
		}//for j1

		cout<<magnitude<<endl;
		magnitude = sqrt(magnitude);

		for (int j= 0; j<3; j++)
		{
			eigenVectors[i][j]=direction*eigenVectors[i][j]/magnitude;
		}//for j2
	
	}//for i

	// write ModeVectors to file. 
	ofstream outputFileW;
	outputFileW.open(protomerMode.c_str());

	for (int i = 0; i<total/3; i++)
	{
		for (int j= 0; j<3; j++)
		{
			outputFileW<<eigenVectors[i][j]<<" ";
		}//for j1

		outputFileW<<endl;
	}//for i 

	outputFileW.close();

	// End cronometer
	std::cout << "Completed at: " << currentDateTime() << std::endl;
	int time_target=(clock() - tStart)/CLOCKS_PER_SEC;

	hour=time_target/ONE_HOUR;
	time_target-=hour*ONE_HOUR;
	min=time_target/ONE_MINUE;
	time_target-=min*ONE_HOUR;
	sec=time_target;
	if (min<10 && sec<10)
	{
		printf("- Total time: %d:0%d:0%d\n",hour,min,sec);
	}
	else if (min<10)
	{
		printf("- Total time: %d:0%d:%d\n",hour,min,sec);
	}
	else if (sec<10)
	{
		printf("- Total time: %d:%d:0%d\n",hour,min,sec);
	}
	else
	{
		printf("- Total time: %d:%d:%d\n",hour,min,sec);
	}

	return 0;
}//main
