//Run individually for each mode
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <math.h> 
#include <stdio.h>

using namespace std;

int main(int argc, char *argv[])
{
	if (argc != 3)
	{
		cout<<"Two parameters expected"<<endl;
		//cout<<"ANM [PDB file]"<<endl;
		return -1;
	}

	string protomerVT = argv[1]; //Protomer3CG_VT.txt
	string protomerMode = argv[2]; //ProtomerMode617.txt
	
	return -1;

	// change these varaibles as required
	int total =  627;
	//int mode1 =  12758;
	int mode1 =  617;
	//int mode1 =  12756;
	//int mode1 =  12755;
	//int mode1 =  12754;

	int direction = 1; //direction of overlap correction

	double eigenVectors[total/3][3];
	int vectorIndex=0;
	ifstream mfile (protomerVT.c_str()); //CHANGE HERE
	string line;
	string element;
	double ve;

	int countLines =0;
	while (! mfile.eof() )
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
	outputFileW.open(protomerMode.c_str()); //CHANGE HERE

	for (int i = 0; i<total/3; i++)
	{
		for (int j= 0; j<3; j++)
		{
			outputFileW<<eigenVectors[i][j]<<" ";
		}//for j1

		outputFileW<<endl;
	}//for i 

	outputFileW.close();

	return 0;
}//main
