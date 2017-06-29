//Performed Normal Mode Analysis using the ANM tool
//The script uses the alglib library and the following must be included for the script to run. These are contained in src folder in a cpp folder of the alglib package. When complying this code use g++ -I path/cpp/src as an option
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <string>
#include <sstream>
#include <limits>
#include <specialfunctions.h>
#include <specialfunctions.cpp>
#include <iostream>
#include <linalg.h> 
#include <linalg.cpp> 
#include <alglibinternal.h>
#include <alglibinternal.cpp>
#include <alglibmisc.h>
#include <alglibmisc.cpp>
#include <ap.h>
#include <ap.cpp>
#include <vector>



// TO GENERALISE THIS WE MUST JUST TAKE A FILE NAME AS A PARAMETER AND HAVE GENERAL OUTPUT FILE NAMES.
using namespace std;
using namespace alglib;

int countAtoms()//counts Carbon Atoms (Beta carbons for all residues but Alpha carbons for Glycine)
{
vector< vector<int> > row;
}//countAtoms


vector< vector<double> > getCoOrds()// gets the x,y,z cordinates for all the carbon atoms in PDB file. Returns 3 X numAtoms Matrix
{	
	vector< vector<double> > C;
	vector<double> atomC;
	
	ifstream mfile ("3VBSPent4_SCA.pdb");// This is the PDB file that will be coarse grained. 
	string li;
	while (! mfile.eof() )
	{
		getline (mfile,li);
		istringstream iss(li);
		string atom, temp, type,res,x,y,z;
		iss>>atom;
		iss>>temp;
		iss>>type;
		iss>>res;
		iss>>temp;
		iss>>temp;
		iss>>x;
		iss>>y;
		iss>>z;
		if (atom=="ATOM")
		{	
			
			if((res=="GLY"&&type=="CA")||type=="CB")// This selects only CA atoms for CB atoms use if((res=="GLY"&&type=="CA")||type=="CB") such that CA are selected in the case of Glycine
			{
				
				atomC.push_back(atof(x.c_str()));
				atomC.push_back(atof(y.c_str()));
				atomC.push_back(atof(z.c_str()));
				C.push_back(atomC);
				atomC.clear();				
			}//if GLY
			
		
		}//if ATOM
	}//while
        mfile.close();
	return C;
}//cords

vector< vector<double> > getForceConstants(vector<double> atom1,vector<double> atom2)// returns a 9x9 vector for the interaction between two nodes
{
	
	vector< vector<double> > dv2k;
	vector<double> dv2kROW;

	//calculate distance squared
	double dist2 =0.0;
	for (int i =0; i<3; i++)
	{
		dist2 = dist2+((atom2[i]-atom1[i])*(atom2[i]-atom1[i]));
	}//for to instantiate distance
	
	double dv2ka;
	
	//if outside cutoff or atom1=atom2			
	if (dist2>576)//parameter 
	{
		for(int i=0;i<3;i++)
		{		
			dv2kROW.push_back(0.0);
			dv2kROW.push_back(0.0);
			dv2kROW.push_back(0.0);
			dv2k.push_back(dv2kROW);

			//reset dv2kROW
			dv2kROW.clear();
		}//for to populate row by row: times 3 rows
							
	}//if
	else 
	{
		for (int i=0; i<3; i++)
		{
			for (int j=0; j<3; j++)
			{
				dv2ka = -((atom2[i]-atom1[i])*(atom2[j]-atom1[j]))/dist2;
				dv2kROW.push_back(dv2ka);  

			}//for
			dv2k.push_back(dv2kROW);
			dv2kROW.clear();
		}//for
	}//else
	return dv2k;

}// getForceConstants

vector< vector<double> > getHessian(vector< vector<double> > C)// this calls the getForceConstants to set up the Hessian Matrix, we manipulate these force constants and the populate the hessian
{
	vector< vector<double> > Hessian; //Full 3Nx3N Hessian
	
	vector<double> rowX;
	vector<double> rowY;
	vector<double> rowZ;
	
	vector< vector<double> > interaction; //9x9 vector for each atom-atom interaction

	vector< vector<double> > diagonal; //holds 9x9 vector of values on a diagonal
	vector<double> diagonalROW; //holds 9x9 vector of values on a diagonal

	vector<double> atom1;
	vector<double> atom2;

	for (int i=0; i < C.size(); i++)
	{
		
		//reset the diagonal
		for (int d = 0; d<3; d++)
		{
			diagonalROW.push_back(0.0);
			diagonalROW.push_back(0.0);
			diagonalROW.push_back(0.0);
			diagonal.push_back(diagonalROW);
			diagonalROW.clear();
		}// set each element = 0
		
		//get x,y,z of atom1
		atom1 = C[i];

		//Nested loop for atom-atom interactions
		for (int j=0; j < C.size(); j++)
		{
			//get x,y,z of interacting atom			
			atom2 = C[j];
			if(i!=j)
			{
				interaction = getForceConstants(atom1,atom2);
				for(int ir = 0; ir<3; ir++)
				{
					rowX.push_back(interaction[0][ir]);
					rowY.push_back(interaction[1][ir]);
					rowZ.push_back(interaction[2][ir]);
				}//itterate RowX,Y,Z

				//itterate summation of digonals
				for (int ir = 0; ir<3; ir++)
				{
					for (int ic = 0; ic<3; ic++)
					{
						diagonal[ir][ic] = diagonal[ir][ic]-interaction[ir][ic];
								
					}// increase each element
				}// increase each element
				interaction.clear();
			}//if i!=j
			else
			{
				for(int ir = 0; ir<3; ir++)
				{
					rowX.push_back(0.0);
					rowY.push_back(0.0);
					rowZ.push_back(0.0);
				}//itterate RowX,Y,Z
			}//set diagonal
		
		}//for atoms 
		//Update diagonal
		for(int ir = 0; ir<3; ir++)
		{
			rowX[i*3+ir] = diagonal[0][ir];
			rowY[i*3+ir] = diagonal[1][ir];
			rowZ[i*3+ir] = diagonal[2][ir];
		}//itterate RowX,Y,Z
		//Reset diagonal
		diagonal.clear();
		
		//Update Hessian by row
		Hessian.push_back(rowX);
		rowX.clear();
		Hessian.push_back(rowY);
		rowY.clear();
		Hessian.push_back(rowZ);
		rowZ.clear();
		
	  	
	}//for atoms
	

	return Hessian;
}// getHessian



int main()
{	
	vector< vector<double> > C = getCoOrds();
	vector< vector<double> > Hessian = getHessian(C);
	
	
	int size = Hessian.size();
	alglib::real_2d_array Hes;
	Hes.setlength(size,size);
	for(int i =size-1; i>=0; i--)
	{
		for(int j =size-1; j>=0; j--)
		{
		Hes[i][j]= Hessian[i][j];
		Hessian[i].erase(Hessian[i].begin()+j);
		}			
	}
	
	Hessian.clear();
	cout<<"Starting Decomposition"<<endl;
	alglib::real_1d_array w;
	alglib::real_2d_array u;
        alglib::real_2d_array vt;
        alglib::rmatrixsvd(Hes,size,size,2, 2,0,w,u,vt);
	

	ofstream outputFileW;
     

        int r = vt.rows();
        int c =vt.cols();
     

        outputFileW.open("4BIP_W.txt");// this is the eigenvalue matrix.
        for (int i=0; i<r; i++)
        {
		double e = w(i);
        	outputFileW<<i<<" "<<e<<endl;
        }
	outputFileW.close();

        ofstream outputFileVT;
        outputFileVT.open("4BIP_VT.txt");

        for (int i=0; i<r; i++)
        {


                for (int j=0; j<c; j++)
                {
                        double e = vt(i,j); //this is the eigenvector matrix - eigenvectors are rows. 
                        outputFileVT<<e<<" ";
                }//for U

                outputFileVT<<endl;
        }// for vt
        outputFileVT.close();

	ofstream outputFileU;
        outputFileVT.open("4BIP_U.txt");

        for (int i=0; i<r; i++)
        {


                for (int j=0; j<c; j++)
                {
                        double e = u(i,j); //this is the eigenvector matrix - eigenvectors are columns, U and VT are square matrics. 
                        outputFileU<<e<<" ";
                }//for U

                outputFileU<<endl;
        }// for vt
        outputFileU.close();


	return 0;
}//main
