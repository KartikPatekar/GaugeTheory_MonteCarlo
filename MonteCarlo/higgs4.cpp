//Include Headers
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <complex.h>
#include <map>
#include <ctime>

using namespace std;
using namespace std::chrono;

#define I1 complex <double> (0., 1.) // imaginary unit
#define PI 3.141592653589793238462643383279502884197169399375105820974944

clock_t starttime=0, endtime=0;

int coolingsweeps=10000;
int usefulsweeps=10000;

//GlobalVariables
long double s=1.0, Kg2=1.0, u0=1.0, u1=1.0;                           //Model Parameters
int L=8, D=3;                         // number of spins in x and y and t (imaginary time)
int x,y,t, mu;                               // loop variables
int steps = 0;							//MC steps counter
double Hvev=1;
vector< vector< vector< vector< vector< double > > > > > H;
vector<  double  > Hret{ 0.0, 0.0, 0.0 };				//For returning 3D Higgs field
vector<  double  > H_old{ 0.0, 0.0, 0.0 };				//For storing old H field
vector<  double  > Hnew{ 0.0, 0.0, 0.0 };				//For new flipped H field


//Random Real Number Function
double randomnumber;
inline double std_rand()
{
    randomnumber=rand() / (RAND_MAX + 1.0);
	return randomnumber;
}

//Dot Product
inline double dotProduct(vector<double> v1, vector<double> v2 )
{
	return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

inline double expmindelS(double S_old, double S_new)
{
	return exp(-(S_new-S_old));
}


//Higgs Amplitude Flip
double randAmp;
inline vector< double > FlipHiggsAmp(vector < double > H, double dA=0.1*Hvev)
{
	double Hamp=sqrt(H[0]*H[0]+H[1]*H[1]+H[2]*H[2]);
	randAmp=Hamp+1.0*dA*(1.0-2.0*std_rand());

	if(Hamp==0.0)
	{	
		Hret[0]=randAmp;
		Hret[1]=randAmp;
		Hret[2]=randAmp;
		return Hret;		
	}

	Hret[0]=randAmp/Hamp*H[0];
	Hret[1]=randAmp/Hamp*H[1];
	Hret[2]=randAmp/Hamp*H[2];

	return Hret;
}

//-----------Higgs rotation updates--------
double randtheta, costheta, sintheta;
inline vector< double > FlipHiggsAnglexy(vector < double > H, double dangle=0.1)
{
	randtheta=dangle*2*PI*std_rand();
	costheta=cos(randtheta);
	sintheta=sin(randtheta);

	Hret[0]=costheta*H[0]-sintheta*H[1];
	Hret[1]=costheta*H[1]+sintheta*H[0];
	Hret[2]=H[2];

	return Hret;
}

inline vector< double > FlipHiggsAnglext(vector < double > H, double dangle=0.1)
{
	randtheta=dangle*2*PI*std_rand();
	costheta=cos(randtheta);
	sintheta=sin(randtheta);

	Hret[0]=costheta*H[0]-sintheta*H[2];
	Hret[2]=costheta*H[2]+sintheta*H[0];
	Hret[1]=H[1];

	return Hret;
}

inline vector< double > FlipHiggsAngleyt(vector < double > H, double dangle=0.1)
{
	randtheta=dangle*2*PI*std_rand();
	costheta=cos(randtheta);
	sintheta=sin(randtheta);

	Hret[1]=costheta*H[1]-sintheta*H[2];
	Hret[2]=costheta*H[2]+sintheta*H[1];
	Hret[0]=H[0];

	return Hret;
}
//--------------------------------------------

//Inititalise Higgs Configurations
void initialiseHiggs()
{
	H.resize(5);
	H[1].resize(L);
	H[2].resize(L);
	H[3].resize(L);
	H[4].resize(L);

	for(x=0; x<L;x++)
	{
		H[1][x].resize(L);
		H[2][x].resize(L);
		H[3][x].resize(L);
		H[4][x].resize(L);

		for(y=0; y<L;y++)
		{
			H[1][x][y].resize(L);
			H[2][x][y].resize(L);
			H[3][x][y].resize(L);
			H[4][x][y].resize(L);

			for(t=0; t<L;t++)
			{
				H[1][x][y][t].resize(D);
				H[2][x][y][t].resize(D);
				H[3][x][y][t].resize(D);
				H[4][x][y][t].resize(D);

				H[1][x][y][t]=FlipHiggsAmp(H[1][x][y][t], 2.0*Hvev);
				H[1][x][y][t]=FlipHiggsAnglexy(H[1][x][y][t], 1.0);
				H[1][x][y][t]=FlipHiggsAnglext(H[1][x][y][t], 1.0);
				H[1][x][y][t]=FlipHiggsAngleyt(H[1][x][y][t], 1.0);				

				H[2][x][y][t]=FlipHiggsAmp(H[2][x][y][t], 2.0*Hvev);
				H[2][x][y][t]=FlipHiggsAnglext(H[2][x][y][t], 1.0);
				H[2][x][y][t]=FlipHiggsAngleyt(H[2][x][y][t], 1.0);				
				H[2][x][y][t]=FlipHiggsAnglexy(H[2][x][y][t], 1.0);

				H[3][x][y][t]=FlipHiggsAmp(H[3][x][y][t], 2.0*Hvev);
				H[3][x][y][t]=FlipHiggsAnglexy(H[3][x][y][t], 1.0);
				H[3][x][y][t]=FlipHiggsAnglext(H[3][x][y][t], 1.0);
				H[3][x][y][t]=FlipHiggsAngleyt(H[3][x][y][t], 1.0);				

				H[4][x][y][t]=FlipHiggsAmp(H[4][x][y][t], 2.0*Hvev);
				H[4][x][y][t]=FlipHiggsAnglexy(H[4][x][y][t], 1.0);
				H[4][x][y][t]=FlipHiggsAnglext(H[4][x][y][t], 1.0);
				H[4][x][y][t]=FlipHiggsAngleyt(H[4][x][y][t], 1.0);				
			}
		}
	}
}

// Calculate Action (Only relevant Part)
vector <int> r{0,0,0}, rplx{0,0,0}, rply{0,0,0}, rplt{0,0,0}, rminx{0,0,0}, rminy{0,0,0}, rmint{0,0,0};
double Seff, S_site, S_derivative;
double calcAction(int xcord, int ycord, int tcord, vector<double> Hfield, int Hnum)
{
	r={xcord, ycord, tcord};
	rplx={(xcord+1)%L, ycord, tcord};
	rply={xcord, (ycord+1)%L, tcord};
	rplt={ xcord, ycord, (tcord+1)%L};
	rminx={(xcord+L-1)%L, ycord, tcord};
	rminy={xcord, (ycord+L-1)%L, tcord};
	rmint={ xcord, ycord, (tcord+L-1)%L};

	H_old=H[Hnum][xcord][ycord][tcord];
	H[Hnum][xcord][ycord][tcord]=Hfield;

	Seff=S_site=S_derivative=0.0;

	for(int l=1; l<=4; l++)
	{
		S_site=S_site+s*dotProduct(H[l][xcord][ycord][tcord],H[l][xcord][ycord][tcord]);

		for(int m=1; m<=4; m++)
		{
			S_site=S_site+u0*dotProduct(H[l][xcord][ycord][tcord],H[l][xcord][ycord][tcord])*dotProduct(H[m][xcord][ycord][tcord],H[m][xcord][ycord][tcord]);
			S_site=S_site+u0*dotProduct(H[l][xcord][ycord][tcord],H[m][xcord][ycord][tcord])*dotProduct(H[l][xcord][ycord][tcord],H[m][xcord][ycord][tcord]);
		
			S_derivative=S_derivative-Kg2*dotProduct(H[l][rplx[0]][rplx[1]][rplx[2]],H[m][rplx[0]][rplx[1]][rplx[2]])*dotProduct(H[l][r[0]][r[1]][r[2]],H[m][r[0]][r[1]][r[2]]);
			S_derivative=S_derivative-Kg2*dotProduct(H[l][rply[0]][rply[1]][rply[2]],H[m][rply[0]][rply[1]][rply[2]])*dotProduct(H[l][r[0]][r[1]][r[2]],H[m][r[0]][r[1]][r[2]]);
			S_derivative=S_derivative-Kg2*dotProduct(H[l][rplt[0]][rplt[1]][rplt[2]],H[m][rplt[0]][rplt[1]][rplt[2]])*dotProduct(H[l][r[0]][r[1]][r[2]],H[m][r[0]][r[1]][r[2]]);
			S_derivative=S_derivative-Kg2*dotProduct(H[l][rminx[0]][rminx[1]][rminx[2]],H[m][rminx[0]][rminx[1]][rminx[2]])*dotProduct(H[l][r[0]][r[1]][r[2]],H[m][r[0]][r[1]][r[2]]);
			S_derivative=S_derivative-Kg2*dotProduct(H[l][rminy[0]][rminy[1]][rminy[2]],H[m][rminy[0]][rminy[1]][rminy[2]])*dotProduct(H[l][r[0]][r[1]][r[2]],H[m][r[0]][r[1]][r[2]]);
			S_derivative=S_derivative-Kg2*dotProduct(H[l][rmint[0]][rmint[1]][rmint[2]],H[m][rmint[0]][rmint[1]][rmint[2]])*dotProduct(H[l][r[0]][r[1]][r[2]],H[m][r[0]][r[1]][r[2]]);
		}
	}

	H[Hnum][xcord][ycord][tcord]=H_old;
	Seff=S_site+S_derivative;
	
	return Seff;
}

double calcAction(int xcord, int ycord, int tcord)
{
	r={xcord, ycord, tcord};
	rplx={(xcord+1)%L, ycord, tcord};
	rply={xcord, (ycord+1)%L, tcord};
	rplt={ xcord, ycord, (tcord+1)%L};
	rminx={(xcord+L-1)%L, ycord, tcord};
	rminy={xcord, (ycord+L-1)%L, tcord};
	rmint={ xcord, ycord, (tcord+L-1)%L};

	Seff=S_site=S_derivative=0.0;

	for(int l=1; l<=4; l++)
	{
		S_site=S_site+s*dotProduct(H[l][xcord][ycord][tcord],H[l][xcord][ycord][tcord]);

		for(int m=1; m<=4; m++)
		{
			S_site=S_site+u0*dotProduct(H[l][xcord][ycord][tcord],H[l][xcord][ycord][tcord])*dotProduct(H[m][xcord][ycord][tcord],H[m][xcord][ycord][tcord]);
			S_site=S_site+u0*dotProduct(H[l][xcord][ycord][tcord],H[m][xcord][ycord][tcord])*dotProduct(H[l][xcord][ycord][tcord],H[m][xcord][ycord][tcord]);
		
			S_derivative=S_derivative-Kg2*dotProduct(H[l][rplx[0]][rplx[1]][rplx[2]],H[m][rplx[0]][rplx[1]][rplx[2]])*dotProduct(H[l][r[0]][r[1]][r[2]],H[m][r[0]][r[1]][r[2]]);
			S_derivative=S_derivative-Kg2*dotProduct(H[l][rply[0]][rply[1]][rply[2]],H[m][rply[0]][rply[1]][rply[2]])*dotProduct(H[l][r[0]][r[1]][r[2]],H[m][r[0]][r[1]][r[2]]);
			S_derivative=S_derivative-Kg2*dotProduct(H[l][rplt[0]][rplt[1]][rplt[2]],H[m][rplt[0]][rplt[1]][rplt[2]])*dotProduct(H[l][r[0]][r[1]][r[2]],H[m][r[0]][r[1]][r[2]]);
			S_derivative=S_derivative-Kg2*dotProduct(H[l][rminx[0]][rminx[1]][rminx[2]],H[m][rminx[0]][rminx[1]][rminx[2]])*dotProduct(H[l][r[0]][r[1]][r[2]],H[m][r[0]][r[1]][r[2]]);
			S_derivative=S_derivative-Kg2*dotProduct(H[l][rminy[0]][rminy[1]][rminy[2]],H[m][rminy[0]][rminy[1]][rminy[2]])*dotProduct(H[l][r[0]][r[1]][r[2]],H[m][r[0]][r[1]][r[2]]);
			S_derivative=S_derivative-Kg2*dotProduct(H[l][rmint[0]][rmint[1]][rmint[2]],H[m][rmint[0]][rmint[1]][rmint[2]])*dotProduct(H[l][r[0]][r[1]][r[2]],H[m][r[0]][r[1]][r[2]]);
		}
	}

	Seff=S_site+S_derivative;
	
	return Seff;
}
//---------------------------------------

//Update Higgs Field at one site
void HiggsUpdate(int xcord, int ycord, int tcord)
{
	double S_old, S_new;
	S_old=calcAction(xcord, ycord, tcord);

	for(int l=1; l<=4; l++)
	{
		Hnew=FlipHiggsAmp(H[l][xcord][ycord][tcord]);
		S_new=calcAction(xcord,ycord,tcord, Hnew, l);
		if(S_old>=S_new)
		{
			H[l][xcord][ycord][tcord]=Hnew;
			S_old=S_new;
		}
		else if(std_rand()<=expmindelS(S_old,S_new))
		{
			H[l][xcord][ycord][tcord]=Hnew;
			S_old=S_new;
		}

		Hnew=FlipHiggsAngleyt(H[l][xcord][ycord][tcord]);
		S_new=calcAction(xcord,ycord,tcord, Hnew, l);
		if(S_old>=S_new)
		{
			H[l][xcord][ycord][tcord]=Hnew;
			S_old=S_new;
		}
		else if(std_rand()<=expmindelS(S_old,S_new))
		{
			H[l][xcord][ycord][tcord]=Hnew;
			S_old=S_new;
		}

		Hnew=FlipHiggsAnglext(H[l][xcord][ycord][tcord]);
		S_new=calcAction(xcord,ycord,tcord, Hnew, l);
		if(S_old>=S_new)
		{
			H[l][xcord][ycord][tcord]=Hnew;
			S_old=S_new;
		}
		else if(std_rand()<=expmindelS(S_old,S_new))
		{
			H[l][xcord][ycord][tcord]=Hnew;
			S_old=S_new;
		}

		Hnew=FlipHiggsAnglexy(H[l][xcord][ycord][tcord]);
		S_new=calcAction(xcord,ycord,tcord, Hnew, l);
		if(S_old>=S_new)
		{
			H[l][xcord][ycord][tcord]=Hnew;
			S_old=S_new;
		}
		else if(std_rand()<=expmindelS(S_old,S_new))
		{
			H[l][xcord][ycord][tcord]=Hnew;
			S_old=S_new;
		}
	}
}

//MC sweep randomly
void MCsweepRandom()
{
	for(int i=0; i<L*L*L; i++)
	{
		x=rand()%L;
		y=rand()%L;
		t=rand()%L;
		HiggsUpdate(x,y,t);
	}
}

//MC sweem all sites
void MCsweepSequential()
{	
	for(x=0; x<L; x++)
	{
		for(y=0; y<L; y++)
		{
			for(t=0; t<L; t++)
			{
				HiggsUpdate(x,y,t);			
			}		
		}
	}
}

// Main function
int main(int argc, char *argv[])
{
	srand (time(NULL));					//Initialise Random Number Generator with seed=computer time

	initialiseHiggs();
	
	starttime=clock();
	for(steps=0; steps<coolingsweeps; steps++)
	{	
		MCsweepSequential();
		if((steps+1)%100==0)
			cout<<steps+1<<" "<<1.0*(steps+1)/coolingsweeps*100<<"%\n";
		
		if(steps<10) continue;
		endtime=clock();
		cout<<double(endtime - starttime) / CLOCKS_PER_SEC<<endl;
		exit(0);
	}

	return 1;
}