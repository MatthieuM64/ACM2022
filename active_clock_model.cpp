/*C++ CODE - MANGEAT MATTHIEU - 2022 */
/*q-STATE ACTIVE CLOCK MODEL*/

//////////////////////
///// LIBRAIRIES /////
//////////////////////

//Public librairies.
#include <cstdlib>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <map>
#include <string>
#include <string.h>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

//Personal libraries.
#include "lib/random.cpp"
#include "lib/special_functions.cpp"

////////////////////////////////////////
///// CLASS FOR ACTIVE CLOCK SPINS /////
////////////////////////////////////////

class particle
{
	public:
	
	double x, y; //position
	double dx, dy; //displacement
	double theta; //orientation
	
	particle(const int &q, const int &LX, const int &LY, const int &init, const double &theta0);
	void hop(const bool &favdir, const int &q, const int &LX, const int &LY);
	void flip(const double &thetap);
};

//Creation of the spin.
particle::particle(const int &q, const int &LX, const int &LY, const int &init, const double &theta0)
{
	//Random initial position + random orientation [0,2Pi].
	if (init==0)
	{
		x=ran()*LX;
		y=ran()*LY;
		theta=(2*M_PI/q)*int(q*ran());
	}
	//Random initial position + ordered orientation.
	else if (init==1)
	{
		x=ran()*LX;
		y=ran()*LY;
		theta=theta0;
	}
	//Initial transverse band.
	else if (init==2)
	{
		x=(2+ran())*0.2*LX;
		y=ran()*LY;
		theta=0.;
	}
	//Initial longitudinal lane.
	else if (init==3)
	{
		x=ran()*LX;
		y=(2+ran())*0.2*LY;
		theta=0.;
	}
	else
	{
		cerr << "BAD INIT VALUE: " << init << endl;
		abort();
	}
}

//Hop in the direction phi and implement the periodicity LX/LY.
void particle::hop(const bool &favdir, const int &q, const int &LX, const int &LY)
{
	double phi;
	//Move to the favoured direction theta if r<epsilon.
	if (favdir)
	{
		phi=theta;
	}
	//Move to a random direction phi if r>epsilon.
	else
	{
		phi=(2*M_PI/q)*int(q*ran());
	}
	//Remark: epsilon=0 the direction is purely random and epsilon=1 the direction is purely the favoured one.
	
	//Update the displacement.
	dx+=cos(phi);
	dy+=sin(phi);
	
	//Update the position (with periodicity).
	x+=cos(phi);
	y+=sin(phi);
	
	while (x<0)
	{
		x+=LX;
	}
	while (x>=LX)
	{
		x-=LX;
	}
	while (y<0)
	{
		y+=LY;
	}
	while (y>=LY)
	{
		y-=LY;
	}
}

//Flip to state thetap.
void particle::flip(const double &thetap)
{
	theta=thetap;
}

/////////////////////////////////////////////////
///// ARRAYS OF PARTICLE INDICES BY SECTORS /////
/////////////////////////////////////////////////

class sectors
{
	vector< vector< vector<int> > > sec;
	
	public:
	
	sectors(const int &LX, const int &LY);
	vector<int> get(const int &x, const int &y) const;
	int density(const int &x, const int &y) const;
	void add(const int &x, const int &y, const int &index);
	void remove(const int &x, const int &y, const int &index);
};

//Creation of the array.
sectors::sectors(const int &LX, const int &LY)
{
	sec=vector< vector< vector<int> > >(LX,vector< vector<int> >(LY, vector<int>(0)));
}

//Return the indices of site (x,y).
vector<int> sectors::get(const int &x, const int &y) const
{
	return sec[x][y];
}

//Return the density of site (x,y).
int sectors::density(const int &x, const int &y) const
{
	return sec[x][y].size();
}

//Add one particle to the site (x,y).
void sectors::add(const int &x, const int &y, const int &index)
{
	sec[x][y].push_back(index);
}

//Remove one particle to the site (x,y).
void sectors::remove(const int &x, const int &y, const int &index)
{
	int remove_pos=0;
	for (int k=0; k<sec[x][y].size(); k++)
	{
		if (sec[x][y][k]==index)
		{
			remove_pos=k;
			break;
		}
	}
	sec[x][y].erase(sec[x][y].begin()+remove_pos);
}

////////////////////////////////
///// AVERAGES IN SUBBOXES /////
////////////////////////////////

class averages
{
	public:
	
	int L0, rx, ry; //parameters of the boxes.
	vector<double> n, n2; //number fluctuations.
	vector<double> m, m2; //magnetization fluctuations.
	int Nav; //number of averages (in time).
	
	averages(const int &LX, const int &LY);
	void update(const int &Npart, const vector<particle> &ACP, const int &LX, const int &LY);
	void exportFile(const int &q, const double &beta, const double &epsilon, const double &rho0, const int &LX, const int &LY, const int &init, const int &RAN);
};

//Creation of the averaged quantities.
averages::averages(const int &LX, const int &LY)
{
	if (LX>=LY)
	{
		L0=LY;
		rx=LX/LY;
		ry=1;
	}
	else
	{
		L0=LX;
		rx=1;
		ry=LY/LX;
	}
	n=vector<double>(L0,0.);
	n2=vector<double>(L0,0.);
	m=vector<double>(L0,0.);
	m2=vector<double>(L0,0.);
	Nav=0;
}

//Update the averaged quantities at time t.
void averages::update(const int &Npart, const vector<particle> &ACP, const int &LX, const int &LY)
{
	//Density and magnetization of each sector at time t.
	vector< vector<int> > RHO(LX,vector<int>(LY,0));
	vector< vector<double> > MX(LX,vector<double>(LY,0.)), MY(LX,vector<double>(LY,0.));
	for (int i=0; i<Npart; i++)
	{
		RHO[int(ACP[i].x)][int(ACP[i].y)]++;
		MX[int(ACP[i].x)][int(ACP[i].y)]+=cos(ACP[i].theta);
		MY[int(ACP[i].x)][int(ACP[i].y)]+=sin(ACP[i].theta);
	}
	
	static const int Nboxes=10;
	//Averages on all sub-boxes.
	for (int l=1; l<L0; l++)
	{
		//Select Nboxes different sub-boxes of size l.
		for (int nbox=0; nbox<Nboxes; nbox++)
		{
			//Random position for the bottom left corner.
			const int x0=int((LX-rx*l)*ran());
			const int y0=int((LY-ry*l)*ran());
			
			//Density and magnetization in this sub-box.
			long unsigned int rho=0;
			double mx=0., my=0.;
			
			for (int x=x0; x<x0+rx*l; x++)
			{
				for (int y=y0; y<y0+ry*l; y++)
				{
					rho+=RHO[x][y];
					mx+=MX[x][y];
					my+=MY[x][y];
				}
			}
			
			const double mag2=mx*mx+my*my;
			
			//Add to n, n2, m, m2.
			n[l]+=double(rho)/Nboxes;
			n2[l]+=double(rho*rho)/Nboxes;
			m[l]+=sqrt(mag2)/Nboxes;
			m2[l]+=mag2/Nboxes;
		}
	}
	//Increase the number of averages.
	Nav++;
}

//Export averages in a file.
void averages::exportFile(const int &q, const double &beta, const double &epsilon, const double &rho0, const int &LX, const int &LY, const int &init, const int &RAN)
{
	static const int returnSystem=system("mkdir -p data_ACM_averages/");
	stringstream ss;	
	ss << "./data_ACM_averages/ACM_fluctuations_q=" << q << "_beta=" << beta << "_epsilon=" << epsilon << "_rho0=" << rho0 << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << ".txt";	
	string nameAV = ss.str();			
	ofstream fileAV(nameAV.c_str(),ios::trunc);	
	fileAV.precision(6);
	
	fileAV << "#The number of averages are: " << Nav << endl;

	for (int l=1; l<L0; l++)
	{
		double N=n[l]/Nav;
		double delN=n2[l]/Nav-N*N;
		double M=m[l]/Nav;
		double delM=m2[l]/Nav-M*M;
		fileAV << l << "\t" << N << "\t" << delN << "\t" << M << "\t" << delM << endl;
	}
	fileAV.close();
}

/////////////////////
///// FUNCTIONS /////
/////////////////////

//Distance between two particles in the periodic domain.
double distance2(const particle &part1, const particle &part2, const int &LX, const int &LY)
{
	const double DX=fabs(part1.x-part2.x);
	const double DY=fabs(part1.y-part2.y);	
	return square(min(DX,LX-DX)) + square(min(DY,LY-DY));
}

//Order parameter (total magnetization).
vector<double> mag(const int &Npart, const vector<particle> &ACP)
{
	double MX=0, MY=0;
	for (int i=0; i<Npart; i++)
	{
		MX+=cos(ACP[i].theta);
		MY+=sin(ACP[i].theta);
	}
	
	vector<double> MAG(2,0.);
	MAG[0]=MX/Npart;
	MAG[1]=MY/Npart;
	
	return MAG;
}

//Mean square displacement.
vector<double> msd(const int &Npart, const vector<particle> &ACP)
{
	double DX=0, DY=0;
	double DX2=0., DY2=0.;
	for (int i=0; i<Npart; i++)
	{
		DX+=ACP[i].dx/Npart;
		DX2+=ACP[i].dx*ACP[i].dx/Npart;
		DY+=ACP[i].dy/Npart;
		DY2+=ACP[i].dy*ACP[i].dy/Npart;
	}	
	vector<double> DR(2,0.);
	DR[0]=DX2+DY2;
	DR[1]=DX*DX+DY*DY;
	return DR;
}

//Modulo in the periodic domain.
int modulo(const int &x, const int &L)
{
	if (x<0)
	{
		return x+L;
	}
	else if (x>=L)
	{
		return x-L;
	}
	else
	{
		return x;
	}
}

//Export the density in a file.
void exportDensity(const int &q, const double &beta, const double &epsilon, const double &rho0, const int &LX, const int &LY, const int &init, const int &RAN, const int &t, const int &Npart, const vector<particle> &ACP)
{
	vector< vector<int> > RHO(LX,vector<int>(LY,0));
	for (int i=0; i<Npart; i++)
	{
		RHO[int(ACP[i].x)][int(ACP[i].y)]++;
	}
	
	static const int returnSystem=system("mkdir -p data_ACM_dynamics/");
	stringstream ssRHO;
	
	ssRHO << "./data_ACM_dynamics/ACM_RHO_q=" << q << "_beta=" << beta << "_epsilon=" << epsilon << "_rho0=" << rho0 << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << "_t=" << t << ".txt";	
	string nameRHO = ssRHO.str();			
	ofstream fileRHO(nameRHO.c_str(),ios::trunc);	
	fileRHO.precision(6);	
	
	for (int Y0=0; Y0<LY; Y0++)
	{
		for (int X0=0; X0<LX; X0++)
		{
			fileRHO << RHO[X0][Y0] << "\t";
		}
		fileRHO << endl;
	}
	
	fileRHO.close();
}

//Export the magnetization in a file.
void exportMagnetization(const int &q, const double &beta, const double &epsilon, const double &rho0, const int &LX, const int &LY, const int &init, const int &RAN, const int &t, const int &Npart, const vector<particle> &ACP)
{
	vector< vector<double> > MX(LX,vector<double>(LY,0.)), MY(LX,vector<double>(LY,0.));
	for (int i=0; i<Npart; i++)
	{
		MX[int(ACP[i].x)][int(ACP[i].y)]+=cos(ACP[i].theta);
		MY[int(ACP[i].x)][int(ACP[i].y)]+=sin(ACP[i].theta);
	}
	
	static const int returnSystem=system("mkdir -p data_ACM_dynamics/");
	stringstream ssMX,ssMY;
	
	ssMX << "./data_ACM_dynamics/ACM_MX_q=" << q << "_beta=" << beta << "_epsilon=" << epsilon << "_rho0=" << rho0 << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << "_t=" << t << ".txt";
	string nameMX = ssMX.str();
	ofstream fileMX(nameMX.c_str(),ios::trunc);
	fileMX.precision(6);
	
	ssMY << "./data_ACM_dynamics/ACM_MY_q=" << q << "_beta=" << beta << "_epsilon=" << epsilon << "_rho0=" << rho0 << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << "_t=" << t << ".txt";
	string nameMY = ssMY.str();
	ofstream fileMY(nameMY.c_str(),ios::trunc);
	fileMY.precision(6);
	
	for (int Y0=0; Y0<LY; Y0++)
	{
		for (int X0=0; X0<LX; X0++)
		{
			fileMX << MX[X0][Y0] << "\t";
			fileMY << MY[X0][Y0] << "\t";
		}
		fileMX << endl;
		fileMY << endl;
	}
	
	fileMX.close();
	fileMY.close();
}

///////////////////////////////////////
///// READ COMMAND LINE ARGUMENTS /////
///////////////////////////////////////

void ReadCommandLine(int argc, char** argv, int &q, double &beta, double &epsilon, double &rho0, int &LX, int &LY, int& init, int &RAN, int &tmax)
{
 	for( int i = 1; i<argc; i++ )
	{
		if (strstr(argv[i], "-q=" ))
		{
			q=atoi(argv[i]+3);
		}
		else if (strstr(argv[i], "-beta=" ))
		{
			beta=atof(argv[i]+6);
		}
		else if (strstr(argv[i], "-epsilon=" ))
		{
			epsilon=atof(argv[i]+9);
		}
		else if (strstr(argv[i], "-rho0=" ))
		{
			rho0=atof(argv[i]+6);
		}
		else if (strstr(argv[i], "-LX=" ))
		{
			LX=atoi(argv[i]+4);
		}
		else if (strstr(argv[i], "-LY=" ))
		{
			LY=atoi(argv[i]+4);
		}
		else if (strstr(argv[i], "-init=" ))
		{
			init=atoi(argv[i]+6);
		}
		else if (strstr(argv[i], "-ran=" ))
		{
			RAN=atoi(argv[i]+5);
		}
		else if (strstr(argv[i], "-tmax=" ))
		{
			tmax=atoi(argv[i]+6);
		}
		else
		{
			cerr << "BAD ARGUMENT : " << argv[i] << endl;
			abort();
		}
	}
}

/////////////////////
///// MAIN CODE /////
/////////////////////

int main(int argc, char *argv[])
{
	//Physical parameters: beta=inverse temperature, epsilon=self-propulsion, rho0=average density, LX*LY=size of the box, q=number of states.
	const double D0=1., J=1.;
	double beta=2., epsilon=0.9, rho0=6;
	int LX=100, LY=100, q=7;
	
	//Numerical parameters: tmax=maximal time, init=initial condition, RAN=index of RNG.
	int tmax=1000000, init=1, RAN=0;

	//Parameters in arguments.
	ReadCommandLine(argc,argv,q,beta,epsilon,rho0,LX,LY,init,RAN,tmax);
	
	//Verify the values of parameters.
	if (init<0 or init>3)
	{
		cerr << "BAD VALUE OF INIT: " << init << endl;
		return 1;
	}

	//Start the random number generator.
	init_gsl_ran();
	cout << "GSL index = " << RAN << "\n";
	gsl_rng_set(GSL_r,RAN);
	
	//Total number of particles.
	const int Npart=int(LX*LY*rho0);
	
	//Number of particles of each state, on each sites.
	sectors SEC(LX,LY);
	
	//Creation of active clock particles.
	vector<particle> ACP;
	const double theta0=(2*M_PI/q)*(RAN%q);
	for (int i=0;i<Npart;i++)
	{
		particle ACM0(q,LX,LY,init,theta0);
		SEC.add(int(ACM0.x),int(ACM0.y),i);
		ACP.push_back(ACM0);
	}
	
	//Creation of file for averages.
	const int returnSystem=system("mkdir -p data_ACM_averages/");
	stringstream strAVERAGES;
	strAVERAGES << "./data_ACM_averages/ACM_AVERAGES_q=" << q << "_beta=" << beta << "_epsilon=" << epsilon << "_rho0=" << rho0 << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << ".txt";
	string nameAVERAGES = strAVERAGES.str();			
	ofstream fileAVERAGES(nameAVERAGES.c_str(),ios::trunc);
	fileAVERAGES.precision(6);
	
	//Number and magnetization fluctuations.
	averages AV(LX,LY);
	double teq=2000;
	
	//Time increment.
	const double DeltaT=1./(D0+exp(2*beta*J));
	
	//Get the probability to hop (constant over the time).
	const double proba_hop=D0*DeltaT;
	
	cout.precision(6);
	
	//Time evolution.
	for(int t=0;t<=tmax;t++)
	{
		if (t%10==0 or t==tmax)
		{
			int rho_average=0;
			for (int Y0=0;Y0<LY;Y0++)
			{
				for (int X0=0; X0<LX; X0++)
				{
					rho_average+=SEC.density(X0,Y0);
				}
			}
			const double RHO0=double(rho_average)/(LX*LY);
			const vector<double> MAG=mag(Npart,ACP);
			const double MX=MAG[0], MY=MAG[1];
			const double VALPHA=sqrt(square(MX)+square(MY)), THETA=atan2(MY,MX);
			const vector<double> MSD=msd(Npart,ACP);
			
			fileAVERAGES << t << "\t" << RHO0 << "\t" << VALPHA << "\t" << THETA << "\t" << MX << "\t" << MY << "\t" << MSD[0] << "\t" << MSD[1] << endl;
			cout << "time=" << t << " -N/V=" << RHO0 << " -M=" << VALPHA << " -THETA=" << THETA << " -MX=" << MX << " -MY=" << MY << " -MSD=" << MSD[0] << " -R2=" << MSD[1] << running_time.TimeRun(" ") << endl;
		}
		if (t>teq)
		{
			AV.update(Npart,ACP,LX,LY);
		}
		if (t%(tmax/100)==0 or t==tmax)
		{
			exportDensity(q,beta,epsilon,rho0,LX,LY,init,RAN,t,Npart,ACP);
			exportMagnetization(q,beta,epsilon,rho0,LX,LY,init,RAN,t,Npart,ACP);
			AV.exportFile(q,beta,epsilon,rho0,LX,LY,init,RAN);
		}
		
		//At each time step update (in average) all the particles.
		for (int i=0;i<Npart;i++)
		{
			//Choose a particle randomly (j), in the sector (X,Y) with spin theta.
			const int j=int(ran()*Npart);
			const int X0=int(ACP[j].x), Y0=int(ACP[j].y);
			const double theta=ACP[j].theta;
			
			//Determination of the new orientation (uniformly).
			double thetap=(2*M_PI/q)*int((q-1)*ran());
			if (thetap>=theta)
			{
				thetap+=2*M_PI/q;
			}
			if (thetap==theta)
			{
				cerr << "BAD VALUE OF THETAP: " << thetap << " THETA=" << theta << endl;
				return 1;
			}
			
			//Calculate the probability to flip in the new orientation.
			int rhoj=1;
			double MX=0., MY=0.;			
			//Take the energy for particles in neighbour sites with a distance smaller than 1.
			for (int XN=X0-1;XN<=X0+1;XN++)
			{
				for (int YN=Y0-1;YN<=Y0+1;YN++)
				{
					const vector<int> neighbours=SEC.get(modulo(XN,LX),modulo(YN,LY)); //Particles in the square XN,YN.
					for (int ll=0; ll<neighbours.size(); ll++)
					{
						const int k=neighbours[ll]; //Index of this particle -> k.
						if (j!=k and distance2(ACP[j],ACP[k],LX,LY)<1)
						{
							MX+=cos(ACP[k].theta);
							MY+=sin(ACP[k].theta);
							rhoj++;
						}
					}
				}
			}
			double DeltaH=MX*(cos(thetap)-cos(theta)) + MY*(sin(thetap)-sin(theta));
			double proba_flip=exp(beta*J*DeltaH/rhoj)*DeltaT;
			
			//Verify that the probability to wait is positive.
			if (proba_hop+proba_flip>1)
			{
				cerr << "THE PROBABILITY TO WAIT IS NEGATIVE: proba_hop=" << proba_hop << " proba_flip=" << proba_flip << endl;
				cerr << "CHANGE THE VALUE OF DELTA_T !" << endl;
				return 1;
			}
			
			double random_number=ran();
			//The particle hops: perform the hopping on the particle and update the population of sectors.
			if (random_number<proba_hop)
			{
				ACP[j].hop(random_number<epsilon*proba_hop,q,LX,LY);
				
				//UPDATE OF SECTOR!
				if (X0!=int(ACP[j].x) or Y0!=int(ACP[j].y))
				{
					SEC.remove(X0,Y0,j);
					SEC.add(int(ACP[j].x),int(ACP[j].y),j);
				}
			}
			//The particle flips: perform the flipping on the particle (on-site, ind. of epsilon).
			else if (random_number<proba_hop+proba_flip)
			{
				ACP[j].flip(thetap);
				//NO UPDATE OF SECTORS! The particle has not moved.
			}
			//Else do nothing (proba_wait).		
		}
	}
	return 0;
}
