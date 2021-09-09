// DigitRecognition.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include<iostream>
#include<vector>
#include <stdio.h>
#include<math.h>
#include <sstream>
#include<fstream>
#include <windows.h>
#include <iomanip>
#include<string>

using namespace std;

#define N 5

# define M 32

#define pi 3.14159265

#define NSamples 320  // no of samples in one frame 

#define p 12  




vector<string>Country;




long double THRESHOLD = pow((long double)10,-30);


int T;

vector< vector<long double> >codebook;


void manipulateB(vector<vector<long double>>&B){
	for(int i=1;i<=N;i++){
		for(int j=1;j<=M;j++){
			if(B[i][j] < THRESHOLD){
				long double dif = THRESHOLD - B[i][j];
				B[i][j] = THRESHOLD;
				long double maxm=0.0;
				int idx=1;
				int k;
				for(k=1;k<=M;k++){
					if(B[i][k] > maxm){
						maxm=B[i][j];
						idx=k;
					}
				}
				B[i][idx] -= dif;
			}
		}
	}
}

vector<long double> read_file(string filename)
{
	fstream file;
	vector<long double>Sample;
	file.open(filename,ios::in);
	int count=0;
	if(file.is_open())
	{
	//	cout<<"hello\n";
		string line; 
		while(getline(file,line))                       // reading line by line
		{
				if(count>=8)
				{

				long double num = stold(line) + 100.0;					// get the double value for string 

				Sample.push_back(num);
				}
				count++;
		}		
		file.close();
	}
	else
	{
		cout<<"File cannot be open\n";
	}
	//cout<<Sample.size()<<"\n";
	return Sample;  // storing values(line by line) present in filename.txt file in vector Sample
}


// creating feed forward model 
void feed_forward_model(vector<vector<long double>>&A,vector<vector<long double>>&B,vector<long double>&PI)
{
	
	for(int i=1;i<=N;i++)
	{
		if(i==1)
		{
			PI[i] = 1;
		}
		else
		{
			PI[i] = 0;
		}
	}

	for(int i=1;i<=N;i++)
	{
		for(int j=1;j<=N;j++)
		{
			if(i==N && j==N)
			{
				A[i][j] = 1;
			}
			else if(i==j)
			{
				A[i][j] = 0.8;
			}
			else if(i+1==j)
			{
				A[i][j] = 0.2;
			}
			else
			{
				A[i][j] = 0;
			}
		}
	}

	for(int i=1;i<=N;i++)
	{
		for(int j=1;j<=M;j++)
		{
			B[i][j] = 1.0/M;
		}
	}
}

long double forward_procedure(vector<vector<long double>>&A,vector<vector<long double>>&B,vector<long double>&PI,vector<vector<long double>>&alpha,vector<int>&obs)
{
	for(int i=1;i<=N;i++)
	{
		alpha[1][i] = PI[i]*B[i][obs[1]];
	}
	// induction

	for(int t = 1;t<=T-1;t++)
	{
		for(int j=1;j<=N;j++)
		{
			
			long double temp = 0;

			for(int i=1;i<=N;i++)
			{
				temp = temp + alpha[t][i]*A[i][j];
			}

			alpha[t+1][j] = temp*B[j][obs[t+1]];

			//cout<<alpha[t+1][j]<<"\n";
		}
	}

	// termination

	long double sum = 0;

	for(int i=1;i<=N;i++)
	{
//		cout<<alpha[T][i];
		sum  = sum + alpha[T][i];
	}
//	cout<<sum<<"?????????\n";
	return sum;
}

void backward_procedure(vector<vector<long double>>&A,vector<vector<long double>>&B,vector<long double>&PI,vector<vector<long double>>&beta,vector<int>obs)
{
	

	
	// initialization

	for(int i=1;i<=N;i++)
	{
		beta[T][i] = 1;
	}

	// induction

	for(int t=T-1;t>=1;t--)
	{
		for(int i=1;i<=N;i++)
		{
			 long double temp = 0;

			 for(int j=1;j<=N;j++)
			 {
				 temp = temp + A[i][j]*B[j][obs[t+1]]*beta[t+1][j];
			 }
			 beta[t][i] = temp;
		}
	}
}

long double viterbi(vector<vector<long double>>&A,vector<vector<long double>>&B,vector<long double>&PI,vector<int>&obs)
{

	vector<vector<long double>>delta(T+1,vector<long double>(N+1));
	vector<vector<int>>psi(T+1,vector<int>(N+1));

	// initialization

	for(int i=1;i<=N;i++)
	{
		delta[1][i] = PI[i]*B[i][obs[1]];
		psi[1][i] = 0;

	}

	// recursion
	
	for(int t=2;t<=T;t++)
	{
		for(int j=1;j<=N;j++)
		{
			long double temp1 = delta[t-1][1]*A[1][j];

			long double temp2 = 1;

			for(int i=2;i<=N;i++)
			{
				if(temp1 < delta[t-1][i]*A[i][j])
				{
					temp1 = delta[t-1][i]*A[i][j];
					temp2 = i;
				}
			}
			delta[t][j] = temp1*B[j][obs[t]];
			psi[t][j] = temp2;


		}
	}
	vector<int>q_star(T+1);

	// termination
	long double p_star=0.0;
	
	for(int i=1;i<=N;i++)
	{
		if(p_star < delta[T][i])
		{
			p_star = delta[T][i];
			q_star[T] = i;
		}
	}
	
	// path backtracking
	

	for(int t=T-1;t>=1;t--)
	{
		q_star[t] = psi[t+1][q_star[t+1]];
	}


	
	return p_star;
}

void restimation(vector<vector<long double>>&A,vector<vector<long double>>&B,vector<long double>&PI,vector<vector<long double>>&alpha,vector<vector<long double>>&beta,vector<int>&obs)
{
	
	vector< vector< vector<long double> > > epsilon(T+1 , vector< vector<long double> > (N+1, vector<long double> (N+1) ) );

	vector<vector<long double>> gamma(T+1,vector<long double>(N+1));

	///cout<<alpha[0].size()<<" "<<beta[0].size()<<"\n";

	for(int t=1;t<=T-1;t++)
	{
		long double temp = 0;

		for(int i=1;i<=N;i++)
		{
			for(int j=1;j<=N;j++)
			{

				temp = temp + alpha[t][i]*A[i][j]*B[j][obs[t+1]]*beta[t+1][j];
			}
		}

		
		
		for(int i=1;i<=N;i++)
		{
			for(int j=1;j<=N;j++)
			{

				epsilon[t][i][j] = alpha[t][i]*A[i][j]*B[j][obs[t+1]]*beta[t+1][j]*1.0/temp;
			}
		}
	
	}


	for(int t=1;t<=T;t++)
	{
		for(int i=1;i<=N;i++)
		{
			
			long double temp = 0;

			for(int j=1;j<=N;j++)
			{
				temp = temp + epsilon[t][i][j];
			}
			
			gamma[t][i] = temp;
		}
	}

	vector<long double>PI_comp(N+1);

	for(int i=1;i<=N;i++)
	{
		PI_comp[i] = gamma[1][i];
	}

	vector<vector<long double>>A_comp(N+1,vector<long double>(N+1));

	for(int i=1;i<=N;i++)
	{
		for(int j=1;j<=N;j++)
		{
			long double temp1=0,temp2=0;

			for(int t=1;t<=T-1;t++)
			{
				temp1 = temp1 + epsilon[t][i][j];
			}

			for(int t=1;t<=T-1;t++)
			{
				temp2 = temp2 + gamma[t][i];
			}

			A_comp[i][j] = temp1*1.0/temp2;
		}
	}

	vector<vector<long double>>B_comp(N+1,vector<long double>(M+1));

	for(int j=1;j<=N;j++)
	{
		for(int k=1;k<=M;k++)
		{

			long double temp1 = 0;
			for(int t=1;t<=T;t++)
			{
				if(obs[t] == k)                    
				{
					temp1 = temp1 + gamma[t][j];
				}
			}

			long double temp2 = 0;

			for(int t=1;t<=T;t++)
			{
				temp2 = temp2 + gamma[t][j];
			}

			B_comp[j][k] = temp1*1.0/temp2;
		}
	}

	PI = PI_comp;

	A = A_comp;

	B = B_comp;
}


void read_model(vector<vector<long double>>&A,vector<vector<long double>>&B,vector<long double>&PI,string filename)
{
	fstream file;

	file.open(filename,ios::in);

	if(file.is_open())
	{
		string line;
		int count1=1;	
		while(getline(file,line))
		{
				
			stringstream s(line); 

			if(count1 == 1)
			{
				int count2 = 1;
				while(getline(s,line,' ')) { 
					
					long double num = stold(line);					// get the double value for string 

					PI[count2] = num;
					count2++;

				}
			}
			else if(count1<=N+1)
			{
				int count2 = 1;
				while(getline(s,line,' ')) { 
					
					long double num = stold(line);					// get the double value for string 

					A[count1-1][count2] = num;
					count2++;

				}
			}
			else
			{
				int count2 = 1;
				while(getline(s,line,' ')) { 
					
					long double num = stold(line);					// get the double value for string 

					B[count1-6][count2] = num;
					count2++;

				}
			}
			count1++;
		}
			
		file.close();
	}
	else
	{
		cout<<"File Cannot be open\n";
	}

}

void print_model(vector<vector<long double>>&A,vector<vector<long double>>&B,vector<long double>&PI)
{
	cout<<"PI : \n";

	for(int i=1;i<PI.size();i++)
	{
		cout<<PI[i]<<" ";
	}
	cout<<"\n";

	cout<<"A : \n";

	for(int i=1;i<A.size();i++)
	{
		for(int j=1;j<A[i].size();j++)
			cout<<A[i][j]<<" ";
		cout<<"\n";
	}

	cout<<"\nB : \n";

	for(int i=1;i<B.size();i++)
	{
		for(int j=1;j<B[i].size();j++)
			cout<<B[i][j]<<" ";
		cout<<"\n";
	}
}



void save_model(vector< vector<long double> >&A,vector< vector<long double> >&B,vector<long double>&PI,int a,int b)
{
	ostringstream strg1;
	strg1<< a;
	string s1 = strg1.str();
	string filename;
	filename += s1;
	filename += "_";
	ostringstream strg2;
	strg2<< b;
	string s2 = strg2.str();
	filename += s2;
	filename += ".txt";
  
	cout<<"created model " <<filename<<"\n";
	ofstream MyFile(filename);

  // Write to the file
 
	for(int i=1;i<PI.size();i++)
	{
		MyFile<<PI[i]<<" ";
	}
	MyFile<<"\n";
	for(int i=1;i<A.size();i++)
	{
		for(int j=1;j<A[i].size();j++)
			MyFile<<A[i][j]<<" ";
		MyFile<<"\n";
		
	}


	for(int i=1;i<B.size();i++)
	{
		for(int j=1;j<B[i].size();j++)
			MyFile<<B[i][j]<<" ";
		MyFile<<"\n";
	}
  // Close the file
  MyFile.close();
}

void avg_model(int a)
{

	vector<vector<long double>>A(N+1,vector<long double>(N+1)); 

	vector<vector<long double>>B(N+1,vector<long double>(M+1));

	vector<long double>PI(N+1);


	for(int i=1;i<PI.size();i++)
	{
		PI[i] = 0;
	}

	for(int i=1;i<A.size();i++)
	{
		for(int j=1;j<A[i].size();j++)
			A[i][j]= 0;
	}

	for(int i=1;i<B.size();i++)
	{
		for(int j=1;j<B[i].size();j++)
			B[i][j]=0;
	}


	for(int i=1;i<=20;i++)
	{


		ostringstream strg1;
		strg1<< a;
		string s1 = strg1.str();
		string filename;
		filename += s1;
		filename += "_";
		ostringstream strg2;
		strg2<< i;
		string s2 = strg2.str();
		filename += s2;
		filename += ".txt";
		fstream file;

		file.open(filename,ios::in);

		if(file.is_open())
		{
			string line;
			int count1=1;	
			while(getline(file,line))
			{
				
				stringstream s(line); 

				if(count1 == 1)
				{
					int count2 = 1;
					while(getline(s,line,' ')) { 
					
						long double num = stold(line);					// get the double value for string 

						PI[count2] += num;
						count2++;

					}
				}
				else if(count1<=N+1)
				{
					int count2 = 1;
					while(getline(s,line,' ')) { 
					
						long double num = stold(line);					// get the double value for string 

						A[count1-1][count2] += num;
						count2++;

					}
				}
				else
				{
					int count2 = 1;
					while(getline(s,line,' ')) { 
					
						long double num = stold(line);					// get the double value for string 

						B[count1-6][count2] += num;
						count2++;

					}
				}
				count1++;
			}
			
			file.close();
		}
		else
		{
			cout<<"File Cannot be open\n";
		}

	}

	
	for(int i=1;i<=N;i++)
	{
		PI[i] = PI[i]/20.0;
	}

	for(int i=1;i<=N;i++)
	{
		for(int j=1;j<=N;j++)
		{
			A[i][j] = A[i][j]/20.0;
		}
	}

	for(int i=1;i<=N;i++)
	{
		for(int j=1;j<=M;j++)
		{
			B[i][j] = B[i][j]/20.0;
		}
	}

	save_model(A,B,PI,a,0);
}

void HMM(vector<int>&obs,int i,int j,int cycle)
{

	
vector<vector<long double>>A(N+1,vector<long double>(N+1)); 

vector<vector<long double>>B(N+1,vector<long double>(M+1));

vector<long double>PI(N+1);
	

	vector<vector<long double> >alpha(T+1,vector<long double>(N+1,0.0));

	vector<vector<long double>>beta(T+1,vector<long double>(N+1,0.0));

	if(cycle == 0)
	{

		feed_forward_model(A,B,PI);
	}
	else
	{
		ostringstream strg1;

		strg1<< i;
		string s1 = strg1.str();
		string filename;
		filename += s1;
		filename += "_0";
		filename += ".txt";
		read_model(A,B,PI,filename);
	}

	manipulateB(B);

	int counter = 0;

	while(counter < 20)
	{
		
		forward_procedure(A,B,PI,alpha,obs);//<<"\n";

		backward_procedure(A,B,PI,beta,obs);
		
		viterbi(A,B,PI,obs);//<<"\n";
		
		restimation(A,B,PI,alpha,beta,obs);
		

		manipulateB(B);

		counter++;
	}

	save_model(A,B,PI,i,j);

	//print_model(A,B,PI);
}
/* DC shift */


void DCShift( vector<long double>&samples )
{
	
	// avg contain the value we have to remove from each sample so to remove DC shift 
	double sum = 0;

	for(int i=0;i<samples.size();i++)
	{
		sum = sum + samples[i];
	}

	double avg = sum*1.0/samples.size();

	// removing avg from each sample

	for(int i=0;i<samples.size();i++)
	{
		samples[i] = samples[i] - avg;
	}

}

/* Normalization */

void Normalization( vector<long double>&samples )
{

	// resizing max value of samples to 10000 
	
	long double max_value = samples[0];

	for(int i=1;i<samples.size();i++)
	{
		if(max_value > samples[i])
		{
			max_value = samples[i];
		}
	}

	
	long double factor = 10000.0/max_value;

	for(int i=0;i<samples.size();i++)
	{
		samples[i] = samples[i]*factor;
	}

}

/* creating frame each of 320 samples */

vector<vector<long double>> framing(vector<long double>&samples)
{
	

	vector< vector<long double> >Frames;

	vector<long double>temp;

	for(int i=0;i<samples.size();i=i+1)
	{
		if(i%320 == 0 && i!=0)
		{
			Frames.push_back(temp);

			temp.clear();
		}

		temp.push_back(samples[i]);
	}
			
	return Frames;
}


/* For applying hamming window on frame */

void hamming_window(vector<double>&s)
{
	for(int i=0;i<s.size();i++)
	{
		s[i] = s[i] * (0.54 - 0.46*cos(2*pi*i/319));
	}
}


/* Correlation() is used to find R's for a frame */

vector<long double> Correlation(vector<long double>&s)
{
	vector<long double>R(13,0);

	for(int i=0;i<=12;i++)
	{
		for(int j=0;j<=320-1-i;j++)
		{
			R[i] += s[j]*s[j+i];
		}
	}

	return R;
}


/* Levinson durbin algorithm -------- used to find a_i's for a frame */

vector<long double> Levinson_durbin(vector<long double>&R)
{
	

	vector<long double>E(13,0);
	
	vector<long double>K(13,0);

	long double alpha[13][13];

	memset(alpha,0,sizeof(alpha));
	
	E[0] = R[0];

	for(int i=1;i<=p;i++)
	{
		
		K[i] = R[i];

		long double sum = 0;

		for(int j=1;j<=i-1;j++)
		{
			
			sum += alpha[i-1][j]*R[i-j];
		}

		K[i] = (K[i] - sum)/E[i-1];

		alpha[i][i] = K[i];


		for(int j=1;j<=i-1;j++)
		{
			alpha[i][j] = alpha[i-1][j] - K[i]*alpha[i-1][i-j];
		}

		E[i] = (1-K[i]*K[i])*E[i-1];
	}

	vector<long double>a;

	for(int i=1;i<=p;i++)
	{
		a.push_back(alpha[p][i]);
	}
	return a;
}


/* It used to find Capestral coefficients for frame */

vector<long double>capestral_coefficients(vector<long double>&A,vector<long double>&R)
{
	vector<long double>C(13,0);

	C[0] = log(R[0]*R[0])/log(long double(2));

	for(int m=1;m<=12;m++)
	{
		C[m] = A[m-1];

		for(int k=1;k<m;k++)
		{
			C[m] += (k*1.0/m)*C[k]*A[m-k-1];
		}
	}

	return C;
}

/* applying raisine sine window on C */
void raisine_sine_window(vector<long double>&C)
{
	int Q = 12;

	vector<long double>w(Q+1,0);

	for(int m=1;m<=Q;m++)
	{
		w[m] = 1 + Q*sin(pi*m*1.0/Q)/2.0;
	}

	for(int i=1;i<C.size();i++)
	{
		C[i] *= w[i];
	}
}

vector<long double> calculate_CR(vector<long double>&Sample)
{


	DCShift(Sample); // applying DC shift

	Normalization(Sample); // applying normalization

	vector<long double> R = Correlation(Sample); // calculate R

	vector<long double> a = Levinson_durbin(R); // calculate a

	vector<long double> C = capestral_coefficients(a,R); // calculate C

	raisine_sine_window(C); // applying raise sine window on C
	
	vector<long double>vec;

	for(int i=1;i<C.size();i++)
	{
		vec.push_back(C[i]);
	}

	return vec;

}

// Tokhuras distance 

long double tokhuras_distance(vector<long double>&A,vector<long double>&B)
{
	long double sum = 0;

	long double w[] = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};

	for(int i=0;i<12;i++)
	{
		sum = sum + w[i] * (A[i] - B[i])*(A[i] - B[i]);
	}

	return sum;
}




vector<int> generate_observation(vector<long double>samples)
{

	vector< vector<long double> > frames = framing(samples);

	vector<int>obs;
	obs.push_back(0);

	for(int j=0;j<frames.size();j++)
	{
		
		vector<long double>vec = calculate_CR(frames[j]);
		long double value;
		int index;
		for(int i=0;i<codebook.size();i++)
		{
			long double temp = tokhuras_distance(codebook[i],vec);
			if(i==0)
			{
				index = 1;
				value = temp;
			}
			else if(temp < value)
			{
				value = temp;
				index = i+1;
			}
		}

		obs.push_back(index);
	}
	return obs;
}

vector<vector<long double> > readDim2(string filename,int c_p){
    ifstream file(filename);
    string str;
    vector<vector< long double> > v;
	vector<long double> temp;
//	v.push_back(temp);
	temp.push_back(0.0);
	int i=0;
    while (file >> str){
      temp.push_back(stold(str));
	  i=(i+1)%c_p;
	  if(i==0){
		  v.push_back(temp);
		  temp.clear();
		  temp.push_back(0.0);
	  }
	}
    return v;
}





bool winner(char matrix[][3],char ch)
{
	if(matrix[0][0] == ch && matrix[0][1] == ch && matrix[0][2] == ch)
	{
		return true;
	}
	if(matrix[1][0] == ch && matrix[1][1] == ch && matrix[1][2] == ch)
	{
		return true;
	}
	if(matrix[2][0] == ch && matrix[2][1] == ch && matrix[2][2] == ch)
	{
		return true;
	}
	if(matrix[0][0] == ch && matrix[1][0] == ch && matrix[2][0] == ch)
	{
		return true;
	}
	if(matrix[0][1] == ch && matrix[1][1] == ch && matrix[2][1] == ch)
	{
		return true;
	}
	if(matrix[0][2] == ch && matrix[1][2] == ch && matrix[2][2] == ch)
	{
		return true;
	}
	if(matrix[0][0] == ch && matrix[1][1] == ch && matrix[2][2] == ch)
	{
		return true;
	}
	if(matrix[0][2] == ch && matrix[1][1] == ch && matrix[2][0] == ch)
	{
		return true;
	}

	return false;
}


int real_time_test()
{
			int index = 1;
			long double prob = 0.0;
			system("Recording_Module.exe 3 input_file.wav filename.txt");
			vector<long double> samples = read_file("filename.txt");
		//	cout<<samples.size()<<"\n";
			if(samples.size() < 100)
			{
				return -1;
			}
			vector<int>obs = generate_observation(samples);
			T = obs.size()-1;
	//		cout<<T<<":::::::\n";

			for(int m=1;m<6;m++)
			{
			//	cout<<"hello\n";
				ostringstream strg1;

				strg1<< m;
				string s1 = strg1.str();
				string filename;
				filename += s1;
				filename += "_0";
				filename += ".txt";
				fstream file;


			
	//			cout<<filename<<"\n";
				vector<vector<long double>>A(N+1,vector<long double>(N+1)); 

				vector<vector<long double>>B(N+1,vector<long double>(M+1));

				vector<long double>PI(N+1);

				vector<vector<long double> >alpha(T+1,vector<long double>(N+1,0.0));
				read_model(A,B,PI,filename);
				
	//			cout<<"hello\n";
			//	cout<<A[0].size()<<" "<<B[0].size()<<" "<<PI.size()<<" ";
				long double temp = forward_procedure(A,B,PI,alpha,obs);
			//	cout<<"hellp\n";
		//		cout<<"hello\n";
		//		cout<<temp<<"\n"; 
			
				if(m==1)
				{
					prob = (long double)temp;
					index = 1;
				}
				else if(temp>prob)
				{
					index = m;
					prob = (long double)temp;
				}

			}

			return index;
}


void Game()
{
	cout<<"Speak Start to begin the game\n";


	Sleep(1000);
	//cout<<real_time_test();//!=4
//	return;

	while(real_time_test()!=4)
	{
		cout<<"Please speak start\n";
	}

	cout<<"Game started......................\n\n";

	cout<<"Player 1 ( 0 ) & Player 2 ( X )"<<"\n\n";

	char matrix[3][3];

	int flag[3][3];

	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			matrix[i][j] = ' ';
			cout<<matrix[i][j]<<" | ";
			
			flag[i][j]=0;
		}
		cout<<"\n------------\n";
	}


	

	

	int count=0;
	while(count<9)
	{

		Sleep(5000);
		
		if(count%2==0)
		{
			cout<<"Player 1 turn : Speak \n";
		}
		else
		{
			cout<<"Player 2 turn : Speak \n";
		}
		
		int index1=0,index2=0;

		while(true)
		{
		

			bool var = true;

			while( var == true)
			{
				cout<<"Speak the row number\n";

				index1 = real_time_test();

				if(index1 > 3)
				{
					cout<<"Please speak again not correctly recognised\n";
				}
				else
				{
					cout<<"Predicted : "<<index1<<"\n";
					cout<<"if correct predicted press(y/n)\n";

					char ch;
					cin>>ch;
					if(ch == 'y')
					{
						var = false;
					}
					else
					{
						cout<<"Please speak again not correctly recognised\n";
					}
				}
			}

			var = true;

			while( var == true)
			{
				cout<<"Speak the column number\n";

				index2 = real_time_test();

				if(index2 > 3)
				{
					cout<<"Please speak again not correctly recognised\n";
				}
				else
				{
					cout<<"Predicted : "<<index2<<"\n";
					cout<<"if correct predicted press(y/n)\n";

					char ch;
					cin>>ch;
					if(ch == 'y')
					{
						var = false;
					}
					else
					{
						cout<<"Please speak again not correctly recognised\n";
					}
				}
			}

			if(flag[index1-1][index2-1] == 1)
			{
				cout<<"already used choose correct location\n";
				continue;
			}

			flag[index1-1][index2-1] = 0;
			if(count % 2 == 0)
			{
			//	cout<<"hello\n";
				matrix[index1-1][index2-1]  = '0';

				cout<<index1-1<<" "<<index2-1<<"\n";

				for(int i=0;i<3;i++)
			{
				for(int j=0;j<3;j++)
				{
					cout<<matrix[i][j]<<" | ";
				}
				cout<<"\n------------\n";
			}
		//		cout<<matrix[index1-1][index2-1];
				
			}
			else
			{
			//	cout<<"hello\n";
				matrix[index1-1][index2-1] = 'X';

				for(int i=0;i<3;i++)
			{
				for(int j=0;j<3;j++)
				{
					cout<<matrix[i][j]<<" | ";
				}
				cout<<"\n------------\n";
			}
				
			}
			

		//	cout<<matrix[index1-1][index2-1]<<"------------------\n";

			


			break;

		}
	//	cout<<matrix[2][2]<<"\n";
		
		if(winner(matrix,'0') == true)
				{
					cout<<"Player 1 wins\n";
					break;
				}
			
		if(winner(matrix,'X') == true)
				{
					cout<<"Player 2 wins\n";
					break;
				}
			count++;

			
		}



		if(count==9)
		{
			cout<<"Match draw\n";
		}
}




int _tmain(int argc, _TCHAR* argv[])

{

	codebook = readDim2("codebook.txt",12);
	Game();
	return 0;
}

