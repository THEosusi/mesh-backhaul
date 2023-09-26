/**
 * C++PhysarumSolver(like C)
 *
 * Please README (including: description of Physarum Solver, Program construction, etc.)
 *
 * Usage
 *
 * Compile... g++ -option main.cpp
 * Execute... ./a.out [number of node] [./topology/filename]
 * here,
 * number of node		: Optional number of node on the simulation topology, unsigned integer.
 * ./topology/filename	: Optional name of simulation topology file. Directory structure is also optional.
 *
 */

#include "pch.h"
#include "MT.h"
#include "ICCG.h"
#include "Pajeck_experience.h"
#include "output.h"
#include "RGAP.h"
#include <cmath>

#define INF 10000.0		// Infinite length representation. Exclude L from the calculation, where L==INF
#define NEG -1.0		// For representation of negative
#define gamma 1.5		// Slope of the sigmoid function 1.5
#define delta_time 0.01 //Δt scale default 0.01
#define plot 1			// Number of loops you want to plot
#define MAX_FLOW 10.0
#define coefficient_tanh 1 //  coefficient of tanh

#define SIG 0 // change of method 0:|Q|, 1:sig
#define NUM 1 // change of method 0:normal, 1:new


// Basis parameter
vector<double> Q_Kirchhoff;				// ΣQ
vector<double> P_tubePressure;			// p_i(t) - p_j(t)
vector<vector<double>> Q_tubeFlow;		// Q_ij(t)
vector<vector<double>> D_tubeThickness; // D_ij(t)
vector<vector<double>> L_tubeLength;	// L_ij(t)


// To calculate parameter
vector<double> Q_Kirchhoff_sinkExcept;				   // ΣQ without sink node
vector<double> P_tubePressure_sinkExcept;			   // p_i(t) without sink node
vector<vector<double>> D_tubeThickness_deltaT;		   // D_ij(Δt)
vector<vector<double>> pressureCoefficient;			   // Coeff. of p(t) derived from the simultaneous equations in Eq.(1) and Eq.(2)
vector<vector<double>> pressureCoefficient_sinkExcept; // without sink node
// Sigmoid function
vector<vector<double>> Q_tubeFlow_sigmoidOutput; // The buffer for storing sigmoid func. solutions

// Pajek
vector<double> red;
vector<double> orange;
vector<double> green;
vector<double> blue;

//SOURCE and DIST
vector<double> Flow_SOURCE;
vector<int> SOURCE;
vector<int> DIST;

// Prototype declaration
void Initialize(int,int);
void VectorFree(int);

// Program Start
int main(int argc, char *argv[])
{
	printf("argc = %d\n", argc);

	for (int i = 0; i < argc; i += 1)
	{
		printf("argv[%d] = %s\n", i, argv[i]);
	}

	// Loop Variables
	static int i, j, k = 0;
	static int a, b = 0;

	// Non-vector parameters initialization
	static double Q_allFlow = 0.0;
	static double degeneracyEffect = 1.0;

	// Node
	vector<double> Flow_SOURCE = {};
	vector<int> SOURCE = {};
	vector<int> DIST = {3};
	static int node = atoi(argv[1]);
	static int node_except = node - DIST.size();

	// flag for DIST
	bool fig_DIST = false;

	// Iteration
	static int num_loop = 0; // Setting
	static int ct = 0;		 // Result

	// Thresholds
	static int test_iter = 10;
	static double eps = 1.0e-10; // Convergence detection

	// The file name is the simulation topology entered at compile.
	const char *NET_file = argv[2];

	// Vectors initialization
	Initialize(node, node_except);

	// Setting for Q_allFlow and num_loop
	int numbers_source;
	printf("Please set  The number of sources ");
	scanf("%d", &numbers_source);

	int sources;
	double flow_source;
	for(int i=0;i<numbers_source;i++){
		printf("Please set  The position of sources ");
		scanf("%d", &sources);
		SOURCE.push_back(sources);
		printf("Please set  The flow of sources ");
		scanf("%lf", &flow_source);
		Flow_SOURCE.push_back(flow_source);
		Q_allFlow += flow_source;
	}
	printf("Please set  The number of iterations: ");
	scanf("%d", &num_loop);

	// Create simulation topology for Pajek
	NodeConfigure(NET_file, node, SOURCE, DIST, D_tubeThickness, L_tubeLength);
    vector<vector<double>> capacity(node, vector<double>(node));
    #if NUM == 0
    capacity = {
    {INF, INF, INF, INF},
    {INF, INF, INF, INF},
    {INF, INF, INF, INF},
    {INF, INF, INF, INF}
    };
    #elif NUM == 1
    capacity = {
    {INF, 30.0, 20.0, INF},
    {30.0, INF, 10.0, 30.0},
    {20.0, 10.0, INF, 20.0},
    {INF, 30.0, 20.0, INF}
    };
    #endif  
	// PhysarumSolver Iteration Start
	while (ct != num_loop)
	{	
		// change situation
		/*if(ct==1500){
			SOURCE = {30,35};
			Flow_SOURCE={100,100};
			Q_allFlow = 0.0;
			for(int i=0;i<SOURCE.size();i++){
			Q_allFlow += flow_source;
			}
		}*/

		// auto start = chrono::system_clock::now();
        Output(ct,Q_tubeFlow);
		// Consider the case where the source and destination change during the loop.
		// The following process* is described in the loop.
		for(int i=0;i<SOURCE.size();i++){
			Q_Kirchhoff[SOURCE[i]] = Flow_SOURCE[i];
		}
		for(int i=0;i<DIST.size();i++){
			Q_Kirchhoff[DIST[i]] = Q_allFlow * NEG * ( 1.0 / DIST.size()); 
		}

		for (i = 0; i < node; i++)
		{
			pressureCoefficient[i][i] = 0.0;
			for(int j=0;j<SOURCE.size();j++){
				for(int k=0;k<DIST.size();k++){
					if (i == SOURCE[j] || i == DIST[k]){
						fig_DIST=true;
					}
				}
			}
			if (!fig_DIST){
				Q_Kirchhoff[i] = 0.0;
			}
			fig_DIST=false;
		}
		// The following process* end
		// Derive pressure gradient for each tube by simultaneous equations
		// Start
		for (i = 0; i < node; i++)
		{
			for (j = 0; j < node; j++)
			{
				if (L_tubeLength[i][j] != INF)
				{
					if (i != j)
					{
						pressureCoefficient[i][j] = D_tubeThickness[i][j] / L_tubeLength[i][j] * NEG;
					}
				}
			}
		}

		k = 0;
		for (i = 0; i < node; i++)
		{
			for (j = 0; j < node; j++)
			{
				if (L_tubeLength[i][j] != INF)
				{
					pressureCoefficient[k][k] += D_tubeThickness[i][j] / L_tubeLength[i][j];
				}
			}
			k++;
		}
		for (i = 0, a = 0; i < node, a < node_except; i++, a++)
		{	
			for(int k=0;k<DIST.size();k++){
				if (i == DIST[k] && DIST[k] != node)
				{
					fig_DIST=true;
				}
			}
			if(fig_DIST){
				i++;
			}
			fig_DIST=false;
			for (j = 0, b = 0; j < node, b < node_except; j++, b++)
			{

				for(int k=0;k<DIST.size();k++){
					if(j==DIST[k]){
						fig_DIST=true;
					}
				}
				if (fig_DIST)
				{
					j++;
					pressureCoefficient_sinkExcept[a][b] = pressureCoefficient[i][j];// it may happen an error especially in case multi sink and a and b.
				}
				else
				{
					pressureCoefficient_sinkExcept[a][b] = pressureCoefficient[i][j];
				}
				fig_DIST=false;
			}
		}
		a = 0;
		for (i = 0; i < node; i++)
		{
			for(int k=0;k<DIST.size();k++){
				if(i==DIST[k]){
					fig_DIST=true;
				}
			}

			if (!fig_DIST)
			{
				Q_Kirchhoff_sinkExcept[a] = Q_Kirchhoff[i];
				a++;
			}
			fig_DIST=false;
		}
		// ICCG:Incomplete Cholesky Conjugate Gradient method
		// This method can solve the sparse matrix linear system of equations inherent to the finite element method with low capacity and high speed.
		if (ICCG(pressureCoefficient_sinkExcept, Q_Kirchhoff_sinkExcept, P_tubePressure_sinkExcept, node_except, test_iter, eps) == 0)
		{
			break;
		}

		a = 0;
		for (i = 0; i < node; i++)
		{
			for(int j=0;j<DIST.size();j++){
				if(i==DIST[j]){
					fig_DIST=true;
				}
			}
			if (fig_DIST)
			{
				P_tubePressure[i] = 0.0;
			}
			else
			{
				P_tubePressure[i] = P_tubePressure_sinkExcept[a];
				a++;
			}
			fig_DIST=false;
		}


		// Derive pressure gradient for each tube by simultaneous equations
		// End

		// Derive the flow of the tube based on D, L, and Pressure gradient
		for (i = 0; i < node; i++)
		{
			for (j = 0; j < node; j++)
			{
				if (L_tubeLength[i][j] != INF)
				{
					Q_tubeFlow[i][j] = (D_tubeThickness[i][j] / L_tubeLength[i][j]) * (P_tubePressure[i] - P_tubePressure[j]);
				}
			}
		}

		// Sigmoid func.
		for (i = 0; i < node; i++)
		{
			for (j = 0; j < node; j++)
			{
				if (L_tubeLength[i][j] != INF)
				{
					Q_tubeFlow_sigmoidOutput[i][j] = (pow(fabs(Q_tubeFlow[i][j]), gamma)) / (1 + pow(fabs(Q_tubeFlow[i][j]), gamma));
				}
			}
		}

		// The thickness of tube after Δt seconds
		// You could even formulate a degeneracy effect
		for (i = 0; i < node; i++)
		{
			for (j = 0; j < node; j++)
			{
				if (L_tubeLength[i][j] != INF)
				{	
                    #if SIG == 0
                    D_tubeThickness_deltaT[i][j] = (fabs(Q_tubeFlow[i][j]) - (degeneracyEffect * D_tubeThickness[i][j])) * delta_time;//Q
                    #elif SIG == 1
                    D_tubeThickness_deltaT[i][j] = (Q_tubeFlow_sigmoidOutput[i][j] - (degeneracyEffect * D_tubeThickness[i][j])) * delta_time;//sig
                    #endif	
				}
			}
		}

		// Reflects tube thickness after Δt seconds
		for (i = 0; i < node; i++)
		{
			for (j = 0; j < node; j++)
			{
				{	
                    #if NUM == 0
                        D_tubeThickness[i][j] = D_tubeThickness[i][j] + (D_tubeThickness_deltaT[i][j]);
                    #elif NUM == 1
					    D_tubeThickness[i][j] = D_tubeThickness[i][j] + (D_tubeThickness_deltaT[i][j]) *  tanh((capacity[i][j] - fabs(Q_tubeFlow[i][j]))*coefficient_tanh);
                    #endif
				}			
			}
		}
		RGAP(ct,Q_tubeFlow);
		cout<<Q_tubeFlow[0][1]<<"aa"<<Q_tubeFlow[0][2]<<"aa"<<Q_tubeFlow[1][2]<<"aa"<<Q_tubeFlow[1][3]<<"aa"<<Q_tubeFlow[2][3]<<endl;
		// Write here the process you want to do every 'plot' times
		if ((ct + 1) % plot == 0)
		{
			SetTopologyColor(node, SOURCE, DIST, Q_tubeFlow, L_tubeLength, eps, Q_allFlow, ct);
		}

		ct++;

	} // while over

	VectorFree(node);
	return 0;
}

// vector initialize
void Initialize(int node, int node_except)
{

	// 1xN matrix
	// Basis parameter
	Q_Kirchhoff.resize(node);
	P_tubePressure.resize(node);

	// To calculate parameter
	Q_Kirchhoff_sinkExcept.resize(node_except);
	P_tubePressure_sinkExcept.resize(node_except);

	// Pajek

	// 2xN matrix
	// Basis parameter
	Q_tubeFlow.resize(node);
	D_tubeThickness.resize(node);
	L_tubeLength.resize(node);

	// To calculate parameter
	D_tubeThickness_deltaT.resize(node);
	pressureCoefficient.resize(node);
	pressureCoefficient_sinkExcept.resize(node_except);

	// Sigmoid function
	Q_tubeFlow_sigmoidOutput.resize(node);


	for (int i = 0; i < node; i++)
	{
		Q_tubeFlow[i].resize(node);
		D_tubeThickness[i].resize(node);
		L_tubeLength[i].resize(node);

		D_tubeThickness_deltaT[i].resize(node);
		pressureCoefficient[i].resize(node);

		Q_tubeFlow_sigmoidOutput[i].resize(node);
	}

	for (int i = 0; i < node_except; i++)
	{
		pressureCoefficient_sinkExcept[i].resize(node_except);
	}
}

// vector swap
void VectorFree(int node)
{

	// 1xN matrix
	// Source and DIST
	vector<double>().swap(Flow_SOURCE);
	vector<int>().swap(SOURCE);
	vector<int>().swap(DIST);

	// Basis parameter
	vector<double>().swap(Q_Kirchhoff);
	vector<double>().swap(P_tubePressure);

	// To calculate parameter
	vector<double>().swap(Q_Kirchhoff_sinkExcept);
	vector<double>().swap(P_tubePressure_sinkExcept);

	// Pajek

	// 2xN matrix
	// Basis parameter
	vector<vector<double>>().swap(Q_tubeFlow);
	vector<vector<double>>().swap(D_tubeThickness);
	vector<vector<double>>().swap(L_tubeLength);

	// To calculate parameter
	vector<vector<double>>().swap(D_tubeThickness_deltaT);
	vector<vector<double>>().swap(pressureCoefficient);
	vector<vector<double>>().swap(pressureCoefficient_sinkExcept);

	// Sigmoid function
	vector<vector<double>>().swap(Q_tubeFlow_sigmoidOutput);
}