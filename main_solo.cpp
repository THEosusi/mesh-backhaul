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
#include "Pajek_solo.h"

#define INF 10000.0		// Infinite length representation. Exclude L from the calculation, where L==INF
#define NEG -1.0		// For representation of negative
#define gamma 1.5		// Slope of the sigmoid function
#define delta_time 0.01 //Δt scale
#define plot 1			// Number of loops you want to plot


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
vector<double> x_coordinate;
vector<double> y_coordinate;
vector<vector<double>> node_distance;
vector<double> red;
vector<double> orange;
vector<double> green;
vector<double> blue;

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

	// Node
	static int SOURCE = 44;
	vector<int> DIST = {18,37,60};
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

	// Non-vector parameters initialization
	static double Q_allFlow = 1.0;
	static double degeneracyEffect = 1.0;

	// The file name is the simulation topology entered at compile.
	const char *NET_file = argv[2];

	// Vectors initialization
	Initialize(node, node_except);

	// Create simulation topology for Pajek
	NodeConfigure(NET_file, node, SOURCE, DIST, x_coordinate, y_coordinate, D_tubeThickness, L_tubeLength, node_distance);

	// Setting for Q_allFlow and num_loop
	printf("Please set  The value to be distributed by PS: ");
	scanf("%lf", &Q_allFlow);
	printf("Please set  The number of iterations: ");
	scanf("%d", &num_loop);

	// PhysarumSolver Iteration Start
	while (ct != num_loop)
	{
		// auto start = chrono::system_clock::now();

		// Consider the case where the source and destination change during the loop.
		// The following process* is described in the loop.
		Q_Kirchhoff[SOURCE] = Q_allFlow;
		for(int i=0;i<DIST.size();i++){
			Q_Kirchhoff[DIST[i]] = Q_allFlow * NEG * ( 1.0 / DIST.size()); // SOUCE LEN is 1 あとでへんこう
		}
		for (i = 0; i < node; i++)
		{
			pressureCoefficient[i][i] = 0.0;
			for(int j=0;j<DIST.size();j++){
				if (i == SOURCE || i == DIST[j]){
					fig_DIST=true;
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

					D_tubeThickness_deltaT[i][j] = (Q_tubeFlow_sigmoidOutput[i][j] - (degeneracyEffect * D_tubeThickness[i][j])) * delta_time;
				}
			}
		}

		// Reflects tube thickness after Δt seconds
		for (i = 0; i < node; i++)
		{
			for (j = 0; j < node; j++)
			{
				if (L_tubeLength[i][j] != INF)
				{
					D_tubeThickness[i][j] = D_tubeThickness[i][j] + D_tubeThickness_deltaT[i][j];
				}
			}
		}

		//cout<<Q_tubeFlow[45][37]<<"aa"<<Q_tubeFlow[36][37]<<"aa"<<Q_tubeFlow[29][37]<<"aa"<<Q_tubeFlow[38][37]<<endl;
		//cout<<Q_tubeFlow[26][18]<<"bb"<<Q_tubeFlow[17][18]<<"bb"<<Q_tubeFlow[10][18]<<"bb"<<Q_tubeFlow[19][18]<<endl;
		/*for(int i=0;i<node;i++){
			cout<<i<<"aa"<<Q_Kirchhoff[i]<<endl;
		}*/
		//cout<<Q_tubeFlow[52][60]<<"aaa"<<Q_tubeFlow[59][60]<<"aaa"<<Q_tubeFlow[61][60]<<endl;
		// Write here the process you want to do every 'plot' times
		if ((ct + 1) % plot == 0)
		{
			SetTopologyColor(node, SOURCE, DIST, x_coordinate, y_coordinate, Q_tubeFlow, L_tubeLength, eps, Q_allFlow, ct, node_distance);
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
	// vector<double> x_coordinate(node);
	x_coordinate.resize(node);
	y_coordinate.resize(node);

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

	node_distance.resize(node);

	for (int i = 0; i < node; i++)
	{
		Q_tubeFlow[i].resize(node);
		D_tubeThickness[i].resize(node);
		L_tubeLength[i].resize(node);

		D_tubeThickness_deltaT[i].resize(node);
		pressureCoefficient[i].resize(node);

		Q_tubeFlow_sigmoidOutput[i].resize(node);

		node_distance[i].resize(node);
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
	// Basis parameter
	vector<double>().swap(Q_Kirchhoff);
	vector<double>().swap(P_tubePressure);

	// To calculate parameter
	vector<double>().swap(Q_Kirchhoff_sinkExcept);
	vector<double>().swap(P_tubePressure_sinkExcept);

	// Pajek
	vector<double>().swap(x_coordinate);
	vector<double>().swap(y_coordinate);

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

	vector<vector<double>>().swap(node_distance);
}