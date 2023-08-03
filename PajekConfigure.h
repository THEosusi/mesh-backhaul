#include "pch.h"

#define INF 10000.0
#define FILENAMESIZE 256

#define INIT_THICKNESS 54.0
#define INIT_LENGTH 1.0

#define THRESHOLD_1 0.5
#define THRESHOLD_2 1.0

void NodeConfigure(const char *NET_file, int node, int SOURCE, int DIST, vector<double> &x_coordinate, vector<double> &y_coordinate, vector<vector<double>> &D_tubeThickness, vector<vector<double>> &L_tubeLength, vector<vector<double>> &node_distance)
{
    int i, j = 0;
    unsigned int x_coordSeed = 101;
    unsigned int y_coordSeed = 201;
    double maxDistance = sqrt(2); //(0,0) ~ (1,1) max distance = sqrt(2)
    FILE *fpN;

    srand(x_coordSeed);
    for (i = 0; i < node; i++)
    {
        x_coordinate[i] = (double)rand() / RAND_MAX;
        printf("x_coordinate[%d] is %lf\n", i, x_coordinate[i]);
    }

    srand(y_coordSeed);
    for (i = 0; i < node; i++)
    {
        y_coordinate[i] = (double)rand() / RAND_MAX;
        printf("y_coordinate[%d] is %lf\n", i, y_coordinate[i]);
    }

    if ((fpN = fopen(NET_file, "w")) == NULL)
    {
        printf("DATA FILE NO EXIST.\n");
        exit(1);
    }

    fprintf(fpN, "*Vertices	%d\n", node);

    for (i = 0; i < node; i++)
    {
        if (i == SOURCE)
        {
            fprintf(fpN, "%d \"%d\" %.4lf %.4lf ic Black\n",
                    i + 1, i + 1, x_coordinate[i], y_coordinate[i]);
        }
        else if (i == DIST)
        {
            fprintf(fpN, "%d \"%d\" %.4lf %.4lf ic Black\n",
                    i + 1, i + 1, x_coordinate[i], y_coordinate[i]);
        }
        else
        {
            fprintf(fpN, "%d \"%d\" %.4lf %.4lf ic White\n",
                    i + 1, i + 1, x_coordinate[i], y_coordinate[i]);
        }
    }

    fprintf(fpN, "*Arcs\n");
    fprintf(fpN, "*Edges\n");
    for (i = 0; i < node; i++)
    {
        for (j = 0; j < node; j++)
        {
            if (i != j)
            {
                node_distance[i][j] = pow((x_coordinate[j] - x_coordinate[i]) * (x_coordinate[j] - x_coordinate[i]) + (y_coordinate[j] - y_coordinate[i]) * (y_coordinate[j] - y_coordinate[i]), 0.5);
            }
        }
    }

    for (i = 0; i < node; i++)
    {
        for (j = 0; j < node; j++)
        {
            if (i == j)
            {
                D_tubeThickness[i][j] = 0;
                L_tubeLength[i][j] = INF;
            }
            else
            {
                if (0.0 < node_distance[i][j] & node_distance[i][j] <= THRESHOLD_1)
                {
                    if (i != j)
                    {
                        fprintf(fpN, "%d %d 1\n", i + 1, j + 1);
                        D_tubeThickness[i][j] = D_tubeThickness[j][i] = INIT_THICKNESS;
                        L_tubeLength[i][j] = L_tubeLength[j][i] = INIT_LENGTH;
                    }
                }
                else if (THRESHOLD_1 < node_distance[i][j] & node_distance[i][j] <= THRESHOLD_2)
                {
                    if (i != j)
                    {
                        fprintf(fpN, "%d %d 1\n", i + 1, j + 1);
                        D_tubeThickness[i][j] = D_tubeThickness[j][i] = INIT_THICKNESS;
                        L_tubeLength[i][j] = L_tubeLength[j][i] = INIT_LENGTH;
                    }
                }
                else if (THRESHOLD_2 < node_distance[i][j] & node_distance[i][j] <= maxDistance)
                {
                    if (i != j)
                    {
                        fprintf(fpN, "%d %d 1\n", i + 1, j + 1);
                        D_tubeThickness[i][j] = D_tubeThickness[j][i] = INIT_THICKNESS;
                        L_tubeLength[i][j] = L_tubeLength[j][i] = INIT_LENGTH;
                    }
                }
                else if (maxDistance < node_distance[i][j])
                {
                    L_tubeLength[i][j] = L_tubeLength[j][i] = INF;
                }
            }
        }
    }

    fclose(fpN);
}

void SetTopologyColor(int node, int SOURCE, int DIST, vector<double> x_coordinate, vector<double> y_coordinate, vector<vector<double>> Q_tubeFlow, vector<vector<double>> L_tubeLength, double eps, double Q_allFlow, int ct, vector<vector<double>> node_distance)
{

    int i, j = 0;
    unsigned int all_path_count = 0;
    unsigned int red_path_count = 0;
    unsigned int blue_path_count = 0;
    unsigned int green_path_count = 0;
    FILE *fpN;

    char filename[FILENAMESIZE];

    sprintf(filename, "./test_topology_%d.net", ct + 1);

    if ((fpN = fopen(filename, "w")) == NULL)
    {
        printf("DATA FILE DOES NOT EXIST.\n");
        exit(1);
    }
    fprintf(fpN, "*Vertices	%d\n", node);
    for (i = 0; i < node; i++)
    {
        if (i == SOURCE)
        {
            fprintf(fpN, "%d \"%d\" %.4lf %.4lf ic Black\n",
                    i + 1, i + 1, x_coordinate[i], y_coordinate[i]);
        }
        else if (i == DIST)
        {
            fprintf(fpN, "%d \"%d\" %.4lf %.4lf ic Black\n",
                    i + 1, i + 1, x_coordinate[i], y_coordinate[i]);
        }
        else
        {
            fprintf(fpN, "%d \"%d\" %.4lf %.4lf ic White\n",
                    i + 1, i + 1, x_coordinate[i], y_coordinate[i]);
        }
    }
    fprintf(fpN, "*Arcs\n");
    fprintf(fpN, "*Edges\n");

    for (i = 0; i < node; i++)
    {
        for (j = 0; j < node; j++)
        {
            if (L_tubeLength[i][j] != INF)
            {
                if (0 < node_distance[i][j] && node_distance[i][j] <= sqrt(2))
                {
                    if (0 < Q_tubeFlow[i][j] && Q_tubeFlow[i][j] <= eps)
                    {   
                        all_path_count++;
                    }
                    else if (eps < Q_tubeFlow[i][j] && Q_tubeFlow[i][j] <= THRESHOLD_1)
                    {
                        fprintf(fpN, "%d %d 1 c Blue\n", i + 1, j + 1);
                        blue_path_count++;
                        all_path_count++;
                    }
                    else if (THRESHOLD_1 < Q_tubeFlow[i][j] && Q_tubeFlow[i][j] <= THRESHOLD_2)
                    {
                        fprintf(fpN, "%d %d 2 c Green\n", i + 1, j + 1);
                        green_path_count++;
                        all_path_count++;
                    }
                    else if (THRESHOLD_2 < Q_tubeFlow[i][j] && Q_tubeFlow[i][j] <= Q_allFlow)
                    {
                        fprintf(fpN, "%d %d 3 c Red\n", i + 1, j + 1);
                        red_path_count++;
                        all_path_count++;
                    }
                }
            }
        }
    }

    fclose(fpN);
}