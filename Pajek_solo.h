#include "pch.h"

#define INF 10000.0
#define FILENAMESIZE 256

#define INIT_THICKNESS 54.0
#define INIT_LENGTH 1.0

#define THRESHOLD_1 0.5
#define THRESHOLD_2 1.0

void NodeConfigure(const char *NET_file, int node, vector<int> &SOURCE, vector<int> &DIST, vector<double> &x_coordinate, vector<double> &y_coordinate, vector<vector<double>> &D_tubeThickness, vector<vector<double>> &L_tubeLength, vector<vector<double>> &node_distance)
{
    int i, j = 0;
    FILE *fpN;
    bool fig_DIST=false;
    bool fig_SOURCE=false;


    if ((fpN = fopen(NET_file, "w")) == NULL)
    {
        printf("DATA FILE NO EXIST.\n");
        exit(1);
    }

    fprintf(fpN, "*Vertices	%d\n", node);

    for (i = 0; i < node; i++)
    {   
        double line=(double)(i%8)/10.0 + 0.15 ;
        double row=(double)(i/8)/10.0 + 0.15 ;
        for(int j=0;j<DIST.size();j++){
            if(i==DIST[j]){
                fig_DIST=true;
            }
        }
        for(int j=0;j<SOURCE.size();j++){
            if(i==SOURCE[j]){
                fig_SOURCE=true;
            }
        }        
        if (fig_SOURCE || fig_DIST)
        {
            fprintf(fpN, "%d \"%d\" %.4lf %.4lf ic Black\n",
                    i + 1, i + 1, line, row );
        }
        else
        {
            fprintf(fpN, "%d \"%d\" %.4lf %.4lf ic White\n",
                    i + 1, i + 1, line, row );
        }
        fig_DIST=false;
        fig_SOURCE=false;
    }


    fprintf(fpN, "*Arcs\n");
    fprintf(fpN, "*Edges\n");
    for (i = 0; i < node; i++)
    {
        for (j = i; j < node; j++)
        {   
            if (i == j)
            {
                D_tubeThickness[i][j] = 0;
                L_tubeLength[i][j] = INF;

            }else if(j-i == 1 && j%8 != 0){
                fprintf(fpN, "%d %d 1\n", i + 1, j + 1);
                D_tubeThickness[i][j] = D_tubeThickness[j][i] = INIT_THICKNESS;
                double PoptoDisti =  abs(((int)DIST[0])/8- i/8) + abs(((int)DIST[0]%8) - i%8) + 1;// due to prevent 0
                double PoptoDistj =  abs(((int)DIST[0])/8- j/8) + abs(((int)DIST[0]%8) - j%8) + 1;                
                for(int k=0;k<DIST.size();k++){
                    double tem_PoptoDisti =  abs(((int)DIST[k])/8- i/8) + abs(((int)DIST[k]%8) - i%8) + 1;
                    double tem_PoptoDistj =  abs(((int)DIST[k])/8- j/8) + abs(((int)DIST[k]%8) - j%8) + 1;
                    if(tem_PoptoDisti<PoptoDisti){
                        PoptoDisti=tem_PoptoDisti;
                    }
                    if(tem_PoptoDistj<PoptoDistj){
                        PoptoDistj=tem_PoptoDistj;
                    }
                }
                L_tubeLength[i][j] =  pow(PoptoDistj,1);
                L_tubeLength[j][i] =  pow(PoptoDisti,1);


            }else if(j-i == 8 && (i<56 || j<56)){
                fprintf(fpN, "%d %d 1\n", i + 1, j + 1);
                double PoptoDisti =  abs(((int)DIST[0])/8- i/8) + abs(((int)DIST[0]%8) - i%8) + 1;// due to prevent 0
                double PoptoDistj =  abs(((int)DIST[0])/8- j/8) + abs(((int)DIST[0]%8) - j%8) + 1;                
                for(int k=0;k<DIST.size();k++){
                    double tem_PoptoDisti =  abs(((int)DIST[k])/8- i/8) + abs(((int)DIST[k]%8) - i%8) + 1;
                    double tem_PoptoDistj =  abs(((int)DIST[k])/8- j/8) + abs(((int)DIST[k]%8) - j%8) + 1;
                    if(tem_PoptoDisti<PoptoDisti){
                        PoptoDisti=tem_PoptoDisti;
                    }
                    if(tem_PoptoDistj<PoptoDistj){
                        PoptoDistj=tem_PoptoDistj;
                    }
                }
                D_tubeThickness[i][j] = D_tubeThickness[j][i] = INIT_THICKNESS;
                L_tubeLength[i][j] =  pow(PoptoDistj,1);
                L_tubeLength[j][i] =  pow(PoptoDistj,1);
            }else{
                L_tubeLength[i][j] = L_tubeLength[j][i] = INF;
            }
        }
    }

    fclose(fpN);
}

void SetTopologyColor(int node, vector<int> &SOURCE, vector<int> &DIST, vector<double> x_coordinate, vector<double> y_coordinate, vector<vector<double>> Q_tubeFlow, vector<vector<double>> L_tubeLength, double eps, double Q_allFlow, int ct, vector<vector<double>> node_distance)
{

    int i, j = 0;
    unsigned int all_path_count = 0;
    unsigned int red_path_count = 0;
    unsigned int blue_path_count = 0;
    unsigned int green_path_count = 0;
    FILE *fpN;
    bool fig_DIST = false;
    bool fig_SOURCE = false;

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
        double line=(double)(i%8)/10.0 + 0.15;
        double row=(double)(i/8)/10.0 + 0.15;
        for(int j=0;j<DIST.size();j++){
            if(i==DIST[j]){
                fig_DIST=true;
            }
        }
        for(int j=0;j<SOURCE.size();j++){
            if(i==SOURCE[j]){
                fig_SOURCE=true;
            }
        }        
        if (fig_SOURCE || fig_DIST )
        {
            fprintf(fpN, "%d \"%d\" %.4lf %.4lf ic Black\n",
                    i + 1, i + 1, line, row );
        }
        else
        {
            fprintf(fpN, "%d \"%d\" %.4lf %.4lf ic White\n",
                    i + 1, i + 1, line, row );
        }
        fig_DIST=false;
        fig_SOURCE=false;
    }
    fprintf(fpN, "*Arcs\n");
    fprintf(fpN, "*Edges\n");

    for (i = 0; i < node; i++)
    {
        for (j = 0; j < node; j++)
        {
            if (L_tubeLength[i][j] != INF)
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

    fclose(fpN);
}