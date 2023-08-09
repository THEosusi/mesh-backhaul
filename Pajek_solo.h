#include "pch.h"

#define INF 10000.0
#define FILENAMESIZE 256

#define INIT_THICKNESS 0.5 //54 default
#define INIT_LENGTH 1.0

#define THRESHOLD_1 0.1
#define THRESHOLD_2 0.3
#define THRESHOLD_3 0.5
#define THRESHOLD_4 0.7
#define THRESHOLD_5 0.9

#define POW 1.0
#define MAX_FLOW 54.0 //depending on the connection

void NodeConfigure(const char *NET_file, int node, vector<int> &SOURCE, vector<int> &DIST, vector<vector<double>> &D_tubeThickness, vector<vector<double>> &L_tubeLength)
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
                double PoptoDistij = (PoptoDisti + PoptoDistj)/2; // case 2 idea            
                for(int k=0;k<DIST.size();k++){
                    double tem_PoptoDisti =  abs(((int)DIST[k])/8- i/8) + abs(((int)DIST[k]%8) - i%8) + 1;
                    double tem_PoptoDistj =  abs(((int)DIST[k])/8- j/8) + abs(((int)DIST[k]%8) - j%8) + 1;
                    double tem_PoptoDistij = (tem_PoptoDisti + tem_PoptoDistj)/2; //case 2
                    if (tem_PoptoDistij<PoptoDistij){
                        PoptoDistij=tem_PoptoDistij;
                    }
                    /*if(tem_PoptoDisti<PoptoDisti){
                        PoptoDisti=tem_PoptoDisti;
                    }
                    if(tem_PoptoDistj<PoptoDistj){
                        PoptoDistj=tem_PoptoDistj;
                    }*/
                }
                //L_tubeLength[i][j] =  pow(PoptoDistj,POW);
                //L_tubeLength[j][i] =  pow(PoptoDisti,POW);
                L_tubeLength[j][i] = L_tubeLength[i][j] = pow(PoptoDistij,POW);

            }else if(j-i == 8 && (i<56 || j<56)){
                fprintf(fpN, "%d %d 1\n", i + 1, j + 1);
                D_tubeThickness[i][j] = D_tubeThickness[j][i] = INIT_THICKNESS;
                double PoptoDisti =  abs(((int)DIST[0])/8- i/8) + abs(((int)DIST[0]%8) - i%8) + 1;// due to prevent 0
                double PoptoDistj =  abs(((int)DIST[0])/8- j/8) + abs(((int)DIST[0]%8) - j%8) + 1;
                double PoptoDistij = (PoptoDisti + PoptoDistj)/2; // case 2 idea               
                for(int k=0;k<DIST.size();k++){
                    double tem_PoptoDisti =  abs(((int)DIST[k])/8- i/8) + abs(((int)DIST[k]%8) - i%8) + 1;
                    double tem_PoptoDistj =  abs(((int)DIST[k])/8- j/8) + abs(((int)DIST[k]%8) - j%8) + 1;
                    double tem_PoptoDistij = (tem_PoptoDisti + tem_PoptoDistj)/2; //case 2
                    if (tem_PoptoDistij<PoptoDistij){
                        PoptoDistij=tem_PoptoDistij;
                    }
                    /*if(tem_PoptoDisti<PoptoDisti){
                        PoptoDisti=tem_PoptoDisti;
                    }
                    if(tem_PoptoDistj<PoptoDistj){
                        PoptoDistj=tem_PoptoDistj;
                    }*/
                }
                //L_tubeLength[i][j] =  pow(PoptoDistj,POW);
                //L_tubeLength[j][i] =  pow(PoptoDisti,POW);
                L_tubeLength[j][i] = L_tubeLength[i][j] = pow(PoptoDistij,POW);
            }else{
                L_tubeLength[i][j] = L_tubeLength[j][i] = INF;
            }
        }
    }

    fclose(fpN);
}

void SetTopologyColor(int node, vector<int> &SOURCE, vector<int> &DIST, vector<vector<double>> Q_tubeFlow, vector<vector<double>> L_tubeLength, double eps, double Q_allFlow, int ct)
{

    int i, j = 0;
    unsigned int all_path_count = 0;
    unsigned int blue_path_count = 0;
    unsigned int cyan_path_count = 0;
    unsigned int yellow_path_count = 0;
    unsigned int orange_path_count = 0;
    unsigned int red_path_count = 0;
    unsigned int redviolet_path_count = 0;
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
                    else if ( Q_tubeFlow[i][j] <= THRESHOLD_1 * MAX_FLOW)
                    {
                        fprintf(fpN, "%d %d 1 c Blue\n", i + 1, j + 1);
                        blue_path_count++;
                        all_path_count++;
                    }
                    else if ( Q_tubeFlow[i][j] <= THRESHOLD_2 * MAX_FLOW)
                    {
                        fprintf(fpN, "%d %d 2 c Cyan\n", i + 1, j + 1);
                        cyan_path_count++;
                        all_path_count++;
                    }
                    else if ( Q_tubeFlow[i][j] <= THRESHOLD_3 * MAX_FLOW)
                    {
                        fprintf(fpN, "%d %d 3 c Yellow\n", i + 1, j + 1);
                        yellow_path_count++;
                        all_path_count++;
                    }
                    else if ( Q_tubeFlow[i][j] <= THRESHOLD_4 * MAX_FLOW)
                    {
                        fprintf(fpN, "%d %d 3 c Orange\n", i + 1, j + 1);
                        orange_path_count++;
                        all_path_count++;
                    }
                    else if ( Q_tubeFlow[i][j] <= THRESHOLD_5 * MAX_FLOW)
                    {
                        fprintf(fpN, "%d %d 3 c Red\n", i + 1, j + 1);
                        red_path_count++;
                        all_path_count++;
                    }
                    else if ( THRESHOLD_5 * MAX_FLOW < Q_tubeFlow[i][j])
                    {
                        fprintf(fpN, "%d %d 3 c RedViolet\n", i + 1, j + 1);
                        redviolet_path_count++;
                        all_path_count++;
                    }                     
            }
        }
    }

    fclose(fpN);
}