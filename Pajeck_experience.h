#include "pch.h"

#define INF 10000.0
#define FILENAMESIZE 256

#define INIT_THICKNESS 0.5 //54 default
#define INIT_LENGTH 1.0
#define THRESHOLD_0 0.01
#define THRESHOLD_1 0.1
#define THRESHOLD_2 0.3
#define THRESHOLD_3 0.5
#define THRESHOLD_4 0.7
#define THRESHOLD_5 0.9

#define POW 1.0
#define MAX_FLOW 10.0

void NodeConfigure(const char *NET_file, int node, vector<int> &SOURCE, vector<int> &DIST, vector<vector<double>> &D_tubeThickness, vector<vector<double>> &L_tubeLength)
{
    int i, j = 0;
    FILE *fpN;
    bool fig_DIST=false;
    bool fig_SOURCE=false;
    vector<double> line = {0.05, 0.35, 0.35, 0.65, 0.65, 0.95};
    vector<double> row = {0.5, 0.75, 0.25, 0.75, 0.25, 0.5};
    vector<vector<double>> thickness = {
    {0, 0.5, 0.5, 0, 0, 0},
    {0.5, 0, 0.5, 0.5, 0.5, 0},
    {0.5, 0.5, 0, 0, 0.5, 0},
    {0, 0.5, 0, 0, 0.5, 0.5},
    {0, 0.5, 0.5, 0.5, 0, 0.5},
    {0, 0, 0, 0.5, 0.5, 0}
    };
    vector<vector<double>> length = {
    {INF, 1.0, 3.0, INF, INF, INF},
    {1.0, INF, 1.0, 2.0, 3.0, INF},
    {3.0, 1.0, INF, INF, 3.0, INF},
    {INF, 2.0, INF, INF, 4.0, 4.0},
    {INF, 3.0, 3.0, 4.0, INF, 1.0},
    {INF, INF, INF, 4.0, 1.0, INF}
    };

    if ((fpN = fopen(NET_file, "w")) == NULL)
    {
        printf("DATA FILE NO EXIST.\n");
        exit(1);
    }

    fprintf(fpN, "*Vertices	%d\n", node);

    for (i = 0; i < node; i++)
    {   
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
                    i + 1, i + 1, line[i], row[i] );
        }
        else
        {
            fprintf(fpN, "%d \"%d\" %.4lf %.4lf ic White\n",
                    i + 1, i + 1, line[i], row[i] );
        }
        fig_DIST=false;
        fig_SOURCE=false;
    }


    fprintf(fpN, "*Arcs\n");
    fprintf(fpN, "*Edges\n");

    fprintf(fpN, "1 2 1\n");
    fprintf(fpN, "1 3 1\n");
    fprintf(fpN, "2 3 1\n");
    fprintf(fpN, "2 4 1\n");
    fprintf(fpN, "2 5 1\n");
    fprintf(fpN, "3 5 1\n");
    fprintf(fpN, "4 5 1\n");
    fprintf(fpN, "4 6 1\n");
    fprintf(fpN, "5 6 1\n");

    for (i = 0; i < node; i++)
    {
        for (j = 0; j < node; j++)
        {   
            D_tubeThickness[i][j] = thickness[i][j];
            L_tubeLength[i][j] = length[i][j];
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

    vector<double> line = {0.05, 0.35, 0.35, 0.65, 0.65, 0.95};
    vector<double> row = {0.5, 0.75, 0.25, 0.75, 0.25, 0.5};
    bool fig_DIST=false;
    bool fig_SOURCE=false;
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
                    i + 1, i + 1, line[i], row[i] );
        }
        else
        {
            fprintf(fpN, "%d \"%d\" %.4lf %.4lf ic White\n",
                    i + 1, i + 1, line[i], row[i] );
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
                    if (0 < Q_tubeFlow[i][j] && Q_tubeFlow[i][j] <=THRESHOLD_0 )//閾値変更可能
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