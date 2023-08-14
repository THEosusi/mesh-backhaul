#include<iostream>
#include<fstream>

void Output(int ct,vector<vector<double>> &Q_tubeFlow)
{   
    ofstream outputfile("test.txt", ios::app);
    if(ct==0){
    outputfile<<"#1.txt"<<endl;
    }
    outputfile<<ct<<" "<<Q_tubeFlow[0][1]<<" "<<Q_tubeFlow[0][2]<<" "<<Q_tubeFlow[1][2]<<" "<<Q_tubeFlow[1][3]<<" "<<Q_tubeFlow[2][3]<<endl;
    outputfile.close();
}