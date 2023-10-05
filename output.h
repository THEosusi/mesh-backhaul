#include<iostream>
#include<fstream>

void Output(int ct,vector<vector<double>> &Q_tubeFlow)
{   
    ofstream outputfile("2dessintest40.txt", ios::app);
    if(ct==0){
    outputfile<<"#1.txt"<<endl;
    }
    outputfile<<ct<<" "<<Q_tubeFlow[0][1]<<" "<<Q_tubeFlow[0][2]<<" "<<Q_tubeFlow[1][2]<<" "<<Q_tubeFlow[1][3]<<" "<<Q_tubeFlow[2][3]<<endl;
    //path0 // path1// path2 //path3 //path 4
    outputfile.close();
}