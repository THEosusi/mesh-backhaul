#include<iostream>
#include<fstream>

void RGAP(int ct,vector<vector<double>> &Q_tubeFlow)
{   //not calculate you need put
    double keisuu=1.0;
    //double ideal=210.0*keisuu;
    double ideal=160.0;                   //multi
    //double ideal=20.0;
    ofstream outputfile("2dessinkRGAP40.txt", ios::app);
    if(ct==0){
    outputfile<<"#1.txt"<<endl;
    }
    double output= abs(1 - ideal/abs((Q_tubeFlow[0][1]*1)+ abs(Q_tubeFlow[0][2]*3)+ abs(Q_tubeFlow[1][2]*1)+abs(Q_tubeFlow[1][3]*5)+abs(Q_tubeFlow[2][3]*2)));
    //double output= abs(1 - (Q_tubeFlow[0][1]*1+ Q_tubeFlow[0][2]*3+ Q_tubeFlow[1][2]*1+Q_tubeFlow[1][3]*5+Q_tubeFlow[2][3]*2)/ideal);
    outputfile<<ct<<" "<<output<<endl;
    outputfile.close();
}