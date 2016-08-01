#include <iostream>
#include <vector>
#include <fstream>
#include "C90-graph.h"

#define MAX_TIMESTEP 1000

using namespace std;


int main()
{
	sexualnetwork Network;
	Network.initiation(1000);				// the size of population 

	string OutputOrderFileName = "order.dat";
	Network.fout_orders(ofstream(OutputOrderFileName));

	int timestep = 0;
	while (timestep < MAX_TIMESTEP)
	{
		Network.timestep();
		timestep++;
	}

	string OutputDotfileName = "test.dot";
	Network.dotoutput(ofstream(OutputDotfileName));

	string OutputFileName = "SIRnum.dat";
	ofstream SIRofs = ofstream(OutputFileName);

	for (int i = 0; i < Network.snapshots_of_SIR_num.size(); i++)
	{
		SIRofs << Network.snapshots_of_SIR_num[i][0] << "\t" << Network.snapshots_of_SIR_num[i][1] << "\t" << Network.snapshots_of_SIR_num[i][2] << endl;
	}

	return 0;
}