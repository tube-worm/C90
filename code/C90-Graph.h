#include <vector>
#include <array>
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;

class node
{
public:
	int id_;
	int order_;
	vector<node*> adj_pointer_vector;

	/* parameter of infection */
	bool changed_at_this_moment;
	int state_of_infection;		// 0-Suceptive, 1-Infected, 2-Recovered
	double rate_of_infected;
	double rate_of_recovery;

	node()
	{
		order_ = 0;
		changed_at_this_moment = false;
		state_of_infection = 0;
	}

	// doubling with no-argument constructer 
	node(int i)
	{
		id_ = i;
		order_ = 0;
		changed_at_this_moment = false;
		state_of_infection = 0;
	}

	void link(node& n)
	{
		adj_pointer_vector.push_back(&n);
		order_++;
	}

	void infectious_targeted()
	{
		if ((changed_at_this_moment == false) && (state_of_infection == 0))
		{
			double r = (double)rand() / RAND_MAX;
			if (r < rate_of_infected)
			{
				changed_at_this_moment = true;
				state_of_infection = 1;
			}
		}
	}

	void infectious_contact()
	{
		for (int i = 0; i < adj_pointer_vector.size(); i++)
		{
			adj_pointer_vector[i]->infectious_targeted();
		}
	}

	void recovery_from_infection()
	{
		double r = (double)rand() / RAND_MAX;
		if (r < rate_of_recovery)
		{
			changed_at_this_moment = true;
			state_of_infection = 2;
		}
	}
};

void annual_link(node& n1, node& n2)
{
	n1.link(n2);
	n2.link(n1);
}

// TODO: initiation
class population
{
public:
	node* array_of_nodes;	// array of the nodes
	int size_;

	population()
	{
		size_ = 0;
	}

	population(int n)
	{
		size_ = n;
		array_of_nodes = new node[n];
	}
};

class scalefreenetwork : public population
{
private:
	const static int initial_num_of_nodes = 2;
	int sum_of_order;

public:
	scalefreenetwork()
	{
		sum_of_order = 0;
	}

	// using Barabasi-Albert algorithm
	void networkconstruction(int NumOfNodes)
	{
		array_of_nodes = new node[NumOfNodes];

		/***** start with complete graph *****/
		for (int i = 0; i < initial_num_of_nodes; i++)
			array_of_nodes[i].order_ = 0;
		for (int i = 0; i < initial_num_of_nodes; i++)
		{
			array_of_nodes[i].id_ = i;
			for (int j = i + 1; j < initial_num_of_nodes; j++)
				annual_link(array_of_nodes[i], array_of_nodes[j]);
		}

		for (int i = 0; i < initial_num_of_nodes; i++)
			sum_of_order += array_of_nodes[i].order_;
		size_ = initial_num_of_nodes;
		/***** end of making complete graph *****/

		/* add the nodes from N0-th to N-th */
		for (int i = initial_num_of_nodes; i < NumOfNodes; i++)
		{
			array_of_nodes[i].order_ = 0;
			array_of_nodes[i].id_ = i;

			int r = rand() % sum_of_order;

			/* put a random value r (0~sum_of_order) in the buckets of order */
			int node_id_to_make_edge_with = -1;
			while (r >= 0)
			{
				node_id_to_make_edge_with++;
				r -= array_of_nodes[node_id_to_make_edge_with].order_;
			}

			annual_link(array_of_nodes[node_id_to_make_edge_with], array_of_nodes[i]);
			sum_of_order += 2;
			size_++;

			// progress report
			if (i % 100 == 0)
				cout << i << "is done!" << endl;
		}
	}
};

class sexualnetwork : scalefreenetwork
{
private:
	static const int max_order = 100;
	static const int num_of_initial_infected = 10;
	int initial_infected_node_id[num_of_initial_infected];	// initiation: start with 10 people infected

	const int threshold_of_order = 4;
	const double infectious_rate_of_normal = 0.2;
	const double infectious_rate_of_hub = 0.01;
	const double recovery_rate = 0.01;

public:
	vector<array<int, 3>> snapshots_of_SIR_num;		// timelapse snapshot vector of the numbers of 0:Suceptive, 1:Infected, and 2:Recovered people

	void initiation(int NumOfNodes)
	{
		networkconstruction(NumOfNodes);

		for (int i = 0; i < size_; i++)
		{
			array_of_nodes[i].rate_of_recovery = recovery_rate;
			if (array_of_nodes[i].order_ > threshold_of_order)
			{
				array_of_nodes[i].rate_of_infected = infectious_rate_of_hub;
			}
			else
			{
				array_of_nodes[i].rate_of_infected = infectious_rate_of_normal;
			}
		}

		for (int i = 0; i < num_of_initial_infected; i++)
		{
			initial_infected_node_id[i] = rand() % size_;
			array_of_nodes[initial_infected_node_id[i]].state_of_infection = 1;
		}

		array<int, 3> array_of_SIR_num = { size_ - num_of_initial_infected, num_of_initial_infected, 0 };

		snapshots_of_SIR_num.push_back(array_of_SIR_num);
	}

	void timestep()
	{
		for (int i = 0; i < size_; i++)
		{
			int state_of_this_node = array_of_nodes[i].state_of_infection;
			if (!array_of_nodes[i].changed_at_this_moment && (state_of_this_node == 1))
			{
				array_of_nodes[i].recovery_from_infection();
				array_of_nodes[i].infectious_contact();
			}
		}

		// end process
		array<int, 3> array_of_SIR_num = {};

		for (int i = 0; i < size_; i++)
		{
			array_of_nodes[i].changed_at_this_moment = false;
			array_of_SIR_num[array_of_nodes[i].state_of_infection] ++;
		}

		snapshots_of_SIR_num.push_back(array_of_SIR_num);

		
	}

	// file output of a histgram of orders
	void fout_orders(ofstream ofs)
	{
		int *array_of_order = new int[size_];
		int *orderhistgram = new int[max_order];

		for (int i = 0; i < size_; i++)
		{
			array_of_order[i] = array_of_nodes[i].order_;
			orderhistgram[array_of_order[i]] ++;
		}
		sort(array_of_order, array_of_order + size_);

		for (int i = 1; i < max_order; i++)
		{
			ofs << orderhistgram[i] << endl;
		}
	}

	// DOT format file output for GraphViz
	void dotoutput(ofstream ofs)
	{
		// header of DOT file
		ofs << "graph graph_name{ \n graph[\n" << endl;
		ofs << "dpi = \"480\", \n layout = sfdp\n ];" << endl;
		ofs << "node[\n shape = point,\n fixedsize = true, \n width = 0.02\n]; " << endl;

		for (int i = 0; i < size_; i++)
		{
			if (array_of_nodes[i].state_of_infection == 0)		// Suceptive: blue
			{
				ofs << i << "[fillcolor = \"#0000FF\",\n color = \"#0000FF\"]\n " << endl;
			}
			else if (array_of_nodes[i].state_of_infection == 1)	// Infected: red
			{
				ofs << i << "[fillcolor = \"#FF0000\",\n color = \"#FF0000\"]\n " << endl;
			}
			else												// Recovered: green
			{	
				ofs << i << "[fillcolor = \"#00FF00\",\n color = \"#00FF00\"]\n " << endl;
			}
			
			for (int j = 0; j < array_of_nodes[i].order_; j++)
			{
				if (i < array_of_nodes[i].adj_pointer_vector[j]->id_)
				{
					ofs << i << "--" << array_of_nodes[i].adj_pointer_vector[j]->id_ << ";" << endl;
				}
			}
		}
		ofs << "}" << endl;
	}
};


