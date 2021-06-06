#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <unistd.h>
#include <cmath>

using namespace std;

int main (int argc, char * argv [])
{
	string static_graph_file_name = argv[1];
	double static_reinsertion_param = atof(argv[2]);
	int static_reinsertion_cap = atoi(argv[3]);
	double general_insertion_param = atof(argv[4]);
	int general_insertion_cap = atoi(argv[5]);
	int num_general_inserts = atoi(argv[6]);
	string stream_file_name = argv[7];

	ifstream static_graph_file{static_graph_file_name};
	long num_nodes, num_edges;
	static_graph_file >> num_nodes >> num_edges;

	ofstream stream_file_out{stream_file_name}; 

	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator(seed);
	uniform_int_distribution<int> rand_char(33, 126);

	int rand_prefix_len = log((double) num_edges * (2 * 1 / static_reinsertion_param + 1) 
			+ num_general_inserts * (2 * 1 / general_insertion_param)) / log(10.0);

	long total_num_updates = 0;	
	char * random_prefix = new char [rand_prefix_len];

	// Static graph geometric inserts/deletes	
	geometric_distribution<int> static_reinsertions(static_reinsertion_param);
	long u, v;
	int num_reinserts, num_updates;
	for (long i = 0; i < num_edges; i++)
	{
		static_graph_file >> u >> v;
		num_reinserts = min(static_reinsertion_cap, static_reinsertions(generator));
		// Want edges to be contained in final state of graph
		num_updates = 2 * num_reinserts + 1;

		for (int j = 0; j < num_updates; j++)
		{	
			for (int k = 0; k < rand_prefix_len; k++)
				random_prefix[k] = (char) rand_char(generator);
		
			stream_file_out << random_prefix << "\t" << u << "\t" << v << "\n";
		}

		total_num_updates += num_updates;
	}	

	// General geometric inserts/deletes	
	geometric_distribution<int> general_insertions(general_insertion_param);
	uniform_int_distribution<long> random_node(0, num_nodes);
	int num_inserts;
	for (long i = 0; i < num_general_inserts; i++)
	{
		long rand_u = random_node(generator);
		long rand_v = random_node(generator);
		
		num_inserts = min(general_insertion_cap, general_insertions(generator));
		// These edges do not persist to the end of the stream
		num_updates = 2 * num_inserts;

		for (int j = 0; j < num_updates; j++)
		{	
			for (int k = 0; k < rand_prefix_len; k++)
				random_prefix[k] = (char) rand_char(generator);
		
			stream_file_out << random_prefix << "\t" << rand_u << "\t" << rand_v << "\n";
		}

		total_num_updates += num_updates;
	}

	delete [] random_prefix;
	stream_file_out.flush();

	// External memory sort by the random prefix to produce 
	// a random permutation

	system((string("sort -k 1,1 ") + stream_file_name + string(" -o ") + stream_file_name).c_str());

	// Remove random prefixes used for permuting
	// Insert prefix denoting insertion or deletion

	vector<bool> edge_present(num_nodes * num_nodes, false);

	stream_file_out.clear();
	stream_file_out.seekp(0);
	// Mantain separate streams to avoid unnecessary flushes
	ifstream stream_file_in{stream_file_name};

	string prefix;
	long index;
	char upd_type;

	long temp_u, temp_v;
	stream_file_in >> prefix >> temp_u >> temp_v;
	stream_file_out << num_nodes << '\t' << total_num_updates << '\n';

	while (!stream_file_in.eof())
	{
		stream_file_in >> prefix >> u >> v;
		
		index = u * num_nodes + v;
		upd_type = edge_present[index] ? '1' : '0';
		edge_present[index] = !edge_present[index];

		stream_file_out << upd_type << '\t' << u << '\t' << v << '\n';
	}

	// Reinsert femporarily removed first update
	index = temp_u * num_nodes + temp_v;
	upd_type = edge_present[index] ? '1' : '0';
	edge_present[index] = !edge_present[index];

	stream_file_out << upd_type << '\t' << temp_u << '\t' << temp_v << '\n';

	// Truncate remainder of file
	
	if (-1 == truncate(stream_file_name.c_str(), stream_file_out.tellp()))
	{
		cout << "Truncation error" << endl;
		return 1;
	}

	return 0;
}
