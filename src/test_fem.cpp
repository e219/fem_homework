#include <iostream>
#include "P1Fem.h"
#include "MG.h"
#include <vector>
#include <fstream>

#define OUTPUT(x) std::cout << x << std::endl;

void write_to_file(const std::string& file_path, const std::vector<double>& h_array, const std::vector<double>& err_array);

void test_p1fem();
void test_mg();

int main() {
	OUTPUT("test P1Fem or MG(0/1): ");
	unsigned int choice;
	while (true) {
		std::cin >> choice;
		if(choice == 0) {
			test_p1fem();
			break;
		}
		else if (choice == 1) {
			test_mg();
			break;
		}
		else {
			OUTPUT("wrong choice!! Choose again!");
		}
	}
}

void write_to_file(const std::string& file_path, const std::vector<double>& h_array, const std::vector<double>& err_array) {
	std::ofstream file(file_path);
	if(file.is_open()) {
		file << "h, err(0)\n";
		for (unsigned int i = 0; i < h_array.size(); ++i) {
			file << h_array[i] << ", " << err_array[i] << "\n";
		}
		std::cout << "success to write to " << file_path << "\n";
		file.close();
	}
	else {
		std::cout << "unable to open file!\n\n";
	}
}

void test_p1fem() {
	std::cout << "input until_level: ";
	unsigned int until_level;
	std::cin >> until_level;

	std::vector<double> h_v(until_level+1, 0);
	std::vector<double> err_v(until_level+1, 0);
	for (unsigned int i = 0; i <= until_level; ++i) {
		unsigned int num_nodes_each_edge = pow(2, i+1) + 1;
		P1Fem temp(num_nodes_each_edge);
		temp.assemble().direct_solve();
		err_v[i] = temp.err();
		h_v[i] = 1.0 / (num_nodes_each_edge-1);
		temp.write_to_file("../data/numerical_result_" + std::to_string(i) + ".csv");
	}

	write_to_file("../data/multi_h_err.csv", h_v, err_v);
}


void test_mg() {
	std::cout << "input the max_level: ";
	unsigned int max_level;
	std::cin >> max_level;

	std::cout << "input the max_iteration: ";
	unsigned int max_iteration;
	std::cin >> max_iteration;

	std::cout << "input the tol: ";
	double tol;
	std::cin >> tol;

	MG my_mg(max_level, max_iteration, tol);

	unsigned int iteration = my_mg.solve();
	std::cout << "l-infinity norm: " << my_mg.err() << "\n";
	std::cout << "iteration times: " << iteration << "\n";

	my_mg.write_to_file("../data/mg_numerical_result.csv");
}