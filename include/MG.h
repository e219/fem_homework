#pragma once

#include <Eigen/Eigen>
#include <vector>
#include "P1Fem.h"
#include <iostream>

using Eigen::seqN, Eigen::Infinity, Eigen::Lower, Eigen::StrictlyUpper;


class MG {
private:
  unsigned int m_lc;
  std::vector<SpMat> m_A_v;
  std::vector<VectorXd> m_b_v;
  std::vector<unsigned int> m_num_nodes_each_edge;
  std::vector<SpMat> m_P_v;
  VectorXd m_numerical_x;
  MatrixXd m_nodes;
  const double PI = 3.14159265358979323846;
  unsigned int m_max_iteration;
  double m_tol;

  void gen_cond();

  VectorXd jacobi(unsigned int l, VectorXd& z0, const VectorXd& b);

  VectorXd gauss_seidel(unsigned int l, VectorXd& z0, const VectorXd& b);

  VectorXd mg(unsigned int l, VectorXd z0, const VectorXd& g);

  SpMat get_prolongation_matrix(unsigned int level);

public:
  MG(unsigned int _lc, unsigned int _max_iteration = 3, double _tol = pow(10, -3)): m_lc(_lc), m_max_iteration(_max_iteration), m_tol(_tol) {
	m_A_v.reserve(m_lc + 1);
	m_b_v.reserve(m_lc + 1);
	m_P_v.reserve(m_lc);
	m_num_nodes_each_edge.reserve(m_lc + 1);
  };

  
  unsigned int solve();

  void printResults() {
	std::cout << "the numerical solution is(mg):\n" << m_numerical_x << '\n';
	std::cout << "the residual is:\n" << (m_b_v[m_lc] - m_A_v[m_lc] * m_numerical_x).array().abs() << '\n';
  }
  
  double err();

  void write_to_file(string file_path);

};