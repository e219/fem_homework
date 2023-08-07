#pragma once

#include <Eigen/Eigen>
#include <fstream>
#include <iostream>

using Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXd, Eigen::Ref, Eigen::RowVectorXd, std::string, Eigen::VectorXi, Eigen::ArrayXd, Eigen::ArrayXXd, Eigen::all, Eigen::Infinity;

using SpMat = Eigen::SparseMatrix<double>;

class P1Fem {
private:
  int m_num_nodes_each_edge;
  int m_num_each_edge;
  int m_num_nodes;
  int m_num_elems;
  MatrixXi m_elem;
  MatrixXi m_rb;
  MatrixXd m_d_lambda;
  VectorXd m_area;
  MatrixXd m_nodes;

  void grad_basis();
  
  void f(const Ref<const VectorXd>& x, const Ref<const VectorXd>& y, Ref<VectorXd> fv);
  
  void calc_fun_g(Ref<ArrayXXd> g_val);
  
  
protected:
  const double PI = 3.14159265358979323846;
  SpMat m_A;
  VectorXd m_numerical_x;
  VectorXd m_b;
  void gen_cond();

public:
  P1Fem(unsigned int _num_nodes_each_edge): 
    m_num_nodes_each_edge(_num_nodes_each_edge),\
    m_num_each_edge(_num_nodes_each_edge - 1),\
	  m_num_nodes(_num_nodes_each_edge * _num_nodes_each_edge),\
	  m_num_elems(m_num_each_edge * m_num_each_edge * 2) {

      m_nodes = MatrixXd::Zero(m_num_nodes, 2);
      m_elem = MatrixXi::Zero(m_num_elems, 3);
      m_rb = MatrixXi::Zero(m_num_each_edge*4, 2);
      m_d_lambda = MatrixXd::Zero(m_num_elems, 6);
      m_area = VectorXd::Zero(m_num_elems);
      m_A.resize(m_num_nodes, m_num_nodes); m_A.setZero();
      m_numerical_x = VectorXd::Zero(m_num_nodes);
      m_b = VectorXd::Zero(m_num_nodes);
    }

  P1Fem& assemble();

  void direct_solve();

  const VectorXd &get_x() { return m_numerical_x; };

  const MatrixXd &get_nodes() { return m_nodes; };

  double err();

  void write_to_file(string file_path);

  SpMat copy_A() { return m_A; };

  SpMat move_A() { return std::move(m_A); };

  VectorXd copy_b() { return m_b; };

  VectorXd move_b() { return std::move(m_b); };

  MatrixXd copy_nodes() { return m_nodes; }  

};