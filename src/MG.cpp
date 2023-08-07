#include "MG.h"
#include "P1Fem.h"
#include <thread>

void MG::gen_cond() {
	// 生成m_A_v, m_b_v, m_num_nodes_each_edge, m_nodes
	for (unsigned int i = 0; i <= m_lc; ++i) {
		m_num_nodes_each_edge[i] = pow(2, i+1) + 1; 
		P1Fem tmp(m_num_nodes_each_edge[i]);
		tmp.assemble();
		m_A_v[i] = tmp.copy_A();
		m_b_v[i] = tmp.copy_b(); 
		if (i == m_lc) {
			m_nodes = tmp.copy_nodes();	
		}
	}
	m_numerical_x.resize(m_num_nodes_each_edge[m_lc]*m_num_nodes_each_edge[m_lc]);
	// 生成各层的延拓矩阵	
	for (unsigned int i = 0; i < m_lc; ++i) {
		m_P_v[i] = get_prolongation_matrix(i);
	}
}


VectorXd MG::mg(unsigned int l, VectorXd z0, const VectorXd& g) {
	if (l == 0) {
  		Eigen::SparseLU<SpMat> solver;
		solver.compute(m_A_v[l]);
		return solver.solve(g);
	}
	else {

		VectorXd pre_result = gauss_seidel(l, z0, g);
		
		// 校正
		  // 将残差限制到粗网格上
		VectorXd residual = g - m_A_v[l] * pre_result;
		VectorXd cu_residual(m_b_v[l-1].size());
		SpMat& prolongation_mat = m_P_v[l-1];
		cu_residual = prolongation_mat.transpose() * residual;

		VectorXd qp = mg(l-1, VectorXd::Zero(m_b_v[l-1].size()), cu_residual);

		// 对粗一级的网格上得到的解进行插值
		z0 = pre_result + prolongation_mat * qp;

		return gauss_seidel(l, z0, g);
		
	}
};

unsigned int MG::solve() {
	gen_cond();
	VectorXd temp_solu(m_numerical_x.size());
	// temp_solu = mg(m_lc, VectorXd::Zero(m_b_v[m_lc].size()), m_b_v[m_lc]);
	double b_norm = m_b_v[m_lc].lpNorm<Infinity>();
	double el; 
	unsigned int iteration = 0;
	do {
		temp_solu = mg(m_lc, temp_solu, m_b_v[m_lc]);
		el = (m_b_v[m_lc] - m_A_v[m_lc]*temp_solu).lpNorm<Infinity>() / b_norm;
		iteration++;
	} while (el > m_tol);

	m_numerical_x = temp_solu;
	
	return iteration;
}

double MG::err() {
	VectorXd exact_x = (2*PI*m_nodes.array().col(0)).sin() * (2*PI*m_nodes.array().col(1)).cos();
	return (m_numerical_x - exact_x).lpNorm<Infinity>();
};

void MG::write_to_file(string file_path) {
  std::ofstream file(file_path);
  if (file.is_open()) {
    double last_x = m_nodes(0, 0);
    int k = 0;
    file << "x,y,z\n";
    for (int i = 0; i < m_nodes.rows(); ++i) {
      // 若两个点的x坐标不同，则先插入一个空行，然后再输出坐标
      if (abs(m_nodes(i, 0) - last_x) > 1e-8) {
        file << "\n";
        last_x = m_nodes(i, 0);
      }
      file << m_nodes(i, 0) << ", " << m_nodes(i, 1) << ", " << m_numerical_x(k++) << "\n";
    }
    file.close();
    std::cout << "The numerical results has been written to " << file_path << "\n";
  }
  else {
    std::cout << "Unable to open file!\n";
  }
}


VectorXd MG::jacobi(unsigned int l, VectorXd& z0, const VectorXd& b) {
	VectorXd x_next(m_b_v[l].size());
	unsigned int iteration = 0;
	VectorXd* curr_ptr = &z0;
	VectorXd* next_ptr = &x_next;
	while (iteration < m_max_iteration) {
		(*next_ptr) = (b - m_A_v[l]*(*curr_ptr)).cwiseQuotient(m_A_v[l].diagonal()) + (*curr_ptr);
		VectorXd* tmp = curr_ptr;
		curr_ptr = next_ptr;
		next_ptr = tmp;
		iteration++;
	}

	return std::move(*curr_ptr);
};

VectorXd MG::gauss_seidel(unsigned int l, VectorXd& z0, const VectorXd& b) {
	VectorXd x_next(m_b_v[l].size());
	unsigned int iteration = 0;
	VectorXd* curr_ptr = &z0;
	VectorXd* next_ptr = &x_next;
	Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
	while (iteration < m_max_iteration) {
		solver.compute(m_A_v[l].triangularView<Lower>());
		*next_ptr = solver.solve(b - m_A_v[l].triangularView<StrictlyUpper>() * (*curr_ptr));
		VectorXd* tmp = curr_ptr;
		curr_ptr = next_ptr;
		next_ptr = tmp;
		iteration++;
	}

	return std::move(*curr_ptr);
};

// 返回从低level层到level+1层的延拓矩阵
SpMat MG::get_prolongation_matrix(unsigned int level) {
	unsigned int num_nodes_xi = m_num_nodes_each_edge[level+1] * m_num_nodes_each_edge[level+1];
	unsigned int num_nodes_cu = m_num_nodes_each_edge[level] * m_num_nodes_each_edge[level];
	SpMat P(num_nodes_xi, num_nodes_cu);
	P.reserve(VectorXi::Constant(P.cols(), 4));
	unsigned int dnl = m_num_nodes_each_edge[level];
	unsigned int dnr = m_num_nodes_each_edge[level+1];
	for (unsigned int i = 0; i < dnl; ++i) {
		for (unsigned j = 0; j < dnl; ++j) {
			P.insert(2*i*dnr+2*j, i*dnl+j) = 1;	
		}
	}
	for (unsigned int i = 0; i < dnl; ++i) {
		for (unsigned int j = 0; j < dnl-1; ++j) {
			P.insert(2*i*dnr+2*j+1, i*dnl+j) = 0.5;
			P.insert(2*i*dnr+2*j+1, i*dnl+j+1) = 0.5;
		}
	}
	for (unsigned int i = 0; i < dnl-1; ++i) {
		for (unsigned int j = 0; j < dnl; ++j) {
			P.insert((2*i+1)*dnr+2*j, i*dnl+j) = 0.5;
			P.insert((2*i+1)*dnr+2*j, (i+1)*dnl+j) = 0.5;
		}
	}

	for (unsigned int i = 0; i < dnl-1; ++i) {
		for (unsigned int j = 1; j < dnl; ++j) {
			P.insert((2*i+1)*dnr+2*j-1, i*dnl+j) = 0.5;
			P.insert((2*i+1)*dnr+2*j-1, (i+1)*dnl+j-1) = 0.5;
		}
	}
	P.makeCompressed();
	return P;
}