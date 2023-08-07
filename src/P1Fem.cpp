#include "P1Fem.h"


void P1Fem::gen_cond() {
  // 下面生成节点的坐标
  // 求解区域边界为[0, 1]x[0, 1], 三角形单元，单元全等
  double dx = 1.0 / (m_num_nodes_each_edge - 1);  
  double dy = dx;
  for (int i = 0; i < m_num_nodes; ++i) {
    int coord_x = i / m_num_nodes_each_edge;
    int coord_y = i % m_num_nodes_each_edge;
    m_nodes(i, 0) = coord_x * dx;
    m_nodes(i, 1) = coord_y * dy;
  }
  // 下面生成边界
  // 左边界
  int kth = 0; 
  for (int i = 0; i < m_num_each_edge; ++i) {
    m_rb(kth, 0) = i;
    m_rb(kth++, 1) = i + 1;
  }
  // 下边界
  for (int i = 0; i < m_num_each_edge; ++i) {
    m_rb(kth, 0) = i * m_num_nodes_each_edge;
    m_rb(kth++, 1) = (i+1) * m_num_nodes_each_edge;
  }
  // 右边界
  for (int i = 0; i < m_num_each_edge; ++i) {
    m_rb(kth, 0) = i + m_num_nodes_each_edge * (m_num_nodes_each_edge - 1);
    m_rb(kth++, 1) = i + 1 + m_num_nodes_each_edge * (m_num_nodes_each_edge - 1);
  }
  // 上边界
  for (int i = 0; i < m_num_each_edge; ++i) {
    m_rb(kth, 0) = i * m_num_nodes_each_edge + m_num_each_edge;
    m_rb(kth++, 1) = (i+1) * m_num_nodes_each_edge + m_num_each_edge;
  }
  // 最后生成单元
  // 三角形单元以右下角的定点为起点，按逆时针方向依次放入elem中
  kth = 0;
  for (int j = 1; j < m_num_nodes_each_edge; ++j) {
    for (int i = 0; i < m_num_each_edge; ++i) {
      // 第一个三角形单元 
      m_elem(kth, 0) = j * m_num_nodes_each_edge + i;
      m_elem(kth, 1) = (j - 1) * m_num_nodes_each_edge + i + 1;
      m_elem(kth++, 2) = (j - 1) * m_num_nodes_each_edge + i;
      // 第二个三角形单元 
      m_elem(kth, 0) = j * m_num_nodes_each_edge + i;
      m_elem(kth, 1) = j * m_num_nodes_each_edge + i + 1;
      m_elem(kth++, 2) = (j - 1) * m_num_nodes_each_edge + i + 1;
    }
  }
}

// 计算基函数的梯度
void P1Fem::grad_basis() {
  MatrixXd ve1 = m_nodes(m_elem.col(1), all) - m_nodes(m_elem.col(2), all);
  MatrixXd ve2 = m_nodes(m_elem.col(2), all) - m_nodes(m_elem.col(0), all);
  MatrixXd ve3 = m_nodes(m_elem.col(0), all) - m_nodes(m_elem.col(1), all);
  m_area = 0.5 * (m_nodes(m_elem.col(0), 0).cwiseProduct(ve1.col(1)) + \
                m_nodes(m_elem.col(1), 0).cwiseProduct(ve2.col(1)) + \
                m_nodes(m_elem.col(2), 0).cwiseProduct(ve3.col(1)) );
  // 三角形单元第一个定点对应的基函数的梯度
  m_d_lambda(all, 0) = ve1.col(1).cwiseQuotient(2 * m_area);
  m_d_lambda(all, 1) = -ve1.col(0).cwiseQuotient(2 * m_area);
  // 三角形单元第二个定点对应的基函数的梯度
  m_d_lambda(all, 2) = ve2.col(1).cwiseQuotient(2 * m_area);
  m_d_lambda(all, 3) = -ve2.col(0).cwiseQuotient(2 * m_area);
  // 三角形单元第三个定点对应的基函数的梯度
  m_d_lambda(all, 4) = ve3.col(1).cwiseQuotient(2 * m_area);
  m_d_lambda(all, 5) = -ve3.col(0).cwiseQuotient(2 * m_area);
}

// 方程右端项函数f = 8 * PI^2 * sin(2*PI*x) * cos(2*PI*y)
inline void P1Fem::f(const Ref<const VectorXd>& x, const Ref<const VectorXd>& y, Ref<VectorXd> fv) {
  fv = 8 * PI * PI * (2*PI*x).array().sin() * (2*PI*y).array().cos();
}


P1Fem& P1Fem::assemble() {
  // 先计算好必要的数据
  gen_cond();
  grad_basis();
  // 组装刚度矩阵
  // 由于刚度矩阵A稀疏，出于效率考虑，给一个每列非零元素数量的估计值
  int num_estimate_nz = std::min(7, int(m_A.rows()));
  m_A.reserve(VectorXi::Constant(m_A.cols(), num_estimate_nz));
  m_A.setZero();
  for (int i = 0; i < m_elem.cols(); ++i) {
    for (int j = i; j < m_elem.cols(); ++j) {
      VectorXd b_ij = m_d_lambda.col(2 * i).array() * m_d_lambda.col(2 * j).array() +\
                      m_d_lambda.col(2 * i + 1).array() * m_d_lambda.col(2 * j + 1).array();
      b_ij = b_ij.array() * m_area.array();
      // 填充稀疏系数矩阵A
      for (int rowth = 0; rowth < m_elem.rows(); ++rowth) {
          m_A.coeffRef(m_elem(rowth, i), m_elem(rowth, j)) += b_ij[rowth];
          if (i != j) {
            m_A.coeffRef(m_elem(rowth, j), m_elem(rowth, i)) += b_ij[rowth];

          }
        }
      }
   }
  for (int i = 0; i < m_rb.rows(); ++i) {
    m_A.coeffRef(m_rb(i, 0), m_rb(i, 1)) += 1.0 / (6 * m_num_each_edge);
    m_A.coeffRef(m_rb(i, 1), m_rb(i, 0)) += 1.0 / (6 * m_num_each_edge);
    m_A.coeffRef(m_rb(i, 0), m_rb(i, 0)) += 1.0 / (3 * m_num_each_edge);
    m_A.coeffRef(m_rb(i, 1), m_rb(i, 1)) += 1.0 / (3 * m_num_each_edge);
  }
  m_A.makeCompressed();
  // 组装右端项b
  // coord_quad是四个求积节点的面积坐标，weight是这四个求积节点的权重
  MatrixXd coord_quad{{1.0 / 3, 1.0 / 3, 1.0 / 3},
                      {0.6, 0.2, 0.2},
                      {0.2, 0.6, 0.2},
                      {0.2, 0.2, 0.6}};
  VectorXd weight{{-27.0 / 48, 25.0 / 48, 25.0 / 48, 25.0 / 48}};
  MatrixXd fv_mat(m_elem.rows(), coord_quad.rows());
  // 将各求积节点的面积坐标转化为笛卡尔坐标
  for (int i = 0; i < coord_quad.rows(); ++i) {
    MatrixXd pxy = coord_quad(i, 0) * m_nodes(m_elem.col(0), all) + coord_quad(i, 1) * m_nodes(m_elem.col(1), all) + coord_quad(i, 2) * m_nodes(m_elem.col(2), all);
    f(pxy.col(0), pxy.col(1), fv_mat.col(i));
    fv_mat.col(i) *= weight(i);
  }
  fv_mat.array().colwise() *= m_area.array();
  m_b(m_elem.reshaped()) += (fv_mat * coord_quad).reshaped();
  // 处理边界条件
  ArrayXXd g_val(m_rb.rows(), 2);
  calc_fun_g(g_val);
  double length_each_edge = 1.0 / m_num_each_edge;
  m_b(m_rb.col(0)) += (g_val.col(0) * (3 - sqrt(3)) + g_val.col(1) * (3 + sqrt(3))).matrix() * length_each_edge / 12;
  m_b(m_rb.col(1)) += (g_val.col(0) * (3 + sqrt(3)) + g_val.col(1) * (3 - sqrt(3))).matrix() * length_each_edge / 12;

  return *this;
};



void P1Fem::calc_fun_g(Ref<ArrayXXd> g_val) {
    ArrayXXd boundary_quad_1 = (m_nodes(m_rb.col(0), all) + m_nodes(m_rb.col(1), all)) / 2.0 + (m_nodes(m_rb.col(1), all) - m_nodes(m_rb.col(0), all)) / (2*sqrt(3));
    ArrayXXd boundary_quad_2 = (m_nodes(m_rb.col(0), all) + m_nodes(m_rb.col(1), all)) / 2.0 - (m_nodes(m_rb.col(1), all) - m_nodes(m_rb.col(0), all)) / (2*sqrt(3));
    // 左边
  // VectorXd exact_x = (2*PI*m_no界
    ArrayXd x1 = boundary_quad_1.topLeftCorner(m_num_each_edge, 1);
    ArrayXd y1 = boundary_quad_1.topRightCorner(m_num_each_edge, 1);
    ArrayXd x2 = boundary_quad_2.topLeftCorner(m_num_each_edge, 1);
    ArrayXd y2 = boundary_quad_2.topRightCorner(m_num_each_edge, 1);
    
    g_val.topLeftCorner(m_num_each_edge, 1) = (2*PI*y1).cos() * ((2*PI*x1).sin() - 2*PI*(2*PI*x1).cos());
    g_val.topRightCorner(m_num_each_edge, 1) = (2*PI*y2).cos() * ((2*PI*x2).sin() - 2*PI*(2*PI*x2).cos());

    // 下边界
    x1 = boundary_quad_1.block(m_num_each_edge, 0, m_num_each_edge, 1);
    y1 = boundary_quad_1.block(m_num_each_edge, 1, m_num_each_edge, 1);
    x2 = boundary_quad_2.block(m_num_each_edge, 0, m_num_each_edge, 1);
    y2 = boundary_quad_2.block(m_num_each_edge, 1, m_num_each_edge, 1);

    g_val.block(m_num_each_edge, 0, m_num_each_edge, 1) = (2*PI*x1).sin() * (2*PI*(2*PI*y1).sin() + (2*PI*y1).cos());
    g_val.block(m_num_each_edge, 1, m_num_each_edge, 1) = (2*PI*x2).sin() * (2*PI*(2*PI*y2).sin() + (2*PI*y2).cos());

    // 右边界
    x1 = boundary_quad_1.block(2*m_num_each_edge, 0, m_num_each_edge, 1);
    y1 = boundary_quad_1.block(2*m_num_each_edge, 1, m_num_each_edge, 1);
    x2 = boundary_quad_2.block(2*m_num_each_edge, 0, m_num_each_edge, 1);
    y2 = boundary_quad_2.block(2*m_num_each_edge, 1, m_num_each_edge, 1);

    g_val.block(2*m_num_each_edge, 0, m_num_each_edge, 1) = (2*PI*y1).cos() * ((2*PI*x1).sin() + 2*PI*(2*PI*x1).cos());
    g_val.block(2*m_num_each_edge, 1, m_num_each_edge, 1) = (2*PI*y2).cos() * ((2*PI*x2).sin() + 2*PI*(2*PI*x2).cos());

    // 上边界
    x1 = boundary_quad_1.bottomLeftCorner(m_num_each_edge, 1);
    y1 = boundary_quad_1.bottomRightCorner(m_num_each_edge, 1);
    x2 = boundary_quad_2.bottomLeftCorner(m_num_each_edge, 1);
    y2 = boundary_quad_2.bottomRightCorner(m_num_each_edge, 1);
    g_val.bottomLeftCorner(m_num_each_edge, 1) = (2*PI*x1).sin() * (-2*PI*(2*PI*y1).sin() + (2*PI*y1).cos());
    g_val.bottomRightCorner(m_num_each_edge, 1) = (2*PI*x2).sin() * (-2*PI*(2*PI*y2).sin() + (2*PI*y2).cos());
};

void P1Fem::write_to_file(string file_path) {
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

void P1Fem::direct_solve() {
  Eigen::SparseLU<SpMat> solver;
  solver.compute(m_A);
  m_numerical_x = solver.solve(m_b);
}

double P1Fem::err() {
  VectorXd exact_x = (2*PI*m_nodes.array().col(0)).sin() * (2*PI*m_nodes.array().col(1)).cos();
  // MatrixXd err_mat(m_elem.rows(), m_elem.cols());
  // err_mat.reshaped() = exact_x(m_elem.reshaped()) - m_numerical_x(m_elem.reshaped());
  // err_mat = err_mat.array().square();
  // VectorXd vl = err_mat.rowwise().sum();
  // return (vl.array() * m_area.array()/3).sum();
  return (m_numerical_x - exact_x).lpNorm<Infinity>();
}