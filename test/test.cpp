#include <iostream>
#include <memory>
#include <string>
#include <thread>
#include <vector>
#include <list>
#include <Eigen/Eigen>

using std::cin, std::cout, std::endl, std::string;

// std::vector<int> g_v {1, 2, 3};

// void myprint(int inum) {
//     cout << "thread id = " << std::this_thread::get_id() << ", gv = " << g_v[0] << g_v[1] << g_v[2] << endl;
// }

// class A {
// private:
//     std::list<int> msg_recv_queue;
// public:
//     void in_msg_recv_queue() {
//         int i = 0;
//         for (unsigned i = 0; i < 10000; ++i) {
//             cout <<  "in_msg_recv_queue is executing, insert a element " << i << endl;
//             msg_recv_queue.push_back(i);
//         }
//     }
    
//     void out_msg_recv_queue() {
//         for (unsigned i = 0; i < 10000; ++i) {
//             if (!msg_recv_queue.empty()) {
//                 msg_recv_queue.pop_front();
//             }
//             else {
//                 cout << "out_msg_recv_queue is executing, but queue is empty now.";
//             }
//         }
//     }        
// };

// void replace_space(string& s) {
//     int count = 0;
//     for (char ch : s) {
//         if (ch == ' ') count++;
//     }
//     int len = s.size();
//     s.resize(len + 2*count);
//     for (int i = len-1, j = s.size()-1; i < j; i--, j--) {
//         if (s[i] != ' ') s[j] = s[i];
//         else {
//            s[j--] = '0' ;
//            s[j--] = '2' ;
//            s[j] = '%' ;
//         }
//     }
// }

int main() {
    Eigen::MatrixXi mat(3, 3);
    mat << 0, 2, 1,
           1, 2, 0,
           0, 1, 2;
    // cout << "initial Matrix is:\n" << mat << endl;
    Eigen::MatrixXi A = Eigen::MatrixXi::Random(3, 3);
    cout << "initial Matrix A is:\n" << A << endl;
    // cout << "A(mat(Eigen::all, {1, 2, 0}), 0)" << A(mat(Eigen::all, {1, 2, 0}), 0) << endl;
    return 0;
}