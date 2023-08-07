#include <iostream>
#include <memory>
#include <string>
#include <thread>
#include <vector>
#include <list>

using std::cout, std::endl, std::string;

std::vector<int> g_v {1, 2, 3};

void myprint(int inum) {
    cout << "thread id = " << std::this_thread::get_id() << ", gv = " << g_v[0] << g_v[1] << g_v[2] << endl;
}

class A {
private:
    std::list<int> msg_recv_queue;
public:
    void in_msg_recv_queue() {
        int i = 0;
        for (unsigned i = 0; i < 10000; ++i) {
            cout <<  "in_msg_recv_queue is executing, insert a element " << i << endl;
            msg_recv_queue.push_back(i);
        }
    }
    
    void out_msg_recv_queue() {
        for (unsigned i = 0; i < 10000; ++i) {
            if (!msg_recv_queue.empty()) {
                msg_recv_queue.pop_front();
            }
            else {
                cout << "out_msg_recv_queue is executing, but queue is empty now.";
            }
        }
    }        
};

int main() {
    A myobj;
    std::thread t1(&A::out_msg_recv_queue, &myobj);
    std::thread t2(&A::in_msg_recv_queue, &myobj);
    t1.join();
    t2.join();

    return 0;
}