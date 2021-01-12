#include "../proto.h"

#define SERVER0_IP "127.0.0.1"
#define SERVER1_IP "127.0.0.1"

int main(int argc, char** argv){
    if(argc != 3){
        return -1;
    }

    const int server_num = atoi(argv[1]);
    const int m = atoi(argv[2]);

    NetIO* io = new NetIO(server_num ? SERVER0_IP : nullptr, 60051);

    gen_boolean_beaver_triples(io,server_num,m);

    return 0;
}