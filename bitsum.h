#include <emp-ot/emp-ot.h>
#include <emp-tool/emp-tool.h>
#include <iostream>
using namespace emp;

template <typename IO, template <typename> class T>
uint64_t bitsum_ot_sender(IO *io,bool *shares, bool *valid, int n){
    PRG prg(fix_key);

    uint64_t b0[n];
    uint64_t b1[n];
    uint64_t r[n];
    uint64_t sum = 0;
    prg.random_data(b0,n*sizeof(uint64_t));

    for(int i = 0; i < n; i++){
        if(shares[i] == true){
            b1[i] = b0[i] - 1;
            r[i] = (UINT64_MAX - b1[i]) + 1;
        }
        else {
            b1[i] = b0[i] + 1;
            r[i] = (UINT64_MAX - b0[i]) + 1;
        }
    }

    block *b0_ot = new block[n];
    block *b1_ot = new block[n];

    for(int i = 0; i < n; i++){
        uint64_t *p = (uint64_t*)&b0_ot[i];
        p[0] = valid[i] ? 0:1;
        p[1] = b0[i];
        uint64_t *q = (uint64_t*)&b1_ot[i];
        q[0] = p[0];
        q[1] = b1[i];
    }

    io->sync();
    auto start = clock_start();
    T<IO> *ot = new T<IO>(io);
    ot->send(b0_ot, b1_ot, n);
    io->flush();
    long long t = time_from(start);

    delete[] b1_ot;
    delete[] b0_ot;

    for(int i = 0; i < n; i++){
        if(valid[i])
            sum += r[i];
    }

    return sum;
}

template <typename IO, template <typename> class T>
uint64_t intsum_ot_sender(IO *io,uint32_t *shares, bool *valid, int n, int num_bits){
    PRG prg(fix_key);

    uint64_t bool_shares[n*num_bits];
    bool bool_valid[n*num_bits];
    uint64_t b0[n*num_bits];
    uint64_t b1[n*num_bits];
    uint32_t r[n*num_bits];
    uint64_t sum = 0;
    prg.random_data(r,n*num_bits*sizeof(uint32_t));
    uint64_t mod = 1 << num_bits;

    for(int i = 0; i < n; i++){
        uint32_t num = shares[i];
        std::cout << "Share : " << num << " num bits " << num_bits << std::endl;
        for(int j = 0; j < num_bits; j++){
            std::cout << num%2 << std::endl;
            bool_shares[i*num_bits + j] = num%2;
            bool_valid[i*num_bits + j] = valid[i];
            r[i*num_bits+j] = r[i*num_bits+j];
            std::cout << "r[" << j << "] : " << r[j] << std::endl;
            num = num >> 1;
        }
    }

    for(int i = 0; i < n; i++){
        for(int j = 0; j < num_bits; j++){
            b0[i*num_bits+j] = ((bool_shares[i*num_bits+j])*(1<<(j)) - r[i*num_bits+j]);
            std::cout << "b0[" << j << "] = " << b0[i*num_bits+j] << std::endl;
            b1[i*num_bits+j] = ((1 - bool_shares[i*num_bits+j])*(1<<(j)) - r[i*num_bits+j]);
            std::cout << "b1[" << j << "] = " << b1[i*num_bits+j] << std::endl;
        }
    }

    block *b0_ot = new block[n*num_bits];
    block *b1_ot = new block[n*num_bits];

    for(int i = 0; i < n*num_bits; i++){
        uint64_t *p = (uint64_t*)&b0_ot[i];
        p[0] = bool_valid[i] ? 0:1;
        p[1] = b0[i];
        uint64_t *q = (uint64_t*)&b1_ot[i];
        q[0] = p[0];
        q[1] = b1[i];
    }

    io->sync();
    auto start = clock_start();
    T<IO> *ot = new T<IO>(io);
    ot->send(b0_ot, b1_ot, n*num_bits);
    io->flush();
    long long t = time_from(start);

    delete[] b1_ot;
    delete[] b0_ot;

    for(int i = 0; i < n; i++){
        if(valid[i]){
            uint32_t local_share = 0;
            for(int j = i*num_bits; j < (i+1)*num_bits; j++)
                local_share += r[j];

            sum += local_share;
        }
            
    }

    return sum;
}

template <typename IO, template <typename> class T>
uint64_t bitsum_ot_receiver(IO *io,bool *shares, int n){

    block *r = new block[n];
    uint64_t ans[n];
    uint64_t sum = 0;
    io->sync();
    auto start = clock_start();
    T<IO> *ot = new T<IO>(io);
    ot->recv(r, shares, n);
    io->flush();
    long long t = time_from(start);

    for(int i = 0; i < n; i++){

        uint64_t *p = (uint64_t*)&r[i];
        ans[i] =(p[0] == 0) ? p[1] : 0;
        if(ans[i])
            sum += ans[i];
    }

    return sum;
}



template <typename IO, template <typename> class T>
uint64_t intsum_ot_receiver(IO *io, uint32_t *shares, int n, int num_bits){

    block *r = new block[n*num_bits];
    bool bool_shares[n*num_bits];

    for(int i = 0; i < n; i++){
        uint32_t num = shares[i];
        std::cout << "Share : " << num << " num bits " << num_bits << std::endl;
        for(int j = 0; j < num_bits; j++){
            bool_shares[i*num_bits + j] = num%2;
            num = num >> 1;
        }
    }

    uint64_t ans[n];
    uint64_t sum = 0;
    io->sync();
    auto start = clock_start();
    T<IO> *ot = new T<IO>(io);
    ot->recv(r, bool_shares, n*num_bits);
    io->flush();
    long long t = time_from(start);

    for(int i = 0; i < n; i++){
        uint64_t valid = ((uint64_t*) &r[i*num_bits])[0];
        if(valid == 0){
            for(int j = 0; j < num_bits; j++){
                uint64_t *p = (uint64_t*)&r[i*num_bits+j];
                sum += p[1];
                std::cout << p[1] << endl;
            }
        }  
    }

    return sum;
}