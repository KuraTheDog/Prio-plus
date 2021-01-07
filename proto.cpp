#include "proto.h"

uint64_t bitsum_ot_sender(NetIO *io,bool *shares, bool *valid, int n){
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
    IKNP<NetIO> *ot = new IKNP<NetIO>(io);
    ot->send(b0_ot, b1_ot, n);
    io->flush();

    delete[] b1_ot;
    delete[] b0_ot;

    for(int i = 0; i < n; i++){
        if(valid[i])
            sum += r[i];
    }

    return sum;
}

uint64_t intsum_ot_sender(NetIO *io,uint32_t *shares, bool *valid, int n, int num_bits){
    PRG prg(fix_key);

    uint64_t bool_shares[n*num_bits];
    bool bool_valid[n*num_bits];
    uint64_t b0[n*num_bits];
    uint64_t b1[n*num_bits];
    uint64_t r[n*num_bits];
    uint64_t sum = 0;
    prg.random_data(r,n*num_bits*sizeof(uint64_t));

    for(int i = 0; i < n; i++){
        uint32_t num = shares[i];
        // std::cout << "Share : " << num << "  valid " << valid[i] << std::endl;
        // std::cout << "Valid : " << valid[i] << " num bits " << num_bits << std::endl;
        for(int j = 0; j < num_bits; j++){
            // std::cout << num%2 << std::endl;
            bool_shares[i*num_bits + j] = num%2;
            bool_valid[i*num_bits + j] = valid[i];
            r[i*num_bits+j] = r[i*num_bits+j];
            // std::cout << "r[" << j << "] : " << r[j] << std::endl;
            num = num >> 1;
        }
    }

    for(int i = 0; i < n; i++){
        for(int j = 0; j < num_bits; j++){
            b0[i*num_bits+j] = ((bool_shares[i*num_bits+j])*(1<<(j)) - r[i*num_bits+j]);
            // std::cout << "b0[" << j << "] = " << b0[i*num_bits+j] << std::endl;
            b1[i*num_bits+j] = ((1 - bool_shares[i*num_bits+j])*(1<<(j)) - r[i*num_bits+j]);
            // std::cout << "b1[" << j << "] = " << b1[i*num_bits+j] << std::endl;
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
    IKNP<NetIO> *ot = new IKNP<NetIO>(io);
    ot->send(b0_ot, b1_ot, n*num_bits);
    io->flush();

    delete[] b1_ot;
    delete[] b0_ot;

    for(int i = 0; i < n; i++){
        if(valid[i]){
            uint64_t local_share = 0;
            for(int j = i*num_bits; j < (i+1)*num_bits; j++)
                local_share += r[j];

            sum += local_share;
        }

    }

    return sum;
}

uint64_t bitsum_ot_receiver(NetIO *io,bool *shares, int n){

    block *r = new block[n];
    uint64_t ans[n];
    uint64_t sum = 0;
    io->sync();
    IKNP<NetIO> *ot = new IKNP<NetIO>(io);
    ot->recv(r, shares, n);
    io->flush();

    for(int i = 0; i < n; i++){

        uint64_t *p = (uint64_t*)&r[i];
        ans[i] =(p[0] == 0) ? p[1] : 0;
        if(ans[i])
            sum += ans[i];
    }

    delete[] r;

    return sum;
}

uint64_t intsum_ot_receiver(NetIO *io, uint32_t *shares, int n, int num_bits){

    block *r = new block[n*num_bits];
    bool bool_shares[n*num_bits];

    for(int i = 0; i < n; i++){
        uint32_t num = shares[i];
        // std::cout << "Share : " << num << " num bits " << num_bits << std::endl;
        for(int j = 0; j < num_bits; j++){
            bool_shares[i*num_bits + j] = num%2;
            num = num >> 1;
        }
    }

    uint64_t sum = 0;
    io->sync();
    IKNP<NetIO> *ot = new IKNP<NetIO>(io);
    ot->recv(r, bool_shares, n*num_bits);
    io->flush();

    for(int i = 0; i < n; i++){
        uint64_t valid = ((uint64_t*) &r[i*num_bits])[0];
        if(valid == 0){
            for(int j = 0; j < num_bits; j++){
                uint64_t *p = (uint64_t*)&r[i*num_bits+j];
                sum += p[1];
            }
        }  
    }

    delete[] r;

    return sum;
}

// Functions to convert XOR shared input values to shares that add up
// Problem : Do they add up with Int_Modulus?
uint64_t xor_to_sum_share_sender(NetIO *io, uint32_t share, int num_bits){
    PRG prg(fix_key);

    uint64_t bool_shares[num_bits];
    uint64_t b0[num_bits];
    uint64_t b1[num_bits];
    uint64_t r[num_bits];
    uint64_t sum = 0;
    prg.random_data(r,num_bits*sizeof(uint64_t));


    uint32_t num = share;
    // std::cout << "Share : " << num << "  valid " << valid[i] << std::endl;
    // std::cout << "Valid : " << valid[i] << " num bits " << num_bits << std::endl;
    for(int j = 0; j < num_bits; j++){
        // std::cout << num%2 << std::endl;
        bool_shares[j] = num%2;

        // std::cout << "r[" << j << "] : " << r[j] << std::endl;
        num = num >> 1;
    }

    for(int j = 0; j < num_bits; j++){
        b0[j] = ((bool_shares[j])*(1<<(j)) - r[j]);
        b1[j] = ((1 - bool_shares[j])*(1<<(j)) - r[j]);
    }


    block *b0_ot = new block[num_bits];
    block *b1_ot = new block[num_bits];

    for(int i = 0; i < num_bits; i++){
        uint64_t *p = (uint64_t*)&b0_ot[i];
        p[0] = 0;
        p[1] = b0[i];
        uint64_t *q = (uint64_t*)&b1_ot[i];
        q[0] = p[0];
        q[1] = b1[i];
    }

    io->sync();
    IKNP<NetIO> *ot = new IKNP<NetIO>(io);
    ot->send(b0_ot, b1_ot, num_bits);
    io->flush();

    delete[] b1_ot;
    delete[] b0_ot;


    uint64_t local_share = 0;
    for(int j = 0; j < num_bits; j++)
        local_share += r[j];

    sum += local_share;

    return sum;
}

uint64_t xor_to_sum_share_receiver(NetIO *io, uint32_t share, int num_bits){
    block *r = new block[num_bits];
    bool bool_shares[num_bits];

    uint32_t num = share;
    // std::cout << "Share : " << num << " num bits " << num_bits << std::endl;
    for(int j = 0; j < num_bits; j++){
        bool_shares[j] = num%2;
        num = num >> 1;
    }

    uint64_t sum = 0;
    io->sync();
    IKNP<NetIO> *ot = new IKNP<NetIO>(io);
    ot->recv(r, bool_shares, num_bits);
    io->flush();

    for(int j = 0; j < num_bits; j++){
        uint64_t *p = (uint64_t*)&r[j];
        sum += p[1];
    }

    return sum;
}