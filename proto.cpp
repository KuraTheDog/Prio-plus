#include "proto.h"

#define SERVER0_IP "127.0.0.1"
#define SERVER1_IP "127.0.0.1"

uint64_t bitsum_ot_sender(NetIO* const io, const bool* const shares, const bool* const valid, const size_t n){
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

    block* const b0_ot = new block[n];
    block* const b1_ot = new block[n];

    for(int i = 0; i < n; i++){
        uint64_t* const p = (uint64_t*)&b0_ot[i];
        p[0] = valid[i] ? 0:1;
        p[1] = b0[i];
        uint64_t* const q = (uint64_t*)&b1_ot[i];
        q[0] = p[0];
        q[1] = b1[i];
    }

    io->sync();
    IKNP<NetIO>* const ot = new IKNP<NetIO>(io);
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

uint64_t intsum_ot_sender(NetIO* const io, const uint32_t* const shares, const bool* const valid, const size_t n, const size_t num_bits){
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

    block* const b0_ot = new block[n*num_bits];
    block* const b1_ot = new block[n*num_bits];

    for(int i = 0; i < n*num_bits; i++){
        uint64_t* const p = (uint64_t*)&b0_ot[i];
        p[0] = bool_valid[i] ? 0:1;
        p[1] = b0[i];
        uint64_t* const q = (uint64_t*)&b1_ot[i];
        q[0] = p[0];
        q[1] = b1[i];
    }

    io->sync();
    IKNP<NetIO>* const ot = new IKNP<NetIO>(io);
    ot->send(b0_ot, b1_ot, n*num_bits);
    io->flush();

    delete[] b1_ot;
    delete[] b0_ot;
    delete ot;

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

uint64_t bitsum_ot_receiver(NetIO* const io, const bool* const shares, const size_t n){

    block* const r = new block[n];
    uint64_t ans[n];
    uint64_t sum = 0;
    io->sync();
    IKNP<NetIO>* const ot = new IKNP<NetIO>(io);
    ot->recv(r, shares, n);
    io->flush();

    for(int i = 0; i < n; i++){
        uint64_t* const p = (uint64_t*)&r[i];
        ans[i] =(p[0] == 0) ? p[1] : 0;
        if(ans[i])
            sum += ans[i];
    }

    delete[] r;
    delete ot;

    return sum;
}

uint64_t intsum_ot_receiver(NetIO* const io, const uint32_t* const shares, const size_t n, const size_t num_bits){

    block* const r = new block[n*num_bits];
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
    IKNP<NetIO>* const ot = new IKNP<NetIO>(io);
    ot->recv(r, bool_shares, n*num_bits);
    io->flush();

    for(int i = 0; i < n; i++){
        uint64_t valid = ((uint64_t*) &r[i*num_bits])[0];
        if(valid == 0){
            for(int j = 0; j < num_bits; j++){
                uint64_t* const p = (uint64_t*)&r[i*num_bits+j];
                sum += p[1];
            }
        }
    }

    delete[] r;
    delete ot;

    return sum;
}

// Functions to convert XOR shared input values to shares that add up
// Problem : Do they add up with Int_Modulus?
uint64_t xor_to_sum_share_sender(NetIO* const io, const uint32_t share, const size_t num_bits){
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


    block* const b0_ot = new block[num_bits];
    block* const b1_ot = new block[num_bits];

    for(int i = 0; i < num_bits; i++){
        uint64_t* const p = (uint64_t*)&b0_ot[i];
        p[0] = 0;
        p[1] = b0[i];
        uint64_t* const q = (uint64_t*)&b1_ot[i];
        q[0] = p[0];
        q[1] = b1[i];
    }

    io->sync();
    IKNP<NetIO>* const ot = new IKNP<NetIO>(io);
    ot->send(b0_ot, b1_ot, num_bits);
    io->flush();

    delete[] b1_ot;
    delete[] b0_ot;
    delete ot;

    uint64_t local_share = 0;
    for(int j = 0; j < num_bits; j++)
        local_share += r[j];

    sum += local_share;

    return sum;
}

uint64_t xor_to_sum_share_receiver(NetIO* const io, const uint32_t share, const size_t num_bits){
    block* const r = new block[num_bits];
    bool bool_shares[num_bits];

    uint32_t num = share;
    // std::cout << "Share : " << num << " num bits " << num_bits << std::endl;
    for(int j = 0; j < num_bits; j++){
        bool_shares[j] = num%2;
        num = num >> 1;
    }

    uint64_t sum = 0;
    io->sync();
    IKNP<NetIO>* const ot = new IKNP<NetIO>(io);
    ot->recv(r, bool_shares, num_bits);
    io->flush();
    delete ot;

    for(int j = 0; j < num_bits; j++){
        uint64_t* const p = (uint64_t*)&r[j];
        sum += p[1];
    }

    return sum;
}

void set_block(block &b, const bool f){
    uint64_t* const p = (uint64_t*) &b;
    p[0] = 0;
    p[1] = (f ? 1 : 0);
}

void block_to_boolean(const block* const B, bool* const b, const int len){
    for(int i = 0; i < len; i++){
        uint64_t* p = (uint64_t*) &(B[i]);
        b[i] = (p[1] == 1) ? true : false;
    }
}

void print_block(const block var) {
    uint64_t v64val[2];
    memcpy(v64val, &var, sizeof(v64val));
    printf("%.16llx %.16llx\n", v64val[1], v64val[0]);
}

// Ref : https://crypto.stackexchange.com/questions/41651/what-are-the-ways-to-generate-beaver-triples-for-multiplication-gate
BooleanBeaverTriple* gen_boolean_beaver_triples(const int server_num, const int m){
    PRG prg;
    BooleanBeaverTriple* ans = new BooleanBeaverTriple[m];
    bool x[m], y[m], z[m], r[m], b[m];
    prg.random_bool(x, m);
    prg.random_bool(y, m);
    prg.random_bool(r, m);

    // B is 16 hex bits? so we can fit 64 bools in one block?
    block b0[m], b1[m], B[m];

    for(int i = 0; i < m; i++){
        set_block(b0[i], r[i]);
        set_block(b1[i], (x[i] != r[i])); // r[i] XOR x[i]
    }

    NetIO* io1 = new NetIO(server_num == 0 ? nullptr : SERVER0_IP, 60051, true);
    NetIO* io2 = new NetIO(server_num == 1 ? nullptr : SERVER1_IP, 60052, true);

    if(server_num == 0){
        io1->sync();
        IKNP<NetIO> ot1(io1);
        // std::cout << "OT1 send" << std::endl;
        ot1.send(b0,b1,m);
        io1->flush();

        io2->sync();
        IKNP<NetIO> ot2(io2);
        // std::cout << "OT2 recv" << std::endl;
        ot2.recv(B,y,m);
        io2->flush();
    }
    else if(server_num == 1){
        io1->sync();
        IKNP<NetIO> ot1(io1);
        // std::cout << "OT1 recv" << std::endl;
        ot1.recv(B,y,m);
        io1->flush();

        io2->sync();
        IKNP<NetIO> ot2(io2);
        // std::cout << "OT2 send" << std::endl;
        ot2.send(b0,b1,m);
        io2->flush();
    }

    delete io1;
    delete io2;

    block_to_boolean(B,b,m);

    for(int i = 0; i < m ; i++){
        // b = r' ^ x'y
        z[i] = b[i] != (r[i] != (x[i] and y[i])); // z[i] = r_A xor x_A.y_A xor r_B xor x_B.y_A
        // std::cout << x[i] << " " << y[i] << " " <<  z[i] << std::endl;
    }

    for(int i = 0 ; i < m ; i++){
        ans[i].set(x[i], y[i], z[i]);
    }

    return ans;
}
