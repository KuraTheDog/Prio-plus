#include "proto.h"

#include "constants.h"
#include "net_share.h"

uint64_t bitsum_ot_sender(NetIO* const io, const bool* const shares, const bool* const valid, const size_t n){
    PRG prg(fix_key);

    uint64_t* b0 = new uint64_t[n];
    uint64_t* b1 = new uint64_t[n];
    uint64_t* r = new uint64_t[n];
    uint64_t sum = 0;
    prg.random_data(b0,n*sizeof(uint64_t));

    for (unsigned int i = 0; i < n; i++){
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

    for (unsigned int i = 0; i < n; i++){
        uint64_t* const p = (uint64_t*)&b0_ot[i];
        p[0] = valid[i] ? 0:1;
        p[1] = b0[i];
        uint64_t* const q = (uint64_t*)&b1_ot[i];
        q[0] = p[0];
        q[1] = b1[i];
    }

    delete[] b0;
    delete[] b1;

    io->sync();
    IKNP<NetIO>* const ot = new IKNP<NetIO>(io);
    ot->send(b0_ot, b1_ot, n);
    io->flush();

    delete[] b1_ot;
    delete[] b0_ot;
    delete ot;

    for (unsigned int i = 0; i < n; i++){
        if(valid[i])
            sum += r[i];
    }

    delete[] r;

    return sum;
}

uint64_t intsum_ot_sender(NetIO* const io, const uint64_t* const shares, const bool* const valid, const size_t n, const size_t num_bits){
    PRG prg(fix_key);

    uint64_t* bool_shares = new uint64_t[n*num_bits];
    bool* bool_valid = new bool[n*num_bits];
    uint64_t* b0 = new uint64_t[n*num_bits];
    uint64_t* b1 = new uint64_t[n*num_bits];
    uint64_t* r = new uint64_t[n*num_bits];
    uint64_t sum = 0;
    prg.random_data(r,n*num_bits*sizeof(uint64_t));

    for (unsigned int i = 0; i < n; i++){
        uint64_t num = shares[i];
        // std::cout << "Share[" << i << "] " << num << "  valid " << valid[i] << std::endl;
        // std::cout << "Valid : " << valid[i] << " num bits " << num_bits << std::endl;
        for (unsigned int j = 0; j < num_bits; j++){
            // std::cout << num%2;
            bool_shares[i*num_bits + j] = num%2;
            bool_valid[i*num_bits + j] = valid[i];
            r[i*num_bits+j] = r[i*num_bits+j];
            // std::cout << "r[" << j << "] : " << r[j] << std::endl;
            num = num >> 1;
        }
        // std::cout << std::endl;
    }

    for (unsigned int i = 0; i < n; i++){
        for (unsigned int j = 0; j < num_bits; j++){
            b0[i*num_bits+j] = ((bool_shares[i*num_bits+j])*(1ULL << j) - r[i*num_bits+j]);
            // std::cout << "b0[" << j << "] = " << b0[i*num_bits+j] << std::endl;
            b1[i*num_bits+j] = ((1 - bool_shares[i*num_bits+j])*(1ULL << j) - r[i*num_bits+j]);
            // std::cout << "b1[" << j << "] = " << b1[i*num_bits+j] << std::endl;
        }
    }

    delete[] bool_shares;

    block* const b0_ot = new block[n*num_bits];
    block* const b1_ot = new block[n*num_bits];

    for (unsigned int i = 0; i < n*num_bits; i++){
        uint64_t* const p = (uint64_t*)&b0_ot[i];
        p[0] = bool_valid[i] ? 0:1;
        p[1] = b0[i];
        uint64_t* const q = (uint64_t*)&b1_ot[i];
        q[0] = p[0];
        q[1] = b1[i];
    }

    delete[] bool_valid;
    delete[] b0;
    delete[] b1;

    io->sync();
    IKNP<NetIO>* const ot = new IKNP<NetIO>(io);
    ot->send(b0_ot, b1_ot, n*num_bits);
    io->flush();

    delete[] b1_ot;
    delete[] b0_ot;
    delete ot;

    for (unsigned int i = 0; i < n; i++){
        if(valid[i]){
            uint64_t local_share = 0;
            for (unsigned int j = i*num_bits; j < (i+1)*num_bits; j++)
                local_share += r[j];

            sum += local_share;
        }
    }
    
    delete[] r;

    return sum;
}

uint64_t bitsum_ot_receiver(NetIO* const io, const bool* const shares, const size_t n){

    block* const r = new block[n];
    uint64_t sum = 0;
    io->sync();
    IKNP<NetIO>* const ot = new IKNP<NetIO>(io);
    ot->recv(r, shares, n);
    io->flush();

    for (unsigned int i = 0; i < n; i++){
        uint64_t* const p = (uint64_t*)&r[i];
        sum += ((p[0] == 0) ? p[1] : 0);
    }

    delete[] r;
    delete ot;

    return sum;
}

uint64_t intsum_ot_receiver(NetIO* const io, const uint64_t* const shares, const size_t n, const size_t num_bits){

    block* const r = new block[n*num_bits];
    bool* bool_shares = new bool[n*num_bits];

    for (unsigned int i = 0; i < n; i++){
        uint64_t num = shares[i];
        // std::cout << "Share[" << i << "] = " << num << " num bits " << num_bits << std::endl;
        for (unsigned int j = 0; j < num_bits; j++){
            // std::cout << num%2;
            bool_shares[i*num_bits + j] = num%2;
            num = num >> 1;
        }
        // std::cout << std::endl;
    }

    uint64_t sum = 0;
    io->sync();
    IKNP<NetIO>* const ot = new IKNP<NetIO>(io);
    ot->recv(r, bool_shares, n*num_bits);
    io->flush();

    delete[] bool_shares;

    for (unsigned int i = 0; i < n; i++){
        uint64_t valid = ((uint64_t*) &r[i*num_bits])[0];
        if(valid == 0){
            for (unsigned int j = 0; j < num_bits; j++){
                uint64_t* const p = (uint64_t*)&r[i*num_bits+j];
                // std::cout << "sum += p[1] at " << j << " = " << p[1] << std::endl;
                sum += p[1];
            }
        }
    }

    delete[] r;
    delete ot;

    return sum;
}

// Functions to convert XOR shared input values to shares that add up
// Add up base 2^64
uint64_t xor_to_sum_share_sender(NetIO* const io, const uint64_t share, const size_t num_bits){
    PRG prg(fix_key);

    uint64_t* bool_shares = new uint64_t[num_bits];
    uint64_t* b0 = new uint64_t[num_bits];
    uint64_t* b1 = new uint64_t[num_bits];
    uint64_t* r = new uint64_t[num_bits];
    uint64_t sum = 0;
    prg.random_data(r,num_bits*sizeof(uint64_t));


    uint64_t num = share;
    // std::cout << "Share : " << num << "  valid " << valid[i] << std::endl;
    // std::cout << "Valid : " << valid[i] << " num bits " << num_bits << std::endl;
    for (unsigned int j = 0; j < num_bits; j++){
        // std::cout << num%2 << std::endl;
        bool_shares[j] = num%2;

        // std::cout << "r[" << j << "] : " << r[j] << std::endl;
        num = num >> 1;
    }

    for (unsigned int j = 0; j < num_bits; j++){
        b0[j] = ((bool_shares[j])*(1ULL << j) - r[j]);
        b1[j] = ((1 - bool_shares[j])*(1ULL << j) - r[j]);
    }

    delete[] bool_shares;

    block* const b0_ot = new block[num_bits];
    block* const b1_ot = new block[num_bits];

    for (unsigned int i = 0; i < num_bits; i++){
        uint64_t* const p = (uint64_t*)&b0_ot[i];
        p[0] = 0;
        p[1] = b0[i];
        uint64_t* const q = (uint64_t*)&b1_ot[i];
        q[0] = p[0];
        q[1] = b1[i];
    }

    delete[] b0;
    delete[] b1;

    io->sync();
    IKNP<NetIO>* const ot = new IKNP<NetIO>(io);
    ot->send(b0_ot, b1_ot, num_bits);
    io->flush();

    delete[] b1_ot;
    delete[] b0_ot;
    delete ot;

    for (unsigned int j = 0; j < num_bits; j++)
        sum += r[j];

    delete[] r;

    return sum;
}

uint64_t xor_to_sum_share_receiver(NetIO* const io, const uint64_t share, const size_t num_bits){
    block* const r = new block[num_bits];
    bool* bool_shares = new bool[num_bits];

    uint64_t num = share;
    // std::cout << "Share : " << num << " num bits " << num_bits << std::endl;
    for (unsigned int j = 0; j < num_bits; j++){
        bool_shares[j] = num%2;
        num = num >> 1;
    }

    uint64_t sum = 0;
    io->sync();
    IKNP<NetIO>* const ot = new IKNP<NetIO>(io);
    ot->recv(r, bool_shares, num_bits);
    io->flush();

    delete[] bool_shares;
    delete ot;

    for (unsigned int j = 0; j < num_bits; j++){
        uint64_t* const p = (uint64_t*)&r[j];
        sum += p[1];
    }

    delete[] r;

    return sum;
}

void set_block(block &b, const bool f) {
    uint64_t* const p = (uint64_t*) &b;
    p[0] = 0;
    p[1] = (f ? 1 : 0);
}

void block_to_boolean(const block* const B, bool* const b, const unsigned int len){
    for (unsigned int i = 0; i < len; i++){
        uint64_t* p = (uint64_t*) &(B[i]);
        b[i] = (p[1] == 1) ? true : false;
    }
}

void print_block(const block var) {
    uint64_t v64val[2];
    memcpy(v64val, &var, sizeof(v64val));
    printf("%.16llux %.16llux\n", (unsigned long long) v64val[1], (unsigned long long) v64val[0]);
}

// Ref : https://crypto.stackexchange.com/questions/41651/what-are-the-ways-to-generate-beaver-triples-for-multiplication-gate
std::queue<BooleanBeaverTriple*> gen_boolean_beaver_triples(const int server_num, const unsigned int m, NetIO* const io0, NetIO* const io1){
    PRG prg;
    std::queue<BooleanBeaverTriple*> ans;
    bool* x = new bool[m];
    bool* y = new bool[m];
    bool* z = new bool[m];
    bool* r = new bool[m];
    bool* b = new bool[m];
    prg.random_bool(x, m);
    prg.random_bool(y, m);
    prg.random_bool(r, m);

    block* b0 = new block[m];
    block* b1 = new block[m];
    block* B = new block[m];

    for (unsigned int i = 0; i < m; i++){
        set_block(b0[i], r[i]);
        set_block(b1[i], (x[i] != r[i])); // r[i] XOR x[i]
    }

    if(server_num == 0){
        io0->sync();
        IKNP<NetIO> ot0(io0);
        ot0.send(b0, b1, m);
        io0->flush();

        io1->sync();
        IKNP<NetIO> ot1(io1);
        ot1.recv(B, y, m);
        io1->flush();
    }
    else if(server_num == 1){
        io0->sync();
        IKNP<NetIO> ot0(io0);
        ot0.recv(B, y, m);
        io0->flush();

        io1->sync();
        IKNP<NetIO> ot1(io1);
        ot1.send(b0, b1, m);
        io1->flush();
    }

    block_to_boolean(B, b, m);

    for (unsigned int i = 0; i < m ; i++){
        // b = r' ^ x'y
        z[i] = b[i] != (r[i] != (x[i] and y[i])); // z[i] = r_A xor x_A.y_A xor r_B xor x_B.y_A
        // std::cout << x[i] << " " << y[i] << " " <<  z[i] << std::endl;
    }

    for (unsigned int i = 0 ; i < m ; i++){
        ans.push(new BooleanBeaverTriple(x[i], y[i], z[i]));
    }

    delete[] b0;
    delete[] b1;
    delete[] B;
    delete[] x;
    delete[] y;
    delete[] z;
    delete[] r;
    delete[] b;

    return ans;
}

// Slow. OT per bit
// TODO: batched creation to reduce rounds.
BeaverTriple* generate_beaver_triple(const int serverfd, const int server_num, NetIO* const io0, NetIO* const io1) {

    // auto start = clock_start();

    BeaverTriple* triple = new BeaverTriple();
    // ceil(log_2 modulus)
    const size_t n = fmpz_clog_ui(Int_Modulus, 2);
    // For n > 128, need multiple blocks to represent a fmpz_t value
    if (n > 128) {
        std::cout << "Int_Modulus is " << n << " bits, which is more than 128, the maximum that the current generate_beaver_triple can handle." << std::endl;
        return triple;
    }
    
    PRG prg(fix_key);

    // Random offset, for debug
    if (server_num == 0) {
        block junk[10];
        prg.random_block(junk, 10);
        fmpz_randm(triple->A, seed, Int_Modulus);
    }

    block r0_block[n], r1_block[n], s_block[n];
    // TODO: replace this with Random OT instead.
    // Note: Random OT breaks. Running 20 varop 8 twice has the second run make bad triples. Future runs line up again. Desync issues?
    prg.random_block(r0_block, n);
    prg.random_block(r1_block, n);

    // Use prg random for this instead?
    fmpz_randm(triple->A, seed, Int_Modulus);
    fmpz_randm(triple->B, seed, Int_Modulus);

    bool a_arr[n];
    for (unsigned int i = 0; i < n; i++)
        a_arr[i] = fmpz_tstbit(triple->A, i);

    // OT random r0/r1 based on bits of a into s
    // s[i] = (ai == 1) ? r1[i] : r0[i]

    // std::cout << "setup : " << (((float)time_from(start))/CLOCKS_PER_SEC) << std::endl;
    // start = clock_start();


    if (server_num == 0) {
        io0->sync();
        IKNP<NetIO> ot0(io0);
        ot0.recv(s_block, a_arr, n);
        io0->flush();

        io1->sync();
        IKNP<NetIO> ot1(io1);
        ot1.send(r0_block, r1_block, n);
        io1->flush();
    } else {
        io0->sync();
        IKNP<NetIO> ot0(io0);
        ot0.send(r0_block, r1_block, n);
        io0->flush();

        io1->sync();
        IKNP<NetIO> ot1(io1);
        ot1.recv(s_block, a_arr, n);
        io1->flush();
    }

    // std::cout << "OT timing: " << (((float)time_from(start))/CLOCKS_PER_SEC) << std::endl;

    fmpz_t t; fmpz_init(t); fmpz_zero(t);  // sum 2^i ti
    fmpz_t q; fmpz_init(q); fmpz_zero(q);  // sum 2^i r0

    fmpz_t r0; fmpz_init(r0);  // r0_block[i]
    fmpz_t r1; fmpz_init(r1);  // r1_block[i]
    fmpz_t d; fmpz_init(d);

    fmpz_t pow; fmpz_init_set_si(pow, 1);  // 2^i
    for (unsigned int i = 0; i < n; i++) {
        fmpz_from_block(r0, r0_block[i], n);
        fmpz_mod(r0, r0, Int_Modulus);
        fmpz_from_block(r1, r1_block[i], n);
        fmpz_mod(r1, r1, Int_Modulus);

        // d = r0 - r1 + b, and swap
        fmpz_sub(d, r0, r1);
        fmpz_add(d, d, triple->B);
        fmpz_mod(d, d, Int_Modulus);
        send_fmpz(serverfd, d);

        fmpz_addmul(q, r0, pow); // r0
        fmpz_mod(q, q, Int_Modulus);

        fmpz_mul_ui(pow, pow, 2);
    }

    fmpz_set_si(pow, 1);  // 2^i
    fmpz_t s; fmpz_init(s);    // s[i]
    fmpz_t ti; fmpz_init(ti);
    for (unsigned int i = 0; i < n; i++) {
        recv_fmpz(serverfd, d);

        // s = r(ai)' = r0' - ai (r0' - r1')
        fmpz_from_block(s, s_block[i], n);
        fmpz_mod(s, s, Int_Modulus);

        // t = s + ai d' = r0' + ai b'
        fmpz_set(ti, s);
        if (a_arr[i]) {
            fmpz_add(ti, ti, d);
            fmpz_mod(ti, ti, Int_Modulus);
        }

        fmpz_addmul(t, ti, pow); // r0' + ai b'
        fmpz_mod(t, t, Int_Modulus);

        fmpz_mul_ui(pow, pow, 2);
    }

    // So t - q' = a b', and vice versa

    // c = a * b + t - q
    fmpz_mul(triple->C, triple->A, triple->B);
    fmpz_add(triple->C, triple->C, t);
    fmpz_sub(triple->C, triple->C, q);
    fmpz_mod(triple->C, triple->C, Int_Modulus);

    return triple;
}

// Simple beaver triple generation. Fast but unsafe.
BeaverTriple* generate_beaver_triple_lazy(const int serverfd, const int server_num) {
    BeaverTriple* triple = new BeaverTriple();

    if (server_num == 0) {
        BeaverTriple* other_triple = new BeaverTriple();
        fmpz_randm(triple->A, seed, Int_Modulus);
        fmpz_randm(triple->B, seed, Int_Modulus);
        fmpz_randm(triple->C, seed, Int_Modulus);
        fmpz_randm(other_triple->A, seed, Int_Modulus);
        fmpz_randm(other_triple->B, seed, Int_Modulus);

        fmpz_add(other_triple->C, triple->A, other_triple->A);
        fmpz_t tmp; fmpz_init(tmp);
        fmpz_add(tmp, triple->B, other_triple->B);
        fmpz_mul(other_triple->C, other_triple->C, tmp);
        fmpz_sub(other_triple->C, other_triple->C, triple->C);
        fmpz_mod(other_triple->C, other_triple->C, Int_Modulus);

        send_BeaverTriple(serverfd, other_triple);

        fmpz_clear(tmp);
        delete other_triple;
    } else {
        recv_BeaverTriple(serverfd, triple);
    }

    return triple;
}
