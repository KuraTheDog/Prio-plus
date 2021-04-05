#include "ot.h"

#include "constants.h"
#include "net_share.h"  // lazy or OT triple, silent 
#include "utils.h"

#if OT_TYPE == EMP_IKNP

OT_Wrapper::OT_Wrapper(const char* address, const int port, const bool is_sender, const int sockfd, const int batch_size)
: is_sender(is_sender)
, io(new emp::NetIO(is_sender ? address : nullptr, port, true))
, ot(new emp::IKNP<emp::NetIO>(io))
{}

OT_Wrapper::~OT_Wrapper() {
    delete ot;
    delete io;
}

void OT_Wrapper::send(const uint64_t* const data0, const uint64_t* const data1,
        const size_t length) {
    if (!is_sender) error_exit("Error: Calling send on non-sender OT Wrapper");

    emp::block* const block0 = new emp::block[length];
    emp::block* const block1 = new emp::block[length];

    for (unsigned int i = 0; i < length; i++) {
        block0[i] = emp::makeBlock(0, data0[i]);
        block1[i] = emp::makeBlock(0, data1[i]);
        // std::cout << "Send[" << i << "] = (" << data0[i] << ", " << data1[i] << ")\n";
    }

    io->sync();
    ot->send(block0, block1, length);
    io->flush();

    delete[] block0;
    delete[] block1;
}

void OT_Wrapper::recv(uint64_t* const data, const bool* b, const size_t length) {
    emp::block* const block = new emp::block[length];
    if (is_sender) error_exit("Error: Calling recv on sender OT Wrapper");

    io->sync();
    ot->recv(block, b, length);
    io->flush();

    for (unsigned int i = 0; i < length; i++) {
        data[i] = *(uint64_t*)&block[i];
        // std::cout << "Recv[" << i << "][" << b[i] << "] = " << data[i] << std::endl;
    }

    delete[] block;
}

#elif OT_TYPE == LIBOTE_IKNP

OT_Wrapper::OT_Wrapper(const char* address, const int port, const bool is_sender, const int sockfd, const int batch_size)
: is_sender(is_sender)
{
    channel = osuCrypto::Session(
        ios, address, port, 
        is_sender ? osuCrypto::SessionMode::Server : osuCrypto::SessionMode::Client
    ).addChannel();

    prng = osuCrypto::PRNG(osuCrypto::sysRandomSeed());
}

OT_Wrapper::~OT_Wrapper() {
}

void OT_Wrapper::send(const uint64_t* const data0, const uint64_t* const data1,
        const size_t length) {
    if (!is_sender) error_exit("Error: Calling send on non-sender OT Wrapper");

    osuCrypto::IknpOtExtSender sender;
    std::vector<std::array<osuCrypto::block, 2>> messages(length);
    for (unsigned int i = 0; i < length; i++) {
        messages[i] = { osuCrypto::toBlock(data0[i]), 
                        osuCrypto::toBlock(data1[i]) };
        // std::cout << "Send[" << i << "] = (" << data0[i] << ", " << data1[i] << ")\n";
    }
    sender.sendChosen(messages, prng, channel);
}

void OT_Wrapper::recv(uint64_t* const data, const bool* b, const size_t length) {
    if (is_sender) error_exit("Error: Calling send on sender OT Wrapper");

    osuCrypto::IknpOtExtReceiver recver;
    osuCrypto::BitVector choices(length);
    for (unsigned int i = 0; i < length; i++)
        choices[i] = b[i];
    std::vector<osuCrypto::block> messages(length);
    recver.receiveChosen(choices, messages, prng, channel);

    for (unsigned int i = 0; i < length; i++) {
        data[i] = *(uint64_t*)&messages[i];
        // std::cout << "Recv[" << i << "][" << b[i] << "] = " << data[i] << std::endl;
    }
}

#elif OT_TYPE == LIBOTE_SILENT

OT_Wrapper::OT_Wrapper(const char* address, const int port, const bool is_sender, const int sockfd, const int batch_size)
: is_sender(is_sender)
, batch_size(batch_size)
, address(address)
, port(port)
, sockfd(sockfd)
{
    std::cout << "Making a " << (is_sender ? "sender" : "receiver") << " on port " << port << " with batch size " << batch_size << std::endl;
    prng = osuCrypto::PRNG(osuCrypto::sysRandomSeed());

    addPrecompute();
}

OT_Wrapper::~OT_Wrapper() {}

void OT_Wrapper::addPrecompute(const size_t n) {

    osuCrypto::IOService ios;
    osuCrypto::Channel channel;

    const int num_to_make = (n > batch_size ? n : batch_size);

    std::cout << "adding precompute of " << num_to_make << std::endl;

    if (is_sender) {
        osuCrypto::SilentOtExtSender sender;
        sender.configure(num_to_make);  // more options? is this needed? todo

        osuCrypto::Channel channel = osuCrypto::Session(
            ios, address, port, osuCrypto::SessionMode::Server
        ).addChannel();

        // build
        std::vector<std::array<osuCrypto::block, 2>> randMessages(num_to_make);
        sender.silentSend(randMessages, prng, channel);

        for (unsigned int i = 0; i < num_to_make; i++) {
            message_cache.push({
                *(uint64_t*)&randMessages[i][0],
                *(uint64_t*)&randMessages[i][1]
            });
        }
    } else {
        osuCrypto::SilentOtExtReceiver recver;
        recver.configure(num_to_make);  // more options? is this needed? todo

        osuCrypto::Channel channel = osuCrypto::Session(
            ios, address, port, osuCrypto::SessionMode::Client
        ).addChannel();

        // build
        osuCrypto::BitVector randChoice(num_to_make);
        std::vector<osuCrypto::block> gotMessages(num_to_make);
        recver.silentReceive(randChoice, gotMessages, prng, channel);

        for (unsigned int i = 0; i < num_to_make; i++) {
            choice_cache.push({
                randChoice[i],
                *(uint64_t*)&gotMessages[i]
            });
        }
    }
    std::cout << "OT Precompute done" << std::endl;
}

void OT_Wrapper::send(const uint64_t* const data0, const uint64_t* const data1,
        const size_t length) {
    if (!is_sender) error_exit("Error: Calling send on non-sender OT Wrapper");
    if (sockfd == -1)
        error_exit("Silent OT Wrapper requires a valid sockfd for online");

    if (message_cache.size() < length)
        addPrecompute(length);

    // work
    bool* const d = new bool[length];
    recv_bool_batch(sockfd, d, length);

    for (unsigned int i = 0; i < length; i++) {
        uint64_t r[2];
        std::tie(r[0], r[1]) = message_cache.front();
        message_cache.pop();

        uint64_t s0 = data0[i] ^ r[d[i]];
        uint64_t s1 = data1[i] ^ r[d[i] ^ 1];
        send_uint64(sockfd, s0);
        send_uint64(sockfd, s1);
    }
    delete[] d;
}

void OT_Wrapper::recv(uint64_t* const data, const bool* b, const size_t length) {
    if (is_sender) error_exit("Error: Calling send on sender OT Wrapper");

    if (sockfd == -1)
        error_exit("Silent OT Wrapper requires a valid sockfd for online");

    if (choice_cache.size() < length)
        addPrecompute(length);

    // work
    bool* const d = new bool[length];
    uint64_t* const rc = new uint64_t[length];
    for (unsigned int i = 0; i < length; i++) {
        bool c;
        std::tie(c, rc[i]) = choice_cache.front();
        choice_cache.pop();

        d[i] = b[i] ^ c;
    }
    send_bool_batch(sockfd, d, length);
    delete[] d;
    for (unsigned int i = 0; i < length; i++) {
        uint64_t s0, s1;
        recv_uint64(sockfd, s0);
        recv_uint64(sockfd, s1);
        data[i] = (b[i] == 0 ? s0 : s1) ^ rc[i];
    }
    delete[] rc;
}

#endif

uint64_t bitsum_ot_sender(OT_Wrapper* const ot, const bool* const shares, const bool* const valid, const size_t n){
    emp::PRG prg(emp::fix_key);

    uint64_t sum = 0;

    uint64_t* const b0 = new uint64_t[n];
    uint64_t* const b1 = new uint64_t[n];

    for (unsigned int i = 0; i < n; i++){

        if (valid[i]) {
            prg.random_data(&b0[i], sizeof(uint64_t));
            b1[i] = b0[i] + 1 - 2 * shares[i];
            sum += (UINT64_MAX - b0[i]) + 1 + shares[i];
        } else {
            b0[i] = 0;
            b1[i] = 0;
        }
    }

    ot->send(b0, b1, n);

    delete[] b1;
    delete[] b0;
    return sum;
}

uint64_t bitsum_ot_receiver(OT_Wrapper* const ot, const bool* const shares, const size_t n){
    uint64_t* const r = new uint64_t[n];
    uint64_t sum = 0;

    ot->recv(r, shares, n);

    for (unsigned int i = 0; i < n; i++)
        sum += r[i];

    delete[] r;

    return sum;
}

uint64_t* intsum_ot_sender(OT_Wrapper* const ot,  const uint64_t* const shares,
                           const bool* const valid, const size_t* const num_bits,
                           const size_t num_shares, const size_t num_values) {
    emp::PRG prg;

    size_t total_bits = 0;
    for (unsigned int j = 0; j < num_values; j++)
        total_bits += num_bits[j];

    uint64_t bool_share, r;
    uint64_t* sum = new uint64_t[num_values];
    memset(sum, 0, num_values * sizeof(uint64_t));

    uint64_t* const b0 = new uint64_t[num_shares * total_bits];
    uint64_t* const b1 = new uint64_t[num_shares * total_bits];

    size_t idx = 0;
    for (unsigned int i = 0; i < num_shares; i++) {
        // std::cout << "valid[" << i << "] = " << valid[i] << std::endl;
        for (unsigned int j = 0; j < num_values; j++) {
            uint64_t num = shares[i * num_values + j];
            // std::cout << "val[" << j << "] = " << num << std::endl;
            for (unsigned int k = 0; k < num_bits[j]; k++) {
                bool_share = num % 2;
                num = num >> 1;

                prg.random_data(&r, sizeof(uint64_t));

                if (valid[i]) {
                    b0[idx] = bool_share * (1ULL << k) - r;
                    b1[idx] = (1 - bool_share) * (1ULL << k) - r;

                    sum[j] += r;
                } else {
                    b0[idx] = 0;
                    b1[idx] = 0;
                }

                idx++;
            }
        }
    }

    ot->send(b0, b1, num_shares * total_bits);

    delete[] b0;
    delete[] b1;

    return sum;
}

uint64_t* intsum_ot_receiver(OT_Wrapper* const ot, const uint64_t* const shares,
                             const size_t* const num_bits,
                             const size_t num_shares, const size_t num_values) {
    size_t total_bits = 0;
    for (unsigned int j = 0; j < num_values; j++)
        total_bits += num_bits[j];

    uint64_t* const r = new uint64_t[num_shares * total_bits];
    bool* const bool_shares = new bool[num_shares * total_bits];
    uint64_t* sum = new uint64_t[num_values];
    memset(sum, 0, num_values * sizeof(uint64_t));

    size_t idx = 0;
    for (unsigned int i = 0; i < num_shares; i++) {
        for (unsigned int j = 0; j < num_values; j++) {
            uint64_t num = shares[i * num_values + j];
            // std::cout << "val[" << j << "] = " << num << std::endl;
            for (unsigned int k = 0; k < num_bits[j]; k++) {
                bool_shares[idx] = num % 2;
                num = num >> 1;
                idx++;
            }
        }
    }

    ot->recv(r, bool_shares, num_shares * total_bits);

    delete[] bool_shares;

    idx = 0;
    for (unsigned int i = 0; i < num_shares; i++) {
        for (unsigned int j = 0; j < num_values; j++) {
            for (unsigned int k = 0; k < num_bits[j]; k++) {
                sum[j] += r[idx];
                idx++;
            }
        }
    }
    delete[] r;
    return sum;
}

// Ref : https://crypto.stackexchange.com/questions/41651/what-are-the-ways-to-generate-beaver-triples-for-multiplication-gate
std::queue<BooleanBeaverTriple*> gen_boolean_beaver_triples(const int server_num, const unsigned int m, OT_Wrapper* const ot0, OT_Wrapper* const ot1){
    emp::PRG prg;
    std::queue<BooleanBeaverTriple*> ans;
    bool* x = new bool[m];
    bool* y = new bool[m];
    bool* z = new bool[m];
    bool* r = new bool[m];
    prg.random_bool(x, m);
    prg.random_bool(y, m);
    prg.random_bool(r, m);

    uint64_t* b0 = new uint64_t[m];
    uint64_t* b1 = new uint64_t[m];
    uint64_t* b = new uint64_t[m];

    for (unsigned int i = 0; i < m; i++){
        b0[i] = r[i];
        b1[i] = (x[i] != r[i]); // r[i] XOR x[i]
    }

    if(server_num == 0){
        ot0->send(b0, b1, m);
        ot1->recv(b, y, m);
    }
    else if(server_num == 1){
        ot0->recv(b, y, m);
        ot1->send(b0, b1, m);
    }

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
    delete[] b;
    delete[] x;
    delete[] y;
    delete[] z;
    delete[] r;

    return ans;
}

// Slow. OT per bit
// Not batched, but we also don't really want to do this
/*
BeaverTriple* generate_beaver_triple(const int serverfd, const int server_num, OT_Wrapper* const ot0, OT_Wrapper* const ot1) {

    // auto start = clock_start();

    BeaverTriple* triple = new BeaverTriple();
    // ceil(log_2 modulus)
    const size_t n = fmpz_clog_ui(Int_Modulus, 2);
    // For n > 128, need multiple blocks to represent a fmpz_t value
    if (n > 128) {
        std::cout << "Int_Modulus is " << n << " bits, which is more than 128, the maximum that the current generate_beaver_triple can handle." << std::endl;
        return triple;
    }

    emp::PRG prg(emp::fix_key);

    uint64_t r0_block[n], r1_block[n], s_block[n];
    // TODO: replace this with Random OT instead.
    // Note: Random OT seems to break. Running 20 varop 8 twice has the second run make bad triples. Future runs line up again. Desync issues?
    prg.random_block(r0_block, n);
    prg.random_block(r1_block, n);

    // Use emp::PRG random for this instead?
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
        ot0->recv(s_block, a_arr, n);
        ot1->send(r0_block, r1_block, n);
    } else {
        ot0->send(r0_block, r1_block, n);
        ot1->recv(s_block, a_arr, n);
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
*/

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
