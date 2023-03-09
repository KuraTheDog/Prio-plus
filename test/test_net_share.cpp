/*
Tests out net_share.cpp

Creates "sender" and "receiver", and sends a bunch of struct.h objects from sender to receiver.
*/

#undef NDEBUG
#include <assert.h>
#include <iostream>

#include "utils_test_connect.h"
#include "../constants.h"
#include "../fmpz_utils.h"
#include "../net_share.h"

void run_sender(int sockfd) {
    int n;
    fmpz_t number; fmpz_init(number);

    bool b = true;
    n = send_bool(sockfd, b);
    std::cout << "send bool \tsize: " << n << " \tval: " << b << std::endl;

    size_t nbool = 15;
    bool b_arr[nbool];
    for (unsigned int i = 0; i < nbool; i++)
        b_arr[i] = i % 2;
    n = send_bool_batch(sockfd, b_arr, nbool);
    std::cout << "send bool arr \tsize: " << n << " \tval: ";
    for (unsigned int i = 0; i < nbool; i++) {
        if (i > 0) std::cout << ", ";
        std::cout << b_arr[i];
    }
    std::cout << std::endl;

    int in = 13579;
    n = send_int(sockfd, in);
    std::cout << "send int, \tsize: " << n << " \tval: " << in << std::endl;

    double d = 3.14159;
    n = send_double(sockfd, d);
    std::cout << "send double, \tsize: " << n << " \tval: " << d << std::endl;

    size_t sz = 92;
    n = send_size(sockfd, sz);
    std::cout << "send size, \tsize: " << n << " \tval: " << sz << std::endl;

    uint32_t i = 1 << 20;
    n = send_uint32(sockfd, i);
    std::cout << "send uint32_t \tsize: " << n << " \tval: " << i << std::endl;

    uint64_t j = 1ul << 50;
    n = send_uint64(sockfd, j);
    std::cout << "send uint64_t \tsize: " << n << " \tval: " << j << std::endl;

    uint64_t batch64[2] = {123, 345};
    n = send_uint64_batch(sockfd, batch64, 2);
    std::cout << "send uint64[] \tsize: " << n << " \tval: " << batch64[0] << ", " << batch64[1] << std::endl;

    std::string s = "Hello World!";
    n = send_string(sockfd, s);
    std::cout << "send string \tsize: " << n << " \tval: " << s << std::endl;

    // Small number
    fmpz_set_d(number, 12345);
    n = send_fmpz(sockfd, number);
    std::cout << "send fmpz \tsize: " << n << " \tval: ";
    fmpz_print(number); std::cout << std::endl;

    fmpz_t* arr; new_fmpz_array(&arr, 2);
    fmpz_set_ui(arr[0], 123454321);
    fmpz_set_ui(arr[1], 987000789);
    n = send_fmpz_batch(sockfd, arr, 2);
    std::cout << "send fmpz[] \tsize: " << n << " \tval: "; fmpz_print(arr[0]);
    std::cout << ", "; fmpz_print(arr[1]); std::cout << std::endl;
    clear_fmpz_array(arr, 2);

    /*  Does not work, currently using fixed fmpz size
    // Large unsigned long
    fmpz_set_ui(number, 12345678900987654321ul);
    n = send_fmpz(sockfd, number);
    std::cout << "send fmpz \tsize: " << n << " \tval: ";
    fmpz_print(number); std::cout << std::endl;

    // Very large multi-limbed number
    fmpz_set_str(number, "3141592653589793238462643383279502884197169399", 10);
    n = send_fmpz(sockfd, number);
    std::cout << "send fmpz \tsize: " << n << " \tval: ";
    fmpz_print(number); std::cout << std::endl;
    */

    BooleanBeaverTriple* btrip = new BooleanBeaverTriple(true, false, true);
    n = send_BooleanBeaverTriple(sockfd, btrip);
    std::cout << "send btriple \tsize: " << n << " \tvals: " << btrip->a << ", " << btrip->b << ", " << btrip->c << std::endl;
    delete btrip;

    const BeaverTriple* const trip = NewBeaverTriple();
    n = send_BeaverTriple(sockfd, trip);
    std::cout << "send triple \tsize: " << n << " \tvals: ";
    fmpz_print(trip->A); std::cout << ", ";
    fmpz_print(trip->B); std::cout << ", ";
    fmpz_print(trip->C); std::cout << std::endl;
    delete trip;

    ClientPacket* const packet = new ClientPacket(1);
    fmpz_set_si(packet->f0_s, 10);
    n = send_ClientPacket(sockfd, packet, 1);
    std::cout << "send packet \tsize: " << n << " \tf0_s: ";
    fmpz_print(packet->f0_s); std::cout << std::endl;
    delete packet;

    ClientPacket* const packet2 = new ClientPacket(2);
    fmpz_set_si(packet2->f0_s, 20);
    n = send_ClientPacket(sockfd, packet2, 2);
    std::cout << "send packet \tsize: " << n << " \tf0_s = ";
    fmpz_print(packet2->f0_s); std::cout << std::endl;
    delete packet2;

    // size_t edasize = 8;
    // EdaBit* b0 = new EdaBit(edasize);
    // EdaBit* b1 = new EdaBit(edasize);
    // makeLocalEdaBit(b0, b1, edasize);
    // n = send_EdaBit(sockfd, b0, edasize);
    // std::cout << "send EdaBit " << edasize << " \tsize: " << n << std::endl;
    // b0->print();
    // delete b0;
    // delete b1;

    // Poly
    // fmpz_set_ui(number, 100);
    // fmpz_mod_poly_t f; fmpz_mod_poly_init(f, number);
    // fmpz_mod_poly_randtest(f, seed, 5);
    // n = send_poly(sockfd, f);
    // std::cout << "send X \tsize: " << n << " \tpoly: ";
    // fmpz_mod_poly_print_pretty(f, "x"); std::cout << std::endl;

    flint_rand_t this_seed; flint_randinit(this_seed);
    n = send_seed(sockfd, this_seed);
    std::cout << "send seed \tsize: " << n << " \tnext random: ";
    fmpz_randm(number, this_seed, Int_Modulus);
    fmpz_print(number); std::cout << std::endl;
    flint_randclear(this_seed);

    // Sanity: sending numbers still works
    fmpz_set_d(number, 54321);
    n = send_fmpz(sockfd, number);
    std::cout << "send fmpz \tsize: " << n << " \tfmpz: ";
    fmpz_print(number); std::cout << std::endl;

    fmpz_clear(number);
}

void run_receiver(int sockfd) {
    int n;
    fmpz_t number; fmpz_init(number);

    bool b;
    n = recv_bool(sockfd, b);
    std::cout << "recv bool \tsize: " << n << " \tval: " << b << std::endl;

    const size_t nbool = 15;
    bool b_arr[nbool];
    n = recv_bool_batch(sockfd, b_arr, nbool);
    std::cout << "recv bool arr \tsize: " << n << " \tval: ";
    for (unsigned int i = 0; i < nbool; i++) {
        if (i > 0) std::cout << ", ";
        std::cout << b_arr[i];
    }
    std::cout << std::endl;

    int in;
    n = recv_int(sockfd, in);
    std::cout << "recv int \tsize: " << n << " \tval: " << in << std::endl;

    double d;
    n = recv_double(sockfd, d);
    std::cout << "recv double \tsize: " << n << " \tval: " << d << std::endl;

    size_t sz;
    n = recv_size(sockfd, sz);
    std::cout << "recv size \tsize: " << n << " \tval: " << sz << std::endl;

    uint32_t i;
    n = recv_uint32(sockfd, i);
    std::cout << "recv uint32_t \tsize: " << n << " \tval: " << i << std::endl;

    uint64_t j;
    n = recv_uint64(sockfd, j);
    std::cout << "recv uint64_t \tsize: " << n << " \tval: " << j << std::endl;

    uint64_t batch64[2];
    n = recv_uint64_batch(sockfd, batch64, 2);
    std::cout << "recv uint64[] \tsize: " << n << " \tval: " << batch64[0] << ", " << batch64[1] << std::endl;

    std::string s;
    n = recv_string(sockfd, s);
    std::cout << "recv string \tsize: " << n << " \tval: " << s << std::endl;

    n = recv_fmpz(sockfd, number);
    std::cout << "recv fmpz \tsize: " << n << " \tval: ";
    fmpz_print(number); std::cout << std::endl;

    fmpz_t* arr; new_fmpz_array(&arr, 2);
    n = recv_fmpz_batch(sockfd, arr, 2);
    std::cout << "recv fmpz[] \tsize: " << n << " \tval: "; fmpz_print(arr[0]);
    std::cout << ", "; fmpz_print(arr[1]); std::cout << std::endl;
    clear_fmpz_array(arr, 2);

    /*
    n = recv_fmpz(sockfd, number);
    std::cout << "recv fmpz \tsize: " << n << " \tval: ";
    fmpz_print(number); std::cout << std::endl;

    n = recv_fmpz(sockfd, number);
    std::cout << "recv fmpz \tsize: " << n << " \tval: ";
    fmpz_print(number); std::cout << std::endl;
    */

    BooleanBeaverTriple* btrip = new BooleanBeaverTriple();
    n = recv_BooleanBeaverTriple(sockfd, btrip);
    std::cout << "send btriple \tsize: " << n << " \tvals: " << btrip->a << ", " << btrip->b << ", " << btrip->c << std::endl;
    delete btrip;

    BeaverTriple* trip = new BeaverTriple();
    n = recv_BeaverTriple(sockfd, trip);
    std::cout << "recv triple \tsize: " << n << " \tvals: ";
    fmpz_print(trip->A); std::cout << ", ";
    fmpz_print(trip->B); std::cout << ", ";
    fmpz_print(trip->C); std::cout << std::endl;
    delete trip;

    ClientPacket* packet = new ClientPacket(1);
    n = recv_ClientPacket(sockfd, packet, 1);
    std::cout << "recv packet \tsize: " << n << " \tf0_s: ";
    fmpz_print(packet->f0_s); std::cout << std::endl;
    delete packet;

    ClientPacket* packet2 = new ClientPacket(2);
    n = recv_ClientPacket(sockfd, packet2, 2);
    std::cout << "recv packet \tsize: " << n << " \tf0_s: ";
    fmpz_print(packet2->f0_s); std::cout << std::endl;
    delete packet2;

    // size_t edasize = 8;
    // EdaBit* b0 = new EdaBit(edasize);
    // n = recv_EdaBit(sockfd, b0, edasize);
    // std::cout << "recv EdaBit " << edasize << " \tsize: " << n << std::endl;
    // b0->print();
    // delete b0;

    // fmpz_set_ui(number, 100);
    // fmpz_mod_poly_t f; fmpz_mod_poly_init(f, number);
    // n = recv_poly(sockfd, f);
    // std::cout << "recv X \tsize: " << n << " \tpoly: ";
    // fmpz_mod_poly_print_pretty(f, "x"); std::cout << std::endl;

    flint_rand_t this_seed; flint_randinit(this_seed);
    // offset from default seed
    fmpz_randm(number, this_seed, Int_Modulus);
    fmpz_randm(number, this_seed, Int_Modulus);
    n = recv_seed(sockfd, this_seed);
    std::cout << "recv seed \tsize: " << n << " \tnext random: ";
    fmpz_randm(number, this_seed, Int_Modulus);
    fmpz_print(number); std::cout << std::endl;
    flint_randclear(this_seed);

    n = recv_fmpz(sockfd, number);
    std::cout << "recv fmpz \tsize: " << n << " \tval: ";
    fmpz_print(number); std::cout << std::endl;

    fmpz_clear(number);
}

void test_swap(int sockfd, int idx, size_t N) {
    std::cout << "Player " << idx << " running swap test size " << N << std::endl;
    bool* buff = new bool[N];
    memset(buff, idx, N * sizeof(bool));

    int bytes = swap_bool_batch(sockfd, buff, N);
    std::cout << "bytes: " << bytes << std::endl;

    assert(buff[0] == 1);
    assert(buff[N - 1] == 1);

    delete[] buff;
}

int main(int argc, char** argv) {
    init_constants();

    int N = 15000000;
    if (argc >= 2) {
        N = atoi(argv[1]);
    }

    std::thread t0([&]() {
        int cli_sockfd = init_sender();

        // sleep();

        run_sender(cli_sockfd);

        test_swap(cli_sockfd, 0, N);

        close(cli_sockfd);
    });

    std::thread t1([&]() {
        int sockfd = init_receiver();
        int newsockfd = accept_receiver(sockfd);

        run_receiver(newsockfd);

        test_swap(newsockfd, 1, N);

        close(newsockfd);
        close(sockfd);
    });

    t0.join();
    t1.join();

    clear_constants();
    return 0;
}
