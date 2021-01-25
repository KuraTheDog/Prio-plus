/*
Tests out net_share.cpp

Creates "sender" and "receiver", and sends a bunch of struct.h objects from sender to receiver.

g++ -std=c++11 -o test_net_share test_net_share.cpp ../fmpz_utils.cpp ../share.cpp -lgmp -lflint -g && ./test_net_share
*/

#include <iostream>

#include "utils_test_connect.h"
#include "../constants.h"
#include "../fmpz_utils.h"
#include "../net_share.h"

void run_sender(int sockfd) {
    int n;
    fmpz_t number;
    fmpz_init(number);

    uint32_t i = 1 << 20;
    n = send_uint32(sockfd, i);
    std::cout << "send: size = " << n << ", uint32_t: " << i << std::endl;

    uint64_t j = 1ul << 50;
    n = send_uint64(sockfd, j);
    std::cout << "send: size = " << n << ", uint64_t: " << j << std::endl;

    // Small number
    fmpz_set_d(number, 12345);
    n = send_fmpz(sockfd, number);
    std::cout << "send: size = " << n << ", fmpz: ";
    fmpz_print(number); std::cout << std::endl;

    // Large unsigned long
    fmpz_set_ui(number, 12345678900987654321ul);
    n = send_fmpz(sockfd, number);
    std::cout << "send: size = " << n << ", fmpz: ";
    fmpz_print(number); std::cout << std::endl;

    // Very large multi-limbed number
    fmpz_set_str(number, "3141592653589793238462643383279502884197169399", 10);
    n = send_fmpz(sockfd, number);
    std::cout << "send: size = " << n << ", fmpz: ";
    fmpz_print(number); std::cout << std::endl;

    BeaverTriple* trip = NewBeaverTriple();
    n = send_BeaverTriple(sockfd, trip);
    std::cout << "send: size = " << n << ", triple: ";
    fmpz_print(trip->A); std::cout << ", ";
    fmpz_print(trip->B); std::cout << ", ";
    fmpz_print(trip->C); std::cout << std::endl;

    ClientPacket packet;
    init_client_packet(packet, 1, 2);
    n = send_ClientPacket(sockfd, packet);
    std::cout << "send: size = " << n << ", packet, N = " << packet->N << ", NWires = " << packet->NWires << std::endl;

    ClientPacket packet2;
    init_client_packet(packet2, 2, 3);
    n = send_ClientPacket(sockfd, packet2);
    std::cout << "send: size = " << n << ", packet, N = " << packet2->N << ", NWires = " << packet2->NWires << std::endl;

    EdaBit* b0 = new EdaBit(4);
    EdaBit* b1 = new EdaBit(4);
    makeLocalEdaBit(b0, b1, 4);
    n = send_EdaBit(sockfd, b0, 4);
    std::cout << "send EdaBit size = " << n << std::endl;
    b0->print();

    // Sanity: sending numbers still works
    fmpz_set_d(number, 54321);
    n = send_fmpz(sockfd, number);
    std::cout << "send: size = " << n << ", fmpz: ";
    fmpz_print(number); std::cout << std::endl;
}

void run_receiver(int newsockfd) {
    int n;
    fmpz_t number;
    fmpz_init(number);

    uint32_t i;
    n = recv_uint32(newsockfd, i);
    std::cout << "recv: size = " << n << ", uint32_t: " << i << std::endl;

    uint64_t j;
    n = recv_uint64(newsockfd, j);
    std::cout << "recv: size = " << n << ", uint64_t: " << j << std::endl;

    n = recv_fmpz(newsockfd, number);
    std::cout << "recv: size = " << n << ", fmpz: ";
    fmpz_print(number); std::cout << std::endl;

    n = recv_fmpz(newsockfd, number);
    std::cout << "recv: size = " << n << ", fmpz: ";
    fmpz_print(number); std::cout << std::endl;

    n = recv_fmpz(newsockfd, number);
    std::cout << "recv: size = " << n << ", fmpz: ";
    fmpz_print(number); std::cout << std::endl;

    BeaverTriple* trip = new BeaverTriple();
    n = recv_BeaverTriple(newsockfd, trip);
    std::cout << "recv: size = " << n << ", triple: ";
    fmpz_print(trip->A); std::cout << ", ";
    fmpz_print(trip->B); std::cout << ", ";
    fmpz_print(trip->C); std::cout << std::endl;

    ClientPacket packet = nullptr;
    n = recv_ClientPacket(newsockfd, packet);
    std::cout << "recv: size = " << n << ", packet, N = " << packet->N << ", NWires = " << packet->NWires << std::endl;

    ClientPacket packet2 = nullptr;
    n = recv_ClientPacket(newsockfd, packet2);
    std::cout << "recv: size = " << n << ", packet, N = " << packet2->N << ". NWires = " << packet2->NWires << std::endl;

    EdaBit* b = new EdaBit(4);
    n = recv_EdaBit(newsockfd, b, 4);
    std::cout << "recv EdaBit size = " << n << std::endl;
    b->print();

    n = recv_fmpz(newsockfd, number);
    std::cout << "recv: size = " << n << ", fmpz: ";
    fmpz_print(number); std::cout << std::endl;
}

int main(int argc, char** argv) {
    /* Initialize constants.h globals constants.
    Mostly for Int_Modulus, for use by share.h. */
    fmpz_init(Int_Modulus);
    fmpz_set_str(Int_Modulus,Int_Modulus_str.c_str(),16);
    flint_randinit(seed);

    /* set up receiver */
    int sockfd = init_receiver();

    /* Launch child to do sending */
    pid_t pid = fork();
    if (pid == 0) {
        int cli_sockfd = init_sender();
        // sleep();
        run_sender(cli_sockfd);
        close(cli_sockfd);
    } else if (pid > 0) {
        /* Do receiver handling */
        int newsockfd = accept_receiver(sockfd);

        run_receiver(newsockfd);

        close(newsockfd);
        close(sockfd);
    } else {
        error_exit("Failed to fork");
    }

    return 0;
}