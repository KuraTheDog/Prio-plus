/*
Tests out net_share.cpp

Creates "sender" and "reciever", and sends a bunch of struct.h objects from sender to reciever.

g++ -std=c++11 -o test_net_share test_net_share.cpp ../fmpz_utils.cpp ../share.cpp -lgmp -lflint -g && ./test_net_share
*/

#include <iostream>

#include <unistd.h>
#include <arpa/inet.h>
#include <netinet/in.h>
#include <sys/socket.h>

#include "../fmpz_utils.h"
#include "../net_share.h"
#include "../prio.h"

int PORT = 8888;

void error(const char *msg){
    perror(msg);
    exit(1);
}

int init_sender() {
    std::cout << "send: start" << std::endl;

    int sockfd;
    struct sockaddr_in addr;

    sockfd = socket(AF_INET, SOCK_STREAM, 0);
    if (sockfd < 0)
        error("send: ERROR opening socket");

    bzero((char *) &addr, sizeof(addr));
    addr.sin_port = htons(PORT);
    addr.sin_family = AF_INET;
    inet_pton(AF_INET, "127.0.0.1", &addr.sin_addr);

    if (connect(sockfd, (sockaddr*) &addr, sizeof(addr)) < 0)
        error("ERROR on connect");

    return sockfd;
}

int init_receiver() {
    std::cout << "recv: start" << std::endl;

    int sockfd;
    struct sockaddr_in rec_addr;

    sockfd = socket(AF_INET, SOCK_STREAM, 0);
    if (sockfd < 0)
        error("recv: ERROR opening socket");

    bzero((char *) &rec_addr, sizeof(rec_addr));
    rec_addr.sin_family = AF_INET;
    rec_addr.sin_addr.s_addr = INADDR_ANY;
    rec_addr.sin_port = htons(PORT);
    if (bind(sockfd, (struct sockaddr *) &rec_addr, sizeof(rec_addr)) < 0)
        error("recv: ERROR on binding");

    if(listen(sockfd, 2) < 0)
        error("recv: ERROR on listen");
    return sockfd;
}

int accept_reciever(int sockfd) {
    int newsockfd;
    socklen_t snd_len;
    struct sockaddr_in snd_addr;
    snd_len = sizeof(snd_addr);
    newsockfd = accept(sockfd, (struct sockaddr *) & snd_addr, &snd_len);
    if (newsockfd < 0)
        error("ERROR on accept");

    printf("recv: got connection from %s port %d\n",
           inet_ntoa(snd_addr.sin_addr), ntohs(snd_addr.sin_port));
    return newsockfd;
}

void run_sender(int sockfd) {
    ShareSender share_sender(sockfd);

    int n;
    fmpz_t number;
    fmpz_init(number);

    // Small number
    fmpz_set_d(number, 12345);
    n = share_sender.fmpz(number);
    std::cout << "send: size = " << n << ", fmpz: ";
    fmpz_print(number); std::cout << std::endl;

    // Large unsigned long
    fmpz_set_ui(number, 12345678900987654321ul);
    n = share_sender.fmpz(number);
    std::cout << "send: size = " << n << ", fmpz: ";
    fmpz_print(number); std::cout << std::endl;

    // Very large multi-limbed number
    fmpz_set_str(number, "3141592653589793238462643383279502884197169399", 10);
    n = share_sender.fmpz(number);
    std::cout << "send: size = " << n << ", fmpz: ";
    fmpz_print(number); std::cout << std::endl;

    BeaverTriple* trip = NewBeaverTriple();
    n = share_sender.BeaverTriple(trip);
    std::cout << "send: size = " << n << ", triple: ";
    fmpz_print(trip->A); std::cout << ", ";
    fmpz_print(trip->B); std::cout << ", ";
    fmpz_print(trip->C); std::cout << std::endl;

    ClientPacket packet;
    init_client_packet(packet, 1, 2);
    n = share_sender.client_packet(packet);
    std::cout << "send: size = " << n << ", packet, N = " << packet->N << ", NWires = " << packet->NWires << std::endl;

    ClientPacket packet2;
    init_client_packet(packet2, 2, 3);
    n = share_sender.client_packet(packet2);
    std::cout << "send: size = " << n << ", packet, N = " << packet2->N << ", NWires = " << packet2->NWires << std::endl;

    // Sanity: sending numbers still works
    fmpz_set_d(number, 54321);
    n = share_sender.fmpz(number);
    std::cout << "send: size = " << n << ", fmpz: ";
    fmpz_print(number); std::cout << std::endl;
}

void run_reciever(int newsockfd) {
    ShareReciever share_reciever(newsockfd);

    int n;
    fmpz_t number;
    fmpz_init(number);

    n = share_reciever.fmpz(number);
    std::cout << "recv: size = " << n << ", fmpz: ";
    fmpz_print(number); std::cout << std::endl;

    n = share_reciever.fmpz(number);
    std::cout << "recv: size = " << n << ", fmpz: ";
    fmpz_print(number); std::cout << std::endl;

    n = share_reciever.fmpz(number);
    std::cout << "recv: size = " << n << ", fmpz: ";
    fmpz_print(number); std::cout << std::endl;

    BeaverTriple* trip = new BeaverTriple();
    n = share_reciever.BeaverTriple(trip);
    std::cout << "recv: size = " << n << ", triple: ";
    fmpz_print(trip->A); std::cout << ", ";
    fmpz_print(trip->B); std::cout << ", ";
    fmpz_print(trip->C); std::cout << std::endl;

    ClientPacket packet = nullptr;
    n = share_reciever.client_packet(packet);
    std::cout << "recv: size = " << n << ", packet, N = " << packet->N << ", NWires = " << packet->NWires << std::endl;

    ClientPacket packet2 = nullptr;
    n = share_reciever.client_packet(packet2);
    std::cout << "recv: size = " << n << ", packet, N = " << packet2->N << ". NWires = " << packet2->NWires << std::endl;
    std::cout << "recv: size = " << n << ", packet" << std::endl;

    n = share_reciever.fmpz(number);
    std::cout << "recv: size = " << n << ", fmpz: ";
    fmpz_print(number); std::cout << std::endl;
}

int main(int argc, char** argv) {
    /* Initialize prio.h globals constants.
    Mostly for Int_Modulus, for use by share.h. */
    fmpz_init(Int_Modulus);
    fmpz_set_str(Int_Modulus,Int_Modulus_str.c_str(),16);
    flint_randinit(seed);

    /* set up reciever */
    int sockfd = init_receiver();

    /* Launch child to do sending */
    pid_t pid = fork();
    if (pid == 0) {
        int cli_sockfd = init_sender();
        // sleep();
        run_sender(cli_sockfd);
        close(cli_sockfd);
    } else if (pid > 0) {
        /* Do reciever handling */
        int newsockfd = accept_reciever(sockfd);

        run_reciever(newsockfd);

        close(newsockfd);
        close(sockfd);
    } else {
        error("Failed to fork");
    }

    return 0;
}