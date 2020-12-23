/*
Tests out net_share.cpp

Creates "sender" and "reciever", and sends a bunch of struct.h objects from sender to reciever.

g++ -std=c++11 -o test_net_share test_net_share.cpp -lgmp -lflint && ./test_net_share

Note: fmpz_print seems to append an extra digit, which seems to the # of digits in the original number.
*/

#include <iostream>

#include <unistd.h>
#include <arpa/inet.h>
#include <netinet/in.h>
#include <sys/socket.h>

# include "../net_share.h"

int PORT = 8888;

void error(const char *msg){
    perror(msg);
    exit(1);
}

int init_sender() {
    std::cout << "snd: start" << std::endl;

    int sockfd;
    struct sockaddr_in addr;

    sockfd = socket(AF_INET, SOCK_STREAM, 0);
    if (sockfd < 0)
        error("snd: ERROR opening socket");

    bzero((char *) &addr, sizeof(addr));
    addr.sin_port = htons(PORT);
    addr.sin_family = AF_INET;
    inet_pton(AF_INET, "127.0.0.1", &addr.sin_addr);

    if (connect(sockfd, (sockaddr*) &addr, sizeof(addr)) < 0)
        error("ERROR on connect");

    return sockfd;
}

int init_receiver() {
    std::cout << "rec: start" << std::endl;

    int sockfd;
    socklen_t snd_len;
    struct sockaddr_in snd_addr, rec_addr;

    sockfd = socket(AF_INET, SOCK_STREAM, 0);
    if (sockfd < 0)
        error("rec: ERROR opening socket");

    bzero((char *) &rec_addr, sizeof(rec_addr));
    rec_addr.sin_family = AF_INET;
    rec_addr.sin_addr.s_addr = INADDR_ANY;
    rec_addr.sin_port = htons(PORT);
    if (bind(sockfd, (struct sockaddr *) &rec_addr, sizeof(rec_addr)) < 0)
        error("rec: ERROR on binding");

    if(listen(sockfd, 2) < 0)
        error("rec: ERROR on listen");
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

    printf("rec: got connection from %s port %d\n",
           inet_ntoa(snd_addr.sin_addr), ntohs(snd_addr.sin_port));
    return newsockfd;
}

void run_sender(int sockfd) {
    ShareSender share_sender(sockfd);

    fmpz_t number;

    fmpz_init(number);
    fmpz_set_d(number, 12345);
    std::cout << "send: fmpz: " << fmpz_print(number) << std::endl;
    share_sender.fmpz(number);

    fmpz_set_str(number, "314159265358979323846264338327950", 10);
    std::cout << "send: fmpz: " << fmpz_print(number) << std::endl;
    share_sender.fmpz(number);

    BeaverTriple* trip = NewBeaverTriple();
    std::cout << "send: triple: " << fmpz_print(trip->A) << ", " << fmpz_print(trip->B) << ", " << fmpz_print(trip->C) << std::endl;
    share_sender.BeaverTriple(trip);

    // TODO: client_packet is currently bugged. Needs more testing.
    // client_packet* packet;
    // init_client_packet(packet, 10);
    // std::cout << "snd: Sending packet,  N = " << packet->N << std::endl;
    // share_sender.client_packet(packet);
}

void run_reciever(int newsockfd) {
    ShareReciever share_reciever(newsockfd);

    fmpz_t number;
    share_reciever.fmpz(number);
    std::cout << "recv: fmpz: " << fmpz_print(number) << std::endl;

    share_reciever.fmpz(number);
    std::cout << "recv: fmpz: " << fmpz_print(number) <<std::endl;

    BeaverTriple* trip = new BeaverTriple();
    share_reciever.BeaverTriple(trip);
    std::cout << "recv: triple: " << fmpz_print(trip->A) << ", " << fmpz_print(trip->B) << ", " << fmpz_print(trip->C) << std::endl;

    // TODO: client_packet is currently bugged. Needs more testing.
    // client_packet* packet = new client_packet();
    // share_reciever.client_packet(packet);
    // std::cout << "rec: Got packet, N = " << packet->N << std::endl;
}

int main(int argc, char** argv) {
    /* Initialize prio.h constants.
    Mostly for Int_Modulus, for use by share.h. */
    fmpz_init(Int_Modulus);
    fmpz_set_str(Int_Modulus,Int_Modulus_str.c_str(),16);
    fmpz_init(Int_Gen);
    fmpz_set_str(Int_Gen,Int_Gen_str.c_str(),16);
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