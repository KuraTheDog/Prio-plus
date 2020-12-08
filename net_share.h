/* 
For sending fmtp_t objects over sockets.

Uses fdopen to use built in fmpz functions to send/recieve via files.

ShareSender: Sends over the socket.
ShareReciever: Recieves over the socket. Takes as argument the variable to set.

E.g. 

client: 
    ShareSender share_sender(socket);
    // produce number and trip somehow.
    share_sender.fmpz(number);
    share_sender.BeaverTriple(trip);
    ...

server: 
    ShareReciever share_reciever(socket);
    ...
    fmpz_t number;  // create new var to set
    share_reciever.fmpz(number);
    BeaverTriple* trip = new BeaverTriple();  // create new var to set
    share_reciever.BeaverTriple(trip);
    // now use number and trip
    ...

See test/test_net_share.cpp for a full example.

All send/recieve return positive on success, negative on failure.
Specifically, on first failure returns, and on success returns last success.

TODO: move this to net_share.cpp? Currently has linker issues when attempted.

TODO: client_packet is currently bugged. Needs more testing.

TODO: maybe hide/delete/comment out ones we don't use.
*/

#ifndef NET_SHARE_H
#define NET_SHARE_H

#include <stdio.h>
#include "share.h"

extern "C" {
    #include "flint/flint.h"
    #include "flint/fmpz.h"
};

class ShareSender {
    int sockfd;
    FILE* file;

    int send_fmpz(fmpz_t x){
        return fmpz_out_raw(file, x);
    }

    int send_int(int x){
        int32_t conv_x = htonl(x);
        return send(sockfd, &conv_x, sizeof(conv_x), 0);
    }

public: 

    ShareSender(int sockfd) {
        this->sockfd = sockfd;
        file = fdopen(sockfd, "w");
    }

    ~ShareSender() {
        fclose(file);
    }

    int fmpz(fmpz_t x) {
        int ret = send_fmpz(x);
        return ret;
    }

    int Cor(Cor *x) {
        int ret = send_fmpz(x->D);
        if (ret < 0) return ret;
        return send_fmpz(x->E);
    }

    int CorShare(CorShare *x) {
        int ret = send_fmpz(x->shareD);
        if (ret < 0) return ret;
        return send_fmpz(x->shareE);
    }

    // Note: This is different from a ClientPacket
    // TODO: client_packet is currently bugged. Needs more testing.
    int client_packet(client_packet *x) {
        int N = x->N;
        int ret = send_int(N);
        if (ret < 0) return ret;

        int i;
        for (i = 0; i < N; i++) {
            ret = send_fmpz(x->WireShares[i]);
            if (ret < 0) return ret;
        }

        ret = send_fmpz(x->f0_s);
        if (ret < 0) return ret;
        ret = send_fmpz(x->g0_s);
        if (ret < 0) return ret;
        ret = send_fmpz(x->h0_s);
        if (ret < 0) return ret;

        for (i = 0; i < N; i++) {
            ret = send_fmpz(x->h_points[i]);
            if (ret < 0) return ret;
        }
        return ret;
    }

    int BeaverTriple(BeaverTriple *x) {
        int ret = send_fmpz(x->A);
        if (ret < 0) return ret;
        ret = send_fmpz(x->B);
        if (ret < 0) return ret;
        return send_fmpz(x->C);
    }

    int BeaverTripleShare(BeaverTripleShare *x) {
        int ret = send_fmpz(x->shareA);
        if (ret < 0) return ret;
        ret = send_fmpz(x->shareB);
        if (ret < 0) return ret;
        return send_fmpz(x->shareC);
    }
};

class ShareReciever {
    int sockfd;
    FILE* file;

    int recv_fmpz(fmpz_t x){
        return fmpz_inp_raw(x, file);
    }

    int recv_int(int& x) {
        int recv_x = 0;
        int ret = recv(sockfd, &recv_x, sizeof(recv_x), 0);
        if (ret < 0) return ret;
        x = ntohl(recv_x);
        return ret;
    }

public: 

    ShareReciever(int sockfd) {
        this->sockfd = sockfd;
        file = fdopen(sockfd, "r");
    }

    ~ShareReciever() {
        fclose(file);
    }

    int fmpz(fmpz_t x) {
        return recv_fmpz(x);
    }

    int Cor(Cor *x) {
        int ret = recv_fmpz(x->D);
        if (ret < 0) return ret;
        return recv_fmpz(x->E);
    }

    int CorShare(CorShare *x) {
        int ret = recv_fmpz(x->shareD);
        if (ret < 0) return ret;
        return recv_fmpz(x->shareE);
    }

    // Note: This is different from a ClientPacket
    // TODO: client_packet is currently bugged. Needs more testing.
    int client_packet(client_packet *x) {
        int i, N;
        int ret = recv_int(N);
        if (ret < 0) return ret;
        x->N = N;

        new_fmpz_array(&x->WireShares, N);
        for (i = 0; i < N; i++) {
            ret = recv_fmpz(x->WireShares[i]);
            if (ret < 0) return ret;
        }

        ret = recv_fmpz(x->f0_s);
        if (ret < 0) return ret;
        ret = recv_fmpz(x->g0_s);
        if (ret < 0) return ret;
        ret = recv_fmpz(x->h0_s);
        if (ret < 0) return ret;

        new_fmpz_array(&x->h_points, N);
        for (i = 0; i < N; i++) {
            ret = recv_fmpz(x->h_points[i]);
            if (ret < 0) return ret;
        }
        return ret;
    }

    int BeaverTriple(BeaverTriple *x) {
        int ret = recv_fmpz(x->A);
        if (ret < 0) return ret;
        ret = recv_fmpz(x->B);
        if (ret < 0) return ret;
        return recv_fmpz(x->C);
    }

    int BeaverTripleShare(BeaverTripleShare *x) {
        int ret = recv_fmpz(x->shareA);
        if (ret < 0) return ret;
        ret = recv_fmpz(x->shareB);
        if (ret < 0) return ret;
        return recv_fmpz(x->shareC);
    }
};

#endif