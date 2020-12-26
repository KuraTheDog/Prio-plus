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

Returns total bytes sent on success, or first fail return of an internal step (typically 0 or negative)

TODO: move this to net_share.cpp? Currently has linker issues when attempted.

TODO: maybe hide/delete/comment out ones we don't use.
*/

#ifndef NET_SHARE_H
#define NET_SHARE_H

#include <stdio.h>
#include "fmpz_utils.h"
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
        // For consistency, just map ints to slong to fmpz, to use the same methedology.
        fmpz_t f;
        fmpz_init(f);
        fmpz_set_si(f, x);
        int ret = send_fmpz(f);
        fmpz_clear(f);
        return ret;
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
        return send_fmpz(x);
    }

    int Cor(Cor *x) {
        int total = 0, ret;
        ret = send_fmpz(x->D);
        if (ret <= 0) return ret; else total += ret;
        ret = send_fmpz(x->E);
        if (ret <= 0) return ret; else total += ret;
        return total;
    }

    int CorShare(CorShare *x) {
        int total = 0, ret;
        ret = send_fmpz(x->shareD);
        if (ret <= 0) return ret; else total += ret;
        ret = send_fmpz(x->shareE);
        if (ret <= 0) return ret; else total += ret;
        return total;
    }

    // Note: This is different from a ClientPacket
    int client_packet(client_packet *x) {
        int N = x->N, total = 0, ret;
        ret = send_int(N);
        if (ret <= 0) return ret; else total += ret;

        int i;
        for (i = 0; i < N; i++) {
            ret = send_fmpz(x->WireShares[i]);
            if (ret <= 0) return ret; else total += ret;
        }

        ret = send_fmpz(x->f0_s);
        if (ret <= 0) return ret; else total += ret;
        ret = send_fmpz(x->g0_s);
        if (ret <= 0) return ret; else total += ret;
        ret = send_fmpz(x->h0_s);
        if (ret <= 0) return ret; else total += ret;

        for (i = 0; i < N; i++) {
            ret = send_fmpz(x->h_points[i]);
            if (ret <= 0) return ret; else total += ret;
        }
        return total;
    }

    int BeaverTriple(BeaverTriple *x) {
        int total = 0, ret;
        ret = send_fmpz(x->A);
        if (ret <= 0) return ret; else total += ret;
        ret = send_fmpz(x->B);
        if (ret <= 0) return ret; else total += ret;
        ret = send_fmpz(x->C);
        if (ret <= 0) return ret; else total += ret;
        return total;
    }

    int BeaverTripleShare(BeaverTripleShare *x) {
        int total = 0, ret;
        ret = send_fmpz(x->shareA);
        if (ret <= 0) return ret; else total += ret;
        ret = send_fmpz(x->shareB);
        if (ret <= 0) return ret; else total += ret;
        ret = send_fmpz(x->shareC);
        return total;
    }
};

class ShareReciever {
    int sockfd;
    FILE* file;

    int recv_fmpz(fmpz_t x){
        return fmpz_inp_raw(x, file);
    }

    int recv_int(int& x) {
        fmpz_t f;
        int ret = recv_fmpz(f);
        x = fmpz_get_si(f);
        fmpz_clear(f);
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
        int total = 0, ret;
        ret = recv_fmpz(x->D);
        if (ret <= 0) return ret; else total += ret;
        ret = recv_fmpz(x->E);
        if (ret <= 0) return ret; else total += ret;
        return total;
    }

    int CorShare(CorShare *x) {
        int total = 0, ret;
        ret = recv_fmpz(x->shareD);
        if (ret <= 0) return ret; else total += ret;
        ret = recv_fmpz(x->shareE);
        if (ret <= 0) return ret; else total += ret;
        return total;
    }

    // Note: This is different from a ClientPacket
    int client_packet(client_packet *x) {
        int i, N, total = 0, ret;
        ret = recv_int(N);
        if (ret <= 0) return ret; else total += ret;
        x->N = N;

        new_fmpz_array(&x->WireShares, N);
        for (i = 0; i < N; i++) {
            ret = recv_fmpz(x->WireShares[i]);
            if (ret <= 0) return ret; else total += ret;
        }

        ret = recv_fmpz(x->f0_s);
        if (ret <= 0) return ret; else total += ret;
        ret = recv_fmpz(x->g0_s);
        if (ret <= 0) return ret; else total += ret;
        ret = recv_fmpz(x->h0_s);
        if (ret <= 0) return ret; else total += ret;

        new_fmpz_array(&x->h_points, N);
        for (i = 0; i < N; i++) {
            ret = recv_fmpz(x->h_points[i]);
            if (ret <= 0) return ret; else total += ret;
        }
        return total;
    }

    int BeaverTriple(BeaverTriple *x) {
        int total = 0, ret;
        ret = recv_fmpz(x->A);
        if (ret <= 0) return ret; else total += ret;
        ret = recv_fmpz(x->B);
        if (ret <= 0) return ret; else total += ret;
        ret = recv_fmpz(x->C);
        if (ret <= 0) return ret; else total += ret;
        return total;
    }

    int BeaverTripleShare(BeaverTripleShare *x) {
        int total = 0, ret;
        ret = recv_fmpz(x->shareA);
        if (ret <= 0) return ret; else total += ret;
        ret = recv_fmpz(x->shareB);
        if (ret <= 0) return ret; else total += ret;
        ret = recv_fmpz(x->shareC);
        if (ret <= 0) return ret; else total += ret;
        return total;
    }
};

#endif