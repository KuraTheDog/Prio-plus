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

#include "fmpz_utils.h"
#include "share.h"

extern "C" {
    #include "flint/flint.h"
    #include "flint/fmpz.h"
};

/* 
Other ideas:
fmpz_in_raw, fmpz_out_raw. Has trouble with using socket for other things.
    can have send_int map to it, but e.g. server has trouble with other general use.
fmpz_sng + fmpz_get_ui_array: best for really large numbers?
fmpz_get_str: best for small numbers. 

ulong array: Always uses ulongs. 32 or 64 bits. "perfect" space efficiency, for really large numbers.
string: best is base 62, so 62/256 ~ 25% space efficiency. So needs ~4x bits compared to numbers. 
*/

class ShareSender {
    int sockfd;

    int send_int(const int x){
        int x_conv = htonl(x);
        const char* data = (const char*) &x_conv;
        int ret = send(sockfd, data, sizeof(int), 0);
        return ret;
    }

    int send_ulong(const ulong x) {
        ulong x_conv = htonll(x);
        const char* data = (const char*) &x_conv;
        int ret = send(sockfd, data, sizeof(ulong), 0);
        return ret;
    }

    int send_fmpz(fmpz_t x){
        int total = 0, ret;
        size_t len = fmpz_size(x);
        ulong arr[len];
        fmpz_get_ui_array(arr, len, x);
        ret = send_int(len);
        if (ret <= 0) return ret; else total += ret;
        for (int i = 0; i < len; i++) {
            ret = send_ulong(arr[i]);
            if (ret <= 0) return ret; else total += ret;
        }
        return total;
    }

public: 

    ShareSender(int sockfd) {
        this->sockfd = sockfd;
    }

    ~ShareSender() {
    }

    int fmpz(fmpz_t x) {
        return send_fmpz(x);
    }

    int integer(int x) {
        return send_int(x);
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

    int client_packet(client_packet *x) {
        int N = x->N, NWires = x->NWires, total = 0, ret;
        ret = send_int(N);
        if (ret <= 0) return ret; else total += ret;
        ret = send_int(NWires);
        if (ret <= 0) return ret; else total += ret;

        int i;
        for (i = 0; i < NWires; i++) {
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

    int recv_int(int& x) {
        int bytes_read = 0, tmp;
        char buf[sizeof(int)];
        while (bytes_read < sizeof(int)) {
            tmp = recv(sockfd, &buf + bytes_read, sizeof(int) - bytes_read, 0);
            if (tmp <= 0) return tmp; else bytes_read += tmp;
        }
        x = ntohl(*((int*)buf));
        return bytes_read;
    }

    int recv_ulong(ulong& x) {
        char buf[sizeof(ulong)];
        int bytes_read = 0, tmp;
        while (bytes_read < sizeof(ulong)) {
            tmp = recv(sockfd, &buf + bytes_read, sizeof(ulong) - bytes_read, 0);
            if (tmp <= 0) return tmp; else bytes_read += tmp;
        }
        x = ntohll(*((ulong*)buf));
        return bytes_read;
    }

    int recv_fmpz(fmpz_t x){
        int total = 0, ret, len;
        ulong tmp;
        ret = recv_int(len);
        if (ret <= 0) return ret; else total += ret;
        if (len == 0) {
            fmpz_set_ui(x, 0);
            return total;
        }
        ulong buf[len];
        for (int i = 0; i < len; i++) {
            ret = recv_ulong(tmp);
            if (ret <= 0) return ret; else total += ret;
            buf[i] = tmp;
        }
        fmpz_set_ui_array(x, buf, len);
        return total;
    }

public: 

    ShareReciever(int sockfd) {
        this->sockfd = sockfd;
    }

    ~ShareReciever() {
    }

    int fmpz(fmpz_t x) {
        return recv_fmpz(x);
    }

    int integer(int& x) {
        return recv_int(x);
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

    int client_packet(client_packet* &x) {
        int N, NWires, total = 0, ret;
        ret = recv_int(N);
        if (ret <= 0) return ret; else total += ret;
        ret = recv_int(NWires);
        if (ret <= 0) return ret; else total += ret;
        init_client_packet(x, N, NWires);

        int i;
        for (i = 0; i < NWires; i++) {
            ret = recv_fmpz(x->WireShares[i]);
            if (ret <= 0) return ret; else total += ret;
        }

        ret = recv_fmpz(x->f0_s);
        if (ret <= 0) return ret; else total += ret;
        ret = recv_fmpz(x->g0_s);
        if (ret <= 0) return ret; else total += ret;
        ret = recv_fmpz(x->h0_s);
        if (ret <= 0) return ret; else total += ret;

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