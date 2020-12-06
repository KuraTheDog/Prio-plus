/* 
For sending fmtp_t objects over sockets.

Send: takes in socket to write to, and object to send. 
Rec: takes in socket to read from, and blank object to overwrite.

In case of success, returns a positive number. In case of failure, returns a non-positive number.
Returns first failing instantly on fail, otherwise on success returns final ret.

TODO: move this to net_share.cpp. Currently has linker issues when attempted.

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

int send_fmpz(int sockfd, fmpz_t x) {
    FILE* file = fdopen(sockfd, "w");
    int ret = fmpz_out_raw(file, x);
    fclose(file);
    return ret;
}

int rec_fmpz(int sockfd, fmpz_t x) {
    FILE* file = fdopen(sockfd, "r");
    int ret = fmpz_inp_raw(x, file);
    fclose(file);
    return ret;
}

int send_Cor(int sockfd, Cor *x) {
    FILE* file = fdopen(sockfd, "w");
    int ret = fmpz_out_raw(file, x->D);
    if (ret < 0) return ret;
    ret = fmpz_out_raw(file, x->E);
    fclose(file);
    return ret;
}

int rec_Cor(int sockfd, Cor *x) {
    FILE* file = fdopen(sockfd, "r");
    int ret = fmpz_inp_raw(x->D, file);
    if (ret < 0) return ret;
    ret = fmpz_inp_raw(x->E, file);
    fclose(file);
    return ret;
}

int send_CorShare(int sockfd, Cor *x) {
    FILE* file = fdopen(sockfd, "w");
    int ret = fmpz_out_raw(file, x->D);
    if (ret < 0) return ret;
    ret = fmpz_out_raw(file, x->E);
    fclose(file);
    return ret;
}

int rec_CorShare(int sockfd, Cor *x) {
    FILE* file = fdopen(sockfd, "r");
    int ret = fmpz_inp_raw(x->D, file);
    if (ret < 0) return ret;
    ret = fmpz_inp_raw(x->E, file);
    fclose(file);
    return ret;
}

// Note: This is different from a ClientPacket
int send_client_packet(int sockfd, client_packet *x) {
    FILE* file = fdopen(sockfd, "w");
    int N = x->N;
    int32_t conv_N = htonl(N);
    int ret = send(sockfd, &conv_N, sizeof(conv_N), 0);
    if (ret < 0) return ret;

    int i;
    for (i = 0; i < N; i++) {
        ret = fmpz_out_raw(file, x->WireShares[i]);
        if (ret < 0) return ret;
    }

    ret = fmpz_out_raw(file, x->f0_s);
    if (ret < 0) return ret;
    ret = fmpz_out_raw(file, x->g0_s);
    if (ret < 0) return ret;
    ret = fmpz_out_raw(file, x->h0_s);
    if (ret < 0) return ret;

    for (i = 0; i < N; i++) {
        ret = fmpz_out_raw(file, x->h_points[i]);
        if (ret < 0) return ret;
    }
    return ret;
}

// Note: This is different from a ClientPacket
int rec_client_packet(int sockfd, client_packet *x) {
    FILE* file = fdopen(sockfd, "r");
    int i, N, recv_N = 0;
    int ret = recv(sockfd, &recv_N, sizeof(recv_N), 0);
    if (ret < 0) return ret;
    N = ntohl(recv_N);
    x->N = N;

    new_fmpz_array(&x->WireShares, N);
    for (i = 0; i < N; i++) {
        ret = fmpz_inp_raw(x->WireShares[i], file);
        if (ret < 0) return ret;
    }

    ret = fmpz_inp_raw(x->f0_s, file);
    if (ret < 0) return ret;
    ret = fmpz_inp_raw(x->g0_s, file);
    if (ret < 0) return ret;
    ret = fmpz_inp_raw(x->h0_s, file);
    if (ret < 0) return ret;

    new_fmpz_array(&x->h_points, N);
    for (i = 0; i < N; i++) {
        ret = fmpz_inp_raw(x->h_points[i], file);
        if (ret < 0) return ret;
    }
    return ret;
}

int send_BeaverTriple(int sockfd, BeaverTriple *x) {
    FILE* file = fdopen(sockfd, "w");
    int ret = fmpz_out_raw(file, x->A);
    if (ret < 0) return ret;
    ret = fmpz_out_raw(file, x->B);
    if (ret < 0) return ret;
    ret = fmpz_out_raw(file, x->C);
    fclose(file);
    // return 1;
    return ret;
}

int rec_BeaverTriple(int sockfd, BeaverTriple *x) {
    FILE* file = fdopen(sockfd, "r");
    int ret = fmpz_inp_raw(x->A, file);
    if (ret < 0) return ret;
    ret = fmpz_inp_raw(x->B, file);
    if (ret < 0) return ret;
    ret = fmpz_inp_raw(x->C, file);
    fclose(file);
    return ret;
    // return 1;
}

int send_BeaverTripleShare(int sockfd, BeaverTripleShare *x) {
    FILE* file = fdopen(sockfd, "w");
    int ret = fmpz_out_raw(file, x->shareA);
    if (ret < 0) return ret;
    ret = fmpz_out_raw(file, x->shareB);
    if (ret < 0) return ret;
    ret = fmpz_out_raw(file, x->shareC);
    fclose(file);
    return ret;
}

int rec_BeaverTripleShare(int sockfd, BeaverTripleShare *x) {
    FILE* file = fdopen(sockfd, "r");
    int ret = fmpz_inp_raw(x->shareA, file);
    if (ret < 0) return ret;
    ret = fmpz_inp_raw(x->shareB, file);
    if (ret < 0) return ret;
    ret = fmpz_inp_raw(x->shareC, file);
    fclose(file);
    return ret;
}

#endif