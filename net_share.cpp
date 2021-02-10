#include "net_share.h"

#include <sys/socket.h>

#include "fmpz_utils.h"

// Ensure this is defined, as it's architecture dependent
#if !defined(htonll) && !defined(ntohll)
#if __BIG_ENDIAN__
# define htonll(x) (x)
# define ntohll(x) (x)
#else
# define htonll(x) (((uint64_t)htonl((x) & 0xFFFFFFFF) << 32) | htonl((x) >> 32))
# define ntohll(x) (((uint64_t)ntohl((x) & 0xFFFFFFFF) << 32) | ntohl((x) >> 32))
#endif
#endif

/* Core functions */

int recv_in(const int sockfd, void* const buf, const size_t len) {
    unsigned int bytes_read = 0, tmp;
    char* bufptr = (char*) buf;
    while (bytes_read < len) {
        tmp = recv(sockfd, bufptr + bytes_read, len - bytes_read, 0);
        if (tmp <= 0) return tmp; else bytes_read += tmp;
    }
    return bytes_read;
}

int send_bool(const int sockfd, const bool x) {
    return send(sockfd, &x, sizeof(bool), 0);
}

int recv_bool(const int sockfd, bool& x) {
    return recv_in(sockfd, &x, sizeof(bool));
}

int send_int(const int sockfd, const int x) {
    int x_conv = htonl(x);
    const char* data = (const char*) &x_conv;
    return send(sockfd, data, sizeof(int), 0);
}

int recv_int(const int sockfd, int& x) {
    int ret = recv_in(sockfd, &x, sizeof(int));
    x = ntohl(x);
    return ret;
}

int send_size(const int sockfd, const size_t x) {
    size_t x_conv = htonl(x);
    const char* data = (const char*) &x_conv;
    return send(sockfd, data, sizeof(size_t), 0);
}

int recv_size(const int sockfd, size_t& x) {
    int ret = recv_in(sockfd, &x, sizeof(size_t));
    x = ntohl(x);
    return ret;
}

int send_uint32(const int sockfd, const uint32_t x) {
    uint32_t x_conv = htonl(x);
    const char* data = (const char*) &x_conv;
    return send(sockfd, data, sizeof(uint32_t), 0);
}

int recv_uint32(const int sockfd, uint32_t& x) {
    int ret = recv_in(sockfd, &x, sizeof(uint32_t));
    x = ntohl(x);
    return ret;
}

int send_uint64(const int sockfd, const uint64_t x) {
    uint64_t x_conv = htonll(x);
    const char* data = (const char*) &x_conv;
    return send(sockfd, data, sizeof(uint64_t), 0);
}

int recv_uint64(const int sockfd, uint64_t& x) {
    int ret = recv_in(sockfd, &x, sizeof(uint64_t));
    x = ntohll(x);
    return ret;
}

int send_ulong(const int sockfd, const ulong x) {
    ulong x_conv = htonll(x);
    const char* data = (const char*) &x_conv;
    return send(sockfd, data, sizeof(ulong), 0);
}

int recv_ulong(const int sockfd, ulong& x) {
    int ret = recv_in(sockfd, &x, sizeof(ulong));
    x = ntohll(x);
    return ret;
}

int send_string(const int sockfd, const std::string x) {
    int total = 0, ret;
    size_t len = x.size();
    ret = send_size(sockfd, len);
    if (ret <= 0) return ret; else total += ret;
    ret = send(sockfd, x.c_str(), len, 0);
    if (ret <= 0) return ret; else total += ret;
    return total;
}

int recv_string(const int sockfd, std::string& x) {
    int total = 0, ret;
    size_t len;
    ret = recv_size(sockfd, len);
    if (ret <= 0) return ret; else total += ret;
    char* buf = new char[len];
    ret = recv_in(sockfd, buf, len);
    if (ret <= 0) return ret; else total += ret;
    x.assign(&buf[0], len);
    delete[] buf;
    return total;
}

int send_fmpz(const int sockfd, const fmpz_t x) {
    int total = 0, ret;
    size_t len = fmpz_size(x);
    ulong arr[len];
    fmpz_get_ui_array(arr, len, x);
    ret = send_size(sockfd, len);
    if (ret <= 0) return ret; else total += ret;
    for (unsigned int i = 0; i < len; i++) {
        ret = send_ulong(sockfd, arr[i]);
        if (ret <= 0) return ret; else total += ret;
    }
    return total;
}

int recv_fmpz(const int sockfd, fmpz_t x) {
    int total = 0, ret;
    size_t len;
    ulong tmp;
    ret = recv_size(sockfd, len);
    if (ret <= 0) return ret; else total += ret;
    if (len == 0) {
        fmpz_set_ui(x, 0);
        return total;
    }
    ulong buf[len];
    for (unsigned int i = 0; i < len; i++) {
        ret = recv_ulong(sockfd, tmp);
        if (ret <= 0) return ret; else total += ret;
        buf[i] = tmp;
    }
    fmpz_set_ui_array(x, buf, len);
    return total;
}

/* Share functions */

int send_Cor(const int sockfd, const Cor* const x) {
    int total = 0, ret;
    ret = send_fmpz(sockfd, x->D);
    if (ret <= 0) return ret; else total += ret;
    ret = send_fmpz(sockfd, x->E);
    if (ret <= 0) return ret; else total += ret;
    return total;
}

int recv_Cor(const int sockfd, Cor* const x) {
    int total = 0, ret;
    ret = recv_fmpz(sockfd, x->D);
    if (ret <= 0) return ret; else total += ret;
    ret = recv_fmpz(sockfd, x->E);
    if (ret <= 0) return ret; else total += ret;
    return total;
}

int send_CorShare(const int sockfd, const CorShare* const x) {
    int total = 0, ret;
    ret = send_fmpz(sockfd, x->shareD);
    if (ret <= 0) return ret; else total += ret;
    ret = send_fmpz(sockfd, x->shareE);
    if (ret <= 0) return ret; else total += ret;
    return total;
}

int recv_CorShare(const int sockfd, CorShare* const x) {
    int total = 0, ret;
    ret = recv_fmpz(sockfd, x->shareD);
    if (ret <= 0) return ret; else total += ret;
    ret = recv_fmpz(sockfd, x->shareE);
    if (ret <= 0) return ret; else total += ret;
    return total;
}

int send_ClientPacket(const int sockfd, const ClientPacket* const x) {
    int total = 0, ret;
    size_t N = x->N, NWires = x->NWires;
    ret = send_size(sockfd, N);
    if (ret <= 0) return ret; else total += ret;
    ret = send_size(sockfd, NWires);
    if (ret <= 0) return ret; else total += ret;

    for (unsigned int i = 0; i < NWires; i++) {
        ret = send_fmpz(sockfd, x->WireShares[i]);
        if (ret <= 0) return ret; else total += ret;
    }

    ret = send_fmpz(sockfd, x->f0_s);
    if (ret <= 0) return ret; else total += ret;
    ret = send_fmpz(sockfd, x->g0_s);
    if (ret <= 0) return ret; else total += ret;
    ret = send_fmpz(sockfd, x->h0_s);
    if (ret <= 0) return ret; else total += ret;

    for (unsigned int i = 0; i < N; i++) {
        ret = send_fmpz(sockfd, x->h_points[i]);
        if (ret <= 0) return ret; else total += ret;
    }

    ret = send_BeaverTripleShare(sockfd, x->triple_share);
    if (ret <= 0) return ret; else total += ret;

    return total;
}

int recv_ClientPacket(const int sockfd, ClientPacket* const x) {
    int total = 0, ret;
    bool is_valid = true;
    size_t N, NWires;
    ret = recv_size(sockfd, N);
    if (ret <= 0) return ret; else total += ret;
    ret = recv_size(sockfd, NWires);
    if (ret <= 0) return ret; else total += ret;

    if (x->N != N or x->NWires != NWires)
        is_valid = false;

    for (unsigned int i = 0; i < NWires; i++) {
        ret = recv_fmpz(sockfd, x->WireShares[i]);
        if (ret <= 0) return ret; else total += ret;
    }

    ret = recv_fmpz(sockfd, x->f0_s);
    if (ret <= 0) return ret; else total += ret;
    ret = recv_fmpz(sockfd, x->g0_s);
    if (ret <= 0) return ret; else total += ret;
    ret = recv_fmpz(sockfd, x->h0_s);
    if (ret <= 0) return ret; else total += ret;

    for (unsigned int i = 0; i < N; i++) {
        ret = recv_fmpz(sockfd, x->h_points[i]);
        if (ret <= 0) return ret; else total += ret;
    }

    ret = recv_BeaverTripleShare(sockfd, x->triple_share);
    if (ret <= 0) return ret; else total += ret;

    return (is_valid ? total : -1);
}

int send_BeaverTriple(const int sockfd, const BeaverTriple* const x) {
    int total = 0, ret;
    ret = send_fmpz(sockfd, x->A);
    if (ret <= 0) return ret; else total += ret;
    ret = send_fmpz(sockfd, x->B);
    if (ret <= 0) return ret; else total += ret;
    ret = send_fmpz(sockfd, x->C);
    if (ret <= 0) return ret; else total += ret;
    return total;
}

int recv_BeaverTriple(const int sockfd, BeaverTriple* const x) {
    int total = 0, ret;
    ret = recv_fmpz(sockfd, x->A);
    if (ret <= 0) return ret; else total += ret;
    ret = recv_fmpz(sockfd, x->B);
    if (ret <= 0) return ret; else total += ret;
    ret = recv_fmpz(sockfd, x->C);
    if (ret <= 0) return ret; else total += ret;
    return total;
}

int send_BeaverTripleShare(const int sockfd, const BeaverTripleShare* const x) {
    int total = 0, ret;
    ret = send_fmpz(sockfd, x->shareA);
    if (ret <= 0) return ret; else total += ret;
    ret = send_fmpz(sockfd, x->shareB);
    if (ret <= 0) return ret; else total += ret;
    ret = send_fmpz(sockfd, x->shareC);
    if (ret <= 0) return ret; else total += ret;
    return total;
}

int recv_BeaverTripleShare(const int sockfd, BeaverTripleShare* const x) {
    int total = 0, ret;
    ret = recv_fmpz(sockfd, x->shareA);
    if (ret <= 0) return ret; else total += ret;
    ret = recv_fmpz(sockfd, x->shareB);
    if (ret <= 0) return ret; else total += ret;
    ret = recv_fmpz(sockfd, x->shareC);
    if (ret <= 0) return ret; else total += ret;
    return total;
}

int send_BooleanBeaverTriple(const int sockfd, const BooleanBeaverTriple* const x) {
    int total = 0, ret;
    ret = send_bool(sockfd, x->a);
    if (ret <= 0) return ret; else total += ret;
    ret = send_bool(sockfd, x->b);
    if (ret <= 0) return ret; else total += ret;
    ret = send_bool(sockfd, x->c);
    if (ret <= 0) return ret; else total += ret;
    return total;
}

int recv_BooleanBeaverTriple(const int sockfd, BooleanBeaverTriple* const x) {
    int total = 0, ret;
    ret = recv_bool(sockfd, x->a);
    if (ret <= 0) return ret; else total += ret;
    ret = recv_bool(sockfd, x->b);
    if (ret <= 0) return ret; else total += ret;
    ret = recv_bool(sockfd, x->c);
    if (ret <= 0) return ret; else total += ret;
    return total;
}

int send_DaBit(const int sockfd, const DaBit* const x) {
    int total = 0, ret;
    ret = send_fmpz(sockfd, x->bp);
    if (ret <= 0) return ret; else total += ret;
    ret = send_bool(sockfd, x->b2);
    if (ret <= 0) return ret; else total += ret;
    return total;
}

int recv_DaBit(const int sockfd, DaBit* const x) {
    int total = 0, ret;
    ret = recv_fmpz(sockfd, x->bp);
    if (ret <= 0) return ret; else total += ret;
    ret = recv_bool(sockfd, x->b2);
    if (ret <= 0) return ret; else total += ret;
    return total;
}

int send_EdaBit(const int sockfd, const EdaBit* const x, const size_t n) {
    int total = 0, ret;
    ret = send_fmpz(sockfd, x->r);
    if (ret <= 0) return ret; else total += ret;
    ret = send(sockfd, &x->b[0], n * sizeof(bool), 0);
    if (ret <= 0) return ret; else total += ret;
    return total;
}

int recv_EdaBit(const int sockfd, EdaBit* const x, const size_t n) {
    int total = 0, ret;
    ret = recv_fmpz(sockfd, x->r);
    if (ret <= 0) return ret; else total += ret;
    ret = recv_in(sockfd, &x->b[0], n * sizeof(bool));
    if (ret <= 0) return ret; else total += ret;
    return total;
}