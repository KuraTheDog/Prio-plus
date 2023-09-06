#include "net_share.h"

#include <sys/socket.h>

#include "constants.h"
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

int recv_in(const int sockfd, void* const buff, const size_t len) {
    unsigned int bytes_read = 0, tmp;
    char* const buffptr = (char*) buff;
    while (bytes_read < len) {
        tmp = recv(sockfd, buffptr + bytes_read, len - bytes_read, 0);
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

int send_bool_batch(const int sockfd, const bool* const x, const size_t n) {
    const size_t len = (n+7) / 8;  // Number of bytes to hold n, aka ceil(n/8)
    char* const buff = new char[len];

    memset(buff, 0, sizeof(char) * len);

    for (unsigned int i = 0; i < n; i++)
        if (x[i])
            buff[i / 8] ^= (1ULL << (i % 8));

    int ret = send(sockfd, buff, len, 0);

    delete[] buff;

    return ret;
}

int recv_bool_batch(const int sockfd, bool* const x, const size_t n) {
    const size_t len = (n+7) / 8;
    char* const buff = new char[len];

    int ret = recv_in(sockfd, buff, len);

    for (unsigned int i = 0; i < n; i++)
        x[i] = (buff[i/8] & (1ULL << (i % 8)));

    delete[] buff;

    return ret;
}

int swap_bool_batch(const int sockfd, const bool* const x, bool* const y, const size_t n) {
    int sent_bytes = 0;
    int recv_bytes = 0;

    std::thread t_send([sockfd, x, n, &sent_bytes]() {
        sent_bytes += send_bool_batch(sockfd, x, n);
    });
    std::thread t_recv([sockfd, &y, n, &recv_bytes]() {
        recv_bytes += recv_bool_batch(sockfd, y, n);
    });
    t_send.join();
    t_recv.join();

    if (sent_bytes != recv_bytes) {
        std::cout << "WARNING: sent " << sent_bytes << " bytes, ";
        std::cout << "but received " << recv_bytes << " bytes" << std::endl;
    }

    return sent_bytes;
}

int reveal_bool_batch(const int sockfd, bool* const x, const size_t n) {
    int sent_bytes = 0;

    bool* buff = new bool[n];

    sent_bytes += swap_bool_batch(sockfd, x, buff, n);

    for (unsigned int i = 0; i < n; i++)
        x[i] ^= buff[i];

    delete[] buff;

    return sent_bytes;
}

int send_int(const int sockfd, const int x) {
    int x_conv = htonl(x);
    const char* const data = (const char*) &x_conv;
    return send(sockfd, data, sizeof(int), 0);
}

int recv_int(const int sockfd, int& x) {
    int ret = recv_in(sockfd, &x, sizeof(int));
    x = ntohl(x);
    return ret;
}

int send_size(const int sockfd, const size_t x) {
    size_t x_conv = htonl(x);
    const char* const data = (const char*) &x_conv;
    return send(sockfd, data, sizeof(size_t), 0);
}

int recv_size(const int sockfd, size_t& x) {
    int ret = recv_in(sockfd, &x, sizeof(size_t));
    x = ntohl(x);
    return ret;
}

int send_double(const int sockfd, const double x) {
    const char* const data = (const char*) &x;
    return send(sockfd, data, sizeof(double), 0);
}

int recv_double(const int sockfd, double& x) {
    int ret = recv_in(sockfd, &x, sizeof(double));
    return ret;
}

int send_uint32(const int sockfd, const uint32_t x) {
    uint32_t x_conv = htonl(x);
    const char* const data = (const char*) &x_conv;
    return send(sockfd, data, sizeof(uint32_t), 0);
}

int recv_uint32(const int sockfd, uint32_t& x) {
    int ret = recv_in(sockfd, &x, sizeof(uint32_t));
    x = ntohl(x);
    return ret;
}

int send_uint64(const int sockfd, const uint64_t x) {
    uint64_t x_conv = htonll(x);
    const char* const data = (const char*) &x_conv;
    return send(sockfd, data, sizeof(uint64_t), 0);
}

int recv_uint64(const int sockfd, uint64_t& x) {
    int ret = recv_in(sockfd, &x, sizeof(uint64_t));
    x = ntohll(x);
    return ret;
}

int send_uint64_batch(const int sockfd, const uint64_t* const x, const size_t n) {
    return send(sockfd, x, n * sizeof(uint64_t), 0);
}

int recv_uint64_batch(const int sockfd, uint64_t* const x, const size_t n) {
    return recv_in(sockfd, x, n * sizeof(uint64_t));
}

int send_ulong(const int sockfd, const ulong x) {
    ulong x_conv = htonll(x);
    const char* const data = (const char*) &x_conv;
    return send(sockfd, data, sizeof(ulong), 0);
}

int recv_ulong(const int sockfd, ulong& x) {
    int ret = recv_in(sockfd, &x, sizeof(ulong));
    x = ntohll(x);
    return ret;
}

int send_ulong_batch(const int sockfd, const ulong* const x, const size_t n) {
    return send(sockfd, x, n * sizeof(ulong), 0);
}

int recv_ulong_batch(const int sockfd, ulong* const x, const size_t n) {
    return recv_in(sockfd, x, n * sizeof(ulong));
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
    char* const buff = new char[len];
    ret = recv_in(sockfd, buff, len);
    if (ret <= 0) return ret; else total += ret;
    x.assign(&buff[0], len);
    delete[] buff;
    return total;
}

int send_fmpz(const int sockfd, const fmpz_t x) {
    int total = 0, ret;
    size_t len;
    if (FIXED_FMPZ_SIZE) {
        len = mod_size;
    } else {
        len = fmpz_size(x);
        ret = send_size(sockfd, len);
        if (ret <= 0) return ret; else total += ret;
    }

    ulong buff[len];
    fmpz_get_ui_array(buff, len, x);
    ret = send_ulong_batch(sockfd, buff, len);
    if (ret <= 0) return ret; else total += ret;

    return total;
}

int recv_fmpz(const int sockfd, fmpz_t x) {
    int total = 0, ret;
    size_t len;

    if (FIXED_FMPZ_SIZE) {
        len = mod_size;
    } else {
        ret = recv_size(sockfd, len);
        if (ret <= 0) return ret; else total += ret;
    }

    if (len == 0) {
        fmpz_set_ui(x, 0);
        return total;
    }
    ulong buff[len];
    ret = recv_ulong_batch(sockfd, buff, len);
    if (ret <= 0) return ret; else total += ret;

    fmpz_set_ui_array(x, buff, len);

    return total;
}

int reveal_fmpz(const int sockfd, fmpz_t x) {
    int sent_bytes = 0;
    int recv_bytes = 0;
    fmpz_t tmp; fmpz_init(tmp);

    std::thread t_send([sockfd, x, &sent_bytes]() {
        sent_bytes += send_fmpz(sockfd, x);
    });
    std::thread t_recv([sockfd, &tmp, &recv_bytes]() {
        recv_bytes += recv_fmpz(sockfd, tmp);
    });
    t_send.join();
    t_recv.join();
    fmpz_mod_add(x, x, tmp, mod_ctx);
    fmpz_clear(tmp);

    return sent_bytes;
}

int send_fmpz_batch(const int sockfd, const fmpz_t* const x, const size_t n) {
    int total = 0, ret;
    if (!FIXED_FMPZ_SIZE) {
        // Lazy version
        for (unsigned int i = 0; i < n; i++) {
            ret += send_fmpz(sockfd, x[i]);
            if (ret <= 0) return ret; else total += ret;
        }
        return total;
    }

    if (n > MAX_FMPZ_BATCH) {
        total += send_fmpz_batch(sockfd, x, MAX_FMPZ_BATCH);
        auto new_x = &x[MAX_FMPZ_BATCH];
        total += send_fmpz_batch(sockfd, new_x, n - MAX_FMPZ_BATCH);
        return total;
    }

    size_t len = mod_size;
    ulong buff[len * n];
    for (unsigned int i = 0; i < n; i++)
        fmpz_get_ui_array(&buff[i * len], len, x[i]);

    ret = send_ulong_batch(sockfd, buff, len * n);
    if (ret <= 0) return ret; else total += ret;
    return total;
}

int recv_fmpz_batch(const int sockfd, fmpz_t* const x, const size_t n) {
    int total = 0, ret;

    if (!FIXED_FMPZ_SIZE) {
        // Lazy version
        for (unsigned int i = 0; i < n; i++) {
            ret += recv_fmpz(sockfd, x[i]);
            if (ret <= 0) return ret; else total += ret;
        }
        return total;
    }

    if (n > MAX_FMPZ_BATCH) {
        total += recv_fmpz_batch(sockfd, x, MAX_FMPZ_BATCH);
        auto new_x = &x[MAX_FMPZ_BATCH];
        total += recv_fmpz_batch(sockfd, new_x, n - MAX_FMPZ_BATCH);
        return total;
    }

    size_t len = mod_size;
    ulong buff[len * n];
    ret = recv_ulong_batch(sockfd, buff, len * n);
    if (ret <= 0) return ret; else total += ret;

    for (unsigned int i = 0; i < n; i++)
        fmpz_set_ui_array(x[i], &buff[i * len], len);
    return total;
}

int swap_fmpz_batch(const int sockfd, const fmpz_t* const x, fmpz_t* const y, const size_t n) {
    int sent_bytes = 0;
    int recv_bytes = 0;

    std::thread t_send([sockfd, x, n, &sent_bytes]() {
        sent_bytes += send_fmpz_batch(sockfd, x, n);
    });
    std::thread t_recv([sockfd, &y, n, &recv_bytes]() {
        recv_bytes += recv_fmpz_batch(sockfd, y, n);
    });
    t_send.join();
    t_recv.join();

    if (sent_bytes != recv_bytes) {
        std::cout << "WARNING: sent " << sent_bytes << " bytes, ";
        std::cout << "but received " << recv_bytes << " bytes" << std::endl;
    }

    return sent_bytes;
}

int reveal_fmpz_batch(const int sockfd, fmpz_t* const x, const size_t n) {
    fmpz_t* buff; new_fmpz_array(&buff, n);

    int sent_bytes = swap_fmpz_batch(sockfd, x, buff, n);

    for (unsigned int i = 0; i < n; i++) {
        fmpz_mod_add(x[i], x[i], buff[i], mod_ctx);
    }

    clear_fmpz_array(buff, n);

    return sent_bytes;
}

int swap_bool_fmpz_batch(const int sockfd,
        const bool* const x, bool* const y, const size_t n,
        const fmpz_t* const xp, fmpz_t* yp, const size_t np) {
    int sent_bytes = 0;
    int recv_bytes = 0;

    std::thread t_send([sockfd, x, n, xp, np, &sent_bytes]() {
        sent_bytes += send_bool_batch(sockfd, x, n);
        sent_bytes += send_fmpz_batch(sockfd, xp, np);
    });
    std::thread t_recv([sockfd, y, n, yp, np, &recv_bytes]() {
        recv_bytes += recv_bool_batch(sockfd, y, n);
        recv_bytes += recv_fmpz_batch(sockfd, yp, np);
    });
    t_send.join();
    t_recv.join();

    if (sent_bytes != recv_bytes) {
        std::cout << "WARNING: sent " << sent_bytes << " bytes, ";
        std::cout << "but received " << recv_bytes << " bytes" << std::endl;
    }

    return sent_bytes;
}

int send_seed(const int sockfd, const flint_rand_t x) {
    const char* const data = (const char*) &x[0];
    return send(sockfd, data, sizeof(x[0]), 0);
}

int recv_seed(const int sockfd, flint_rand_t x) {
    return recv_in(sockfd, &x[0], sizeof(x[0]));
}

/* Share functions */

int send_Cor(const int sockfd, const Cor* const x) {
    int total = 0, ret;
    fmpz_t* buff; new_fmpz_array(&buff, 2);
    fmpz_set(buff[0], x->D);
    fmpz_set(buff[1], x->E);
    ret = send_fmpz_batch(sockfd, buff, 2);
    if (ret <= 0) return ret; else total += ret;
    clear_fmpz_array(buff, 2);
    return total;
}

int recv_Cor(const int sockfd, Cor* const x) {
    int total = 0, ret;
    fmpz_t* buff; new_fmpz_array(&buff, 2);
    ret = recv_fmpz_batch(sockfd, buff, 2);
    if (ret <= 0) return ret; else total += ret;
    fmpz_set(x->D, buff[0]);
    fmpz_set(x->E, buff[1]);
    clear_fmpz_array(buff, 2);
    return total;
}

int send_Cor_batch(const int sockfd, const Cor* const * const x, const size_t n) {
    int total = 0, ret;
    fmpz_t* buff; new_fmpz_array(&buff, 2 * n);

    for (unsigned int i = 0; i < n; i++) {
        fmpz_set(buff[i], x[i]->D);
        fmpz_set(buff[i + n], x[i]->E);
    }
    ret = send_fmpz_batch(sockfd, buff, 2 * n);
    if (ret <= 0) return ret; else total += ret;

    clear_fmpz_array(buff, 2 * n);
    return total;
}

int recv_Cor_batch(const int sockfd, Cor* const * const x, const size_t n) {
    fmpz_t* buff; new_fmpz_array(&buff, 2 * n);

    for (unsigned int i = 0; i < n; i++) {
        fmpz_set(x[i]->D, buff[i]);
        fmpz_set(x[i]->E, buff[i + n]);
    }

    int ret = recv_fmpz_batch(sockfd, buff, 2 * n);

    clear_fmpz_array(buff, 2 * n);
    return ret;
}

int reveal_Cor_batch(const int sockfd, Cor* const * const x, const size_t n) {
    fmpz_t* buff; new_fmpz_array(&buff, 2 * n);
    for (unsigned int i = 0; i < n; i++) {
        fmpz_set(buff[i], x[i]->D);
        fmpz_set(buff[i + n], x[i]->E);
    }
    int sent_bytes = reveal_fmpz_batch(sockfd, buff, 2 * n);

    for (unsigned int i = 0; i < n; i++) {
        fmpz_set(x[i]->D, buff[i]);
        fmpz_set(x[i]->E, buff[i + n]);
    }

    clear_fmpz_array(buff, 2 * n);
    return sent_bytes;
}

int send_ClientPacket(const int sockfd, const ClientPacket* const x,
                      const size_t NMul) {
    int total = 0, ret;
    const size_t N = NextPowerOfTwo(NMul);
    const size_t len = NMul + N + 3 + 3;
    fmpz_t* buff; new_fmpz_array(&buff, len);

    for (unsigned int i = 0; i < NMul; i++)
        fmpz_set(buff[i], x->MulShares[i]);
    for (unsigned int i = 0; i < N; i++)
        fmpz_set(buff[NMul + i], x->h_points[i]);
    fmpz_set(buff[NMul+N], x->f0_s);
    fmpz_set(buff[NMul+N+1], x->g0_s);
    fmpz_set(buff[NMul+N+2], x->h0_s);
    fmpz_set(buff[NMul+N+3], x->triple->A);
    fmpz_set(buff[NMul+N+4], x->triple->B);
    fmpz_set(buff[NMul+N+5], x->triple->C);

    ret = send_fmpz_batch(sockfd, buff, len);
    if (ret <= 0) return ret; else total += ret;

    clear_fmpz_array(buff, len);
    return total;
}

int recv_ClientPacket(const int sockfd, ClientPacket* const x,
                      const size_t NMul) {
    int total = 0, ret;
    const size_t N = NextPowerOfTwo(NMul);
    const size_t len = NMul + N + 3 + 3;
    fmpz_t* buff; new_fmpz_array(&buff, len);

    ret = recv_fmpz_batch(sockfd, buff, len);
    if (ret <= 0) return ret; else total += ret;

    for (unsigned int i = 0; i < NMul; i++)
        fmpz_set(x->MulShares[i], buff[i]);
    for (unsigned int i = 0; i < N; i++)
        fmpz_set(x->h_points[i], buff[NMul + i]);
    fmpz_set(x->f0_s, buff[NMul+N]);
    fmpz_set(x->g0_s, buff[NMul+N+1]);
    fmpz_set(x->h0_s, buff[NMul+N+2]);
    fmpz_set(x->triple->A, buff[NMul+N+3]);
    fmpz_set(x->triple->B, buff[NMul+N+4]);
    fmpz_set(x->triple->C, buff[NMul+N+5]);

    clear_fmpz_array(buff, len);
    return total;
}

int send_BeaverTriple(const int sockfd, const BeaverTriple* const x) {
    int total = 0, ret;
    fmpz_t* buff; new_fmpz_array(&buff, 3);
    fmpz_set(buff[0], x->A);
    fmpz_set(buff[1], x->B);
    fmpz_set(buff[2], x->C);
    ret = send_fmpz_batch(sockfd, buff, 3);
    if (ret <= 0) return ret; else total += ret;
    clear_fmpz_array(buff, 3);
    return total;
}

int recv_BeaverTriple(const int sockfd, BeaverTriple* const x) {
    int total = 0, ret;
    fmpz_t* buff; new_fmpz_array(&buff, 3);
    ret = recv_fmpz_batch(sockfd, buff, 3);
    if (ret <= 0) return ret; else total += ret;
    fmpz_set(x->A, buff[0]);
    fmpz_set(x->B, buff[1]);
    fmpz_set(x->C, buff[2]);
    clear_fmpz_array(buff, 3);
    return total;
}

int send_BoolTriple(const int sockfd, const BooleanBeaverTriple* const x) {
    bool buff[3] = {x->a, x->b, x->c};
    return send_bool_batch(sockfd, buff, 3);
}

int recv_BoolTriple(const int sockfd, BooleanBeaverTriple* const x) {
    bool buff[3];
    int ret = recv_bool_batch(sockfd, buff, 3);
    x->a = buff[0];
    x->b = buff[1];
    x->c = buff[2];
    return ret;
}

int send_BoolTriple_batch(const int sockfd, const BooleanBeaverTriple* const * const x, const size_t n) {
    bool buff[3 * n];
    for (unsigned int i = 0; i < n; i++) {
        buff[3*i] = x[i]->a;
        buff[3*i+1] = x[i]->b;
        buff[3*i+2] = x[i]->c;
    }
    return send_bool_batch(sockfd, buff, 3 * n);
}

int recv_BoolTriple_batch(const int sockfd, BooleanBeaverTriple* const * const x, const size_t n) {
    bool buff[3 * n];
    int ret = recv_bool_batch(sockfd, buff, 3 * n);
    for (unsigned int i = 0; i < n; i++) {
        x[i]->a = buff[3*i];
        x[i]->b = buff[3*i+1];
        x[i]->c = buff[3*i+2];
    }
    return ret;
}

int send_AltTriple(const int sockfd, const AltTriple* const x) {
    int total = 0, ret;
    fmpz_t* buff; new_fmpz_array(&buff, 2);
    fmpz_set(buff[0], x->AB);
    fmpz_set(buff[1], x->C);
    ret = send_fmpz_batch(sockfd, buff, 2);
    if (ret <= 0) return ret; else total += ret;
    clear_fmpz_array(buff, 2);
    return total;
}

int recv_AltTriple(const int sockfd, AltTriple* const x) {
    int total = 0, ret;
    fmpz_t* buff; new_fmpz_array(&buff, 2);
    ret = recv_fmpz_batch(sockfd, buff, 2);
    if (ret <= 0) return ret; else total += ret;
    fmpz_set(x->AB, buff[0]);
    fmpz_set(x->C, buff[1]);
    clear_fmpz_array(buff, 2);
    return total;
}

int send_AltTriple_batch(const int sockfd, const AltTriple* const * const x, const size_t n) {
    int total = 0, ret;
    fmpz_t* buff; new_fmpz_array(&buff, 2 * n);

    for (unsigned int i = 0; i < n; i++) {
        fmpz_set(buff[2*i], x[i]->AB);
        fmpz_set(buff[2*i+1], x[i]->C);
    }

    ret = send_fmpz_batch(sockfd, buff, 2 * n);
    if (ret <= 0) return ret; else total += ret;

    clear_fmpz_array(buff, 2 * n);
    return total;
}

int recv_AltTriple_batch(const int sockfd, AltTriple* const * const x, const size_t n) {
    int total = 0, ret;
    fmpz_t* buff; new_fmpz_array(&buff, 2 * n);

    ret = recv_fmpz_batch(sockfd, buff, 2 * n);
    if (ret <= 0) return ret; else total += ret;

    for (unsigned int i = 0; i < n; i++) {
        fmpz_set(x[i]->AB, buff[2*i]);
        fmpz_set(x[i]->C, buff[2*i+1]);
    }

    clear_fmpz_array(buff, 2 * n);
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

int send_DaBit_batch(const int sockfd, const DaBit* const * const x, const size_t n) {
    int total = 0, ret;

    if (n > MAX_DABIT_BATCH) {
        total += send_DaBit_batch(sockfd, x, MAX_DABIT_BATCH);
        auto new_x = &x[MAX_DABIT_BATCH];
        total += send_DaBit_batch(sockfd, new_x, n - MAX_DABIT_BATCH);
        return total;
    }

    fmpz_t* bp; new_fmpz_array(&bp, n);
    bool* const b2 = new bool[n];

    for (unsigned int i = 0; i < n; i++) {
        fmpz_set(bp[i], x[i]->bp);
        b2[i] = x[i]->b2;
    }

    ret = send_fmpz_batch(sockfd, bp, n);
    if (ret <= 0) return ret; else total += ret;

    ret = send_bool_batch(sockfd, b2, n);
    if (ret <= 0) return ret; else total += ret;

    clear_fmpz_array(bp, n);
    delete[] b2;

    return total;
}

int recv_DaBit_batch(const int sockfd, DaBit* const * const x, const size_t n) {
    int total = 0, ret;

    if (n > MAX_DABIT_BATCH) {
        total += recv_DaBit_batch(sockfd, x, MAX_DABIT_BATCH);
        auto new_x = &x[MAX_DABIT_BATCH];
        total += recv_DaBit_batch(sockfd, new_x, n - MAX_DABIT_BATCH);
        return total;
    }

    fmpz_t* bp; new_fmpz_array(&bp, n);
    bool* const b2 = new bool[n];

    ret = recv_fmpz_batch(sockfd, bp, n);
    if (ret <= 0) return ret; else total += ret;
    ret = recv_bool_batch(sockfd, b2, n);
    if (ret <= 0) return ret; else total += ret;

    for (unsigned int i = 0; i < n; i++) {
        fmpz_set(x[i]->bp, bp[i]);
        x[i]->b2 = b2[i];
    }

    clear_fmpz_array(bp, n);
    delete[] b2;

    return total;
}

/*
int send_EdaBit(const int sockfd, const EdaBit* const x, const size_t nbits) {
    int total = 0, ret;
    ret = send_fmpz(sockfd, x->r);
    if (ret <= 0) return ret; else total += ret;
    ret = send_bool_batch(sockfd, &x->b[0], nbits);
    if (ret <= 0) return ret; else total += ret;
    return total;
}

int recv_EdaBit(const int sockfd, EdaBit* const x, const size_t nbits) {
    int total = 0, ret;
    ret = recv_fmpz(sockfd, x->r);
    if (ret <= 0) return ret; else total += ret;
    ret = recv_bool_batch(sockfd, &x->b[0], nbits);
    if (ret <= 0) return ret; else total += ret;
    return total;
}

int send_EdaBit_batch(const int sockfd, const EdaBit* const * const x, const size_t nbits, const size_t n) {
    int total = 0, ret;
    fmpz_t* r; new_fmpz_array(&r, n);
    bool* const b = new bool[n * nbits];

    for (unsigned int i = 0; i < n; i++) {
        fmpz_set(r[i], x[i]->r);
        memcpy(&b[i * nbits], x[i]->b, nbits * sizeof(bool));
    }

    ret = send_fmpz_batch(sockfd, r, n);
    if (ret <= 0) return ret; else total += ret;
    ret = send_bool_batch(sockfd, b, n * nbits);
    if (ret <= 0) return ret; else total += ret;

    clear_fmpz_array(r, n);
    delete[] b;

    return total;
}

int recv_EdaBit_batch(const int sockfd, EdaBit* const * const x, const size_t nbits, const size_t n) {
    int total = 0, ret;
    fmpz_t* r; new_fmpz_array(&r, n);
    bool* const b = new bool[n * nbits];

    ret = recv_fmpz_batch(sockfd, r, n);
    if (ret <= 0) return ret; else total += ret;
    ret = recv_bool_batch(sockfd, b, n * nbits);
    if (ret <= 0) return ret; else total += ret;

    for (unsigned int i = 0; i < n; i++) {
        fmpz_set(x[i]->r, r[i]);
        memcpy(x[i]->b, &b[i * nbits], nbits * sizeof(bool));
    }

    clear_fmpz_array(r, n);
    delete[] b;

    return total;
}
*/
