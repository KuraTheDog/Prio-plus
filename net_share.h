/*
For sending fmtp_t objects over sockets.

Also has htonl/ntohl wrappers for sending various int-like types.

See test/test_net_share.cpp for a full example.

Returns total bytes sent on success, or first fail return of an internal step (typically 0 or negative)
*/

#ifndef NET_SHARE_H
#define NET_SHARE_H

#include <string>

#include "share.h"
#include "types.h"

extern "C" {
  #include "flint/flint.h"
  #include "flint/fmpz.h"
};

/*
Other ideas for fmpz sending:
fmpz_in_raw, fmpz_out_raw. Has trouble with using socket for other things.
    can have send_int map to it, but e.g. server has trouble with other general use.
fmpz_sng + fmpz_get_ui_array: best for really large numbers?
fmpz_get_str: best for small numbers.

ulong array: Always uses ulongs. 32 or 64 bits. "perfect" space efficiency, for really large numbers.
string: best it can do is base 62, so 62/256 ~ 25% space efficiency. So needs ~4x bits compared to numbers.
*/

// Since fmpz is bound by Int modulus, we can leave off sending size?
// Reduces number of bits to send fmpz
#define FIXED_FMPZ_SIZE true

// We batch things together into single send/recieves, to reduce overhead (mainly recv wrapper I think). However, it segfaults if given too large batches. So this makes sure batches are capped
#define MAX_FMPZ_BATCH  1000000
#define MAX_DABIT_BATCH 320000

/* Core functions */

// Send trivial
int recv_in(const int sockfd, void* const buf, const size_t len);

int send_bool(const int sockfd, const bool x);
int recv_bool(const int sockfd, bool& x);

// Batch send bools compacted into char array
// Versus 1 bool at a time takes up a whole byte per bool
int send_bool_batch(const int sockfd, const bool* const x, const size_t n);
int recv_bool_batch(const int sockfd, bool* const x, const size_t n);

// Unused
int send_int(const int sockfd, const int x);
int recv_int(const int sockfd, int& x);

int send_size(const int sockfd, const size_t x);
int recv_size(const int sockfd, size_t& x);

int send_double(const int sockfd, const double x);
int recv_double(const int sockfd, double& x);

// Unused
int send_uint32(const int sockfd, const uint32_t x);
int recv_uint32(const int sockfd, uint32_t& x);

int send_uint64(const int sockfd, const uint64_t x);
int recv_uint64(const int sockfd, uint64_t& x);

int send_uint64_batch(const int sockfd, const uint64_t* const x, const size_t n);
int recv_uint64_batch(const int sockfd, uint64_t* const x, const size_t n);

int send_ulong(const int sockfd, const ulong x);
int recv_ulong(const int sockfd, ulong& x);

int send_ulong_batch(const int sockfd, const ulong* const x, const size_t n);
int recv_ulong_batch(const int sockfd, ulong* const x, const size_t n);

int send_string(const int sockfd, const std::string x);
int recv_string(const int sockfd, std::string& x);

// If FIXED_FMPZ_SIZE then uses less bytes
int send_fmpz(const int sockfd, const fmpz_t x);
int recv_fmpz(const int sockfd, fmpz_t x);

int send_fmpz_batch(const int sockfd, const fmpz_t* const x, const size_t n);
int recv_fmpz_batch(const int sockfd, fmpz_t* const x, const size_t n);

int send_seed(const int sockfd, const flint_rand_t x);
int recv_seed(const int sockfd, flint_rand_t x);

/* Share functions */

// Unused
int send_Cor(const int sockfd, const Cor* const x);
int recv_Cor(const int sockfd, Cor* const x);

int send_CorShare(const int sockfd, const CorShare* const x);
int recv_CorShare(const int sockfd, CorShare* const x);

int send_CorShare_batch(const int sockfd, const CorShare* const * const x, const size_t n);
int recv_CorShare_batch(const int sockfd, CorShare* const * const x, const size_t n);

int send_ClientPacket(const int sockfd, const ClientPacket* const x,
                      const size_t NMul);
int recv_ClientPacket(const int sockfd, ClientPacket* const x,
                      const size_t NMul);

// Only used in making lazy triples. Also batch?
int send_BeaverTriple(const int sockfd, const BeaverTriple* const x);
int recv_BeaverTriple(const int sockfd, BeaverTriple* const x);

// Unused
int send_BooleanBeaverTriple(const int sockfd, const BooleanBeaverTriple* const x);
int recv_BooleanBeaverTriple(const int sockfd, BooleanBeaverTriple* const x);

// Bits stuff
int send_DaBit(const int sockfd, const DaBit* const x);
int recv_DaBit(const int sockfd, DaBit* const x);

int send_DaBit_batch(const int sockfd, const DaBit* const * const x, const size_t n);
int recv_DaBit_batch(const int sockfd, DaBit* const * const x, const size_t n);

// Assumes n is already known.
int send_EdaBit(const int sockfd, const EdaBit* const x, const size_t nbits);
int recv_EdaBit(const int sockfd, EdaBit* const x, const size_t nbits);

int send_EdaBit_batch(const int sockfd, const EdaBit* const * const x, const size_t nbits, const size_t n);
int recv_EdaBit_batch(const int sockfd, EdaBit* const * const x, const size_t nbits, const size_t n);

#endif