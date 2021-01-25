/*
For sending fmtp_t objects over sockets.

Also has htonl/ntohl wrappers for sending various int-like types.

See test/test_net_share.cpp for a full example.

Returns total bytes sent on success, or first fail return of an internal step (typically 0 or negative)

TODO: move this to net_share.cpp? Currently has linker issues when attempted.

TODO: maybe hide/delete/comment out ones we don't use.
*/

#ifndef NET_SHARE_H
#define NET_SHARE_H

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

/* Core functions */

// Send trivial. 
int recv_in(const int sockfd, void* const buf, const size_t len);

int send_bool(const int sockfd, const bool x);
int recv_bool(const int sockfd, bool& x);

int send_int(const int sockfd, const int x);
int recv_int(const int sockfd, int& x);

int send_size(const int sockfd, const size_t x);
int recv_size(const int sockfd, size_t& x);

int send_uint32(const int sockfd, const uint32_t x);
int recv_uint32(const int sockfd, uint32_t& x);

int send_uint64(const int sockfd, const uint64_t x);
int recv_uint64(const int sockfd, uint64_t& x);

int send_ulong(const int sockfd, const ulong x);
int recv_ulong(const int sockfd, ulong& x);

int send_fmpz(const int sockfd, const fmpz_t x);
int recv_fmpz(const int sockfd, fmpz_t x);

/* Share functions */

int send_Cor(const int sockfd, const Cor* const x);
int recv_Cor(const int sockfd, Cor* const x);

int send_CorShare(const int sockfd, const CorShare* const x);
int recv_CorShare(const int sockfd, CorShare* const x);

// ClientPacket = client_packet*
int send_ClientPacket(const int sockfd, const ClientPacket x);
int recv_ClientPacket(const int sockfd, ClientPacket &x);

int send_BeaverTriple(const int sockfd, const BeaverTriple* const x);
int recv_BeaverTriple(const int sockfd, BeaverTriple* const x);

int send_BeaverTripleShare(const int sockfd, const BeaverTripleShare* const x);
int recv_BeaverTripleShare(const int sockfd, BeaverTripleShare* const x);

int send_BooleanBeaverTriple(const int sockfd, const BooleanBeaverTriple x);
int recv_BooleanBeaverTriple(const int sockfd, BooleanBeaverTriple& x);

// Bits stuff
int send_DaBit(const int sockfd, const DaBit* const x);
int recv_DaBit(const int sockfd, DaBit* const x);

// Assumes n is already known.
int send_EdaBit(const int sockfd, const EdaBit* const x, const size_t n);
int recv_EdaBit(const int sockfd, EdaBit* const x, const size_t n);

#endif