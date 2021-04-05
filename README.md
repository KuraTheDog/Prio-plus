# Background

This is the software prototype that accompanies the research paper

> Prio+ (TODO: link)
By Eli Jaffe, Kevin Garbe, Surya Addanki, Rafail Ostrovsky, and Antigoni Polychroniadou

This is an improvement of the original Prio code by Henry Corrigan-Gibs and Dan Boneh, which can be found [here](https://github.com/henrycg/prio), detailed in ["Prio: Private, Robust, and Scalable Computation of Aggregate Statistics"](https://crypto.stanford.edu/prio/paper.pdf)

Some of this code, such as the fast polynomial operations, is directly based on their code. 


## Dependencies

1. [Flint 2.7.0+](https://flintlib.org)
2. [emp-ot](https://github.com/emp-toolkit/emp-ot)
3. [PALISADE](https://gitlab.com/palisade/palisade-release)
4. [libOTe](https://github.com/osu-crypto/libOTe/tree/master/libOTe/TwoChooseOne)

# Getting Started

## Install dependencies

Follow the links above to install the corresponding packages

### PALISADE Instructions

For full install instructions, see [here](https://gitlab.com/palisade/palisade-release/-/wikis/Build-instructions).  
Currently, Palisade is built with default args (`cmake ..`).

PALISADE encryption is used to generate beaver triples.

### libOTe Instructions

Build libOTe with `cmake . -DENABLE_RELIC=ON -DENABLE_NP=ON -DENABLE_IKNP=ON -DENABLE_SILENTOT=ON`.
Becaue installing is not complete yet, it's currently hard coded into this `CMakeLists.txt`.
Change `libOTe_Dirs` as needed based on where libOTe ended up.

Used for alternative IKNP OTs, and silent OTs

Swapped between them based on OT_TYPE define in `ot.h`

## Build Prio+

1. `mkdir build`
2. `cd build`
3. `cmake ..`
4. `make`

## Run

There are two relevant binaries:
* server: Runs a server instance
* client: Mimics a cluster of individual clients

By default, they are configured to connect to localhost.
For testing across multiple machines, configure the constants on top of `client.cpp` and `server.cpp` to the IP addresses of the two servers.

The code runs two servers, 0 and 1, each of which needs to be started separately. 
Server 0 needs to be started before server 1.

* Server arguments are `server_num client_listen_port server0_port max_bits`
* Client arguments are `num_inputs server0_port server1_port operand max_bits linreg_degree`

* Ports and max bits need to be consistent across runs and both servers and the client.
* `max_bits` is used for int based summations, and must match the server value in this case.  
  * For MAXOP, client `max_bits` instead determines the max value (e.g. 7 -> 128), and does not have to match the servers.  
* For server communication, `server0_port` tells Server 0 which port to open, and server 1 which port of server 0 is open.

### Usage example

0. Run `cd build` and `make`
1. Run `./bin/server 0 8800 8888 8` to start the first server
2. In another window, run `./bin/server 1 8801 8888 8` to start the second server
3. In another window, run `./bin/client 10 8800 8801 VAROP 8` to run a meta-client that sends out client messages
3. Repeat step 4 as desired with different parameters

## Supported protocols

* BITSUM: 1 bit integers sum
* INTSUM: `max_bits` -bit integer sum
* ANDOP / OROP: Boolean and/or
* MAXOP / MINOP: Max/min, with values between 0 and `2^max_bits`
* VAROP / STDDEVOP: Variance / Standard Deviation of `max_bits`-bit integers
* LINREGOP: `linreg_degree` degree linear regression on `max_bits`-bit integers

# Code flow outline

0. Servers connect to each other
1. Servers do initial communication and precomputation
2. Client produces shares
3. Client connects to servers, sends corresponding shares
4. Servers recieve all shares
5. For each share, servers check if public keys line up
6. If relevant, servers use EdaBit share conversion to get wire shares for SNIPS, or OT share conversion if SNIPS are not needed
7. If relevent, servers run SNIPS validations
8. Servers combine their aggregates, taking into account validitiy of shares
9. Server 0 computes the final answer using the combined aggregates
