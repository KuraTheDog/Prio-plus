# Background

This is the software prototype based on Prio+ (below), which adds functionality for identifying the most frequent heavy hitters.

## Original background

This is the software prototype that accompanies the research paper

Prio+ https://eprint.iacr.org/2021/576.pdf
By Surya Addanki, Kevin Garbe, Eli Jaffe, Rafail Ostrovsky, and Antigoni Polychroniadou

This is an improvement of the original Prio code by Henry Corrigan-Gibs and Dan Boneh, which can be found [here](https://github.com/henrycg/prio), detailed in ["Prio: Private, Robust, and Scalable Computation of Aggregate Statistics"](https://crypto.stanford.edu/prio/paper.pdf)
Some of this code, such as the fast polynomial operations, is directly based on their code.


## Dependencies

1. [Flint 2.7.0+](https://flintlib.org)
2. [emp toolkit](https://github.com/emp-toolkit)
3. [PALISADE](https://gitlab.com/palisade/palisade-release)

# Getting Started

## Install dependencies

Follow the links above to install the corresponding packages

### EMP Instructions

Uses EMP OT and SH2PC.
Install via `python install.py --deps --tool --ot --sh2pc` to install all relevant packages.
(May need to run `install_name_tool -id '/usr/local/lib/libemp-tool.dylib' /usr/local/lib/libemp-tool.dylib` to get linkages to work, at least on MAC).

### PALISADE Instructions

For full install instructions, see [here](https://gitlab.com/palisade/palisade-release/-/wikis/Build-instructions).
Currently, Palisade is built with default args (`cmake ..`).

PALISADE encryption is used to generate arithmetic Beaver triples.

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

* Server arguments are `server_num client_listen_port server0_port`
* Client arguments are `./bin/client num_submissions server0_port server1_port OPERATION num_bits (extra params)`
  * `num_bits` required for ops using integers.
  * `extra params` depends on protocols

* Ports and max bits need to be consistent across runs and both servers and the client.
* `max_bits` is used for int based summations, and must match the server value in this case.
  * For MAX, client `max_bits` instead determines the max value (e.g. 7 -> 128), and does not have to match the servers.
* For server communication, `server0_port` tells Server 0 which port to open, and server 1 which port of server 0 is open.

### Usage example

0. Run `cd build` and `make`
1. For the first server, run `./bin/server 0 8800 8888` to start it
2. For the second server, run `./bin/server 1 8801 8888` to start it (in another window / instance / computer)
3. Use e.g. `./bin/client 10 8800 8801 VAR 8` to run a meta-client that sends out client messages (in another window / instance / computer)
4. Repeat step 3 as desired with different parameters

## Supported protocols

* BITSUM: 1 bit integers sum
* INTSUM: `max_bits` -bit integer sum
* AND / OR: Boolean and/or
* MAX / MIN: Max/min, with values between 0 and `2^max_bits`
* VAR / STDDEV: Variance / Standard Deviation of `max_bits`-bit integers
* LINREG: `degree` degree linear regression on `max_bits`-bit integers
* FREQ: Standard frequency counting
* HEAVY: CountMin example. Queries all values, returns those that are heavy. As example, uses uniform stream of `max_bits`-bit integers, except with one value occuring a fraction `f` of the stream.
* MULTIHEAVY: ZipF distribution of `max_bits` values. Tries to return with accuracy at least `1 - delta` be within `eps` range for the top `K`.
* TOPK: ZipF distribution of `max_bits` values. Like MultiHeavy, but with more repetition and halving, to work better on more even distributions.

# Code flow outline example

0. Servers connect to each other
1. Servers do initial communication and precomputation
2. Client produces shares
3. Client connects to servers, sends corresponding shares
4. Servers receive all shares
5. Servers align shares using tags (for if in different order)
6. If relevant, servers use daBit or OT share conversion to convert the values to arithmetic shares, for aggregation and SNIPS
7. If relevent, servers run validations on shares, including SNIPS
8. Servers combine their aggregates, taking into account validitiy of shares. This may take multiple rounds. This may also happen in parallel with validation, depending on protocol.
9. Each server locally computes the final answer using the aggregates, either by just revealing them to each other, or performing additional rounds for evaluation.
