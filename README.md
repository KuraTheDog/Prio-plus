# Background

This is the software prototype for our paper on Private Computation on Streaming Data.

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

## Build

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

* Server arguments are `server_num client_listen_port server_communicate_port`
* Client arguments are `./bin/client num_submissions server0_listen_port server1_listen_port TOPK num_bits delta K`

The client samples from a ZipF distribution of `max_bits`-bit values.
With probability at least `1 - delta`, the approximately top K most frequent values are returned.

### Usage example

0. Run `cd build` and `make`
1. For the first server, run `./bin/server 0 8800 8888` to start it
2. For the second server, run `./bin/server 1 8801 8888` to start it (in another window / instance / computer)
3. Use e.g. `./bin/client 10 8800 8801 TOPK 8 0.05 5` to run a meta-client that sends out client messages (in another window / instance / computer)
4. Eventually the server will output, and wait for a fresh set of inputs (`waiting for connection...`)
5. Repeat step 3 as desired with different parameters

# Code flow outline

0. Servers connect to each other
1. Servers do initial communication and precomputation
2. Client produces shares
3. Client connects to servers, sends corresponding shares
4. Servers receive all shares
5. Servers align shares using tags (for if in different order)
6. The servers compute on the shares, through communication rounds, and updates their accumulators with values from valid shares.
7. Each server locally computes the final answer using the aggregates, by performing additional rounds for evaluation.
