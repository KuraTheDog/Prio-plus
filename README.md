# Installation

1. Install [emp-ot](https://github.com/emp-toolkit/emp-ot)
2. Install Boost
3. Install [PALISADE](https://gitlab.com/palisade/palisade-release)
3. `cd AggrProject && mkdir build`
4. `cd build && cmake ..`
5. `make`

## PALISADE Instructions

For full install instructions, see [here](https://gitlab.com/palisade/palisade-release/-/wikis/Build-instructions).  
Currently, Palisade is built with default args (`cmake ..`).  
Also make sure to e.g. install with `make install` after `make`.

PALISADE encryption is used to generate beaver triples.  
It only supports up to 60 bit moduli, so otherwise we use slower methods for now.

# Usage example

1. make
2. Run `./bin/server 0 8800 8888 8`
3. In another window, run `./bin/server 1 8801 8888 8`
4. In another window, run `./bin/client 10 8800 8801 VAROP 8`
5. Repeat step 4 as desired with different num_inputs and operands. 

* Server arguments are `server_num Client_listen_port server_0_port max_bits`
* Client arguments are `num_inputs server0_port server1_port operand max_bits do_batch include_invalid`
  * `do_batch` (default `true`)
    * If `true`, the client make all shares, then send them all
    * If `false`, the client makes and sends shares one at a time
  * `include_invalid` (default false)
    * For testing/debugging. If true, does a modified run where some clients intentionally send bad data.

Ports and max bits need to be consistent across runs and both servers and the client.
Max bits is only used for int based summations.
Server 0 port tells Server 0 which port to open, and server 1 which port of server 0 is open.

# Outline

0. Servers connect to each other
1. Servers do initial communication and precomputation
2. Client produces shares
3. Client connects to servers, sends corresponding shares
4. Servers recieve all shares
5. For each share, servers check if public keys line up
6. If relevent, servers also do SNIPS validation for each share
7. If doing SNIPS, servers also ensure shares match the SNIPS using share conversion.
8. If relevant, Server 1 uses an OT Reciever with all of it's shares.
9. Server 1 sends it's aggregate value to server 0. 
10. If relevant, Server 0 uses an OT sender with it's shares, and the shared information about share validity.
11. Server 0 computes it's own aggregate value, and combines it with the one from server 1 to produce a final answer.

block - 128 bit number
In the ot protocols, the first 64 bits to encode validity(done at server0) and second 64 bits contain the actual digit.
