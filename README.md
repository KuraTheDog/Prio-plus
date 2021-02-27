# Installation

1. Install [emp-ot](https://github.com/emp-toolkit/emp-ot)
2. Install Boost
3. `cd AggrProject && mkdir build`
4. `cd build && cmake ..`
5. `make`

# Usage example

1. make
2. Run `./bin/server 0 8800 8888 8`
3. In another window, run `./bin/server 1 8801 8888 8`
4. In another window, run `./bin/client 10 8800 8801 VAROP 8`
5. Repeat step 4 as desired with different num_inputs and operands. 

* Server arguments are `server_num Client_listen_port server_0_port max_bits`
* Client arguments are `num_inputs server0_port server1_port operand max_bits linreg_degree`

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
6. If relevant, servers use share conversion to get wire shares for SNIPS
7. If relevent, servers run SNIPS validations
8. If relevant, Server 1 does OT share conversion and aggregates.
9. Server 1 sends it's aggregate values to server 0. 
10. If relevant, Server 0 does OT share conversion and aggregates, accounting for validity
11. Server 0 computes the final answer using the combined aggregates
