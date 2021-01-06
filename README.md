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

* Client arguments are `server_num Client_listen_port server_0_port max_bits`
* Client arguments are `num_inputs server0_port server1_port operand max_bits`

Ports and max bits need to be consistent across runs.
Max bits is only used for int based summations.
Server 0 port tells Server 0 which port to open, and server 1 which port of server 0 is open.

# Outline

bit/int sum protocol

client ---shares0--> server0, server1
server1 ----Init_bit/int_sum msg ----> server0
server1 ---pks of all submissions----> server0
server1 starts a ot_receiver with choice bit as its bitshare
server0 filters the pks sent by server1, creates b0, b1 (r, 1-r arrays) ot arrays and uses ot library to send them over
each server gets sum of the values sent in the OT and exchange final results and sum it up.

block - 128 bit number
In the ot protocols, I use first 64 bits to encode validity(done at server0) and second 64 bits contain the actual digit.
