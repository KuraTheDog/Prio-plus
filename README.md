# Installation

1. Install [emp-ot](https://github.com/emp-toolkit/emp-ot)
2. Install Boost
3. cd AggrProject && mkdir build
4. cd build && cmake ..
5. make

bit/int sum protocol

client ---shares0--> server0, server1
server1 ----Init_bit/int_sum msg ----> server0
server1 ---pks of all submissions----> server0
server1 starts a ot_receiver with choice bit as its bitshare
server0 filters the pks sent by server1, creates b0, b1 (r, 1-r arrays) ot arrays and uses ot library to send them over
each server gets sum of the values sent in the OT and exchange final results and sum it up.

block - 128 bit number
In the ot protocols, I use first 64 bits to encode validity(done at server0) and second 64 bits contain the actual digit.