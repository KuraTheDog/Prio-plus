#include <iostream>

#include "cryptoTools/Common/Defines.h"   // block
#include "cryptoTools/Network/IOService.h"  // IOService
#include "libOTe/TwoChooseOne/IknpOtExtReceiver.h"
#include "libOTe/TwoChooseOne/IknpOtExtSender.h"

#define ADDR "127.0.0.1"
#define PORT 60050

void test() {
  // Setup networking. See cryptoTools\frontend_cryptoTools\Tutorials\Network.cpp
  osuCrypto::IOService ios;
  osuCrypto::IOService ios2;
  osuCrypto::Channel senderChl = osuCrypto::Session(ios, ADDR, PORT, osuCrypto::SessionMode::Server).addChannel();
  osuCrypto::Channel recverChl = osuCrypto::Session(ios2, ADDR, PORT, osuCrypto::SessionMode::Client).addChannel();

  // The number of OTs.
  int n = 1;

  // The code to be run by the OT receiver.
  auto recverThread = std::thread([&]() {
    osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());
    osuCrypto::IknpOtExtReceiver recver;

    // Choose which messages should be received.
    osuCrypto::BitVector choices(n);
    choices[0] = 1;
    //...

    // Receive the messages
    std::vector<osuCrypto::block> messages(n);
    recver.receiveChosen(choices, messages, prng, recverChl);

    // messages[i] = sendMessages[i][choices[i]];
    std::cout << "Selecting 0 on " << choices[0] << std::endl;
    std::cout << "Got 0 = " << messages[0] << std::endl;
  });

  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());
  osuCrypto::IknpOtExtSender sender;

  // Choose which messages should be sent.
  std::vector<std::array<osuCrypto::block, 2>> sendMessages(n);
  sendMessages[0] = { osuCrypto::toBlock(54), osuCrypto::toBlock(33) };
  //...

  std::cout << "Send 0, 0 = " << sendMessages[0][0] << std::endl;
  std::cout << "Send 0, 1 = " << sendMessages[0][1] << std::endl;

  // Send the messages.
  sender.sendChosen(sendMessages, prng, senderChl);
  recverThread.join();
}

void test_send() {
  osuCrypto::IOService ios;
  osuCrypto::Channel channel = osuCrypto::Session(ios, ADDR, PORT, osuCrypto::SessionMode::Server).addChannel();
  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());
  int n = 1;

  osuCrypto::IknpOtExtSender sender;
  std::vector<std::array<osuCrypto::block, 2>> sendMessages(n);
  sendMessages[0] = { osuCrypto::toBlock(54), osuCrypto::toBlock(33) };
  //...

  std::cout << "Send 0, 0 = " << sendMessages[0][0] << std::endl;
  std::cout << "Send 0, 1 = " << sendMessages[0][1] << std::endl;

  sender.sendChosen(sendMessages, prng, channel);
}

void test_recv() {
  osuCrypto::IOService ios;
  osuCrypto::Channel channel = osuCrypto::Session(ios, ADDR, PORT, osuCrypto::SessionMode::Client).addChannel();
  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());
  int n = 1;

  osuCrypto::IknpOtExtReceiver recver;

  osuCrypto::BitVector choices(n);
  choices[0] = 1;
  //...

  // Receive the messages
  std::vector<osuCrypto::block> messages(n);
  recver.receiveChosen(choices, messages, prng, channel);

  // messages[i] = sendMessages[i][choices[i]];
  std::cout << "Selecting 0 on " << choices[0] << std::endl;
  std::cout << "Got 0 = " << messages[0] << std::endl;
}

int main(int argc, char* argv[]) {
  int server_num = -1;
  if(argc >= 2){
    server_num = atoi(argv[1]);
  }

  if (server_num == -1) {
    test();
  } else if (server_num == 0) {
    test_send();
  } else if (server_num == 1) {
    test_recv();
  }
}