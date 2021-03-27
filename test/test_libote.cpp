#include <iostream>

#include "cryptoTools/Common/Defines.h"   // block
#include "cryptoTools/Network/IOService.h"  // IOService
// #include "libOTe/Base/SimplestOT.h"
#include "libOTe/TwoChooseOne/IknpOtExtReceiver.h"
#include "libOTe/TwoChooseOne/IknpOtExtSender.h"

int main(int argc, char* argv[]) {

  // Setup networking. See cryptoTools\frontend_cryptoTools\Tutorials\Network.cpp
  osuCrypto::IOService ios;
  osuCrypto::Channel senderChl = osuCrypto::Session(ios, "localhost", 1212, osuCrypto::SessionMode::Server).addChannel();
  osuCrypto::Channel recverChl = osuCrypto::Session(ios, "localhost:1212", osuCrypto::SessionMode::Client).addChannel();

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