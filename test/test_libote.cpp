#include <iostream>

#include "cryptoTools/Common/Defines.h"   // block
#include "cryptoTools/Network/IOService.h"  // IOService
#include "libOTe/TwoChooseOne/IknpOtExtReceiver.h"
#include "libOTe/TwoChooseOne/IknpOtExtSender.h"
#include "libOTe/TwoChooseOne/SilentOtExtReceiver.h"
#include "libOTe/TwoChooseOne/SilentOtExtSender.h"

#define ADDR "127.0.0.1"
#define PORT 60050

/*
silentSend and silentReceive
- configures at start w/ messges.size, 2, 128, chls.size (default anyways)
- if mGen.hasBaseOts = false,
  - if mIknpSender.hasBase = false, genBase
  - then gen silent base
*/

void test(bool silent) {
  // Setup networking. See cryptoTools\frontend_cryptoTools\Tutorials\Network.cpp
  osuCrypto::IOService ios;
  osuCrypto::IOService ios2;
  osuCrypto::Channel senderChl = osuCrypto::Session(ios, ADDR, PORT, osuCrypto::SessionMode::Server).addChannel();
  osuCrypto::Channel recverChl = osuCrypto::Session(ios2, ADDR, PORT, osuCrypto::SessionMode::Client).addChannel();

  std::cout << (silent ? "Silent test" : "IKNP test") << std::endl;

  /* The number of OTs.
  Silent: Requires >= 888... whyyyy
  ./frontend/frontend_libOTe -u 63 -v -n 127
  requires at least 128, which makes more sense
  */
  int n = 888;

  // The code to be run by the OT receiver.
  auto recverThread = std::thread([&]() {
    osuCrypto::PRNG prng(osuCrypto::toBlock(0, 1));

    // Choose which messages should be received.
    osuCrypto::BitVector choices(n);
    for (int i = 0; i < n; i++)
      choices[i] = i % 2;
    //...

    // Receive the messages
    std::vector<osuCrypto::block> messages(n);

    if(silent) {
      osuCrypto::SilentOtExtReceiver recver;

      // recver.mDebug = true;
      // recver.configure(n);
      // auto count = recver.silentBaseOtCount();
      // std::cout << "receiver base ot count: " << count << std::endl;
      // // fake base ots
      // osuCrypto::BitVector basechoices = recver.sampleBaseChoiceBits(prng);
      // std::vector<osuCrypto::block> msg(basechoices.size());
      // for (unsigned int i = 0; i < msg.size(); ++i)
      //   msg[i] = osuCrypto::toBlock(i, basechoices[i]);
      // recver.setSlientBaseOts(msg);

      recver.silentReceive(choices, messages, prng, recverChl);
    } else {
      osuCrypto::IknpOtExtReceiver recver;
      recver.receiveChosen(choices, messages, prng, recverChl);
    }

    // messages[i] = sendMessages[i][choices[i]];
    std::cout << "Got 0 [" << choices[0] << "] \t= " << messages[0] << std::endl;
    std::cout << "Got " << n-1 << " [" << choices[n-1] << "] \t= " << messages[n-1] << std::endl;
  });

  // osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());
  osuCrypto::PRNG prng(osuCrypto::toBlock(0, 0));

  // Choose which messages should be sent.
  std::vector<std::array<osuCrypto::block, 2>> sendMessages(n);

  // Send the messages.
  if (silent) {
    osuCrypto::SilentOtExtSender sender;

    // sender.mDebug = true;
    // sender.configure(n);
    // auto count = sender.silentBaseOtCount();
    // std::cout << "sender base ot count: " << count << std::endl;
    // // fake base ots
    // std::vector<std::array<osuCrypto::block, 2>> msg(count);
    // for (unsigned int i = 0; i < msg.size(); ++i)
    // {
    //   msg[i][0] = osuCrypto::toBlock(i, 0);
    //   msg[i][1] = osuCrypto::toBlock(i, 1);
    // }
    // sender.setSlientBaseOts(msg);

    sender.silentSend(sendMessages, prng, senderChl);

    std::cout << "Send 0, 0 \t= " << sendMessages[0][0] << std::endl;
    std::cout << "Send 0, 1 \t= " << sendMessages[0][1] << std::endl;
    std::cout << "Send " << n-1 << ", 0 \t= " << sendMessages[n-1][0] << std::endl;
    std::cout << "Send " << n-1 << ", 1 \t= " << sendMessages[n-1][1] << std::endl;
  } else {
    osuCrypto::IknpOtExtSender sender;
    sender.sendChosen(sendMessages, prng, senderChl);
  }
  recverThread.join();
}

void test_send(bool silent) {
  std::cout << (silent ? "Silent test" : "IKNP test") << std::endl;
  osuCrypto::IOService ios;
  osuCrypto::Channel channel = osuCrypto::Session(ios, ADDR, PORT, osuCrypto::SessionMode::Server).addChannel();
  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());
  int n = 1;

  // osuCrypto::IknpOtExtSender sender;
  osuCrypto::SilentOtExtSender sender;

  std::vector<std::array<osuCrypto::block, 2>> sendMessages(n);
  sendMessages[0] = { osuCrypto::toBlock(54), osuCrypto::toBlock(33) };
  //...

  std::cout << "Send 0, 0 = " << sendMessages[0][0] << std::endl;
  std::cout << "Send 0, 1 = " << sendMessages[0][1] << std::endl;

  if (silent) {
    osuCrypto::SilentOtExtSender sender;
    sender.silentSend(sendMessages, prng, channel);
  } else {
    osuCrypto::IknpOtExtSender sender;
    sender.sendChosen(sendMessages, prng, channel);
  }
}

void test_recv(bool silent) {
  std::cout << (silent ? "Silent test" : "IKNP test") << std::endl;
  osuCrypto::IOService ios;
  osuCrypto::Channel channel = osuCrypto::Session(ios, ADDR, PORT, osuCrypto::SessionMode::Client).addChannel();
  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());
  int n = 1;

  osuCrypto::BitVector choices(n);
  choices[0] = 1;
  //...

  // Receive the messages
  std::vector<osuCrypto::block> messages(n);
  if(silent) {
    osuCrypto::SilentOtExtReceiver recver;
    recver.silentReceive(choices, messages, prng, channel);
  } else {
    osuCrypto::IknpOtExtReceiver recver;
    recver.receiveChosen(choices, messages, prng, channel);
  }

  // messages[i] = sendMessages[i][choices[i]];
  std::cout << "Selecting 0 on " << choices[0] << std::endl;
  std::cout << "Got 0 = " << messages[0] << std::endl;
}

void test2() {
  std::cout << "Test 2: silent test mimic" << std::endl;
  int threads = 1;
  int n = 888;

  // ??? defaults? effects minimum?
  int s = 4;
  int sec = 80;

  osuCrypto::IOService ios;
  osuCrypto::Session s0(ios, "localhost:1212", osuCrypto::SessionMode::Server);
  osuCrypto::Session s1(ios, "localhost:1212", osuCrypto::SessionMode::Client);
  // std::vector<osuCrypto::Channel> chls0(threads), chls1(threads);
  // for (unsigned int i = 0; i < threads; ++i){
  //   chls0[i] = s0.addChannel();
  //   chls1[i] = s1.addChannel();
  // }
  osuCrypto::Channel chls0 = s0.addChannel();
  osuCrypto::Channel chls1 = s1.addChannel();
  osuCrypto::PRNG prng(osuCrypto::toBlock(0, 0));

  osuCrypto::SilentOtExtSender sender;
  osuCrypto::SilentOtExtReceiver recver;

  // fake base ots
  {
    recver.configure(n, s, sec, threads);
    // recver.configure(n);
    // osuCrypto::BitVector choices = recver.sampleBaseChoiceBits(prng);
    // std::vector<osuCrypto::block> msg(choices.size());
    // for (unsigned int i = 0; i < msg.size(); ++i)
    //   msg[i] = osuCrypto::toBlock(i, choices[i]);
    // recver.setSlientBaseOts(msg);
  }

  {
    sender.configure(n, s, sec, threads);
    // sender.configure(n);
    // auto count = sender.silentBaseOtCount();
    // std::cout << "sender silentBaseOtCount = " << count << std::endl;
    // std::vector<std::array<osuCrypto::block, 2>> msg(count);
    // osuCrypto::PRNG prngz(osuCrypto::ZeroBlock);
    // for (unsigned int i = 0; i < msg.size(); ++i) {
    //   msg[i][0] = osuCrypto::toBlock(i, 0);
    //   msg[i][1] = osuCrypto::toBlock(i, 1); 
    // }
    // sender.setSlientBaseOts(msg);
  }

  /*
  Sender: 
    has m0, m1
    OT gives (r0, r1)
    get d
    send s0 = (m0 ^ r(d)), s1 = (m1 ^ r(!d))
  Receiver:
    has c
    OT gives (rb, b)
    send d = c ^ b
    get both si = (mi ^ r(d ^ i))
    find mc = sc ^ rb (= (mc ^ r(c ^ b ^ c)) ^ rb = mc ^ rb ^ rb)
  */

  osuCrypto::BitVector choice(n);
  std::vector<std::array<osuCrypto::block, 2>> sendMessages(n);

  for (int i = 0; i < n; i++) {
    sendMessages[i] = { osuCrypto::toBlock(i, 0), osuCrypto::toBlock(i, 1) };
    choice[i] = i % 2;
  }

  osuCrypto::BitVector randChoice(n);
  std::vector<std::array<osuCrypto::block, 2>> randMessages(n);
  std::vector<osuCrypto::block> messages2(n);

  std::cout << "sending..." << std::endl;
  sender.silentSend(randMessages, prng, chls0);
  std::cout << "sent. Receiving..." << std::endl;
  recver.silentReceive(randChoice, messages2, prng, chls1);
  std::cout << "received" << std::endl;
  // sender.silentSend(randMessages, prng, chls0);
  // recver.silentReceive(randChoice, messages2, prng, chls1);
  // sender.send(sendMessages, prng, chls0);
  // recver.receive(choice, messages2, prng, chls1);

  std::vector<osuCrypto::block> gotMessages(n);

  for (int i = 0; i < n; i++) {
    bool d = choice[i] ^ randChoice[i];
    // r->s: d
    osuCrypto::block s0 = sendMessages[i][0] ^ randMessages[i][d];
    osuCrypto::block s1 = sendMessages[i][1] ^ randMessages[i][d^1];
    // s->r: s0, s1
    osuCrypto::block mb = (choice[i] == 0 ? s0 : s1) ^ messages2[i];
    gotMessages[i] = mb;

    if (i % 101 == 0) {
      std::cout << std::endl << "i: " << i << std::endl;
      std::cout << "set m0  = " << sendMessages[i][0] << std::endl;
      std::cout << "set m1  = " << sendMessages[i][1] << std::endl;
      // std::cout << "rand0   = " << randMessages[i][0] << std::endl;
      // std::cout << "rand1   = " << randMessages[i][1] << std::endl;
      // std::cout << "rand c  = " << randChoice[i] << std::endl;
      // std::cout << "rand[c] = " << messages2[i] << std::endl;
      // std::cout << "d = b^c = " << d << std::endl;
      // std::cout << "s0      = " << s0 << std::endl;
      // std::cout << "s1      = " << s1 << std::endl;
      std::cout << "Got " << choice[i] << "   = " << gotMessages[i] << std::endl;
    }
  }
}

int main(int argc, char* argv[]) {
  int server_num = -1;
  if(argc >= 2){
    server_num = atoi(argv[1]);
  }

  if (server_num == -1) {
    // test(false);
    test(true);
    // test2();
  } else if (server_num == 0) {
    // test_send(false);
    test_send(true);
  } else if (server_num == 1) {
    // test_recv(false);
    test_send(true);
  }
}