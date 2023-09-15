#ifndef UTILSTESTCONNECT_H
#define UTILSTESTCONNECT_H

#include <arpa/inet.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <unistd.h>

#include <cstring>
#include <iostream>

#include "utils.h"

const int TEST_PORT = 8887;

int init_sender(const char* ip = "127.0.0.1") {
  std::cout << "send: start" << std::endl;

  int sockfd;
  struct sockaddr_in addr;

  sockfd = socket(AF_INET, SOCK_STREAM, 0);
  if (sockfd < 0) error_exit("send: ERROR opening socket");

  int reuse = 1;
  if (setsockopt(sockfd, SOL_SOCKET, SO_REUSEADDR, &reuse, sizeof(reuse)))
    error_exit("Sockopt failed");
  if (setsockopt(sockfd, SOL_SOCKET, SO_REUSEPORT, &reuse, sizeof(reuse)))
    error_exit("Sockopt failed");

  bzero((char *) &addr, sizeof(addr));
  addr.sin_port = htons(TEST_PORT);
  addr.sin_family = AF_INET;
  inet_pton(AF_INET, ip, &addr.sin_addr);

  if (connect(sockfd, (sockaddr*) &addr, sizeof(addr)) < 0)
      error_exit("ERROR on connect");

  return sockfd;
}

int init_receiver() {
  std::cout << "recv: start" << std::endl;

  int sockfd;
  struct sockaddr_in rec_addr;

  sockfd = socket(AF_INET, SOCK_STREAM, 0);
  if (sockfd < 0) error_exit("recv: ERROR opening socket");

  int reuse = 1;
  if (setsockopt(sockfd, SOL_SOCKET, SO_REUSEADDR, &reuse, sizeof(reuse)))
    error_exit("Sockopt failed");
  if (setsockopt(sockfd, SOL_SOCKET, SO_REUSEPORT, &reuse, sizeof(reuse)))
    error_exit("Sockopt failed");

  bzero((char *) &rec_addr, sizeof(rec_addr));
  rec_addr.sin_family = AF_INET;
  rec_addr.sin_addr.s_addr = INADDR_ANY;
  rec_addr.sin_port = htons(TEST_PORT);

  if (bind(sockfd, (struct sockaddr *) &rec_addr, sizeof(rec_addr)) < 0)
    error_exit("recv: ERROR on binding");

  if (listen(sockfd, 2) < 0) error_exit("recv: ERROR on listen");
  return sockfd;
}

int accept_receiver(int sockfd) {
  int newsockfd;
  socklen_t snd_len;
  struct sockaddr_in snd_addr;
  snd_len = sizeof(snd_addr);
  newsockfd = accept(sockfd, (struct sockaddr *) & snd_addr, &snd_len);
  if (newsockfd < 0) error_exit("ERROR on accept");

  printf("recv: got connection from %s port %d\n",
         inet_ntoa(snd_addr.sin_addr), ntohs(snd_addr.sin_port));
  return newsockfd;
}

#endif
