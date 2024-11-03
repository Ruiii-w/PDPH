
// Copyright (c) Simon Fraser University & The Chinese University of Hong Kong. All rights reserved.
// Licensed under the MIT license.
#ifndef HASH_INTERFACE_H_
#define HASH_INTERFACE_H_

// #include "./pair.h"
#include "util.h"

/*
* Parent function of all hash indexes
* Used to define the interface of the hash indexes
*/

template <class T>
class Hash {
 public:
  Hash(void) = default;
  ~Hash(void) = default;
  virtual int Insert(T, ValueType) = 0;
  virtual int Insert(T, ValueType, bool) = 0;
  virtual ValueType Get(T) = 0;
  virtual ValueType Get(T key, bool is_in_epoch) = 0;
  virtual void getNumber() = 0;
};

#endif  // _HASH_INTERFACE_H_
