/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#ifndef __SETUPAIDE
#define __SETUPAIDE

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "nekrsSys.hpp"
#include <map>

class setupAide
{
private:
public:
  setupAide(){};
  ~setupAide() = default;

  setupAide(const setupAide &) = default;
  setupAide &operator=(const setupAide &) = default;

  std::map<std::string, std::string> keyWordToDataMap;

  std::string getArgs(std::string) const;

  void removeArgs(std::string key);

  void setArgs(std::string key, std::string value);

  template <class T> int getArgs(std::string, T &) const;

  template <class T> int getArgs(std::string, std::vector<T> &) const;

  int getArgs(std::string, std::vector<std::string> &, std::string) const;

  int compareArgs(std::string key, std::string token) const;

  auto begin() const
  {
    return keyWordToDataMap.begin();
  }

  auto end() const
  {
    return keyWordToDataMap.end();
  }

  friend std::ostream &operator<<(std::ostream &out, const setupAide &aide);
};

#include <setupAide.tpp>

#endif

// Combined header/source needs explicit instantiation
