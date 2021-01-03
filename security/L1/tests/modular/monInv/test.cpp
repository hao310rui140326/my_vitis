/*
 * Copyright 2019 Xilinx, Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <hls_stream.h>
#define AP_INT_MAX_W 4097
#include <ap_int.h>

#include <cstdio>
#include <cstring>
#include <iomanip>
#include <iostream>
using namespace std;
#include <sstream>
#include <string>
#include <vector>

#include "xf_security/modular.hpp"

void test(ap_uint<2048> a, ap_uint<2048> m, ap_uint<2048>& inv) {
    inv = xf::security::internal::monInv<2048>(a, m);
}

/*
void test(ap_uint<5> a, ap_uint<5> m, ap_uint<5>& res) {
    res = xf::security::internal::modularInv<5>(a, m);
}
*/
/*
void test(ap_uint<32> a, ap_uint<32> b, ap_uint<32> m, ap_uint<32>& ab) {
    ab = xf::security::internal::monProduct<32>(a, b, m);
}
*/
#ifndef __SYNTHESIS__
int main() {
    /*
    ap_uint<5> a = 10;
    ap_uint<5> m = 17;
    ap_uint<5> res;
    test(a, m, res);
    std::cout << "res: " << res << std::endl;
    ap_uint<10> tmp;
    tmp = res * a;
    tmp %= m;
    ap_uint<10> tmp2;
    tmp2 = 1;
    tmp2 <<= 5;
    tmp2 %= m;
    std::cout << "tmp:" << tmp << std::endl;
    std::cout << "tmp2:" << tmp2 << std::endl;
    */

    ap_uint<2048> m = ap_uint<2048>(
        "0x9d41cd0d38339220ebd110e8c31feb279c5fae3c23090045a0886301588d4c8114fa5cdde708ea77ba0f527e6f6ea8f5634acf517f04"
        "ca6399e188d5c2d7f03cc90e04dbf7d5d0056ee1b14b8baaf90ef78f5142ddce9ba2eff84c0295f656c29aecaae80ddd5c7127ddc60215"
        "9458f272316100f726a71362516223f26ddeafa425d3eb2c7f61de7e8586e77d475037563425d931885f03693618bb885ab9b58de74f60"
        "4a86f28e494dcd819bd8c0bb42f699596969b84f680819e4c9fc0ba687558775f770a302d5b266905defe47bc53c98ce261523b49db624"
        "1567f4b48c661482ef9c453750c6d420a0b1a3bd4d3d05b060c026ce8efd9bb9456dfe2f5d");

    ap_uint<2048> a;
    unsigned char rawm[256] =
        "RSA TEST FILE : "
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuv"
        "wxyz1234567890ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefgh"
        "ijklmnopqrstuvwxyz.";
    for (int i = 0; i < 256; i++) {
        a.range(i * 8 + 7, i * 8) = (unsigned int)rawm[255 - i];
    }

    ap_uint<2048> result;

    test(a, m, result);

    ap_uint<4096> tmp = a * result;
    tmp %= m;

    ap_uint<4096> tmp2 = 1;
    tmp2 <<= 2048;
    tmp2 %= m;

    std::cout << "tmp: " << std::hex << tmp << std::endl;
    std::cout << "tmp2: " << std::hex << tmp2 << std::endl;

    if (tmp != tmp2) {
        return 1;
    } else {
        return 0;
    }
}

#endif
