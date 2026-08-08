/*****************************************************************************

 AUTHORS

   Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com> and other contributors

 LICENSE AND COPYRIGHT

   Copyright (c) 2008-2010 Michigan State University
   Copyright (c) 2011-2025 The University of Texas at Austin
   Copyright (c) 2025-     Michigan State University

   breseq is free software; you can redistribute it and/or modify it under the
   terms of the GNU General Public License as published by the Free Software
   Foundation; either version 2, or (at your option) any later version.

   SPDX-License-Identifier: GPL-2.0-or-later

*****************************************************************************/

#include "portable_random.h"

using namespace std;

namespace breseq {

  cPortableRandom::cPortableRandom(uint32_t seed, uint32_t stream_index)
  {
    // std::seed_seq::generate is specified exactly by the standard, so this mixing step gives the
    // same engine state everywhere. Seeding the engine directly with (seed + stream_index) would
    // also be portable, but adjacent seeds give mt19937 highly correlated initial output; seed_seq
    // exists precisely to avoid that.
    std::seed_seq seq{static_cast<uint32_t>(seed), static_cast<uint32_t>(stream_index)};
    m_engine.seed(seq);
  }

  uint32_t cPortableRandom::next_below(uint32_t n)
  {
    ASSERT(n > 0, "cPortableRandom::next_below requires n > 0");

    // mt19937 emits exactly 2^32 equally likely values. If n does not divide 2^32 there is a
    // ragged bucket at the top; discarding it makes every residue class equally likely, so the
    // modulo below is unbiased. Expected discards are under one draw for any n.
    const uint64_t kEngineOutcomes = 4294967296ULL;   // 2^32
    const uint64_t limit = (kEngineOutcomes / static_cast<uint64_t>(n)) * static_cast<uint64_t>(n);

    uint64_t x;
    do {
      x = static_cast<uint64_t>(m_engine());
    } while (x >= limit);

    return static_cast<uint32_t>(x % static_cast<uint64_t>(n));
  }

  uint32_t cPortableRandom::next_in_range_1(uint32_t lo_1, uint32_t hi_1)
  {
    ASSERT(lo_1 <= hi_1, "cPortableRandom::next_in_range_1 requires lo_1 <= hi_1");
    return lo_1 + next_below(hi_1 - lo_1 + 1);
  }

  int32_t cPortableRandom::next_gaussian(int32_t mean, int32_t stdev)
  {
    if (stdev == 0) return mean;

    // Irwin-Hall: the sum of 12 independent uniforms on {0..65535} has mean 12*65535/2 = 393210 and
    // variance 12*(65536^2-1)/12 = 65536^2-1, i.e. a standard deviation of essentially 65536. So
    // centering the sum and dividing by 65536 yields a unit-variance deviate, and multiplying by
    // stdev scales it -- all without a square root, because the 12 was chosen to make the variance
    // come out to exactly the range squared. The result is truncated at +/- 6 sigma by construction.
    //
    // next_below(65536) never rejects, since 65536 divides 2^32 exactly, so this costs exactly 12
    // engine draws every time.
    int64_t sum = 0;
    for (int i = 0; i < 12; i++) sum += static_cast<int64_t>(next_below(65536u));
    sum -= 393210;

    // Round-half-up divide by 65536, for negative sums as well as positive ones. The bias term is a
    // multiple of 65536, so adding it before the division and subtracting it after is exact; going
    // through it avoids relying on the sign behavior of >> for negative operands, which was
    // implementation-defined before C++20.
    const int64_t kBias = 1LL << 40;
    int64_t scaled = sum * static_cast<int64_t>(stdev) + 32768 + kBias;
    int64_t offset = (scaled / 65536) - (kBias / 65536);

    return mean + static_cast<int32_t>(offset);
  }

} // namespace breseq
