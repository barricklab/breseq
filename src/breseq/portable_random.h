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

#ifndef _BRESEQ_PORTABLE_RANDOM_H_
#define _BRESEQ_PORTABLE_RANDOM_H_

#include "common.h"

#include <random>

namespace breseq {

  //! A random number generator whose output is IDENTICAL on every platform and compiler.
  //
  //  WHY THIS EXISTS
  //  ---------------
  //  breseq's read simulator generates the input FASTQ for several regression tests at test time,
  //  instead of those reads being committed to the repository. The tests then compare a committed
  //  expected.gd, and CI runs them on BOTH ubuntu-latest and macos-latest
  //  (.github/workflows/build-and-test.yml). So a given seed has to produce byte-identical reads
  //  under glibc/libstdc++ and under Apple libc/libc++, or the suite cannot pass on both.
  //
  //  Three things in ordinary C++ break that, and all three are easy to reintroduce by accident:
  //
  //    1. rand()/srand() are not specified beyond their interface. glibc and Apple's libc use
  //       different algorithms and return different sequences from the same seed.
  //
  //    2. The <random> DISTRIBUTIONS (uniform_int_distribution, normal_distribution, ...) are not
  //       specified either. The standard pins the ENGINES exactly -- [rand.eng.mers] even fixes
  //       mt19937's 10000th output -- but leaves each distribution's mapping from engine output to
  //       result up to the implementation, and libstdc++ and libc++ genuinely differ. std::seed_seq
  //       IS fully specified, so seeding through it is safe.
  //
  //    3. libm transcendentals (pow, log, sqrt, cos, exp) are not required to be correctly rounded.
  //       Results differ in the last ulp between glibc and Apple's libm, which is enough to flip a
  //       comparison and take a different branch. (Plain IEEE + - * / and comparisons ARE correctly
  //       rounded and would be safe, but the cheapest way to be certain is to use no floating point
  //       at all -- which is what this class does.)
  //
  //  THE RULES FOR THIS FILE AND FOR EVERY CALLER THAT GENERATES SIMULATED DATA
  //  --------------------------------------------------------------------------
  //  ALLOWED:   std::mt19937, std::seed_seq, fixed-width integer arithmetic, table lookups.
  //  FORBIDDEN: rand()/srand(), everything in <random> except the engine and seed_seq, every libm
  //             function, and any double whose rounded value feeds a branch.
  //
  //  If you are tempted to "simplify" next_below() into a std::uniform_int_distribution, or
  //  next_gaussian() into a std::normal_distribution or a Box-Muller transform, please don't --
  //  that is precisely the change that makes the test suite pass on your machine and fail in CI.
  class cPortableRandom {
  public:

    //! Construct an independent stream for (seed, stream_index).
    //
    //  Giving each simulated read pair its own stream (rather than drawing them all from one
    //  sequence) means a pair's reads depend only on its own index, so changing how many draws
    //  earlier pairs make -- adding an option, reordering a call -- does not perturb every read
    //  after it. Seeding runs through std::seed_seq, whose generate() is specified exactly.
    cPortableRandom(uint32_t seed, uint32_t stream_index);

    //! Raw engine output, uniform over the full 32-bit range.
    uint32_t next_uint32() { return m_engine(); }

    //! Uniform in [0, n). Unbiased: the top partial bucket is rejected rather than folded in,
    //  which a bare modulo would do (making small residues slightly more likely).
    uint32_t next_below(uint32_t n);

    //! Uniform in [lo_1, hi_1], inclusive, 1-based-friendly.
    uint32_t next_in_range_1(uint32_t lo_1, uint32_t hi_1);

    //! Uniform true/false.
    bool next_coin() { return next_below(2) == 1; }

    //! True with probability rate/1e9. Integer arithmetic, so an error rate expressed this way is
    //  exactly reproducible where a double comparison against pow(10, -q/10.0) would not be.
    bool next_event_per_billion(uint32_t rate) { return next_below(1000000000u) < rate; }

    //! Approximately normal, rounded to an integer, truncated near +/- 6 sigma.
    //
    //  Uses the Irwin-Hall construction (sum of 12 uniforms) rather than Box-Muller, because the
    //  sum of 12 uniforms over a power-of-two range has variance almost exactly the square of that
    //  range -- so converting to a unit-variance deviate is a shift, with no sqrt/log/cos anywhere.
    //  Accuracy is far more than a fragment-length or quality-score distribution needs.
    int32_t next_gaussian(int32_t mean, int32_t stdev);

  private:
    std::mt19937 m_engine;
  };

} // namespace breseq

#endif
