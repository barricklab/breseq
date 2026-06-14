// ============================================================================
// pgzstream, parallel C++ iostream classes wrapping the zlib compression
// library.
//
// Drop-in replacement for gzstream's igzstream/ogzstream/iogzstream that
// compresses on output using a pool of worker threads (writing a valid
// multi-member gzip stream, RFC 1952) and that decompresses on input using a
// background thread so decompression overlaps with the consumer's
// processing.
// ============================================================================

#ifndef PGZSTREAM_H
#define PGZSTREAM_H 1

#include <iostream>
#include <fstream>
#include <string>
#include <deque>
#include <vector>
#include <future>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <zlib.h>

#include "ctpl_stl.h"

// ----------------------------------------------------------------------------
// Internal class implementing the parallel gzip streambuf.
// ----------------------------------------------------------------------------

class pgzstreambuf : public std::streambuf {
public:
  pgzstreambuf();
  ~pgzstreambuf();

  int is_open() const { return opened; }

  pgzstreambuf* open(const char* name, int open_mode, unsigned int num_threads);
  pgzstreambuf* close();

protected:
  virtual int underflow();
  virtual int overflow(int c = EOF);
  virtual int sync();

private:
  static const size_t CHUNK_SIZE = 1 << 20; // 1 MiB uncompressed chunks

  static std::string compress_chunk(const std::string& data);

  // write side
  void flush_chunk(bool force_even_if_empty);
  void drain_one();

  // read side
  void reader_thread_func();
  bool pop_decompressed(std::string& out);
  void push_decompressed(std::string&& chunk);

  char opened;
  int mode;
  unsigned int num_threads;

  std::fstream m_file;

  // write side state
  std::vector<char> m_obuf;
  std::unique_ptr<ctpl::thread_pool> m_pool;
  std::deque<std::future<std::string>> m_pending;
  size_t m_max_pending;
  bool m_any_chunk_written;

  // read side state
  std::thread m_reader_thread;
  std::deque<std::string> m_decompressed_queue;
  std::mutex m_queue_mutex;
  std::condition_variable m_queue_cv;
  size_t m_max_queue_size;
  bool m_reader_done;
  bool m_stop_requested;
  std::string m_current_chunk;
};

class pgzstreambase : virtual public std::ios {
protected:
  pgzstreambuf buf;
public:
  pgzstreambase() { init(&buf); }
  pgzstreambase(const char* name, int mode, unsigned int num_threads) {
    init(&buf);
    open(name, mode, num_threads);
  }
  ~pgzstreambase() {
    buf.close();
  }
  void open(const char* name, int open_mode, unsigned int num_threads) {
    if (!buf.open(name, open_mode, num_threads))
      clear(rdstate() | std::ios::badbit);
  }
  void close() {
    if (buf.is_open())
      if (!buf.close())
        clear(rdstate() | std::ios::badbit);
  }
  pgzstreambuf* rdbuf() { return &buf; }
};

// ----------------------------------------------------------------------------
// User classes. Use ipgzstream and opgzstream analogously to igzstream and
// ogzstream respectively.
// ----------------------------------------------------------------------------

class ipgzstream : public pgzstreambase, public std::istream {
public:
  ipgzstream() : std::istream(&buf) {}
  ipgzstream(const char* name, int open_mode = std::ios::in, unsigned int num_threads = 1)
    : pgzstreambase(name, open_mode, num_threads), std::istream(&buf) {}
  pgzstreambuf* rdbuf() { return pgzstreambase::rdbuf(); }
  void open(const char* name, int open_mode = std::ios::in, unsigned int num_threads = 1) {
    pgzstreambase::open(name, open_mode, num_threads);
  }
};

class opgzstream : public pgzstreambase, public std::ostream {
public:
  opgzstream() : std::ostream(&buf) {}
  opgzstream(const char* name, int open_mode = std::ios::out, unsigned int num_threads = 1)
    : pgzstreambase(name, open_mode, num_threads), std::ostream(&buf) {}
  pgzstreambuf* rdbuf() { return pgzstreambase::rdbuf(); }
  void open(const char* name, int open_mode = std::ios::out, unsigned int num_threads = 1) {
    pgzstreambase::open(name, open_mode, num_threads);
  }
};

class iopgzstream : public pgzstreambase, public std::iostream {
public:
  iopgzstream() : std::iostream(&buf) {}
  iopgzstream(const char* name, int open_mode = std::ios::out, unsigned int num_threads = 1)
    : pgzstreambase(name, open_mode, num_threads), std::iostream(&buf) {}
  pgzstreambuf* rdbuf() { return pgzstreambase::rdbuf(); }
  void open(const char* name, int open_mode = std::ios::out, unsigned int num_threads = 1) {
    pgzstreambase::open(name, open_mode, num_threads);
  }
};

#endif // PGZSTREAM_H
