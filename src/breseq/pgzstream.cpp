// ============================================================================
// pgzstream, parallel C++ iostream classes wrapping the zlib compression
// library. See pgzstream.h for an overview.
// ============================================================================

#include "pgzstream.h"

#include <cstring>

pgzstreambuf::pgzstreambuf()
: opened(0)
, mode(0)
, num_threads(1)
, m_max_pending(1)
, m_any_chunk_written(false)
, m_max_queue_size(4)
, m_reader_done(false)
, m_stop_requested(false)
{
}

pgzstreambuf::~pgzstreambuf()
{
  close();
}

pgzstreambuf* pgzstreambuf::open(const char* name, int open_mode, unsigned int _num_threads)
{
  if (is_open())
    return 0;

  mode = open_mode;
  num_threads = std::max<unsigned int>(1, _num_threads);

  // no append nor read/write mode
  if ((mode & std::ios::ate) || (mode & std::ios::app)
      || ((mode & std::ios::in) && (mode & std::ios::out)))
    return 0;

  if (mode & std::ios::in) {

    m_file.open(name, std::ios::in | std::ios::binary);
    if (!m_file.is_open())
      return 0;

    m_reader_done = false;
    m_stop_requested = false;
    m_reader_thread = std::thread(&pgzstreambuf::reader_thread_func, this);

  } else if (mode & std::ios::out) {

    m_file.open(name, std::ios::out | std::ios::trunc | std::ios::binary);
    if (!m_file.is_open())
      return 0;

    m_obuf.resize(CHUNK_SIZE);
    setp(&m_obuf[0], &m_obuf[0] + m_obuf.size());

    m_max_pending = 2 * num_threads;
    m_pool.reset(new ctpl::thread_pool(static_cast<int>(num_threads)));
    m_any_chunk_written = false;

  } else {
    return 0;
  }

  opened = 1;
  return this;
}

pgzstreambuf* pgzstreambuf::close()
{
  if (!is_open())
    return 0;

  if (mode & std::ios::out) {

    // Flush any remaining buffered bytes as the final chunk. If nothing was
    // ever written, this writes a single empty gzip member so the output is
    // still a valid (empty) gzip file.
    flush_chunk(true);
    while (!m_pending.empty())
      drain_one();

    m_pool.reset();
    m_file.close();

  } else if (mode & std::ios::in) {

    {
      std::unique_lock<std::mutex> lock(m_queue_mutex);
      m_stop_requested = true;
      m_queue_cv.notify_all();
    }
    if (m_reader_thread.joinable())
      m_reader_thread.join();

    m_file.close();
  }

  opened = 0;
  return this;
}

int pgzstreambuf::sync()
{
  // Intentionally a no-op: the internal chunk buffer is not flushed early so
  // that std::endl / std::flush on the wrapping stream do not fragment the
  // gzip output into many tiny members. Final flushing happens in close().
  return 0;
}

int pgzstreambuf::overflow(int c)
{
  if (!(mode & std::ios::out) || !opened)
    return EOF;

  flush_chunk(false);

  if (c != EOF) {
    *pptr() = static_cast<char>(c);
    pbump(1);
  }
  return c;
}

void pgzstreambuf::flush_chunk(bool force_even_if_empty)
{
  size_t n = pptr() - pbase();

  if (n == 0) {
    if (!force_even_if_empty || m_any_chunk_written)
      return;
  } else {
    setp(&m_obuf[0], &m_obuf[0] + m_obuf.size());
  }

  std::string data(pbase(), n);
  m_any_chunk_written = true;

  std::future<std::string> fut = m_pool->push([data](int) {
    return pgzstreambuf::compress_chunk(data);
  });
  m_pending.push_back(std::move(fut));

  while (m_pending.size() > m_max_pending)
    drain_one();
}

void pgzstreambuf::drain_one()
{
  std::string compressed = m_pending.front().get();
  m_pending.pop_front();
  m_file.write(compressed.data(), static_cast<std::streamsize>(compressed.size()));
}

std::string pgzstreambuf::compress_chunk(const std::string& data)
{
  z_stream strm;
  memset(&strm, 0, sizeof(strm));

  // windowBits = 15 + 16 selects gzip framing (header + CRC32 trailer), so
  // each compressed chunk is itself a complete, independent gzip member.
  // Concatenating such members produces a valid multi-member gzip stream.
  deflateInit2(&strm, Z_DEFAULT_COMPRESSION, Z_DEFLATED, 15 + 16, 8, Z_DEFAULT_STRATEGY);

  std::string out;
  out.resize(deflateBound(&strm, static_cast<uLong>(data.size())));

  strm.next_in = reinterpret_cast<Bytef*>(const_cast<char*>(data.data()));
  strm.avail_in = static_cast<uInt>(data.size());
  strm.next_out = reinterpret_cast<Bytef*>(&out[0]);
  strm.avail_out = static_cast<uInt>(out.size());

  deflate(&strm, Z_FINISH);
  out.resize(strm.total_out);

  deflateEnd(&strm);
  return out;
}

int pgzstreambuf::underflow()
{
  if (gptr() && (gptr() < egptr()))
    return *reinterpret_cast<unsigned char*>(gptr());

  if (!(mode & std::ios::in) || !opened)
    return EOF;

  if (!pop_decompressed(m_current_chunk))
    return EOF;

  char* base = &m_current_chunk[0];
  setg(base, base, base + m_current_chunk.size());
  return *reinterpret_cast<unsigned char*>(gptr());
}

bool pgzstreambuf::pop_decompressed(std::string& out)
{
  std::unique_lock<std::mutex> lock(m_queue_mutex);
  m_queue_cv.wait(lock, [this] {
    return !m_decompressed_queue.empty() || m_reader_done;
  });

  if (m_decompressed_queue.empty())
    return false;

  out = std::move(m_decompressed_queue.front());
  m_decompressed_queue.pop_front();
  m_queue_cv.notify_all();
  return true;
}

void pgzstreambuf::push_decompressed(std::string&& chunk)
{
  std::unique_lock<std::mutex> lock(m_queue_mutex);
  m_queue_cv.wait(lock, [this] {
    return (m_decompressed_queue.size() < m_max_queue_size) || m_stop_requested;
  });

  if (m_stop_requested)
    return;

  m_decompressed_queue.push_back(std::move(chunk));
  m_queue_cv.notify_all();
}

void pgzstreambuf::reader_thread_func()
{
  z_stream strm;
  memset(&strm, 0, sizeof(strm));

  // windowBits = 15 + 32 enables automatic gzip/zlib header detection.
  // After Z_STREAM_END, inflateReset() preserves this setting, so
  // concatenated (multi-member) gzip streams are decoded correctly.
  inflateInit2(&strm, 15 + 32);

  std::vector<char> inbuf(CHUNK_SIZE);
  std::vector<char> outbuf(CHUNK_SIZE);
  bool input_eof = false;

  while (true) {
    {
      std::unique_lock<std::mutex> lock(m_queue_mutex);
      if (m_stop_requested)
        break;
    }

    if (strm.avail_in == 0 && !input_eof) {
      m_file.read(inbuf.data(), static_cast<std::streamsize>(inbuf.size()));
      std::streamsize got = m_file.gcount();
      if (got > 0) {
        strm.avail_in = static_cast<uInt>(got);
        strm.next_in = reinterpret_cast<Bytef*>(inbuf.data());
      } else {
        input_eof = true;
      }
    }

    if (strm.avail_in == 0 && input_eof)
      break;

    strm.avail_out = static_cast<uInt>(outbuf.size());
    strm.next_out = reinterpret_cast<Bytef*>(outbuf.data());

    int ret = inflate(&strm, Z_NO_FLUSH);

    size_t produced = outbuf.size() - strm.avail_out;
    if (produced > 0)
      push_decompressed(std::string(outbuf.data(), produced));

    if (ret == Z_STREAM_END) {
      inflateReset(&strm);
    } else if (ret != Z_OK && ret != Z_BUF_ERROR) {
      break;
    }
  }

  inflateEnd(&strm);

  std::unique_lock<std::mutex> lock(m_queue_mutex);
  m_reader_done = true;
  m_queue_cv.notify_all();
}
