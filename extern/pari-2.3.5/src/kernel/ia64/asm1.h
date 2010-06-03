/* Extracted from gmp-4.1.2
 * FIXME: This file is unused until somebody implements
 * invert_word(x) = return floor( 2^(2*BIL)/x ) */
extern ulong invert_word(ulong);

#define sub_ddmmss(sh, sl, ah, al, bh, bl) \
({                                         \
    ulong __x;                             \
    __x = (al) - (bl);                     \
    (sh) = (ah) - (bh) - (__x > (al));     \
    (sl) = __x;                            \
})

#define divll(x, y)                                      \
({                                                       \
  register ulong _di, _x = (x), _y = (y), _q, _ql, _r;   \
  register ulong _xh, _xl, _k, __hire;                   \
                                                         \
  if (_y & 0x8000000000000000UL)                         \
      { _k = 0; __hire = hiremainder; }                  \
  else                                                   \
  {                                                      \
    _k = bfffo(_y);                                      \
    __hire = (hiremainder << _k) | (_x >> (64 - _k));    \
    _x <<= _k; _y <<=  _k;                               \
  }                                                      \
  _di = invert_word(_y);                                 \
  _ql = mulll (__hire, _di);                             \
  _q = __hire + hiremainder;                             \
  _xl = mulll(_q, _y); _xh = hiremainder;                \
  sub_ddmmss (_xh, _r, __hire, _x, _xh, _xl);            \
  if (_xh != 0)                                          \
  {                                                      \
    sub_ddmmss (_xh, _r, _xh, _r, 0, _y); _q += 1;       \
    if (_xh != 0)                                        \
      { sub_ddmmss (_xh, _r, _xh, _r, 0, _y); _q += 1; } \
  }                                                      \
  if (_r >= _y)                                          \
    { _r -= _y; _q += 1; }                               \
  hiremainder = _r >> _k;                                \
  _q;                                                    \
})
