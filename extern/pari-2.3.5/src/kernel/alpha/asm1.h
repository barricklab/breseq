/* $Id: asm1.h 7697 2006-02-16 17:34:46Z kb $

Copyright (C) 2000  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */

/* This file is a slight adaptation of source code extracted from gmp-3.1.1
  (from T. Granlund), files longlong.h and gmp-impl.h

  Copyright (C) 2000 Free Software Foundation, Inc.

 * FIXME: This file is unused until somebody implements
 * invert_word(x) = return floor( 2^(2*BIL)/x ) */

extern ulong invert_word(ulong);

#define sub_ddmmss(sh, sl, ah, al, bh, bl) \
  do {                                     \
    ulong __x;                             \
    __x = (al) - (bl);                     \
    (sh) = (ah) - (bh) - (__x > (al));     \
    (sl) = __x;                            \
  } while (0)

#ifdef __GNUC__

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

#else /* __GNUC__ */

static ulong
divll(ulong x, ulong y)
{
  register ulong _di, _x = (x), _y = (y), _q, _ql, _r;
  register ulong _xh, _xl, _k, __hire;

  if (_y & 0x8000000000000000UL)
      { _k = 0; __hire = hiremainder; }
  else
  {
    _k = bfffo(_y);
    __hire = (hiremainder << _k) | (_x >> (64 - _k));
    _x <<= _k; _y <<=  _k;
  }
  _di = invert_word(_y);
  _ql = mulll (__hire, _di);
  _q = __hire + hiremainder;
  _xl = mulll(_q, _y); _xh = hiremainder;
  sub_ddmmss (_xh, _r, __hire, _x, _xh, _xl);
  if (_xh != 0)
  {
    sub_ddmmss (_xh, _r, _xh, _r, 0, _y); _q += 1;
    if (_xh != 0)
      { sub_ddmmss (_xh, _r, _xh, _r, 0, _y); _q += 1; }
  }
  if (_r >= _y)
    { _r -= _y; _q += 1; }
  hiremainder = _r >> _k;
  return _q;
}

#endif /* __GNUC__ */
