#ifndef _PYRODIGAL_IMPL_NEON_H
#define _PYRODIGAL_IMPL_NEON_H

#ifdef WIN32
#define export __declspec( dllexport )
#else
#define export extern
#endif

#include <stdint.h>

export void skippable_neon(const int8_t*, const uint8_t*, const uint8_t*, const int, const int, uint8_t*);

#endif
