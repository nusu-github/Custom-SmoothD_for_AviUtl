/*******************************************************************
 
 *******************************************************************/

#pragma once

#define IDCT_CLIP_TABLE_OFFSET 512

#ifdef __cplusplus
extern "C" {
#endif

#ifndef IDCT_CLIP_TABLE_C
extern const int idct_clip_table[1024];
#endif

#ifdef __cplusplus
}
#endif
