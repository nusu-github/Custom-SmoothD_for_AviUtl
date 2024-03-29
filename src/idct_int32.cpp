/*******************************************************************
  IDCT module used LLM algorithm (based IJG idct_int.c)
 *******************************************************************/

#include "idct_int32.hpp"

#include "idct_clip_table.hpp"

enum {
  FIX_0_298631336 = 2446,
  FIX_0_390180644 = 3196,
  FIX_0_541196100 = 4433,
  FIX_0_765366865 = 6270,
  FIX_0_899976223 = 7373,
  FIX_1_175875602 = 9633,
  FIX_1_501321110 = 12299,
  FIX_1_847759065 = 15137,
  FIX_1_961570560 = 16069,
  FIX_2_053119869 = 16819,
  FIX_2_562915447 = 20995,
  FIX_3_072711026 = 25172
};

void idct_int32(int32_t &block) {

  int32_t w0, w1, w2, w3, w4, w5, w6, w7;
  int32_t z1, z2, z3, z4, z5;
  int32_t *w, *d, *s;

  int32_t work[64]{};

  s = &block;
  w = work;

  for (uint8_t i = 0; i < 8; ++i) {
    if ((s[1] | s[2] | s[3] | s[4] | s[5] | s[6] | s[7]) == 0) {
      w[0] = w[1] = w[2] = w[3] = w[4] = w[5] = w[6] = w[7] = s[0] << 4;
      s += 8;
      w += 8;
      continue;
    }

    z2 = s[2];
    z3 = s[6];

    z1 = (z2 + z3) * FIX_0_541196100;
    w2 = z1 + z3 * -FIX_1_847759065;
    w3 = z1 + z2 * FIX_0_765366865;

    w0 = (s[0] + s[4]) << 13;
    w1 = (s[0] - s[4]) << 13;

    w4 = w0 + w3;
    w5 = w1 + w2;
    w6 = w1 - w2;
    w7 = w0 - w3;

    w0 = s[7];
    w1 = s[5];
    w2 = s[3];
    w3 = s[1];

    z1 = w0 + w3;
    z2 = w1 + w2;
    z3 = w0 + w2;
    z4 = w1 + w3;
    z5 = (z3 + z4) * FIX_1_175875602;

    w0 *= FIX_0_298631336;
    w1 *= FIX_2_053119869;
    w2 *= FIX_3_072711026;
    w3 *= FIX_1_501321110;
    z1 *= -FIX_0_899976223;
    z2 *= -FIX_2_562915447;
    z3 *= -FIX_1_961570560;
    z4 *= -FIX_0_390180644;

    z3 += z5;
    z4 += z5;

    w0 += z1 + z3;
    w1 += z2 + z4;
    w2 += z2 + z3;
    w3 += z1 + z4;

    w[0] = (w4 + w3) >> 9;
    w[1] = (w5 + w2) >> 9;
    w[2] = (w6 + w1) >> 9;
    w[3] = (w7 + w0) >> 9;
    w[4] = (w7 - w0) >> 9;
    w[5] = (w6 - w1) >> 9;
    w[6] = (w5 - w2) >> 9;
    w[7] = (w4 - w3) >> 9;

    s += 8;
    w += 8;
  }

  w = work;
  d = &block;

  for (uint8_t i = 0; i < 8; ++i) {
    if ((w[1 * 8] | w[2 * 8] | w[3 * 8] | w[4 * 8] | w[5 * 8] | w[6 * 8] | w[7 * 8]) == 0) {
      d[0 * 8] = d[1 * 8] = d[2 * 8] = d[3 * 8] = d[4 * 8] = d[5 * 8] = d[6 * 8] = d[7 * 8] =
          idct_clip_table[IDCT_CLIP_TABLE_OFFSET + ((w[0] + 64) >> 7)];
      w++;
      d++;
      continue;
    }

    z2 = w[2 * 8];
    z3 = w[6 * 8];

    z1 = (z2 + z3) * FIX_0_541196100;
    w2 = z1 + z3 * -FIX_1_847759065;
    w3 = z1 + z2 * FIX_0_765366865;

    w0 = (w[0 * 8] + w[4 * 8]) << 13;
    w1 = (w[0 * 8] - w[4 * 8]) << 13;

    w4 = w0 + w3;
    w7 = w0 - w3;
    w5 = w1 + w2;
    w6 = w1 - w2;

    w0 = w[7 * 8];
    w1 = w[5 * 8];
    w2 = w[3 * 8];
    w3 = w[1 * 8];

    z1 = w0 + w3;
    z2 = w1 + w2;
    z3 = w0 + w2;
    z4 = w1 + w3;
    z5 = (z3 + z4) * FIX_1_175875602;

    w0 *= FIX_0_298631336;
    w1 *= FIX_2_053119869;
    w2 *= FIX_3_072711026;
    w3 *= FIX_1_501321110;
    z1 *= -FIX_0_899976223;
    z2 *= -FIX_2_562915447;
    z3 *= -FIX_1_961570560;
    z4 *= -FIX_0_390180644;

    z3 += z5;
    z4 += z5;

    w0 += z1 + z3;
    w1 += z2 + z4;
    w2 += z2 + z3;
    w3 += z1 + z4;

    d[0 * 8] = idct_clip_table[IDCT_CLIP_TABLE_OFFSET + ((w4 + w3 + 524288) >> 20)];
    d[1 * 8] = idct_clip_table[IDCT_CLIP_TABLE_OFFSET + ((w5 + w2 + 524288) >> 20)];
    d[2 * 8] = idct_clip_table[IDCT_CLIP_TABLE_OFFSET + ((w6 + w1 + 524288) >> 20)];
    d[3 * 8] = idct_clip_table[IDCT_CLIP_TABLE_OFFSET + ((w7 + w0 + 524288) >> 20)];
    d[4 * 8] = idct_clip_table[IDCT_CLIP_TABLE_OFFSET + ((w7 - w0 + 524288) >> 20)];
    d[5 * 8] = idct_clip_table[IDCT_CLIP_TABLE_OFFSET + ((w6 - w1 + 524288) >> 20)];
    d[6 * 8] = idct_clip_table[IDCT_CLIP_TABLE_OFFSET + ((w5 - w2 + 524288) >> 20)];
    d[7 * 8] = idct_clip_table[IDCT_CLIP_TABLE_OFFSET + ((w4 - w3 + 524288) >> 20)];

    w++;
    d++;
  }
}
