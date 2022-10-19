#pragma once

// Aviutl SDK
#include <aviutl/FilterPlugin.hpp>
#include <aviutl/FilterProcInfo.hpp>

// DCT-iDCTを始める開始位置
// 最初の八個の位置は本家SmoothDから。残りは適当に手入力しました。
// 使ってみた印象では八個以上はほとんど変化が無い印象でしたが、試しに全ての位置を網羅しました。
constexpr int_fast8_t shift[] = {
    3, 1, 1, 4, 4, 7, 7, 3, 1, 1, 7, 1, 1, 7, 7, 7, 5, 2, 7, 6, 0, 1, 1, 7, 6, 0, 4, 6, 0, 3, 4, 1,
    7, 3, 4, 3, 3, 7, 3, 0, 1, 4, 6, 4, 0, 5, 5, 1, 5, 7, 0, 2, 7, 5, 0, 7, 4, 2, 2, 1, 7, 7, 1, 2,
    4, 5, 1, 6, 6, 2, 2, 3, 6, 6, 0, 4, 6, 1, 2, 6, 3, 1, 4, 7, 1, 5, 5, 3, 1, 3, 5, 0, 3, 5, 6, 5,
    7, 0, 3, 2, 5, 6, 1, 1, 3, 6, 7, 4, 2, 0, 5, 4, 6, 3, 2, 7, 3, 3, 4, 0, 2, 4, 4, 4, 7, 2, 0, 0,
};

struct smoothdfa_struct {
  int32_t *threshold;    // 閾値
  int32_t *quantization; // 量子化係数
  int32_t n_shift;       // DCT-iDCTを実行する回数
  int32_t *zero_weight;  // ノイズ除去後の画像に混ぜる元の画像の割合
  bool all_quantization; // すべて量子化するか
  bool dct_on;           // 離散コサイン変換を使うか
};

struct main_smoothdfa_struct {
  int32_t pictur_width;  // 画像の幅
  int32_t pictur_height; // 画像の高さ
  int32_t pictur_size;   // 画像のサイズ
  int32_t mem_size;      // メモリのサイズ
};

auto func_init(AviUtl::FilterPlugin *fp) -> BOOL;
auto func_exit(AviUtl::FilterPlugin *fp) -> BOOL;
auto func_update(AviUtl::FilterPlugin *fp, AviUtl::FilterPluginDLL::UpdateStatus status) -> BOOL;
auto func_proc(AviUtl::FilterPlugin *fp, AviUtl::FilterProcInfo *fpip) -> BOOL;

void get_multi_thread(int thread_id, int thread_num, void *param1, void * /*param2*/);

void shift_data(int thread_id, int thread_num, void *param1, void *param2);
void Loop(int thread_id, int thread_num, void *param1, void *param2);
void copy_pix(AviUtl::FilterProcInfo *fpip, AviUtl::PixelYC *wsp, const int32_t &shiftx, const int32_t &shifty);

void transform_quantization(const int32_t &threshold, const int32_t &quantization, const bool &all_quantization,
                            std::array<int32_t, 64> &block, std::array<int32_t, 64> &blockCb,
                            std::array<int32_t, 64> &blockCr);

void Loop_dct(AviUtl::FilterProcInfo *fpip, const int32_t &threshold, const int32_t &quantization,
              const bool &all_quantization, const int &thread_id, const int16_t &shiftx, const int16_t &shifty,
              const int16_t &xblock, const int16_t &yblock);
void Loop_fwht(AviUtl::FilterProcInfo *fpip, const int32_t &threshold, const int32_t &quantization,
               const bool &all_quantization, const int &thread_id, const int16_t &shiftx, const int16_t &shifty,
               const int16_t &xblock, const int16_t &yblock);