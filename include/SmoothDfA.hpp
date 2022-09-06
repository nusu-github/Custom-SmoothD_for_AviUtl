#pragma once

/*
 * 以下のグローバル変数ですが。
 * Aviutl内臓のマルチスレッド機能は二つの値しか関数で渡せないので、やむを得ずグローバル変数にしました。
 * ここまで増えたら構造体を使うべきなのでしょうか……?
 */
struct smoothdfa_struct {
  int threshold{};         // 閾値
  int quantization{};      // 量子化係数
  int n_shift{};           // DCT-iDCTを実行する回数
  int zero_weight{};       // ノイズ除去後の画像に混ぜる元の画像の割合
  bool rounding_on{};      // 8ビットの計算を丸めるか
  bool all_quantization{}; // すべて量子化するか
  bool dct_on{};           // 離散コサイン変換を使うか
};

struct Soa_PixelYC {
  int16_t *y;
  int16_t *cb;
  int16_t *cr;
};

void shift_data(AviUtl::FilterProcInfo *fpip);
void Loop(int thread_id, int thread_num, void *param1, void *param2);
auto get_multi_thread(AviUtl::FilterPlugin *fp) -> int;
void copy_pix(AviUtl::FilterProcInfo *fpip, AviUtl::PixelYC *wsp, int shiftx, int shifty);
auto func_proc(AviUtl::FilterPlugin *fp, AviUtl::FilterProcInfo *fpip) -> BOOL;
auto func_init(AviUtl::FilterPlugin *fp) -> BOOL;
auto func_update(AviUtl::FilterPlugin *fp, AviUtl::FilterPluginDLL::UpdateStatus status) -> BOOL;