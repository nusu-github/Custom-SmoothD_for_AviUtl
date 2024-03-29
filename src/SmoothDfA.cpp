//----------------------------------------------------------------------------------
//  SmoothD for AviUtl ver.0.01
//----------------------------------------------------------------------------------
/*
avisynthのフィルタのSmoothDが面白いと感じたので(使えるのかは、微妙ですが)Aviutlに移植しました。
SmoothDのアルゴリズムですが、SmoothDに同胞されているindex.htmlに書かれているように、以下で紹介されているものを使用しています。

http://www.utdallas.edu/~aria/mcl/post/

画像の0,0からずらした場所をスタート地点として、8*8の離散コサイン変換の量子化によるぼかしをおこないます。何度かスタート地点を変えてぼかしを行い、その画像を重ねることで、ブロックノイズを消すことができます。なんでもffdshowのSSP
deblockingもこのアルゴリズムを使用しているようです。 ただし、このSmoothD for
AviUtlでは、速度と作者の趣味のため、直交変換にアダマール変換を使用しています。とは言え。本家のSmoothDはMMX化されたxvidのDCT-iDCTを使用しているのでそれなりに高速です。
紹介したサイトで書いているように、このアルゴリズムの利点は、DCT-iDCTというありふれた(?)技術を使用していることです。このフィルタも、直交変換の関数に他のものを使用することで、処理を高速にすることが可能となるはずです。

//色差の処理はしていません。処理時間のこともありますが、どうしても色が変わってしまうので。16で除算して16で乗算するだけで色がたいそう変わるのはなぜなのでしょう?//
//色差を処理するようにしました。12ビットから8ビットまでデータ量を落とさずとも飽和しないはずなので、10ビットで処理することにし、8で除算しています。これならば色がほとんど変わりません。
また、輝度は9ビットに落とすようにしました。飽和しないのならば、できるだけデータ量を削減しないに越したことは無いですから。//
画素データを8ビットに落とす処理に戻しました。茂木氏のDCT-iDCTが8ビットを前提に作られていたこと、また、直交変換の関数に他のものを使うよう変更する際にも、8ビットで処理していた方が無難でしょうから。変色するのは今後の課題と言うことで。

マルチスレッド処理で画像にゴミが出るため、処理した画素データをwork_spaceに一時退避する処理に変更しました。メモリの同じ場所に同時に書き込もうとすると問題が出るようです。当然、なのかな?
コンパイラが何とかしてくれるのではないかと、ちょっと希望は持っていたのですが。ただし。マルチスレッドで処理しない場合、まったく無駄な処理が一つ増えただけです。ごめんなさい。

//量子化係数を使用した量子化ではなく、直交変換後の値が閾値以下なら0とする処理に変更しました。これにより輪郭がぼけにくくなりましたが、ブロックノイズが消えづらくなりました。閾値処理と量子化を併用する処理も試したのですが、ブロックノイズ除去には効果があったものの、処理時間がかなり増加してしまったため閾値処理のみとしました。//
//輝度のみ量子化と輪郭強調を行う処理にしたところ、それほど遅くならなかったので、量子化の処理を併用することにしました。ただし、デフォルトは量子化係数が"1"なので、実質おこなっていません。
閾値と量子化係数を上手く調整することで、2DNRとしても、デブロッキングとしても、機能するはずです。//
色々と試してみたところ、"閾値以下は量子化"の処理が最も効果的と判明しました。そのように処理するようにしています。

//閾値処理後の値を増加or減算する処理を追加しました。輪郭強調orぼかし処理ということで。はたして使える機能なのかは謎ですが、周波数で画像を処理するならば(それが例え足し算引き算のみのアダマール変換であったとしても)、この機能を付けたくなってしまうものなのです。//
この機能は削除しました。どう考えても使用しない機能なので、処理時間が遅くなるだけですから。

注意! 画像サイズが8の倍数専用です。判定はしていません。



2008/9/13
・正式公開。ver番号をつける。
・DCT-iDCTに使う配列をshort型にしたため、茂木和洋氏作のDCT-iDCTを使う部分をやむなく削除。
2008/9/14 ver0.002a
・マルチスレッド処理で画像にゴミが出るため、work_spaceに一時退避する処理に変更(処理時間が大幅増加する……)。
・輝度と色差の内部処理のビット数を増やす。
・色差を処理するようにしてコンパイル。
2008/9/15 ver0.003
・直交変換後に量子化ではなく閾値以下を削る処理に変更。
・輪郭強調(ぼかし)の処理を追加。
2008/9/16 ver0.003b
・輪郭強調(ぼかし)は輝度のみとした。
・輝度のみ量子化の処理を併用するようにした。
2008/9/22 ver0.01
・輝度は閾値以下を量子化する処理に変更。
・輪郭強調(ぼかし)の処理を削除。
・輝度と色差の内部処理のビット数を8ビットに。
・茂木氏のDCT-iDCTの部分を復活。
・トラックのデフォルト値を変更。
・ソースの整理。12ビット→8ビットの処理をマルチスレッドに。
*/
// C++ stl
#include <array>
#include <cmath>
#include <span>

#include <eve/module/core.hpp>
#include <eve/wide.hpp>

// Main Library
#include "dct_int32.hpp"
#include "ht.hpp"
#include "idct_int32.hpp"

#include "SmoothDfA.hpp"

// 処理した画像を一時的に代入するための変数
AviUtl::PixelYC *work_space;

// グローバル変数
// スレッド数
int32_t MT = 1;
// パラメーター系
smoothdfa_struct smoothdfa;

main_smoothdfa_struct main_smoothdfa;

//---------------------------------------------------------------------
//  フィルタ構造体定義
//---------------------------------------------------------------------
constexpr int TRACK_N = 4; // トラックバーの数
const char **track_name = new const char *[TRACK_N] { "閾値", "量子化係数", "品質", "重み" }; // トラックバーの名前
int track_default[TRACK_N] = {24, 8, 5, 0};     // トラックバーの初期値
int track_s[TRACK_N]       = {0, 1, 1, 0};      // トラックバーの下限値
int track_e[TRACK_N]       = {255, 255, 6, 64}; // トラックバーの上限値

constexpr int CHECK_N      = 3; // 　チェックボックスの数
const char **check_name = new const char *[CHECK_N] { "DCT", "変色防止", "全量子化" }; // チェックボックスの名前
int check_default[CHECK_N] = {1, 1, 1}; // チェックボックスの初期値

#ifdef _DEBUG
const char *filter_name = "SmoothD for AviUtl [DEBUG]"; // フィルタ名
#else
const char *filter_name = "SmoothD for AviUtl"; // フィルタ名
#endif

AviUtl::FilterPluginDLL filter{
    .flag          = AviUtl::FilterPlugin::Flag::ExInformation, // フィルタのフラグ
    .name          = filter_name,                               // フィルタの名前
    .track_n       = TRACK_N,                                   // トラックバーの数
    .track_name    = track_name,                                // トラックバーの名前
    .track_default = track_default,                             // トラックバーの初期値
    .track_s       = track_s,                                   // トラックバーの下限値
    .track_e       = track_e,                                   // トラックバーの上限値
    .check_n       = CHECK_N,                                   // チェックボックスの数
    .check_name    = check_name,                                // チェックボックスの名前
    .check_default = check_default,                             // チェックボックスの初期値 (値は0か1)
    .func_proc     = func_proc,                                 // フィルタの処理関数
    .func_init     = func_init,                                 // 初期化関数
    .func_exit     = func_exit,                                 // 終了関数
    .func_update   = func_update,                               // トラックバーの変更時に呼ばれる関数
};

//---------------------------------------------------------------------
//  フィルタ構造体のポインタを渡す関数
//---------------------------------------------------------------------
extern "C" [[maybe_unused]] AviUtl::FilterPluginDLL *GetFilterTable(void) { return &filter; }

//---------------------------------------------------------------------
//  Aviutl初期化関数
//---------------------------------------------------------------------
BOOL func_init(AviUtl::FilterPlugin *fp) {

  smoothdfa.threshold = &fp->track[0]; // 閾値。直交変換後の周波数の値がこの値以下ならば量子化を行う。
  smoothdfa.quantization = &fp->track[1]; // 量子化係数。
  smoothdfa.zero_weight  = &fp->track[3]; // ノイズ除去後の画像に混ぜる元の画像の割合

  return TRUE;
}

//---------------------------------------------------------------------
//  Aviutl終了処理関数
//---------------------------------------------------------------------
BOOL func_exit(AviUtl::FilterPlugin * /*fp*/) {
  _aligned_free(work_space);

  main_smoothdfa.pictur_width  = 0;
  main_smoothdfa.pictur_height = 0;
  main_smoothdfa.pictur_size   = 0;
  main_smoothdfa.mem_size      = 0;

  return TRUE;
}

//---------------------------------------------------------------------
//  トラックバーの変更時に呼ばれる関数
//---------------------------------------------------------------------
BOOL func_update(AviUtl::FilterPlugin *fp, AviUtl::FilterPluginDLL::UpdateStatus /*status*/) {

  // マルチスレッド数を取得
  fp->exfunc->exec_multi_thread_func(get_multi_thread, &MT, nullptr);

  // トラック番号を間違えないためトラックの値を変数に代入しておく
  smoothdfa.n_shift = static_cast<int32_t>(std::ceil(std::pow(2, fp->track[2]))); // DCT-iDCTを実行する回数
  smoothdfa.dct_on  = static_cast<bool>(fp->check[0]);          // 茂木氏作のDCT-iDCTを使用するかどうか
  smoothdfa.all_quantization = static_cast<bool>(fp->check[2]); // 全て量子化するかどうか

  return TRUE;
}

//---------------------------------------------------------------------
//  フィルタ処理関数
//---------------------------------------------------------------------
BOOL func_proc(AviUtl::FilterPlugin *fp, AviUtl::FilterProcInfo *fpip) {

  // 画像の幅と高さが変わっていたら初期化
  if (main_smoothdfa.pictur_width != fpip->w || main_smoothdfa.pictur_height != fpip->h) {
    _aligned_free(work_space);

    main_smoothdfa.pictur_width  = fpip->w;
    main_smoothdfa.pictur_height = fpip->h;
    main_smoothdfa.pictur_size   = fpip->w * fpip->h;
    main_smoothdfa.mem_size      = sizeof(AviUtl::PixelYC) * main_smoothdfa.pictur_size * MT;
    // マルチスレッド数だけ画像サイズのメモリを確保して0で埋める。
    // アライメントされたメモリを確保するためにmallocを使用。
    work_space = (AviUtl::PixelYC *)_aligned_malloc(main_smoothdfa.mem_size, alignof(AviUtl::PixelYC));
  } else {
    // メモリを0で埋める
    ZeroMemory(work_space, main_smoothdfa.mem_size);
  }

  const int32_t &zero_weight = *smoothdfa.zero_weight;
  const int32_t &n_shift     = smoothdfa.n_shift;

  // 画素データを8ビット化する。マルチスレッドにしてみました。
  fp->exfunc->exec_multi_thread_func(shift_data, fpip, &fp->check[1]);

  /*
  以下でぼかした画像を作業領域に加算する処理を繰り返す
  マルチスレッドで関数を呼び出して並列で処理しています。なのでn_shiftの値に関係なく4スレッドならば最低4回、そして4の倍数ずつ処理することになります。
  */
  int32_t count = zero_weight;

  for (int32_t i = 0; i < n_shift; i += MT) {
    fp->exfunc->exec_multi_thread_func(Loop, fpip, &i);
    count += MT;
  }

  /*
  作業領域に加算された値を変数 temp
  に合計、元画像をzero_weightだけ乗算して、12ビット化してから count で除算。
  */
  for (int32_t y = 0; y < fpip->h; ++y) {
    const std::span span_ycp(static_cast<AviUtl::PixelYC *>(fpip->ycp_edit) + y * fpip->max_w, fpip->w);
    for (int32_t x = 0; x < fpip->w; ++x) {
      AviUtl::PixelYC &ycp = span_ycp[x];
      AviUtl::PixelYC temp = {};

      for (int32_t i = 0; i < MT; ++i) {
        const AviUtl::PixelYC &wsp = work_space[y * fpip->w + x + main_smoothdfa.pictur_size * i];
        temp.y += wsp.y;
        temp.cb += wsp.cb;
        temp.cr += wsp.cr;
      }

      ycp.y  = ((ycp.y * zero_weight + temp.y) << 4) / count;
      ycp.cb = ((ycp.cb * zero_weight + temp.cb) << 4) / count;
      ycp.cr = ((ycp.cr * zero_weight + temp.cr) << 4) / count;
    }
  }

  return TRUE;
}

//---------------------------------------------------------------------
//  画素データ8ビット化関数
//---------------------------------------------------------------------
void shift_data(int thread_id, int thread_num, void *param1, void *param2) {

  if ((thread_id & 0x1) != 0)
    return;
  thread_id  = thread_id >> 1;
  thread_num = (thread_num + 1) >> 1;

  // マルチスレッド対応サンプルフィルタの方法そのままでマルチスレッド処理をしています
  const auto &fpip        = static_cast<AviUtl::FilterProcInfo *>(param1);

  const auto &rounding_on = static_cast<bool *>(param2); // 8bitにする時に四捨五入するかどうか
  //	スレッド毎の画像を処理する場所を計算する
  const int32_t &y_start = fpip->h * thread_id / thread_num;
  const int32_t &y_end   = fpip->h * (thread_id + 1) / thread_num;

  /*
   *	12ビットは一般的でないことと、加算していくと飽和してしまうかもしれないので、画素データのビット数を落とす。
   *	輝度と色差、どちらも8ビットに。
   */

  const auto x_end = fpip->w * 3;

  if (*rounding_on) {
    for (int32_t y = y_start; y < y_end; ++y) {
      constexpr auto loop_count = eve::expected_cardinal<int32_t>() * 2; // コンパイル時にSIMDのレーン数を取得する
      auto *ycp = static_cast<int16_t *>(fpip->ycp_edit) + (y * fpip->max_w) * 3;
      for (int32_t x = 0; x < x_end; x += loop_count) {
        const eve::wide<int16_t, eve::fixed<loop_count>> _ycp{&ycp[x]};
        const auto _ycp2 = eve::if_else(_ycp > 0, ((eve::int32(_ycp) * 100 / 16 + 45) / 100),
                                        ((eve::int32(_ycp) * 100 / 16 - 45) / 100));
        eve::store(eve::int16(_ycp2), ycp + x);
      }
    }
  } else {
    for (int32_t y = y_start; y < y_end; ++y) {
      constexpr auto loop_count = eve::expected_cardinal<int16_t>() * 2; // コンパイル時にSIMDのレーン数を取得する
      auto *ycp = static_cast<int16_t *>(fpip->ycp_edit) + (y * fpip->max_w) * 3;
      for (int32_t x = 0; x < x_end; x += loop_count) {
        eve::wide<int16_t, eve::fixed<loop_count>> _ycp{&ycp[x]};
        _ycp >>= 4;
        eve::store(_ycp, ycp + x);
      }
    }
  }
}

//---------------------------------------------------------------------
//  ぼかし処理関数
//---------------------------------------------------------------------
void Loop(int thread_id, int /*thread_num*/, void *param1, void *param2) {

  if ((thread_id & 0x1) != 0)
    return;
  thread_id                    = thread_id >> 1;

  auto *fpip                   = static_cast<AviUtl::FilterProcInfo *>(param1);

  const int32_t &threshold     = *smoothdfa.threshold;
  const int32_t &quantization  = *smoothdfa.quantization;
  const bool &all_quantization = smoothdfa.all_quantization;
  const bool &dct_on           = smoothdfa.dct_on;

  // DCT-iDCTの開始位置
  const int16_t &shiftx = shift[(*static_cast<int32_t *>(param2) + thread_id) * 2];
  const int16_t &shifty = shift[(*static_cast<int32_t *>(param2) + thread_id) * 2 + 1];
  // 0,0からずらしてDCT-iDCTをおこなうので画像の端まで処理はしない
  const int16_t h = fpip->h - 8;
  const int16_t w = fpip->w - 8;

  if (dct_on) {
    for (int16_t yblock = shifty; yblock < h; yblock += 8) {
      for (int16_t xblock = shiftx; xblock < w; xblock += 8) {
        Loop_dct(fpip, threshold, quantization, all_quantization, thread_id, shiftx, shifty, xblock, yblock);
      }
    }
  } else {
    for (int16_t yblock = shifty; yblock < h; yblock += 8) {
      for (int16_t xblock = shiftx; xblock < w; xblock += 8) {
        Loop_fwht(fpip, threshold, quantization, all_quantization, thread_id, shiftx, shifty, xblock, yblock);
      }
    }
  }

  // 画像四辺の処理をしていない部分にデータを代入
  copy_pix(fpip, work_space + fpip->w * fpip->h * thread_id, shiftx, shifty);
}

//---------------------------------------------------------------------
//  スレッド数取得関数
//---------------------------------------------------------------------
void get_multi_thread(int thread_id, int thread_num, void *param1, void * /*param2*/) {
  if (thread_id == 0) {
    *static_cast<int32_t *>(param1) = thread_num / 2;
  }
}

//---------------------------------------------------------------------
//  四辺未処理部分コピー関数
//---------------------------------------------------------------------
/*
長くなるので関数にして隔離。
*/
void copy_pix(AviUtl::FilterProcInfo *fpip, AviUtl::PixelYC *wsp, const int32_t &shiftx, const int32_t &shifty) {

  for (int32_t y = 0; y < shifty; ++y) {
    std::span span_ycp(static_cast<AviUtl::PixelYC *>(fpip->ycp_edit) + y * fpip->max_w, fpip->w);
    std::span span_ycp2(wsp + y * fpip->w, fpip->w);
    for (int32_t x = 0; x < fpip->w; ++x) {
      const AviUtl::PixelYC &ycp = span_ycp[x];
      AviUtl::PixelYC &ycp2      = span_ycp2[x];
      ycp2.y += ycp.y;
      ycp2.cb += ycp.cb;
      ycp2.cr += ycp.cr;
    }
  }

  for (int32_t y = fpip->h - (8 - shifty); y < fpip->h; ++y) {
    std::span span_ycp(static_cast<AviUtl::PixelYC *>(fpip->ycp_edit) + y * fpip->max_w, fpip->w);
    std::span span_ycp2(wsp + y * fpip->w, fpip->w);
    for (int32_t x = 0; x < fpip->w; ++x) {
      const AviUtl::PixelYC &ycp = span_ycp[x];
      AviUtl::PixelYC &ycp2      = span_ycp2[x];
      ycp2.y += ycp.y;
      ycp2.cb += ycp.cb;
      ycp2.cr += ycp.cr;
    }
  }

  for (int32_t y = shifty; y < fpip->h - (8 - shifty); ++y) {
    std::span span_ycp(static_cast<AviUtl::PixelYC *>(fpip->ycp_edit) + y * fpip->max_w, shiftx);
    std::span span_ycp2(wsp + y * fpip->w, shiftx);
    for (int32_t x = 0; x < shiftx; ++x) {
      const AviUtl::PixelYC &ycp = span_ycp[x];
      AviUtl::PixelYC &ycp2      = span_ycp2[x];
      ycp2.y += ycp.y;
      ycp2.cb += ycp.cb;
      ycp2.cr += ycp.cr;
    }
  }

  for (int32_t y = shifty; y < fpip->h - (8 - shifty); ++y) {
    std::span span_ycp(static_cast<AviUtl::PixelYC *>(fpip->ycp_edit) + y * fpip->max_w + fpip->w - (8 - shiftx),
                       fpip->w);
    std::span span_ycp2(wsp + y * fpip->w + fpip->w - (8 - shiftx), fpip->w);
    for (int32_t x = fpip->w - (8 - shiftx); x < fpip->w; ++x) {
      const AviUtl::PixelYC &ycp = span_ycp[x - (fpip->w - (8 - shiftx))];
      AviUtl::PixelYC &ycp2      = span_ycp2[x - (fpip->w - (8 - shiftx))];
      ycp2.y += ycp.y;
      ycp2.cb += ycp.cb;
      ycp2.cr += ycp.cr;
    }
  }
}

//---------------------------------------------------------------------
//  直交変換後量子化関数
//---------------------------------------------------------------------
/*
 * 直交変換をし閾値以下の値に量子化を行います。試してみたところ、n_shiftの値を少なくすると、閾値以下の値を0"の処理では輪郭に粗ができました。閾値を少し高くし量子化を行うことで、少しぼけるものの、n_shiftの値が小さくとも輪郭の粗は目立たなくなりました。
 * 高速化のため、色差は量子化ではなく0にします。見た目は違いが分かりません。
 * なお、圧縮するわけではないので、直流成分はそのままにするべく、iは1からスタートしています。
 */
void transform_quantization(const int32_t &threshold, const int32_t &quantization, const bool &all_quantization,
                            std::array<int32_t, 64> &block, std::array<int32_t, 64> &blockCb,
                            std::array<int32_t, 64> &blockCr) {
  constexpr auto loop_count = eve::expected_cardinal<int32_t>(); // コンパイル時にSIMDのレーン数を取得する

  int32_t first_block   = block[0];
  int32_t first_blockCb = blockCb[0];
  int32_t first_blockCr = blockCr[0];
  for (int32_t i = 0; i < 64; i += loop_count) {
    eve::wide<int32_t, eve::fixed<loop_count>> _block{&block[i]};
    eve::wide<int32_t, eve::fixed<loop_count>> _blockCb{&blockCb[i]};
    eve::wide<int32_t, eve::fixed<loop_count>> _blockCr{&blockCr[i]};
    if (all_quantization) {
      _block   = eve::if_else(eve::abs(_block) < threshold, _block / quantization * quantization, _block);
      _blockCb = eve::if_else(eve::abs(_blockCb) < threshold, _blockCb / quantization * quantization, _blockCb);
      _blockCr = eve::if_else(eve::abs(_blockCr) < threshold, _blockCr / quantization * quantization, _blockCr);
    } else {
      _block   = eve::if_else(eve::abs(_block) < threshold, 0, _block);
      _blockCb = eve::if_else(eve::abs(_blockCb) < threshold, 0, _blockCb);
      _blockCr = eve::if_else(eve::abs(_blockCr) < threshold, 0, _blockCr);
    }
    eve::store(_block, &block[i]);
    eve::store(_blockCb, &blockCb[i]);
    eve::store(_blockCr, &blockCr[i]);
  }
  block[0]   = first_block;
  blockCb[0] = first_blockCb;
  blockCr[0] = first_blockCr;
}

void Loop_dct(AviUtl::FilterProcInfo *fpip, const int32_t &threshold, const int32_t &quantization,
              const bool &all_quantization, const int &thread_id, const int16_t &shiftx, const int16_t &shifty,
              const int16_t &xblock, const int16_t &yblock) {
  std::array<int32_t, 64> block{}, blockCb{}, blockCr{};

  // 8*8のブロックに画像を代入する
  for (int16_t y = 0; y < 8; ++y) {
    std::span span_ycp(static_cast<AviUtl::PixelYC *>(fpip->ycp_edit) + (y + yblock) * fpip->max_w + xblock, 8);
    for (int16_t x = 0; x < 8; x++) {
      const auto &ycp    = span_ycp[x];
      block[x + y * 8]   = ycp.y;
      blockCb[x + y * 8] = ycp.cb;
      blockCr[x + y * 8] = ycp.cr;
    }
  }

  dct_int32(*block.data());   // 離散コサイン変換
  dct_int32(*blockCb.data()); // 離散コサイン変換
  dct_int32(*blockCr.data()); // 離散コサイン変換

  transform_quantization(threshold, quantization, all_quantization, block, blockCb, blockCr);

  idct_int32(*block.data());   // 離散コサイン変換
  idct_int32(*blockCb.data()); // 離散コサイン変換
  idct_int32(*blockCr.data()); // 離散コサイン変換

  // 直接テンポラリ領域に書き込むと画像にゴミが出るため、一度work_spaceにぼかしたブロックを代入する。
  // おそらくマルチスレッドで同時に同じメモリの場所に書き込むことからの不具合と思われます。
  for (int16_t y = 0; y < 8; ++y) {
    std::span span_ycp(work_space + (y + yblock) * fpip->w + xblock + fpip->w * fpip->h * thread_id, 8);
    for (int16_t x = 0; x < 8; x++) {
      auto &ycp = span_ycp[x];
      ycp.y += block[x + y * 8];
      ycp.cb += blockCb[x + y * 8];
      ycp.cr += blockCr[x + y * 8];
    }
  }
}

void Loop_fwht(AviUtl::FilterProcInfo *fpip, const int32_t &threshold, const int32_t &quantization,
               const bool &all_quantization, const int &thread_id, const int16_t &shiftx, const int16_t &shifty,
               const int16_t &xblock, const int16_t &yblock) {
  std::array<int32_t, 64> block{}, blockCb{}, blockCr{};

  // 8*8のブロックに画像を代入する
  for (int16_t y = 0; y < 8; ++y) {
    std::span span_ycp(static_cast<AviUtl::PixelYC *>(fpip->ycp_edit) + (y + yblock) * fpip->max_w + xblock, 8);
    for (int16_t x = 0; x < 8; x++) {
      const auto &ycp    = span_ycp[x];
      block[x + y * 8]   = ycp.y;
      blockCb[x + y * 8] = ycp.cb;
      blockCr[x + y * 8] = ycp.cr;
    }
  }

  fwht(*block.data());   // アダマール変換
  fwht(*blockCb.data()); // アダマール変換
  fwht(*blockCr.data()); // アダマール変換

  transform_quantization(threshold, quantization, all_quantization, block, blockCb, blockCr);

  fwht(*block.data());   // アダマール変換
  fwht(*blockCb.data()); // アダマール変換
  fwht(*blockCr.data()); // アダマール変換

  // 直接テンポラリ領域に書き込むと画像にゴミが出るため、一度work_spaceにぼかしたブロックを代入する。
  // おそらくマルチスレッドで同時に同じメモリの場所に書き込むことからの不具合と思われます。
  for (int16_t y = 0; y < 8; ++y) {
    std::span span_ycp(work_space + (y + yblock) * fpip->w + xblock + fpip->w * fpip->h * thread_id, 8);
    for (int16_t x = 0; x < 8; x++) {
      auto &ycp = span_ycp[x];
      ycp.y += block[x + y * 8];
      ycp.cb += blockCb[x + y * 8];
      ycp.cr += blockCr[x + y * 8];
    }
  }
}