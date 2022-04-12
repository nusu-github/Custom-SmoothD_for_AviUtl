//----------------------------------------------------------------------------------
//		SmoothD for AviUtl ver.0.01
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

#include <Windows.h>
#include <cmath>

#include "dct_int32.c"
#include "filter.h"
#include "ht.cpp"

/*
 *	以下のグローバル変数ですが。
 *	Aviutl内臓のマルチスレッド機能は二つの値しか関数で渡せないので、やむを得ずグローバル変数にしました。
 *	ここまで増えたら構造体を使うべきなのでしょうか……?
 */
struct smoothdfa_struct {
  int  threshold{};        // 閾値
  int  quantization{};     // 量子化係数
  bool rounding_on{};      // 8ビットの計算を丸めるか
  bool all_quantization{}; // すべて量子化するか
  bool dct_on{};           // 離散コサイン変換を使うか
};

smoothdfa_struct smoothdfa;

// DCT-iDCTを始める開始位置
//最初の八個の位置は本家SmoothDから。残りは適当に手入力しました。
//使ってみた印象では八個以上はほとんど変化が無い印象でしたが、試しに全ての位置を網羅しました。
int shift[] = {
    3, 1, 1, 4, 4, 7, 7, 3, 1, 1, 7, 1, 1, 7, 7, 7, 5, 2, 7, 6, 0, 1,
    1, 7, 6, 0, 4, 6, 0, 3, 4, 1, 7, 3, 4, 3, 3, 7, 3, 0, 1, 4, 6, 4,
    0, 5, 5, 1, 5, 7, 0, 2, 7, 5, 0, 7, 4, 2, 2, 1, 7, 7, 1, 2, 4, 5,
    1, 6, 6, 2, 2, 3, 6, 6, 0, 4, 6, 1, 2, 6, 3, 1, 4, 7, 1, 5, 5, 3,
    1, 3, 5, 0, 3, 5, 6, 5, 7, 0, 3, 2, 5, 6, 1, 1, 3, 6, 7, 4, 2, 0,
    5, 4, 6, 3, 2, 7, 3, 3, 4, 0, 2, 4, 4, 4, 7, 2, 0, 0,
};

//処理した画像を一時的に代入するための変数
PIXEL_YC *work_space;

//---------------------------------------------------------------------
//		フィルタ構造体定義
//---------------------------------------------------------------------
#define TRACK_N 4 //	トラックバーの数
TCHAR *track_name[] = {
    "threshold",
    "quant",
    "n_shift",
    "zero_weight",
}; //	トラックバーの名前
int track_default[] = {
    24,
    8,
    5,
    0,
}; //	トラックバーの初期値
int track_s[] = {
    0,
    1,
    1,
    0,
}; //	トラックバーの下限値
int track_e[] = {
    255,
    255,
    6,
    64,
}; //	トラックバーの上限値

#define CHECK_N 3 //　チェックボックスの数
TCHAR *check_name[] = {
    "Rounding",
    "All_Quantization",
    "DCT",
}; //　チェックボックスの名前
int check_default[] = {
    0,
    0,
    0,
}; //　チェックボックスの初期値 (値は0か1)

FILTER_DLL filter = {
    FILTER_FLAG_EX_INFORMATION, //	フィルタのフラグ
    0, 0,
    "SmoothD for AviUtl", //	フィルタの名前
    TRACK_N, //	トラックバーの数 (0なら名前初期値等もNULLでよい)
    track_name,    //	トラックバーの名前郡へのポインタ
    track_default, //	トラックバーの初期値郡へのポインタ
    track_s,
    track_e, //	トラックバーの数値の下限上限 (NULLなら全て0～256)
    CHECK_N, //	チェックボックスの数 (0なら名前初期値等もNULLでよい)
    check_name,    //	チェックボックスの名前郡へのポインタ
    check_default, //	チェックボックスの初期値郡へのポインタ
    func_proc, //	フィルタ処理関数へのポインタ (NULLなら呼ばれません)
    nullptr, //	開始時に呼ばれる関数へのポインタ (NULLなら呼ばれません)
    nullptr, //	終了時に呼ばれる関数へのポインタ (NULLなら呼ばれません)
    nullptr, //	設定が変更されたときに呼ばれる関数へのポインタ
             //(NULLなら呼ばれません)
    nullptr, //	設定ウィンドウにウィンドウメッセージが来た時に呼ばれる関数へのポインタ
             //(NULLなら呼ばれません)
    nullptr,
    nullptr, //	システムで使いますので使用しないでください
    nullptr, //  拡張データ領域へのポインタ
             //  (FILTER_FLAG_EX_DATAが立っている時に有効)
    NULL, //  拡張データサイズ (FILTER_FLAG_EX_DATAが立っている時に有効)
    "SmoothD for AviUtl", //  フィルタ情報へのポインタ
                          //  (FILTER_FLAG_EX_INFORMATIONが立っている時に有効)
    nullptr, //	セーブが開始される直前に呼ばれる関数へのポインタ
             //(NULLなら呼ばれません)
    nullptr, //	セーブが終了した直前に呼ばれる関数へのポインタ
             //(NULLなら呼ばれません)
};

//---------------------------------------------------------------------
//		フィルタ構造体のポインタを渡す関数
//---------------------------------------------------------------------
EXTERN_C FILTER_DLL __declspec(dllexport) * __stdcall GetFilterTable(void) {
  return &filter;
}

//---------------------------------------------------------------------
//		関数定義
//---------------------------------------------------------------------
void shift_data(int thread_id, int thread_num, void *param1, void *param2);
void Loop(int thread_id, int thread_num, void *param1, void *param2);
auto get_multi_thread(FILTER *fp) -> int;
void copy_pix(FILTER_PROC_INFO *fpip, PIXEL_YC *wsp, int shiftx, int shifty);

//---------------------------------------------------------------------
//		フィルタ処理関数
//---------------------------------------------------------------------
auto func_proc(FILTER *fp, FILTER_PROC_INFO *fpip) -> BOOL {
  int       i;
  const int pictur_size = fpip->w * fpip->h;
  PIXEL_YC  temp;

  //トラック番号を間違えないためトラックの値を変数に代入しておく
  smoothdfa.threshold =
      fp->track
          [0]; //閾値。直交変換後の周波数の値がこの値以下ならば量子化を行う。
  smoothdfa.quantization = fp->track[1]; //量子化係数。
  const int n_shift      = fp->track[2]; // DCT-iDCTを実行する回数
  const int zero_weight =
      fp->track[3]; //ノイズ除去後の画像に混ぜる元の画像の割合
  smoothdfa.rounding_on      = static_cast<bool>(fp->check[0]);
  smoothdfa.all_quantization = static_cast<bool>(fp->check[1]);
  smoothdfa.dct_on           = static_cast<bool>(fp->check[2]);

  //マルチスレッド数を取得
  const int MT = get_multi_thread(fp);

  //マルチスレッド数だけ画像サイズのメモリを確保して0で埋める。
  work_space = new PIXEL_YC[sizeof(PIXEL_YC) * pictur_size * MT];
#ifdef DEBUG
  ZeroMemory(work_space, (sizeof(PIXEL_YC) * pictur_size * MT));
#endif

  //"0"で埋めなくとも問題なく動作します。デフォルトコンストラクタが"0"ということ、なのでしょうか?
  //この処理をしないとかなり処理速度が上がるため、とりあえず行わないことに。

  //画素データを8ビット化する。マルチスレッドにしてみました。
  fp->exfunc->exec_multi_thread_func(shift_data, static_cast<void *>(fpip),
                                     nullptr);

  /*
  以下でぼかした画像を作業領域に加算する処理を繰り返す
  マルチスレッドで関数を呼び出して並列で処理しています。なのでn_shiftの値に関係なく4スレッドならば最低4回、そして4の倍数ずつ処理することになります。
  */
  int count = zero_weight;

  for (i = 0; i < std::pow(2, n_shift); i += MT) {
    fp->exfunc->exec_multi_thread_func(Loop, static_cast<void *>(fpip),
                                       static_cast<void *>(&i));

    count += MT;
  }

  /*
  作業領域に加算された値を変数 temp
  に合計、元画像をzero_weightだけ乗算して、12ビット化してから count で除算。
  */
  PIXEL_YC *wsp = work_space;
  for (int y = 0; y < fpip->h; y++) {
    PIXEL_YC *ycp = fpip->ycp_edit + y * fpip->max_w;
    for (int x = 0; x < fpip->w; x++) {
      PIXEL_YC *ycp2 = wsp;
      ZeroMemory(&temp,
                 6); // short型が三つなので、6バイト決め打ちでtempの値をクリア
      for (i = 0; i < MT; i++) {
        temp.y += ycp2->y;
        temp.cb += ycp2->cb;
        temp.cr += ycp2->cr;
        ycp2 += pictur_size;
      }

      ycp->y =
          static_cast<short>(((ycp->y * zero_weight + temp.y) << 4) / count);
      ycp->cb =
          static_cast<short>(((ycp->cb * zero_weight + temp.cb) << 4) / count);
      ycp->cr =
          static_cast<short>(((ycp->cr * zero_weight + temp.cr) << 4) / count);
      ycp++;
      wsp++;
    }
  }

  delete[] work_space;

  return TRUE;
}

//---------------------------------------------------------------------
//		画素データ8ビット化関数
//---------------------------------------------------------------------
/*
マルチスレッドで処理するべく関数にしてみました
*/
void shift_data(int thread_id, int thread_num, void *param1,
                void * /*param2*/) {
  //マルチスレッド対応サンプルフィルタの方法そのままでマルチスレッド処理をしています
  const FILTER_PROC_INFO *fpip = static_cast<FILTER_PROC_INFO *>(param1);

  const bool rounding_on = smoothdfa.rounding_on;
  //	スレッド毎の画像を処理する場所を計算する
  int y_start = fpip->h * thread_id / thread_num;
  int y_end   = fpip->h * (thread_id + 1) / thread_num;

  /*
   *	12ビットは一般的でないことと、加算していくと飽和してしまうかもしれないので、画素データのビット数を落とす。
   *	輝度と色差、どちらも8ビットに。
   */
  for (int y = y_start; y < y_end; y++) {
    PIXEL_YC *ycp = fpip->ycp_edit + y * fpip->max_w;
    for (int x = 0; x < fpip->w; x++) {
      if (rounding_on) {
        ycp->y  = (ycp->y * 100 / 16 + 45) / 100;
        ycp->cb = (ycp->cb > 0) ? (ycp->cb * 100 / 16 + 45) / 100
                                : (ycp->cb * 100 / 16 - 45) / 100;
        ycp->cr = (ycp->cr > 0) ? (ycp->cr * 100 / 16 + 45) / 100
                                : (ycp->cr * 100 / 16 - 45) / 100;
      } else {
        ycp->y >>= 4;
        ycp->cb >>= 4;
        ycp->cr >>= 4;
      }
      ycp++;
    }
  }
}

//---------------------------------------------------------------------
//		ぼかし処理関数
//---------------------------------------------------------------------
void Loop(int thread_id, int /*thread_num*/, void *param1, void *param2) {
  auto *fpip = static_cast<FILTER_PROC_INFO *>(param1);

  const int  threshold        = smoothdfa.threshold;
  const int  quantization     = smoothdfa.quantization;
  const bool all_quantization = smoothdfa.all_quantization;
  const bool dct_on           = smoothdfa.dct_on;

  int       x;
  int       y;
  PIXEL_YC *ycp;

  /*
   *	茂木氏のDCT-iDCTがint型で処理するようになっているので、宣言が二通りに。
   *	茂木氏のソースを書き換えても良いのですが、他の方のソースに手を入れるのもどうかな、と。
   */
  int block[64];
  int blockCb[64];
  int blockCr[64];

  // DCT-iDCTの開始位置
  const int shiftx = shift[(*static_cast<int *>(param2) + thread_id) * 2];
  const int shifty = shift[(*static_cast<int *>(param2) + thread_id) * 2 + 1];
  // 0,0からずらしてDCT-iDCTをおこなうので画像の端まで処理はしない
  const int h = fpip->h - 8;
  const int w = fpip->w - 8;

  for (int yblock = shifty; yblock < h; yblock += 8) {
    for (int xblock = shiftx; xblock < w; xblock += 8) {
      // 8*8のブロックに画像を代入する
      for (y = 0; y < 8; y++) {
        ycp = fpip->ycp_edit + (y + yblock) * fpip->max_w + xblock;
        for (x = 0; x < 8; x++) {
          block[x + y * 8]   = ycp->y;
          blockCb[x + y * 8] = ycp->cb;
          blockCr[x + y * 8] = ycp->cr;
          ycp++;
        }
      }

      /*
       *	直交変換をし閾値以下の値に量子化を行います。試してみたところ、n_shiftの値を少なくすると、"閾値以下の値を0"の処理では輪郭に粗ができました。閾値を少し高くし量子化を行うことで、少しぼけるものの、n_shiftの値が小さくとも輪郭の粗は目立たなくなりました。
       *	高速化のため、色差は量子化ではなく0にします。見た目は違いが分かりません。
       *	なお、圧縮するわけではないので、直流成分はそのままにするべく、iは1からスタートしています。
       */
      if (dct_on) {
        dct_int32(block);   //離散コサイン変換
        dct_int32(blockCb); //離散コサイン変換
        dct_int32(blockCr); //離散コサイン変換
      } else {
        fwht(block);   //アダマール変換
        fwht(blockCb); //アダマール変換
        fwht(blockCr); //アダマール変換
      }
      for (int i = 1; i < 64; i++) {
        if (all_quantization) {
          if (abs(block[i]) < threshold)
            block[i] = (block[i] / quantization) * quantization;
          if (abs(blockCb[i]) < threshold)
            blockCb[i] = (blockCb[i] / quantization) * quantization;
          if (abs(blockCr[i]) < threshold)
            blockCr[i] = (blockCr[i] / quantization) * quantization;
        } else {
          if (abs(block[i]) < threshold)
            block[i] = block[i] / quantization * quantization;
          if (abs(blockCb[i]) < threshold)
            blockCb[i] = 0;
          if (abs(blockCr[i]) < threshold)
            blockCr[i] = 0;
        }
      }
      if (dct_on) {
        idct_int32(block);   //離散コサイン変換
        idct_int32(blockCb); //離散コサイン変換
        idct_int32(blockCr); //離散コサイン変換
      } else {
        fwht(block);   //アダマール変換
        fwht(blockCb); //アダマール変換
        fwht(blockCr); //アダマール変換
      }
      //直接テンポラリ領域に書き込むと画像にゴミが出るため、一度work_spaceにぼかしたブロックを代入する。
      //おそらくマルチスレッドで同時に同じメモリの場所に書き込むことからの不具合と思われます。
      for (y = 0; y < 8; y++) {
        ycp = work_space + (y + yblock) * fpip->w + xblock +
              fpip->w * fpip->h * thread_id;
        for (x = 0; x < 8; x++) {
          ycp->y += static_cast<short>(block[x + y * 8]);
          ycp->cb += static_cast<short>(blockCb[x + y * 8]);
          ycp->cr += static_cast<short>(blockCr[x + y * 8]);
          ycp++;
        }
      }
    }
  }

  //画像四辺の処理をしていない部分にデータを代入
  copy_pix(fpip, work_space + fpip->w * fpip->h * thread_id, shiftx, shifty);
}

//---------------------------------------------------------------------
//		スレッド数取得関数
//---------------------------------------------------------------------
/*
マルチスレッド数の値がシステムインフォメーション構造体に入っていないので、関数を作って自前で取得しています。
もっと簡単にできる方法があるのでしょうか?
*/
void multi_thread_count(int /*thread_id*/, int thread_num, void *param1,
                        void * /*param2*/) {
  auto *MT_temp = static_cast<int *>(param1);
  MT_temp[0]    = thread_num;
}

auto get_multi_thread(FILTER *fp) -> int {
  int MT_temp[1] = {0};

  fp->exfunc->exec_multi_thread_func(multi_thread_count, (void *)MT_temp,
                                     static_cast<void *>(nullptr));

  int temp = MT_temp[0];

  return temp;
}

//---------------------------------------------------------------------
//		四辺未処理部分コピー関数
//---------------------------------------------------------------------
/*
長くなるので関数にして隔離。
*/
void copy_pix(FILTER_PROC_INFO *fpip, PIXEL_YC *wsp, int shiftx, int shifty) {
  int       x;
  int       y;
  PIXEL_YC *ycp;
  PIXEL_YC *ycp2;

  for (y = 0; y < shifty; y++) {
    ycp  = fpip->ycp_edit + y * fpip->max_w;
    ycp2 = wsp + y * fpip->w;
    for (x = 0; x < fpip->w; x++) {
      ycp2->y += ycp->y;
      ycp2->cb += ycp->cb;
      ycp2->cr += ycp->cr;
      ycp++;
      ycp2++;
    }
  }

  for (y = fpip->h - (8 - shifty); y < fpip->h; y++) {
    ycp  = fpip->ycp_edit + y * fpip->max_w;
    ycp2 = wsp + y * fpip->w;
    for (x = 0; x < fpip->w; x++) {
      ycp2->y += ycp->y;
      ycp2->cb += ycp->cb;
      ycp2->cr += ycp->cr;
      ycp++;
      ycp2++;
    }
  }

  for (y = shifty; y < fpip->h - (8 - shifty); y++) {
    ycp  = fpip->ycp_edit + y * fpip->max_w;
    ycp2 = wsp + y * fpip->w;
    for (x = 0; x < shiftx; x++) {
      ycp2->y += ycp->y;
      ycp2->cb += ycp->cb;
      ycp2->cr += ycp->cr;
      ycp++;
      ycp2++;
    }
  }

  for (y = shifty; y < fpip->h - (8 - shifty); y++) {
    ycp  = fpip->ycp_edit + y * fpip->max_w + fpip->w - (8 - shiftx);
    ycp2 = wsp + y * fpip->w + fpip->w - (8 - shiftx);
    for (x = fpip->w - (8 - shiftx); x < fpip->w; x++) {
      ycp2->y += ycp->y;
      ycp2->cb += ycp->cb;
      ycp2->cr += ycp->cr;
      ycp++;
      ycp2++;
    }
  }
}
