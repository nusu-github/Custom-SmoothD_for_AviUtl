# Custom-SmoothD_for_AviUtl

avisynthのフィルタのSmoothDをAviutlに移植した[SmoothD for AviUtl](https://w.atwiki.jp/aviutl41991/pages/41.html)の改造

## 要求

- Windows10
- Aviutl 1.10

## インストール方法

Custom-SmoothDfA.aufをpluginsフォルダーに入れてください。

## 説明

- threshold 閾値。直交変換後の周波数の値がこの値以下ならば量子化を行う。
- quantization 量子化係数。
- n_shift 品質。1~6の範囲で調整可能。Aviutlのコア数設定で最低値が変わるためどの環境でも1が最低品質ではない。
- zero_weight ノイズ除去後の画像に混ぜる元の画像の割合
- Rounding 8ビットにするときに値を丸めるか。
- All_Quantization すべて量子化するか。
- DCT 茂木氏作の DCT-iDCT を使い、離散コサイン変換を行うか。オフの時は二次元高速ウォルシュ-アダマール変換を行う。

## 改造箇所

- 8スレッド以上のCPUで動作しない問題を修正。
- 処理方法別で分かれていたSmoothDfA_DCT.aufとSmoothDfA.aufを切り替えられるようにして一本化
- 12ビットから8ビットにする際に変色を減らすオプションを追加
- 色差も量子化するオプションを追加

## 参考元

- [AviUtl オリジナルプラグイン公開サイト @ wiki](https://w.atwiki.jp/aviutl41991/)
- [SmoothD for AviUtl](https://w.atwiki.jp/aviutl41991/pages/41.html)
- [まるも製作所](https://www.marumo.ne.jp/)
- [DCT-iDCT AviUtl フィルタプラグイン](https://www.marumo.ne.jp/auf/)

## Licence

- [AviUtl オリジナルプラグイン公開サイト @ wiki](https://w.atwiki.jp/aviutl41991/)から[SmoothD for AviUtl](https://w.atwiki.jp/aviutl41991/pages/41.html)をダウンロードして改造・使用しています。

- DCT 版は[まるも製作所(茂木氏)](https://www.marumo.ne.jp/)の[DCT-iDCT AviUtl フィルタプラグイン](https://www.marumo.ne.jp/auf/)のソースコードを改造・使用しています。

本リポジトリをフォークする場合は以上のライセンスを明記してください。
