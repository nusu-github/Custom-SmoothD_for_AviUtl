//----------------------------------------------------------------------------------
//  SmoothD for AviUtl ver.0.01
//----------------------------------------------------------------------------------
/*
avisynth�̃t�B���^��SmoothD���ʔ����Ɗ������̂�(�g����̂��́A�����ł���)Aviutl�ɈڐA���܂����B
SmoothD�̃A���S���Y���ł����ASmoothD�ɓ��E����Ă���index.html�ɏ�����Ă���悤�ɁA�ȉ��ŏЉ��Ă�����̂��g�p���Ă��܂��B

http://www.utdallas.edu/~aria/mcl/post/

�摜��0,0���炸�炵���ꏊ���X�^�[�g�n�_�Ƃ��āA8*8�̗��U�R�T�C���ϊ��̗ʎq���ɂ��ڂ����������Ȃ��܂��B���x���X�^�[�g�n�_��ς��Ăڂ������s���A���̉摜���d�˂邱�ƂŁA�u���b�N�m�C�Y���������Ƃ��ł��܂��B�Ȃ�ł�ffdshow��SSP
deblocking�����̃A���S���Y�����g�p���Ă���悤�ł��B �������A����SmoothD for
AviUtl�ł́A���x�ƍ�҂̎�̂��߁A����ϊ��ɃA�_�}�[���ϊ����g�p���Ă��܂��B�Ƃ͌����B�{�Ƃ�SmoothD��MMX�����ꂽxvid��DCT-iDCT���g�p���Ă���̂ł���Ȃ�ɍ����ł��B
�Љ���T�C�g�ŏ����Ă���悤�ɁA���̃A���S���Y���̗��_�́ADCT-iDCT�Ƃ�������ӂꂽ(?)�Z�p���g�p���Ă��邱�Ƃł��B���̃t�B���^���A����ϊ��̊֐��ɑ��̂��̂��g�p���邱�ƂŁA�����������ɂ��邱�Ƃ��\�ƂȂ�͂��ł��B

//�F���̏����͂��Ă��܂���B�������Ԃ̂��Ƃ�����܂����A�ǂ����Ă��F���ς���Ă��܂��̂ŁB16�ŏ��Z����16�ŏ�Z���邾���ŐF�����������ς��̂͂Ȃ��Ȃ̂ł��傤?//
//�F������������悤�ɂ��܂����B12�r�b�g����8�r�b�g�܂Ńf�[�^�ʂ𗎂Ƃ����Ƃ��O�a���Ȃ��͂��Ȃ̂ŁA10�r�b�g�ŏ������邱�Ƃɂ��A8�ŏ��Z���Ă��܂��B����Ȃ�ΐF���قƂ�Ǖς��܂���B
�܂��A�P�x��9�r�b�g�ɗ��Ƃ��悤�ɂ��܂����B�O�a���Ȃ��̂Ȃ�΁A�ł��邾���f�[�^�ʂ��팸���Ȃ��ɉz�������Ƃ͖����ł�����B//
��f�f�[�^��8�r�b�g�ɗ��Ƃ������ɖ߂��܂����B�Ζ؎���DCT-iDCT��8�r�b�g��O��ɍ���Ă������ƁA�܂��A����ϊ��̊֐��ɑ��̂��̂��g���悤�ύX����ۂɂ��A8�r�b�g�ŏ������Ă�����������ł��傤����B�ϐF����͍̂���̉ۑ�ƌ������ƂŁB

�}���`�X���b�h�����ŉ摜�ɃS�~���o�邽�߁A����������f�f�[�^��work_space�Ɉꎞ�ޔ����鏈���ɕύX���܂����B�������̓����ꏊ�ɓ����ɏ����������Ƃ���Ɩ�肪�o��悤�ł��B���R�A�Ȃ̂���?
�R���p�C�������Ƃ����Ă����̂ł͂Ȃ����ƁA������Ɗ�]�͎����Ă����̂ł����B�������B�}���`�X���b�h�ŏ������Ȃ��ꍇ�A�܂��������ʂȏ�����������������ł��B���߂�Ȃ����B

//�ʎq���W�����g�p�����ʎq���ł͂Ȃ��A����ϊ���̒l��臒l�ȉ��Ȃ�0�Ƃ��鏈���ɕύX���܂����B����ɂ��֊s���ڂ��ɂ����Ȃ�܂������A�u���b�N�m�C�Y�������Â炭�Ȃ�܂����B臒l�����Ɨʎq���𕹗p���鏈�����������̂ł����A�u���b�N�m�C�Y�����ɂ͌��ʂ����������̂́A�������Ԃ����Ȃ葝�����Ă��܂�������臒l�����݂̂Ƃ��܂����B//
//�P�x�̂ݗʎq���Ɨ֊s�������s�������ɂ����Ƃ���A����قǒx���Ȃ�Ȃ������̂ŁA�ʎq���̏����𕹗p���邱�Ƃɂ��܂����B�������A�f�t�H���g�͗ʎq���W����"1"�Ȃ̂ŁA���������Ȃ��Ă��܂���B
臒l�Ɨʎq���W������肭�������邱�ƂŁA2DNR�Ƃ��Ă��A�f�u���b�L���O�Ƃ��Ă��A�@�\����͂��ł��B//
�F�X�Ǝ����Ă݂��Ƃ���A"臒l�ȉ��͗ʎq��"�̏������ł����ʓI�Ɣ������܂����B���̂悤�ɏ�������悤�ɂ��Ă��܂��B

//臒l������̒l�𑝉�or���Z���鏈����ǉ����܂����B�֊s����or�ڂ��������Ƃ������ƂŁB�͂����Ďg����@�\�Ȃ̂��͓�ł����A���g���ŉ摜����������Ȃ��(���ꂪ�Ⴆ�����Z�����Z�݂̂̃A�_�}�[���ϊ��ł������Ƃ��Ă�)�A���̋@�\��t�������Ȃ��Ă��܂����̂Ȃ̂ł��B//
���̋@�\�͍폜���܂����B�ǂ��l���Ă��g�p���Ȃ��@�\�Ȃ̂ŁA�������Ԃ��x���Ȃ邾���ł�����B

����! �摜�T�C�Y��8�̔{����p�ł��B����͂��Ă��܂���B



2008/9/13
�E�������J�Bver�ԍ�������B
�EDCT-iDCT�Ɏg���z���short�^�ɂ������߁A�Ζؘa�m�����DCT-iDCT���g����������ނȂ��폜�B
2008/9/14 ver0.002a
�E�}���`�X���b�h�����ŉ摜�ɃS�~���o�邽�߁Awork_space�Ɉꎞ�ޔ����鏈���ɕύX(�������Ԃ��啝��������c�c)�B
�E�P�x�ƐF���̓��������̃r�b�g���𑝂₷�B
�E�F������������悤�ɂ��ăR���p�C���B
2008/9/15 ver0.003
�E����ϊ���ɗʎq���ł͂Ȃ�臒l�ȉ�����鏈���ɕύX�B
�E�֊s����(�ڂ���)�̏�����ǉ��B
2008/9/16 ver0.003b
�E�֊s����(�ڂ���)�͋P�x�݂̂Ƃ����B
�E�P�x�̂ݗʎq���̏����𕹗p����悤�ɂ����B
2008/9/22 ver0.01
�E�P�x��臒l�ȉ���ʎq�����鏈���ɕύX�B
�E�֊s����(�ڂ���)�̏������폜�B
�E�P�x�ƐF���̓��������̃r�b�g����8�r�b�g�ɁB
�E�Ζ؎���DCT-iDCT�̕����𕜊��B
�E�g���b�N�̃f�t�H���g�l��ύX�B
�E�\�[�X�̐����B12�r�b�g��8�r�b�g�̏������}���`�X���b�h�ɁB
*/
// C++ stl
#include <algorithm>
#include <cmath>
#include <span>

// Main Library
#include "dct_int32.hpp"
#include "ht.hpp"
#include "idct_int32.hpp"

#include "SmoothDfA.hpp"

// ���������摜���ꎞ�I�ɑ�����邽�߂̕ϐ�
AviUtl::PixelYC *work_space;

// �O���[�o���ϐ�
// �X���b�h��
int_fast32_t MT = 1;
// �p�����[�^�[�n
smoothdfa_struct smoothdfa;

main_smoothdfa_struct main_smoothdfa;

//---------------------------------------------------------------------
//  �t�B���^�\���̒�`
//---------------------------------------------------------------------
constexpr int TRACK_N = 4; // �g���b�N�o�[�̐�
const char **track_name = new const char *[TRACK_N] { "臒l", "�ʎq���W��", "�i��", "�d��" }; // �g���b�N�o�[�̖��O
int track_default[TRACK_N] = {24, 8, 5, 0};     // �g���b�N�o�[�̏����l
int track_s[TRACK_N]       = {0, 1, 1, 0};      // �g���b�N�o�[�̉����l
int track_e[TRACK_N]       = {255, 255, 6, 64}; // �g���b�N�o�[�̏���l

constexpr int CHECK_N      = 3; // �@�`�F�b�N�{�b�N�X�̐�
const char **check_name = new const char *[CHECK_N] { "DCT", "�ϐF�h�~", "�S�ʎq��" }; // �`�F�b�N�{�b�N�X�̖��O
int check_default[CHECK_N] = {1, 1, 1}; // �`�F�b�N�{�b�N�X�̏����l

#ifdef _DEBUG
const char *filter_name = "SmoothD for AviUtl [DEBUG]"; // �t�B���^��
#else
const char *filter_name = "SmoothD for AviUtl"; // �t�B���^��
#endif

AviUtl::FilterPluginDLL filter{
    .flag          = AviUtl::FilterPlugin::Flag::ExInformation, // �t�B���^�̃t���O
    .name          = filter_name,                               // �t�B���^�̖��O
    .track_n       = TRACK_N,                                   // �g���b�N�o�[�̐�
    .track_name    = track_name,                                // �g���b�N�o�[�̖��O
    .track_default = track_default,                             // �g���b�N�o�[�̏����l
    .track_s       = track_s,                                   // �g���b�N�o�[�̉����l
    .track_e       = track_e,                                   // �g���b�N�o�[�̏���l
    .check_n       = CHECK_N,                                   // �`�F�b�N�{�b�N�X�̐�
    .check_name    = check_name,                                // �`�F�b�N�{�b�N�X�̖��O
    .check_default = check_default,                             // �`�F�b�N�{�b�N�X�̏����l (�l��0��1)
    .func_proc     = func_proc,                                 // �t�B���^�̏����֐�
    .func_init     = func_init,                                 // �������֐�
    .func_exit     = func_exit,                                 // �I���֐�
    .func_update   = func_update,                               // �g���b�N�o�[�̕ύX���ɌĂ΂��֐�
};

//---------------------------------------------------------------------
//  �t�B���^�\���̂̃|�C���^��n���֐�
//---------------------------------------------------------------------
extern "C" auto GetFilterTable(void) -> AviUtl::FilterPluginDLL * { return &filter; }

//---------------------------------------------------------------------
//  Aviutl�������֐�
//---------------------------------------------------------------------
auto func_init(AviUtl::FilterPlugin *fp) -> BOOL {
  // �}���`�X���b�h�����擾
  fp->exfunc->exec_multi_thread_func(get_multi_thread, &MT, nullptr);

  return TRUE;
}

//---------------------------------------------------------------------
//  Aviutl�������֐�
//---------------------------------------------------------------------
auto func_exit(AviUtl::FilterPlugin * /*fp*/) -> BOOL {
  delete[] work_space;
  main_smoothdfa.pictur_width  = 0;
  main_smoothdfa.pictur_height = 0;
  main_smoothdfa.pictur_size   = 0;
  main_smoothdfa.mem_size      = 0;

  return TRUE;
}

//---------------------------------------------------------------------
//  �g���b�N�o�[�̕ύX���ɌĂ΂��֐�
//---------------------------------------------------------------------
auto func_update(AviUtl::FilterPlugin *fp, AviUtl::FilterPluginDLL::UpdateStatus /*status*/) -> BOOL {
  // �g���b�N�ԍ����ԈႦ�Ȃ����߃g���b�N�̒l��ϐ��ɑ�����Ă���
  smoothdfa.threshold = fp->track[0]; // 臒l�B����ϊ���̎��g���̒l�����̒l�ȉ��Ȃ�Ηʎq�����s���B
  smoothdfa.quantization = fp->track[1];                                               // �ʎq���W���B
  smoothdfa.n_shift = static_cast<int_fast16_t>(std::ceil(std::pow(2, fp->track[2]))); // DCT-iDCT�����s�����
  smoothdfa.zero_weight      = fp->track[3]; // �m�C�Y������̉摜�ɍ����錳�̉摜�̊���
  smoothdfa.dct_on           = static_cast<bool>(fp->check[0]); // �Ζ؎����DCT-iDCT���g�p���邩�ǂ���
  smoothdfa.rounding_on      = static_cast<bool>(fp->check[1]); // 8bit�ɂ��鎞�Ɏl�̌ܓ����邩�ǂ���
  smoothdfa.all_quantization = static_cast<bool>(fp->check[2]); // �S�ėʎq�����邩�ǂ���
  return TRUE;
}

//---------------------------------------------------------------------
//  �t�B���^�����֐�
//---------------------------------------------------------------------
auto func_proc(AviUtl::FilterPlugin *fp, AviUtl::FilterProcInfo *fpip) -> BOOL {

  // �摜�̕��ƍ������ς���Ă����珉����
  if (main_smoothdfa.pictur_width != fpip->w || main_smoothdfa.pictur_height != fpip->h) {
    delete[] work_space;

    main_smoothdfa.pictur_width  = fpip->w;
    main_smoothdfa.pictur_height = fpip->h;
    main_smoothdfa.pictur_size   = fpip->w * fpip->h;
    main_smoothdfa.mem_size      = sizeof(AviUtl::PixelYC) * main_smoothdfa.pictur_size * MT;
    // �}���`�X���b�h�������摜�T�C�Y�̃��������m�ۂ���0�Ŗ��߂�B
    work_space = new AviUtl::PixelYC[main_smoothdfa.mem_size]{};
  }

  const int_fast16_t zero_weight = smoothdfa.zero_weight;
  const int_fast16_t n_shift     = smoothdfa.n_shift;

  // ��f�f�[�^��8�r�b�g������B�}���`�X���b�h�ɂ��Ă݂܂����B
  fp->exfunc->exec_multi_thread_func(shift_data, fpip, nullptr);

  /*
  �ȉ��łڂ������摜����Ɨ̈�ɉ��Z���鏈�����J��Ԃ�
  �}���`�X���b�h�Ŋ֐����Ăяo���ĕ���ŏ������Ă��܂��B�Ȃ̂�n_shift�̒l�Ɋ֌W�Ȃ�4�X���b�h�Ȃ�΍Œ�4��A������4�̔{�����������邱�ƂɂȂ�܂��B
  */
  int_fast16_t count = zero_weight;

  for (int_fast16_t i = 0; i < n_shift; i += MT) {
    // Loop(omp_get_thread_num(), omp_get_max_threads(), fpip, &i);
    fp->exfunc->exec_multi_thread_func(Loop, fpip, &i);
    count += MT;
  }

  /*
  ��Ɨ̈�ɉ��Z���ꂽ�l��ϐ� temp
  �ɍ��v�A���摜��zero_weight������Z���āA12�r�b�g�����Ă��� count �ŏ��Z�B
  */

  for (int_fast16_t y = 0; y < fpip->h; y++) {
    const std::span span_ycp(static_cast<AviUtl::PixelYC *>(fpip->ycp_edit) + y * fpip->max_w, fpip->w);
    for (int_fast16_t x = 0; x < fpip->w; x++) {
      AviUtl::PixelYC &ycp = span_ycp[x];
      AviUtl::PixelYC temp = {};

      for (int_fast16_t i = 0; i < MT; i++) {
        AviUtl::PixelYC &wsp = work_space[y * fpip->w + x + main_smoothdfa.pictur_size * i];
        temp.y += wsp.y;
        temp.cb += wsp.cb;
        temp.cr += wsp.cr;
        wsp = {};
      }

      ycp.y  = ((ycp.y * zero_weight + temp.y) << 4) / count;
      ycp.cb = ((ycp.cb * zero_weight + temp.cb) << 4) / count;
      ycp.cr = ((ycp.cr * zero_weight + temp.cr) << 4) / count;
    }
  }

  return TRUE;
}

//---------------------------------------------------------------------
//  ��f�f�[�^8�r�b�g���֐�
//---------------------------------------------------------------------
void shift_data(int thread_id, int thread_num, void *param1, void * /*param2*/) {

  // �}���`�X���b�h�Ή��T���v���t�B���^�̕��@���̂܂܂Ń}���`�X���b�h���������Ă��܂�
  const auto *fpip       = static_cast<AviUtl::FilterProcInfo *>(param1);

  const auto rounding_on = smoothdfa.rounding_on;
  //	�X���b�h���̉摜����������ꏊ���v�Z����
  const int_fast16_t y_start = fpip->h * thread_id / thread_num;
  const int_fast16_t y_end   = fpip->h * (thread_id + 1) / thread_num;

  /*
   *	12�r�b�g�͈�ʓI�łȂ����ƂƁA���Z���Ă����ƖO�a���Ă��܂���������Ȃ��̂ŁA��f�f�[�^�̃r�b�g���𗎂Ƃ��B
   *	�P�x�ƐF���A�ǂ����8�r�b�g�ɁB
   */

  if (rounding_on) {
    for (int_fast16_t y = y_start; y < y_end; y++) {
      const std::span span_ycp(static_cast<AviUtl::PixelYC *>(fpip->ycp_edit) + y * fpip->max_w, fpip->w);
      std::for_each(span_ycp.begin(), span_ycp.end(), [&](AviUtl::PixelYC &ycp) {
        ycp.y  = (ycp.y * 100 / 16 + 45) / 100;
        ycp.cb = ycp.cb > 0 ? (ycp.cb * 100 / 16 + 45) / 100 : (ycp.cb * 100 / 16 - 45) / 100;
        ycp.cr = ycp.cr > 0 ? (ycp.cr * 100 / 16 + 45) / 100 : (ycp.cr * 100 / 16 - 45) / 100;
      });
    }
  } else {
    for (int_fast16_t y = y_start; y < y_end; y++) {
      const std::span span_ycp(static_cast<AviUtl::PixelYC *>(fpip->ycp_edit) + y * fpip->max_w, fpip->w);
      std::for_each(span_ycp.begin(), span_ycp.end(), [&](AviUtl::PixelYC &ycp) {
        ycp.y >>= 4;
        ycp.cb >>= 4;
        ycp.cr >>= 4;
      });
    }
  }
}

//---------------------------------------------------------------------
//  �ڂ��������֐�
//---------------------------------------------------------------------
void Loop(int thread_id, int /*thread_num*/, void *param1, void *param2) {

  auto *fpip                      = static_cast<AviUtl::FilterProcInfo *>(param1);

  const int_fast16_t threshold    = smoothdfa.threshold;
  const int_fast16_t quantization = smoothdfa.quantization;
  const bool all_quantization     = smoothdfa.all_quantization;
  const bool dct_on               = smoothdfa.dct_on;

  // DCT-iDCT�̊J�n�ʒu
  const int_fast16_t shiftx = shift[(*static_cast<int_fast16_t *>(param2) + thread_id) * 2];
  const int_fast16_t shifty = shift[(*static_cast<int_fast16_t *>(param2) + thread_id) * 2 + 1];
  // 0,0���炸�炵��DCT-iDCT�������Ȃ��̂ŉ摜�̒[�܂ŏ����͂��Ȃ�
  const int_fast16_t h = fpip->h - 8;
  const int_fast16_t w = fpip->w - 8;

  for (int_fast16_t yblock = shifty; yblock < h; yblock += 8) {
    for (int_fast16_t xblock = shiftx; xblock < w; xblock += 8) {

      int_fast32_t block[64]{};
      int_fast32_t blockCb[64]{};
      int_fast32_t blockCr[64]{};

      // 8*8�̃u���b�N�ɉ摜��������
      for (int_fast16_t y = 0; y < 8; y++) {
        std::span span_ycp(static_cast<AviUtl::PixelYC *>(fpip->ycp_edit) + (y + yblock) * fpip->max_w + xblock, 8);
        for (int_fast16_t x = 0; x < 8; x++) {
          const auto &ycp    = span_ycp[x];
          block[x + y * 8]   = ycp.y;
          blockCb[x + y * 8] = ycp.cb;
          blockCr[x + y * 8] = ycp.cr;
        }
      }

      /*
       * ����ϊ�����臒l�ȉ��̒l�ɗʎq�����s���܂��B�����Ă݂��Ƃ���An_shift�̒l�����Ȃ�����ƁA"臒l�ȉ��̒l��0"�̏����ł͗֊s�ɑe���ł��܂����B臒l�������������ʎq�����s�����ƂŁA�����ڂ�����̂́An_shift�̒l���������Ƃ��֊s�̑e�͖ڗ����Ȃ��Ȃ�܂����B
       * �������̂��߁A�F���͗ʎq���ł͂Ȃ�0�ɂ��܂��B�����ڂ͈Ⴂ��������܂���B
       * �Ȃ��A���k����킯�ł͂Ȃ��̂ŁA���������͂��̂܂܂ɂ���ׂ��Ai��1����X�^�[�g���Ă��܂��B
       */
      if (dct_on) {
        dct_int32(*block);   // ���U�R�T�C���ϊ�
        dct_int32(*blockCb); // ���U�R�T�C���ϊ�
        dct_int32(*blockCr); // ���U�R�T�C���ϊ�
      } else {
        fwht(*block);   // �A�_�}�[���ϊ�
        fwht(*blockCb); // �A�_�}�[���ϊ�
        fwht(*blockCr); // �A�_�}�[���ϊ�
      }

#pragma omp simd
      for (int_fast16_t i = 1; i < 64; i++) {
        block[i]   = std::abs(block[i]) < threshold ? block[i] / quantization * quantization : block[i];
        blockCb[i] = std::abs(blockCb[i]) < threshold ? all_quantization ? blockCb[i] / quantization * quantization : 0
                                                      : blockCb[i];
        blockCr[i] = std::abs(blockCr[i]) < threshold ? all_quantization ? blockCr[i] / quantization * quantization : 0
                                                      : blockCr[i];
      }

      if (dct_on) {
        idct_int32(*block);   // ���U�R�T�C���ϊ�
        idct_int32(*blockCb); // ���U�R�T�C���ϊ�
        idct_int32(*blockCr); // ���U�R�T�C���ϊ�
      } else {
        fwht(*block);   // �A�_�}�[���ϊ�
        fwht(*blockCb); // �A�_�}�[���ϊ�
        fwht(*blockCr); // �A�_�}�[���ϊ�
      }

      // ���ڃe���|�����̈�ɏ������ނƉ摜�ɃS�~���o�邽�߁A��xwork_space�ɂڂ������u���b�N��������B
      // �����炭�}���`�X���b�h�œ����ɓ����������̏ꏊ�ɏ������ނ��Ƃ���̕s��Ǝv���܂��B
      for (int_fast16_t y = 0; y < 8; y++) {
        std::span span_ycp(work_space + (y + yblock) * fpip->w + xblock + fpip->w * fpip->h * thread_id, 8);
        for (int_fast16_t x = 0; x < 8; x++) {
          auto &ycp = span_ycp[x];
          ycp.y += block[x + y * 8];
          ycp.cb += blockCb[x + y * 8];
          ycp.cr += blockCr[x + y * 8];
        }
      }
    }
  }

  // �摜�l�ӂ̏��������Ă��Ȃ������Ƀf�[�^����
  copy_pix(fpip, work_space + fpip->w * fpip->h * thread_id, shiftx, shifty);
}

//---------------------------------------------------------------------
//  �X���b�h���擾�֐�
//---------------------------------------------------------------------
void get_multi_thread(int thread_id, int thread_num, void *param1, void * /*param2*/) {
  if (thread_id == 0) {
    *static_cast<int_fast16_t *>(param1) = thread_num;
  }
}

//---------------------------------------------------------------------
//  �l�Ӗ����������R�s�[�֐�
//---------------------------------------------------------------------
/*
�����Ȃ�̂Ŋ֐��ɂ��Ċu���B
*/
void copy_pix(AviUtl::FilterProcInfo *fpip, AviUtl::PixelYC *wsp, int shiftx, int shifty) {

  for (int_fast16_t y = 0; y < shifty; y++) {
    std::span span_ycp(static_cast<AviUtl::PixelYC *>(fpip->ycp_edit) + y * fpip->max_w, fpip->w);
    std::span span_ycp2(wsp + y * fpip->w, fpip->w);
    for (int_fast16_t x = 0; x < fpip->w; x++) {
      const AviUtl::PixelYC &ycp = span_ycp[x];
      AviUtl::PixelYC &ycp2      = span_ycp2[x];
      ycp2.y += ycp.y;
      ycp2.cb += ycp.cb;
      ycp2.cr += ycp.cr;
    }
  }

  for (int_fast16_t y = fpip->h - (8 - shifty); y < fpip->h; y++) {
    std::span span_ycp(static_cast<AviUtl::PixelYC *>(fpip->ycp_edit) + y * fpip->max_w, fpip->w);
    std::span span_ycp2(wsp + y * fpip->w, fpip->w);
    for (int_fast16_t x = 0; x < fpip->w; x++) {
      const AviUtl::PixelYC &ycp = span_ycp[x];
      AviUtl::PixelYC &ycp2      = span_ycp2[x];
      ycp2.y += ycp.y;
      ycp2.cb += ycp.cb;
      ycp2.cr += ycp.cr;
    }
  }

  for (int_fast16_t y = shifty; y < fpip->h - (8 - shifty); y++) {
    std::span span_ycp(static_cast<AviUtl::PixelYC *>(fpip->ycp_edit) + y * fpip->max_w, shiftx);
    std::span span_ycp2(wsp + y * fpip->w, shiftx);
    for (int_fast16_t x = 0; x < shiftx; x++) {
      const AviUtl::PixelYC &ycp = span_ycp[x];
      AviUtl::PixelYC &ycp2      = span_ycp2[x];
      ycp2.y += ycp.y;
      ycp2.cb += ycp.cb;
      ycp2.cr += ycp.cr;
    }
  }

  for (int_fast16_t y = shifty; y < fpip->h - (8 - shifty); y++) {
    std::span span_ycp(static_cast<AviUtl::PixelYC *>(fpip->ycp_edit) + y * fpip->max_w + fpip->w - (8 - shiftx),
                       fpip->w);
    std::span span_ycp2(wsp + y * fpip->w + fpip->w - (8 - shiftx), fpip->w);
    for (int_fast16_t x = fpip->w - (8 - shiftx); x < fpip->w; x++) {
      const AviUtl::PixelYC &ycp = span_ycp[x - (fpip->w - (8 - shiftx))];
      AviUtl::PixelYC &ycp2      = span_ycp2[x - (fpip->w - (8 - shiftx))];
      ycp2.y += ycp.y;
      ycp2.cb += ycp.cb;
      ycp2.cr += ycp.cr;
    }
  }
}
