//----------------------------------------------------------------------------------
//		SmoothD for AviUtl ver.0.01
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

#include <Windows.h>
#include <cmath>

#include "dct_int32.c"
#include "filter.h"
#include "ht.cpp"

/*
 *	�ȉ��̃O���[�o���ϐ��ł����B
 *	Aviutl�����̃}���`�X���b�h�@�\�͓�̒l�����֐��œn���Ȃ��̂ŁA��ނ𓾂��O���[�o���ϐ��ɂ��܂����B
 *	�����܂ő�������\���̂��g���ׂ��Ȃ̂ł��傤���c�c?
 */
struct smoothdfa_struct {
  int  threshold{};        // 臒l
  int  quantization{};     // �ʎq���W��
  bool rounding_on{};      // 8�r�b�g�̌v�Z���ۂ߂邩
  bool all_quantization{}; // ���ׂėʎq�����邩
  bool dct_on{};           // ���U�R�T�C���ϊ����g����
};

smoothdfa_struct smoothdfa;

// DCT-iDCT���n�߂�J�n�ʒu
//�ŏ��̔��̈ʒu�͖{��SmoothD����B�c��͓K���Ɏ���͂��܂����B
//�g���Ă݂���ۂł͔��ȏ�͂قƂ�Ǖω���������ۂł������A�����ɑS�Ă̈ʒu��ԗ����܂����B
int shift[] = {
    3, 1, 1, 4, 4, 7, 7, 3, 1, 1, 7, 1, 1, 7, 7, 7, 5, 2, 7, 6, 0, 1,
    1, 7, 6, 0, 4, 6, 0, 3, 4, 1, 7, 3, 4, 3, 3, 7, 3, 0, 1, 4, 6, 4,
    0, 5, 5, 1, 5, 7, 0, 2, 7, 5, 0, 7, 4, 2, 2, 1, 7, 7, 1, 2, 4, 5,
    1, 6, 6, 2, 2, 3, 6, 6, 0, 4, 6, 1, 2, 6, 3, 1, 4, 7, 1, 5, 5, 3,
    1, 3, 5, 0, 3, 5, 6, 5, 7, 0, 3, 2, 5, 6, 1, 1, 3, 6, 7, 4, 2, 0,
    5, 4, 6, 3, 2, 7, 3, 3, 4, 0, 2, 4, 4, 4, 7, 2, 0, 0,
};

//���������摜���ꎞ�I�ɑ�����邽�߂̕ϐ�
PIXEL_YC *work_space;

//---------------------------------------------------------------------
//		�t�B���^�\���̒�`
//---------------------------------------------------------------------
#define TRACK_N 4 //	�g���b�N�o�[�̐�
TCHAR *track_name[] = {
    "threshold",
    "quant",
    "n_shift",
    "zero_weight",
}; //	�g���b�N�o�[�̖��O
int track_default[] = {
    24,
    8,
    5,
    0,
}; //	�g���b�N�o�[�̏����l
int track_s[] = {
    0,
    1,
    1,
    0,
}; //	�g���b�N�o�[�̉����l
int track_e[] = {
    255,
    255,
    6,
    64,
}; //	�g���b�N�o�[�̏���l

#define CHECK_N 3 //�@�`�F�b�N�{�b�N�X�̐�
TCHAR *check_name[] = {
    "Rounding",
    "All_Quantization",
    "DCT",
}; //�@�`�F�b�N�{�b�N�X�̖��O
int check_default[] = {
    0,
    0,
    0,
}; //�@�`�F�b�N�{�b�N�X�̏����l (�l��0��1)

FILTER_DLL filter = {
    FILTER_FLAG_EX_INFORMATION, //	�t�B���^�̃t���O
    0, 0,
    "SmoothD for AviUtl", //	�t�B���^�̖��O
    TRACK_N, //	�g���b�N�o�[�̐� (0�Ȃ疼�O�����l����NULL�ł悢)
    track_name,    //	�g���b�N�o�[�̖��O�S�ւ̃|�C���^
    track_default, //	�g���b�N�o�[�̏����l�S�ւ̃|�C���^
    track_s,
    track_e, //	�g���b�N�o�[�̐��l�̉������ (NULL�Ȃ�S��0�`256)
    CHECK_N, //	�`�F�b�N�{�b�N�X�̐� (0�Ȃ疼�O�����l����NULL�ł悢)
    check_name,    //	�`�F�b�N�{�b�N�X�̖��O�S�ւ̃|�C���^
    check_default, //	�`�F�b�N�{�b�N�X�̏����l�S�ւ̃|�C���^
    func_proc, //	�t�B���^�����֐��ւ̃|�C���^ (NULL�Ȃ�Ă΂�܂���)
    nullptr, //	�J�n���ɌĂ΂��֐��ւ̃|�C���^ (NULL�Ȃ�Ă΂�܂���)
    nullptr, //	�I�����ɌĂ΂��֐��ւ̃|�C���^ (NULL�Ȃ�Ă΂�܂���)
    nullptr, //	�ݒ肪�ύX���ꂽ�Ƃ��ɌĂ΂��֐��ւ̃|�C���^
             //(NULL�Ȃ�Ă΂�܂���)
    nullptr, //	�ݒ�E�B���h�E�ɃE�B���h�E���b�Z�[�W���������ɌĂ΂��֐��ւ̃|�C���^
             //(NULL�Ȃ�Ă΂�܂���)
    nullptr,
    nullptr, //	�V�X�e���Ŏg���܂��̂Ŏg�p���Ȃ��ł�������
    nullptr, //  �g���f�[�^�̈�ւ̃|�C���^
             //  (FILTER_FLAG_EX_DATA�������Ă��鎞�ɗL��)
    NULL, //  �g���f�[�^�T�C�Y (FILTER_FLAG_EX_DATA�������Ă��鎞�ɗL��)
    "SmoothD for AviUtl", //  �t�B���^���ւ̃|�C���^
                          //  (FILTER_FLAG_EX_INFORMATION�������Ă��鎞�ɗL��)
    nullptr, //	�Z�[�u���J�n����钼�O�ɌĂ΂��֐��ւ̃|�C���^
             //(NULL�Ȃ�Ă΂�܂���)
    nullptr, //	�Z�[�u���I���������O�ɌĂ΂��֐��ւ̃|�C���^
             //(NULL�Ȃ�Ă΂�܂���)
};

//---------------------------------------------------------------------
//		�t�B���^�\���̂̃|�C���^��n���֐�
//---------------------------------------------------------------------
EXTERN_C FILTER_DLL __declspec(dllexport) * __stdcall GetFilterTable(void) {
  return &filter;
}

//---------------------------------------------------------------------
//		�֐���`
//---------------------------------------------------------------------
void shift_data(int thread_id, int thread_num, void *param1, void *param2);
void Loop(int thread_id, int thread_num, void *param1, void *param2);
auto get_multi_thread(FILTER *fp) -> int;
void copy_pix(FILTER_PROC_INFO *fpip, PIXEL_YC *wsp, int shiftx, int shifty);

//---------------------------------------------------------------------
//		�t�B���^�����֐�
//---------------------------------------------------------------------
auto func_proc(FILTER *fp, FILTER_PROC_INFO *fpip) -> BOOL {
  int       i;
  const int pictur_size = fpip->w * fpip->h;
  PIXEL_YC  temp;

  //�g���b�N�ԍ����ԈႦ�Ȃ����߃g���b�N�̒l��ϐ��ɑ�����Ă���
  smoothdfa.threshold =
      fp->track
          [0]; //臒l�B����ϊ���̎��g���̒l�����̒l�ȉ��Ȃ�Ηʎq�����s���B
  smoothdfa.quantization = fp->track[1]; //�ʎq���W���B
  const int n_shift      = fp->track[2]; // DCT-iDCT�����s�����
  const int zero_weight =
      fp->track[3]; //�m�C�Y������̉摜�ɍ����錳�̉摜�̊���
  smoothdfa.rounding_on      = static_cast<bool>(fp->check[0]);
  smoothdfa.all_quantization = static_cast<bool>(fp->check[1]);
  smoothdfa.dct_on           = static_cast<bool>(fp->check[2]);

  //�}���`�X���b�h�����擾
  const int MT = get_multi_thread(fp);

  //�}���`�X���b�h�������摜�T�C�Y�̃��������m�ۂ���0�Ŗ��߂�B
  work_space = new PIXEL_YC[sizeof(PIXEL_YC) * pictur_size * MT];
#ifdef DEBUG
  ZeroMemory(work_space, (sizeof(PIXEL_YC) * pictur_size * MT));
#endif

  //"0"�Ŗ��߂Ȃ��Ƃ����Ȃ����삵�܂��B�f�t�H���g�R���X�g���N�^��"0"�Ƃ������ƁA�Ȃ̂ł��傤��?
  //���̏��������Ȃ��Ƃ��Ȃ菈�����x���オ�邽�߁A�Ƃ肠�����s��Ȃ����ƂɁB

  //��f�f�[�^��8�r�b�g������B�}���`�X���b�h�ɂ��Ă݂܂����B
  fp->exfunc->exec_multi_thread_func(shift_data, static_cast<void *>(fpip),
                                     nullptr);

  /*
  �ȉ��łڂ������摜����Ɨ̈�ɉ��Z���鏈�����J��Ԃ�
  �}���`�X���b�h�Ŋ֐����Ăяo���ĕ���ŏ������Ă��܂��B�Ȃ̂�n_shift�̒l�Ɋ֌W�Ȃ�4�X���b�h�Ȃ�΍Œ�4��A������4�̔{�����������邱�ƂɂȂ�܂��B
  */
  int count = zero_weight;

  for (i = 0; i < std::pow(2, n_shift); i += MT) {
    fp->exfunc->exec_multi_thread_func(Loop, static_cast<void *>(fpip),
                                       static_cast<void *>(&i));

    count += MT;
  }

  /*
  ��Ɨ̈�ɉ��Z���ꂽ�l��ϐ� temp
  �ɍ��v�A���摜��zero_weight������Z���āA12�r�b�g�����Ă��� count �ŏ��Z�B
  */
  PIXEL_YC *wsp = work_space;
  for (int y = 0; y < fpip->h; y++) {
    PIXEL_YC *ycp = fpip->ycp_edit + y * fpip->max_w;
    for (int x = 0; x < fpip->w; x++) {
      PIXEL_YC *ycp2 = wsp;
      ZeroMemory(&temp,
                 6); // short�^���O�Ȃ̂ŁA6�o�C�g���ߑł���temp�̒l���N���A
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
//		��f�f�[�^8�r�b�g���֐�
//---------------------------------------------------------------------
/*
�}���`�X���b�h�ŏ�������ׂ��֐��ɂ��Ă݂܂���
*/
void shift_data(int thread_id, int thread_num, void *param1,
                void * /*param2*/) {
  //�}���`�X���b�h�Ή��T���v���t�B���^�̕��@���̂܂܂Ń}���`�X���b�h���������Ă��܂�
  const FILTER_PROC_INFO *fpip = static_cast<FILTER_PROC_INFO *>(param1);

  const bool rounding_on = smoothdfa.rounding_on;
  //	�X���b�h���̉摜����������ꏊ���v�Z����
  int y_start = fpip->h * thread_id / thread_num;
  int y_end   = fpip->h * (thread_id + 1) / thread_num;

  /*
   *	12�r�b�g�͈�ʓI�łȂ����ƂƁA���Z���Ă����ƖO�a���Ă��܂���������Ȃ��̂ŁA��f�f�[�^�̃r�b�g���𗎂Ƃ��B
   *	�P�x�ƐF���A�ǂ����8�r�b�g�ɁB
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
//		�ڂ��������֐�
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
   *	�Ζ؎���DCT-iDCT��int�^�ŏ�������悤�ɂȂ��Ă���̂ŁA�錾����ʂ�ɁB
   *	�Ζ؎��̃\�[�X�����������Ă��ǂ��̂ł����A���̕��̃\�[�X�Ɏ������̂��ǂ����ȁA�ƁB
   */
  int block[64];
  int blockCb[64];
  int blockCr[64];

  // DCT-iDCT�̊J�n�ʒu
  const int shiftx = shift[(*static_cast<int *>(param2) + thread_id) * 2];
  const int shifty = shift[(*static_cast<int *>(param2) + thread_id) * 2 + 1];
  // 0,0���炸�炵��DCT-iDCT�������Ȃ��̂ŉ摜�̒[�܂ŏ����͂��Ȃ�
  const int h = fpip->h - 8;
  const int w = fpip->w - 8;

  for (int yblock = shifty; yblock < h; yblock += 8) {
    for (int xblock = shiftx; xblock < w; xblock += 8) {
      // 8*8�̃u���b�N�ɉ摜��������
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
       *	����ϊ�����臒l�ȉ��̒l�ɗʎq�����s���܂��B�����Ă݂��Ƃ���An_shift�̒l�����Ȃ�����ƁA"臒l�ȉ��̒l��0"�̏����ł͗֊s�ɑe���ł��܂����B臒l�������������ʎq�����s�����ƂŁA�����ڂ�����̂́An_shift�̒l���������Ƃ��֊s�̑e�͖ڗ����Ȃ��Ȃ�܂����B
       *	�������̂��߁A�F���͗ʎq���ł͂Ȃ�0�ɂ��܂��B�����ڂ͈Ⴂ��������܂���B
       *	�Ȃ��A���k����킯�ł͂Ȃ��̂ŁA���������͂��̂܂܂ɂ���ׂ��Ai��1����X�^�[�g���Ă��܂��B
       */
      if (dct_on) {
        dct_int32(block);   //���U�R�T�C���ϊ�
        dct_int32(blockCb); //���U�R�T�C���ϊ�
        dct_int32(blockCr); //���U�R�T�C���ϊ�
      } else {
        fwht(block);   //�A�_�}�[���ϊ�
        fwht(blockCb); //�A�_�}�[���ϊ�
        fwht(blockCr); //�A�_�}�[���ϊ�
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
        idct_int32(block);   //���U�R�T�C���ϊ�
        idct_int32(blockCb); //���U�R�T�C���ϊ�
        idct_int32(blockCr); //���U�R�T�C���ϊ�
      } else {
        fwht(block);   //�A�_�}�[���ϊ�
        fwht(blockCb); //�A�_�}�[���ϊ�
        fwht(blockCr); //�A�_�}�[���ϊ�
      }
      //���ڃe���|�����̈�ɏ������ނƉ摜�ɃS�~���o�邽�߁A��xwork_space�ɂڂ������u���b�N��������B
      //�����炭�}���`�X���b�h�œ����ɓ����������̏ꏊ�ɏ������ނ��Ƃ���̕s��Ǝv���܂��B
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

  //�摜�l�ӂ̏��������Ă��Ȃ������Ƀf�[�^����
  copy_pix(fpip, work_space + fpip->w * fpip->h * thread_id, shiftx, shifty);
}

//---------------------------------------------------------------------
//		�X���b�h���擾�֐�
//---------------------------------------------------------------------
/*
�}���`�X���b�h���̒l���V�X�e���C���t�H���[�V�����\���̂ɓ����Ă��Ȃ��̂ŁA�֐�������Ď��O�Ŏ擾���Ă��܂��B
�����ƊȒP�ɂł�����@������̂ł��傤��?
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
//		�l�Ӗ����������R�s�[�֐�
//---------------------------------------------------------------------
/*
�����Ȃ�̂Ŋ֐��ɂ��Ċu���B
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
