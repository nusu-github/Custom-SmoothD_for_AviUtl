#pragma once

// Aviutl SDK
#include <aviutl/FilterPlugin.hpp>
#include <aviutl/FilterProcInfo.hpp>

// DCT-iDCT���n�߂�J�n�ʒu
// �ŏ��̔��̈ʒu�͖{��SmoothD����B�c��͓K���Ɏ���͂��܂����B
// �g���Ă݂���ۂł͔��ȏ�͂قƂ�Ǖω���������ۂł������A�����ɑS�Ă̈ʒu��ԗ����܂����B
constexpr int_fast8_t shift[] = {
    3, 1, 1, 4, 4, 7, 7, 3, 1, 1, 7, 1, 1, 7, 7, 7, 5, 2, 7, 6, 0, 1, 1, 7, 6, 0, 4, 6, 0, 3, 4, 1,
    7, 3, 4, 3, 3, 7, 3, 0, 1, 4, 6, 4, 0, 5, 5, 1, 5, 7, 0, 2, 7, 5, 0, 7, 4, 2, 2, 1, 7, 7, 1, 2,
    4, 5, 1, 6, 6, 2, 2, 3, 6, 6, 0, 4, 6, 1, 2, 6, 3, 1, 4, 7, 1, 5, 5, 3, 1, 3, 5, 0, 3, 5, 6, 5,
    7, 0, 3, 2, 5, 6, 1, 1, 3, 6, 7, 4, 2, 0, 5, 4, 6, 3, 2, 7, 3, 3, 4, 0, 2, 4, 4, 4, 7, 2, 0, 0,
};

struct smoothdfa_struct {
  int_fast16_t threshold;    // 臒l
  int_fast16_t quantization; // �ʎq���W��
  int_fast16_t n_shift;      // DCT-iDCT�����s�����
  int_fast16_t zero_weight;  // �m�C�Y������̉摜�ɍ����錳�̉摜�̊���
  bool rounding_on;          // 8�r�b�g�̌v�Z���ۂ߂邩
  bool all_quantization;     // ���ׂėʎq�����邩
  bool dct_on;               // ���U�R�T�C���ϊ����g����
};

struct main_smoothdfa_struct {
  int_fast16_t pictur_width;  // �摜�̕�
  int_fast16_t pictur_height; // �摜�̍���
  int_fast32_t pictur_size;   // �摜�̃T�C�Y
  int_fast32_t mem_size;      // �������̃T�C�Y
};

auto func_proc(AviUtl::FilterPlugin *fp, AviUtl::FilterProcInfo *fpip) -> BOOL;
auto func_init(AviUtl::FilterPlugin *fp) -> BOOL;
auto func_exit(AviUtl::FilterPlugin *fp) -> BOOL;
auto func_update(AviUtl::FilterPlugin *fp, AviUtl::FilterPluginDLL::UpdateStatus status) -> BOOL;

void get_multi_thread(int thread_id, int thread_num, void *param1, void * /*param2*/);

void shift_data(int thread_id, int thread_num, void *param1, void * /*param2*/);
void Loop(int thread_id, int thread_num, void *param1, void *param2);
void copy_pix(AviUtl::FilterProcInfo *fpip, AviUtl::PixelYC *wsp, int shiftx, int shifty);