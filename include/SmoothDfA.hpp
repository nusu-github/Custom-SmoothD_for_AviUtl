#pragma once

/*
 * �ȉ��̃O���[�o���ϐ��ł����B
 * Aviutl�����̃}���`�X���b�h�@�\�͓�̒l�����֐��œn���Ȃ��̂ŁA��ނ𓾂��O���[�o���ϐ��ɂ��܂����B
 * �����܂ő�������\���̂��g���ׂ��Ȃ̂ł��傤���c�c?
 */
struct smoothdfa_struct {
  int threshold{};         // 臒l
  int quantization{};      // �ʎq���W��
  int n_shift{};           // DCT-iDCT�����s�����
  int zero_weight{};       // �m�C�Y������̉摜�ɍ����錳�̉摜�̊���
  bool rounding_on{};      // 8�r�b�g�̌v�Z���ۂ߂邩
  bool all_quantization{}; // ���ׂėʎq�����邩
  bool dct_on{};           // ���U�R�T�C���ϊ����g����
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