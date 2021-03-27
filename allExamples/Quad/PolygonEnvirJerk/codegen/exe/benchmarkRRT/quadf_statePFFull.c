/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * quadf_statePFFull.c
 *
 * Code generation for function 'quadf_statePFFull'
 *
 */

/* Include files */
#include <math.h>
#include "rt_nonfinite.h"
#include "benchmarkRRT.h"
#include "quadf_statePFFull.h"
#include "quadf_cost.h"
#include "benchmarkRRT_rtwutil.h"

/* Function Definitions */
void quadf_statePFFull(const creal_T t, const creal_T t_s, double x01, double
  x02, double x03, double x04, double x05, double x06, double x07, double x08,
  double x09, double x11, double x12, double x13, double states[9])
{
  double y_re;
  double y_im;
  double x_re;
  double x_im;
  double b_y_re;
  double b_y_im;
  double c_y_re;
  double c_y_im;
  double r;
  double d_y_re;
  double d_y_im;
  double e_y_re;
  double e_y_im;
  double f_y_re;
  double f_y_im;
  double g_y_re;
  double g_y_im;
  double h_y_re;
  double h_y_im;
  double i_y_re;
  double i_y_im;
  double j_y_re;
  double j_y_im;
  double k_y_re;
  double k_y_im;
  double l_y_re;
  double l_y_im;
  double m_y_re;
  double m_y_im;
  double n_y_re;
  double n_y_im;
  double o_y_re;
  double o_y_im;
  double p_y_re;
  double p_y_im;
  double q_y_re;
  double q_y_im;
  double r_y_re;
  double r_y_im;
  double s_y_re;
  double s_y_im;
  double t_y_re;
  double t_y_im;
  double u_y_re;
  double u_y_im;
  double v_y_re;
  double v_y_im;
  double w_y_re;
  double w_y_im;
  double x_y_re;
  double x_y_im;
  double y_y_re;
  double y_y_im;
  double ab_y_re;
  double ab_y_im;
  double bb_y_re;
  double bb_y_im;
  double cb_y_re;
  double cb_y_im;
  double db_y_re;
  double db_y_im;
  double eb_y_re;
  double eb_y_im;
  double fb_y_re;
  double fb_y_im;
  double gb_y_re;
  double gb_y_im;
  double hb_y_re;
  double hb_y_im;
  double ib_y_re;
  double ib_y_im;
  double jb_y_re;
  double jb_y_im;
  double kb_y_re;
  double kb_y_im;
  double lb_y_re;
  double lb_y_im;
  double mb_y_re;
  double mb_y_im;
  double nb_y_re;
  double nb_y_im;
  double ob_y_re;
  double ob_y_im;
  double pb_y_re;
  double pb_y_im;
  double qb_y_re;
  double qb_y_im;
  double rb_y_re;
  double rb_y_im;
  double sb_y_re;
  double sb_y_im;
  double tb_y_re;
  double tb_y_im;
  double ub_y_re;
  double ub_y_im;
  double vb_y_re;
  double vb_y_im;
  double wb_y_re;
  double wb_y_im;
  double xb_y_re;
  double xb_y_im;
  double yb_y_re;
  double yb_y_im;
  double ac_y_re;
  double ac_y_im;
  double bc_y_re;
  double bc_y_im;
  double cc_y_re;
  double cc_y_im;
  double dc_y_re;
  double dc_y_im;
  double ec_y_re;
  double ec_y_im;
  double fc_y_re;
  double fc_y_im;
  double gc_y_re;
  double gc_y_im;
  double hc_y_re;
  double hc_y_im;
  double ic_y_re;
  double ic_y_im;
  double jc_y_re;
  double jc_y_im;
  double kc_y_re;
  double kc_y_im;
  double lc_y_re;
  double lc_y_im;
  double mc_y_re;
  double mc_y_im;
  double nc_y_re;
  double nc_y_im;
  double oc_y_re;
  double oc_y_im;
  double pc_y_re;
  double pc_y_im;
  double qc_y_re;
  double qc_y_im;
  double rc_y_re;
  double rc_y_im;
  double sc_y_re;
  double sc_y_im;
  double tc_y_re;
  double tc_y_im;
  double uc_y_re;
  double uc_y_im;
  double vc_y_re;
  double vc_y_im;
  double wc_y_re;
  double wc_y_im;
  double xc_y_re;
  double xc_y_im;
  double yc_y_re;
  double yc_y_im;
  double ad_y_re;
  double ad_y_im;
  double bd_y_re;
  double bd_y_im;
  double cd_y_re;
  double cd_y_im;
  double dd_y_re;
  double dd_y_im;
  double ed_y_re;
  double ed_y_im;
  double fd_y_re;
  double fd_y_im;
  double gd_y_re;
  double gd_y_im;
  double hd_y_re;
  double hd_y_im;
  double id_y_re;
  double id_y_im;
  double jd_y_re;
  double jd_y_im;
  double kd_y_re;
  double kd_y_im;
  double ld_y_re;
  double ld_y_im;
  double md_y_re;
  double md_y_im;
  double nd_y_re;
  double nd_y_im;
  double od_y_re;
  double od_y_im;
  double pd_y_re;
  double pd_y_im;
  double qd_y_re;
  double qd_y_im;
  double rd_y_re;
  double rd_y_im;
  double sd_y_re;
  double sd_y_im;
  double td_y_re;
  double td_y_im;
  double ud_y_re;
  double ud_y_im;
  double vd_y_re;
  double vd_y_im;
  double wd_y_re;
  double wd_y_im;
  double xd_y_re;
  double xd_y_im;
  double yd_y_re;
  double yd_y_im;
  double ae_y_re;
  double ae_y_im;
  double be_y_re;
  double be_y_im;
  double ce_y_re;
  double ce_y_im;
  double de_y_re;
  double de_y_im;
  double ee_y_re;
  double ee_y_im;
  double fe_y_re;
  double fe_y_im;
  double ge_y_re;
  double ge_y_im;
  double he_y_re;
  double he_y_im;
  double ie_y_re;
  double ie_y_im;
  double je_y_re;
  double je_y_im;
  double ke_y_re;
  double ke_y_im;
  double le_y_re;
  double le_y_im;
  double me_y_re;
  double me_y_im;
  double ne_y_re;
  double ne_y_im;
  double oe_y_re;
  double oe_y_im;
  double pe_y_re;
  double pe_y_im;
  double qe_y_re;
  double qe_y_im;
  double re_y_re;
  double re_y_im;
  double se_y_re;
  double se_y_im;
  double te_y_re;
  double te_y_im;
  double ue_y_re;
  double ue_y_im;
  double ve_y_re;
  double ve_y_im;
  double we_y_re;
  double we_y_im;
  double xe_y_re;
  double xe_y_im;
  double ye_y_re;
  double ye_y_im;
  double af_y_re;
  double af_y_im;
  double bf_y_re;
  double bf_y_im;
  double cf_y_re;
  double cf_y_im;
  double df_y_re;
  double df_y_im;
  double ef_y_re;
  double ef_y_im;
  double ff_y_re;
  double ff_y_im;
  double gf_y_re;
  double gf_y_im;
  double hf_y_re;
  double hf_y_im;
  double if_y_re;
  double if_y_im;
  double jf_y_re;
  double jf_y_im;
  double kf_y_re;
  double kf_y_im;
  double lf_y_re;
  double lf_y_im;
  double mf_y_re;
  double mf_y_im;
  double nf_y_re;
  double nf_y_im;
  double of_y_re;
  double of_y_im;
  double pf_y_re;
  double pf_y_im;
  double qf_y_re;
  double qf_y_im;
  double rf_y_re;
  double rf_y_im;
  double sf_y_re;
  double sf_y_im;
  double tf_y_re;
  double tf_y_im;
  double uf_y_re;
  double uf_y_im;
  double vf_y_re;
  double vf_y_im;
  double wf_y_re;
  double wf_y_im;
  double xf_y_re;
  double xf_y_im;
  double yf_y_re;
  double yf_y_im;
  double ag_y_re;
  double ag_y_im;
  double bg_y_re;
  double bg_y_im;
  double cg_y_re;
  double cg_y_im;
  double dg_y_re;
  double dg_y_im;
  double eg_y_re;
  double eg_y_im;
  double fg_y_re;
  double fg_y_im;
  double gg_y_re;
  double gg_y_im;
  double hg_y_re;
  double hg_y_im;
  double ig_y_re;
  double ig_y_im;
  double jg_y_re;
  double jg_y_im;
  double kg_y_re;
  double kg_y_im;
  double lg_y_re;
  double lg_y_im;
  double mg_y_re;
  double mg_y_im;
  double ng_y_re;
  double ng_y_im;
  double og_y_re;
  double og_y_im;
  double pg_y_re;
  double pg_y_im;
  double qg_y_re;
  double qg_y_im;
  double rg_y_re;
  double rg_y_im;
  double sg_y_re;
  double sg_y_im;
  double tg_y_re;
  double tg_y_im;
  double ug_y_re;
  double ug_y_im;
  double vg_y_re;
  double vg_y_im;
  double wg_y_re;
  double wg_y_im;
  double xg_y_re;
  double xg_y_im;
  double yg_y_re;
  double yg_y_im;
  double ah_y_re;
  double ah_y_im;
  double bh_y_re;
  double bh_y_im;
  double ch_y_re;
  double ch_y_im;
  double dh_y_re;
  double dh_y_im;
  double eh_y_re;
  double eh_y_im;
  double fh_y_re;
  double fh_y_im;
  double gh_y_re;
  double gh_y_im;
  double hh_y_re;
  double hh_y_im;
  double ih_y_re;
  double ih_y_im;
  double jh_y_re;
  double jh_y_im;
  double kh_y_re;
  double kh_y_im;
  double lh_y_re;
  double lh_y_im;
  double mh_y_re;
  double mh_y_im;
  double nh_y_re;
  double nh_y_im;
  double oh_y_re;
  double oh_y_im;
  double ph_y_re;
  double ph_y_im;
  double qh_y_re;
  double qh_y_im;
  double rh_y_re;
  double rh_y_im;
  double sh_y_re;
  double sh_y_im;
  double th_y_re;
  double th_y_im;
  double uh_y_re;
  double uh_y_im;
  double vh_y_re;
  double vh_y_im;
  double wh_y_re;
  double wh_y_im;
  double xh_y_re;
  double xh_y_im;
  double yh_y_re;
  double yh_y_im;
  double ai_y_re;
  double ai_y_im;
  double bi_y_re;
  double bi_y_im;
  double ci_y_re;
  double ci_y_im;
  double di_y_re;
  double di_y_im;
  double ei_y_re;
  double ei_y_im;
  double fi_y_re;
  double fi_y_im;
  double gi_y_re;
  double gi_y_im;
  double hi_y_re;
  double hi_y_im;
  double ii_y_re;
  double ii_y_im;
  double ji_y_re;
  double ji_y_im;
  double ki_y_re;
  double ki_y_im;
  double li_y_re;
  double li_y_im;
  double mi_y_re;
  double mi_y_im;
  double ni_y_re;
  double ni_y_im;
  double oi_y_re;
  double oi_y_im;
  double pi_y_re;
  double pi_y_im;
  double qi_y_re;
  double qi_y_im;
  double ri_y_re;
  double ri_y_im;
  double si_y_re;
  double si_y_im;
  double ti_y_re;
  double ti_y_im;
  double ui_y_re;
  double ui_y_im;
  double vi_y_re;
  double vi_y_im;
  double wi_y_re;
  double wi_y_im;
  double xi_y_re;
  double xi_y_im;
  double yi_y_re;
  double yi_y_im;
  double aj_y_re;
  double aj_y_im;
  double bj_y_re;
  double bj_y_im;
  double cj_y_re;
  double cj_y_im;
  double dj_y_re;
  double dj_y_im;
  double ej_y_re;
  double ej_y_im;
  double fj_y_re;
  double fj_y_im;
  double gj_y_re;
  double gj_y_im;
  double hj_y_re;
  double hj_y_im;
  double ij_y_re;
  double ij_y_im;
  double jj_y_re;
  double jj_y_im;
  double kj_y_re;
  double kj_y_im;
  double lj_y_re;
  double lj_y_im;
  double mj_y_re;
  double mj_y_im;
  double nj_y_re;
  double nj_y_im;
  double oj_y_re;
  double oj_y_im;
  double pj_y_re;
  double pj_y_im;
  double qj_y_re;
  double qj_y_im;
  double rj_y_re;
  double rj_y_im;
  double sj_y_re;
  double sj_y_im;
  double tj_y_re;
  double tj_y_im;
  double uj_y_re;
  double uj_y_im;
  double vj_y_re;
  double vj_y_im;
  double wj_y_re;
  double wj_y_im;
  double xj_y_re;
  double xj_y_im;
  double yj_y_re;
  double yj_y_im;
  double ak_y_re;
  double ak_y_im;
  double bk_y_re;
  double bk_y_im;
  double ck_y_re;
  double ck_y_im;
  double dk_y_re;
  double dk_y_im;
  double ek_y_re;
  double ek_y_im;
  double fk_y_re;
  double fk_y_im;
  double gk_y_re;
  double gk_y_im;
  double hk_y_re;
  double hk_y_im;
  double ik_y_re;
  double ik_y_im;
  double jk_y_re;
  double jk_y_im;
  double kk_y_re;
  double kk_y_im;
  double lk_y_re;
  double lk_y_im;
  double mk_y_re;
  double mk_y_im;
  double nk_y_re;
  double nk_y_im;
  double ok_y_re;
  double ok_y_im;
  double pk_y_re;
  double pk_y_im;
  double qk_y_re;
  double qk_y_im;
  double rk_y_re;
  double rk_y_im;
  double sk_y_re;
  double sk_y_im;
  double tk_y_re;
  double tk_y_im;
  double uk_y_re;
  double uk_y_im;
  double vk_y_re;
  double vk_y_im;
  double t_re_tmp;
  double t_im_tmp;
  double ar;
  double t_re;
  double t_s_re_tmp;
  double t_s_im_tmp;
  double ai;
  double brm;
  double re;
  double bim;
  double wk_y_re;
  double b_r;
  double b_re;
  double xk_y_re;
  double c_r;
  double c_re;
  double yk_y_re;
  double d_r;
  double d_re;
  double al_y_re;
  double e_r;
  double bl_y_re;
  double f_r;
  double cl_y_re;
  double g_r;
  double dl_y_re;
  double h_r;
  double el_y_re;
  double i_r;
  double fl_y_re;
  double j_r;
  double gl_y_re;
  double k_r;
  double hl_y_re;
  double l_r;
  double il_y_re;
  double m_r;
  double jl_y_re;
  double n_r;
  double kl_y_re;
  double o_r;
  double ll_y_re;
  double p_r;
  double ml_y_re;
  double q_r;
  double nl_y_re;
  double r_r;
  double ol_y_re;
  double s_r;
  double pl_y_re;
  double t_r;
  double ql_y_re;
  double u_r;
  double rl_y_re;
  double v_r;
  double sl_y_re;
  double w_r;
  double tl_y_re;
  double x_r;
  double ul_y_re;
  double y_r;
  double vl_y_re;
  double ab_r;
  double wl_y_re;
  double bb_r;
  double xl_y_re;
  double cb_r;
  double yl_y_re;
  double db_r;
  double am_y_re;
  double eb_r;
  double bm_y_re;
  double fb_r;
  double cm_y_re;
  double gb_r;
  double dm_y_re;
  double hb_r;
  double em_y_re;
  double ib_r;
  double fm_y_re;
  double jb_r;
  double gm_y_re;
  double kb_r;
  double hm_y_re;
  double lb_r;
  double im_y_re;
  double mb_r;
  double jm_y_re;
  double nb_r;
  double km_y_re;
  double ob_r;
  double lm_y_re;
  double pb_r;
  double mm_y_re;
  double qb_r;
  double nm_y_re;
  double rb_r;
  double om_y_re;
  double sb_r;
  double pm_y_re;
  double tb_r;
  double qm_y_re;
  double ub_r;
  double rm_y_re;
  double vb_r;
  double sm_y_re;
  double wb_r;
  double tm_y_re;
  double xb_r;
  double um_y_re;
  double yb_r;
  double vm_y_re;
  double ac_r;
  double wm_y_re;
  double bc_r;
  double xm_y_re;
  double cc_r;
  double ym_y_re;
  double dc_r;
  double an_y_re;
  double ec_r;
  double bn_y_re;
  double fc_r;
  double cn_y_re;
  double gc_r;
  double dn_y_re;
  double hc_r;
  double en_y_re;
  double ic_r;
  double fn_y_re;
  double jc_r;
  double gn_y_re;
  double kc_r;
  double hn_y_re;
  double lc_r;
  double in_y_re;
  double mc_r;
  double jn_y_re;
  double nc_r;
  double kn_y_re;
  double oc_r;
  double ln_y_re;
  double pc_r;
  double mn_y_re;
  double qc_r;
  double nn_y_re;
  double rc_r;
  double on_y_re;
  double sc_r;
  double pn_y_re;
  double tc_r;
  double qn_y_re;
  double uc_r;
  double rn_y_re;
  double vc_r;
  double sn_y_re;
  double wc_r;
  double tn_y_re;
  double xc_r;
  double un_y_re;
  double yc_r;
  double vn_y_re;
  double ad_r;
  double wn_y_re;
  double bd_r;
  double xn_y_re;
  double cd_r;
  double yn_y_re;
  double dd_r;
  double ao_y_re;
  double ed_r;
  double bo_y_re;
  double fd_r;
  double co_y_re;
  double gd_r;
  double do_y_re;
  double hd_r;
  double eo_y_re;
  double id_r;
  double fo_y_re;
  double jd_r;
  double go_y_re;
  double wk_y_im;
  double ho_y_re;
  double xk_y_im;
  double io_y_re;
  double yk_y_im;
  double jo_y_re;
  double al_y_im;
  double ko_y_re;
  double bl_y_im;
  double lo_y_re;
  double b_x_im;
  if ((t.im == 0.0) && (t.re >= 0.0)) {
    y_re = rt_powd_snf(t.re, 3.0);
    y_im = 0.0;
  } else if (t.re == 0.0) {
    y_re = 0.0;
    y_im = -rt_powd_snf(t.im, 3.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    y_re = 3.0 * x_re;
    y_im = 3.0 * x_im;
    if (y_im == 0.0) {
      y_re = exp(y_re);
      y_im = 0.0;
    } else if (rtIsInf(y_im) && rtIsInf(y_re) && (y_re < 0.0)) {
      y_re = 0.0;
      y_im = 0.0;
    } else {
      r = exp(y_re / 2.0);
      y_re = r * (r * cos(y_im));
      y_im = r * (r * sin(y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    b_y_re = rt_powd_snf(t_s.re, 3.0);
    b_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    b_y_re = 0.0;
    b_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    b_y_re = 3.0 * x_re;
    b_y_im = 3.0 * x_im;
    if (b_y_im == 0.0) {
      b_y_re = exp(b_y_re);
      b_y_im = 0.0;
    } else if (rtIsInf(b_y_im) && rtIsInf(b_y_re) && (b_y_re < 0.0)) {
      b_y_re = 0.0;
      b_y_im = 0.0;
    } else {
      r = exp(b_y_re / 2.0);
      b_y_re = r * (r * cos(b_y_im));
      b_y_im = r * (r * sin(b_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    c_y_re = rt_powd_snf(t.re, 3.0);
    c_y_im = 0.0;
  } else if (t.re == 0.0) {
    c_y_re = 0.0;
    c_y_im = -rt_powd_snf(t.im, 3.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    c_y_re = 3.0 * x_re;
    c_y_im = 3.0 * x_im;
    if (c_y_im == 0.0) {
      c_y_re = exp(c_y_re);
      c_y_im = 0.0;
    } else if (rtIsInf(c_y_im) && rtIsInf(c_y_re) && (c_y_re < 0.0)) {
      c_y_re = 0.0;
      c_y_im = 0.0;
    } else {
      r = exp(c_y_re / 2.0);
      c_y_re = r * (r * cos(c_y_im));
      c_y_im = r * (r * sin(c_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    d_y_re = rt_powd_snf(t_s.re, 4.0);
    d_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    d_y_re = rt_powd_snf(t_s.im, 4.0);
    d_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    d_y_re = 4.0 * x_re;
    d_y_im = 4.0 * x_im;
    if (d_y_im == 0.0) {
      d_y_re = exp(d_y_re);
      d_y_im = 0.0;
    } else if (rtIsInf(d_y_im) && rtIsInf(d_y_re) && (d_y_re < 0.0)) {
      d_y_re = 0.0;
      d_y_im = 0.0;
    } else {
      r = exp(d_y_re / 2.0);
      d_y_re = r * (r * cos(d_y_im));
      d_y_im = r * (r * sin(d_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    e_y_re = rt_powd_snf(t_s.re, 3.0);
    e_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    e_y_re = 0.0;
    e_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    e_y_re = 3.0 * x_re;
    e_y_im = 3.0 * x_im;
    if (e_y_im == 0.0) {
      e_y_re = exp(e_y_re);
      e_y_im = 0.0;
    } else if (rtIsInf(e_y_im) && rtIsInf(e_y_re) && (e_y_re < 0.0)) {
      e_y_re = 0.0;
      e_y_im = 0.0;
    } else {
      r = exp(e_y_re / 2.0);
      e_y_re = r * (r * cos(e_y_im));
      e_y_im = r * (r * sin(e_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    f_y_re = rt_powd_snf(t.re, 4.0);
    f_y_im = 0.0;
  } else if (t.re == 0.0) {
    f_y_re = rt_powd_snf(t.im, 4.0);
    f_y_im = 0.0;
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    f_y_re = 4.0 * x_re;
    f_y_im = 4.0 * x_im;
    if (f_y_im == 0.0) {
      f_y_re = exp(f_y_re);
      f_y_im = 0.0;
    } else if (rtIsInf(f_y_im) && rtIsInf(f_y_re) && (f_y_re < 0.0)) {
      f_y_re = 0.0;
      f_y_im = 0.0;
    } else {
      r = exp(f_y_re / 2.0);
      f_y_re = r * (r * cos(f_y_im));
      f_y_im = r * (r * sin(f_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    g_y_re = rt_powd_snf(t_s.re, 5.0);
    g_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    g_y_re = 0.0;
    g_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    g_y_re = 5.0 * x_re;
    g_y_im = 5.0 * x_im;
    if (g_y_im == 0.0) {
      g_y_re = exp(g_y_re);
      g_y_im = 0.0;
    } else if (rtIsInf(g_y_im) && rtIsInf(g_y_re) && (g_y_re < 0.0)) {
      g_y_re = 0.0;
      g_y_im = 0.0;
    } else {
      r = exp(g_y_re / 2.0);
      g_y_re = r * (r * cos(g_y_im));
      g_y_im = r * (r * sin(g_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    h_y_re = rt_powd_snf(t_s.re, 4.0);
    h_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    h_y_re = rt_powd_snf(t_s.im, 4.0);
    h_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    h_y_re = 4.0 * x_re;
    h_y_im = 4.0 * x_im;
    if (h_y_im == 0.0) {
      h_y_re = exp(h_y_re);
      h_y_im = 0.0;
    } else if (rtIsInf(h_y_im) && rtIsInf(h_y_re) && (h_y_re < 0.0)) {
      h_y_re = 0.0;
      h_y_im = 0.0;
    } else {
      r = exp(h_y_re / 2.0);
      h_y_re = r * (r * cos(h_y_im));
      h_y_im = r * (r * sin(h_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    i_y_re = rt_powd_snf(t.re, 5.0);
    i_y_im = 0.0;
  } else if (t.re == 0.0) {
    i_y_re = 0.0;
    i_y_im = rt_powd_snf(t.im, 5.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    i_y_re = 5.0 * x_re;
    i_y_im = 5.0 * x_im;
    if (i_y_im == 0.0) {
      i_y_re = exp(i_y_re);
      i_y_im = 0.0;
    } else if (rtIsInf(i_y_im) && rtIsInf(i_y_re) && (i_y_re < 0.0)) {
      i_y_re = 0.0;
      i_y_im = 0.0;
    } else {
      r = exp(i_y_re / 2.0);
      i_y_re = r * (r * cos(i_y_im));
      i_y_im = r * (r * sin(i_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    j_y_re = rt_powd_snf(t_s.re, 6.0);
    j_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    j_y_re = -rt_powd_snf(t_s.im, 6.0);
    j_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    j_y_re = 6.0 * x_re;
    j_y_im = 6.0 * x_im;
    if (j_y_im == 0.0) {
      j_y_re = exp(j_y_re);
      j_y_im = 0.0;
    } else if (rtIsInf(j_y_im) && rtIsInf(j_y_re) && (j_y_re < 0.0)) {
      j_y_re = 0.0;
      j_y_im = 0.0;
    } else {
      r = exp(j_y_re / 2.0);
      j_y_re = r * (r * cos(j_y_im));
      j_y_im = r * (r * sin(j_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    k_y_re = rt_powd_snf(t_s.re, 5.0);
    k_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    k_y_re = 0.0;
    k_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    k_y_re = 5.0 * x_re;
    k_y_im = 5.0 * x_im;
    if (k_y_im == 0.0) {
      k_y_re = exp(k_y_re);
      k_y_im = 0.0;
    } else if (rtIsInf(k_y_im) && rtIsInf(k_y_re) && (k_y_re < 0.0)) {
      k_y_re = 0.0;
      k_y_im = 0.0;
    } else {
      r = exp(k_y_re / 2.0);
      k_y_re = r * (r * cos(k_y_im));
      k_y_im = r * (r * sin(k_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    l_y_re = rt_powd_snf(t.re, 3.0);
    l_y_im = 0.0;
  } else if (t.re == 0.0) {
    l_y_re = 0.0;
    l_y_im = -rt_powd_snf(t.im, 3.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    l_y_re = 3.0 * x_re;
    l_y_im = 3.0 * x_im;
    if (l_y_im == 0.0) {
      l_y_re = exp(l_y_re);
      l_y_im = 0.0;
    } else if (rtIsInf(l_y_im) && rtIsInf(l_y_re) && (l_y_re < 0.0)) {
      l_y_re = 0.0;
      l_y_im = 0.0;
    } else {
      r = exp(l_y_re / 2.0);
      l_y_re = r * (r * cos(l_y_im));
      l_y_im = r * (r * sin(l_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    m_y_re = rt_powd_snf(t_s.re, 3.0);
    m_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    m_y_re = 0.0;
    m_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    m_y_re = 3.0 * x_re;
    m_y_im = 3.0 * x_im;
    if (m_y_im == 0.0) {
      m_y_re = exp(m_y_re);
      m_y_im = 0.0;
    } else if (rtIsInf(m_y_im) && rtIsInf(m_y_re) && (m_y_re < 0.0)) {
      m_y_re = 0.0;
      m_y_im = 0.0;
    } else {
      r = exp(m_y_re / 2.0);
      m_y_re = r * (r * cos(m_y_im));
      m_y_im = r * (r * sin(m_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    n_y_re = rt_powd_snf(t.re, 3.0);
    n_y_im = 0.0;
  } else if (t.re == 0.0) {
    n_y_re = 0.0;
    n_y_im = -rt_powd_snf(t.im, 3.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    n_y_re = 3.0 * x_re;
    n_y_im = 3.0 * x_im;
    if (n_y_im == 0.0) {
      n_y_re = exp(n_y_re);
      n_y_im = 0.0;
    } else if (rtIsInf(n_y_im) && rtIsInf(n_y_re) && (n_y_re < 0.0)) {
      n_y_re = 0.0;
      n_y_im = 0.0;
    } else {
      r = exp(n_y_re / 2.0);
      n_y_re = r * (r * cos(n_y_im));
      n_y_im = r * (r * sin(n_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    o_y_re = rt_powd_snf(t_s.re, 4.0);
    o_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    o_y_re = rt_powd_snf(t_s.im, 4.0);
    o_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    o_y_re = 4.0 * x_re;
    o_y_im = 4.0 * x_im;
    if (o_y_im == 0.0) {
      o_y_re = exp(o_y_re);
      o_y_im = 0.0;
    } else if (rtIsInf(o_y_im) && rtIsInf(o_y_re) && (o_y_re < 0.0)) {
      o_y_re = 0.0;
      o_y_im = 0.0;
    } else {
      r = exp(o_y_re / 2.0);
      o_y_re = r * (r * cos(o_y_im));
      o_y_im = r * (r * sin(o_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    p_y_re = rt_powd_snf(t_s.re, 3.0);
    p_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    p_y_re = 0.0;
    p_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    p_y_re = 3.0 * x_re;
    p_y_im = 3.0 * x_im;
    if (p_y_im == 0.0) {
      p_y_re = exp(p_y_re);
      p_y_im = 0.0;
    } else if (rtIsInf(p_y_im) && rtIsInf(p_y_re) && (p_y_re < 0.0)) {
      p_y_re = 0.0;
      p_y_im = 0.0;
    } else {
      r = exp(p_y_re / 2.0);
      p_y_re = r * (r * cos(p_y_im));
      p_y_im = r * (r * sin(p_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    q_y_re = rt_powd_snf(t.re, 4.0);
    q_y_im = 0.0;
  } else if (t.re == 0.0) {
    q_y_re = rt_powd_snf(t.im, 4.0);
    q_y_im = 0.0;
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    q_y_re = 4.0 * x_re;
    q_y_im = 4.0 * x_im;
    if (q_y_im == 0.0) {
      q_y_re = exp(q_y_re);
      q_y_im = 0.0;
    } else if (rtIsInf(q_y_im) && rtIsInf(q_y_re) && (q_y_re < 0.0)) {
      q_y_re = 0.0;
      q_y_im = 0.0;
    } else {
      r = exp(q_y_re / 2.0);
      q_y_re = r * (r * cos(q_y_im));
      q_y_im = r * (r * sin(q_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    r_y_re = rt_powd_snf(t_s.re, 5.0);
    r_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    r_y_re = 0.0;
    r_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    r_y_re = 5.0 * x_re;
    r_y_im = 5.0 * x_im;
    if (r_y_im == 0.0) {
      r_y_re = exp(r_y_re);
      r_y_im = 0.0;
    } else if (rtIsInf(r_y_im) && rtIsInf(r_y_re) && (r_y_re < 0.0)) {
      r_y_re = 0.0;
      r_y_im = 0.0;
    } else {
      r = exp(r_y_re / 2.0);
      r_y_re = r * (r * cos(r_y_im));
      r_y_im = r * (r * sin(r_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    s_y_re = rt_powd_snf(t_s.re, 4.0);
    s_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    s_y_re = rt_powd_snf(t_s.im, 4.0);
    s_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    s_y_re = 4.0 * x_re;
    s_y_im = 4.0 * x_im;
    if (s_y_im == 0.0) {
      s_y_re = exp(s_y_re);
      s_y_im = 0.0;
    } else if (rtIsInf(s_y_im) && rtIsInf(s_y_re) && (s_y_re < 0.0)) {
      s_y_re = 0.0;
      s_y_im = 0.0;
    } else {
      r = exp(s_y_re / 2.0);
      s_y_re = r * (r * cos(s_y_im));
      s_y_im = r * (r * sin(s_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    t_y_re = rt_powd_snf(t.re, 5.0);
    t_y_im = 0.0;
  } else if (t.re == 0.0) {
    t_y_re = 0.0;
    t_y_im = rt_powd_snf(t.im, 5.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    t_y_re = 5.0 * x_re;
    t_y_im = 5.0 * x_im;
    if (t_y_im == 0.0) {
      t_y_re = exp(t_y_re);
      t_y_im = 0.0;
    } else if (rtIsInf(t_y_im) && rtIsInf(t_y_re) && (t_y_re < 0.0)) {
      t_y_re = 0.0;
      t_y_im = 0.0;
    } else {
      r = exp(t_y_re / 2.0);
      t_y_re = r * (r * cos(t_y_im));
      t_y_im = r * (r * sin(t_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    u_y_re = rt_powd_snf(t_s.re, 6.0);
    u_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    u_y_re = -rt_powd_snf(t_s.im, 6.0);
    u_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    u_y_re = 6.0 * x_re;
    u_y_im = 6.0 * x_im;
    if (u_y_im == 0.0) {
      u_y_re = exp(u_y_re);
      u_y_im = 0.0;
    } else if (rtIsInf(u_y_im) && rtIsInf(u_y_re) && (u_y_re < 0.0)) {
      u_y_re = 0.0;
      u_y_im = 0.0;
    } else {
      r = exp(u_y_re / 2.0);
      u_y_re = r * (r * cos(u_y_im));
      u_y_im = r * (r * sin(u_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    v_y_re = rt_powd_snf(t_s.re, 5.0);
    v_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    v_y_re = 0.0;
    v_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    v_y_re = 5.0 * x_re;
    v_y_im = 5.0 * x_im;
    if (v_y_im == 0.0) {
      v_y_re = exp(v_y_re);
      v_y_im = 0.0;
    } else if (rtIsInf(v_y_im) && rtIsInf(v_y_re) && (v_y_re < 0.0)) {
      v_y_re = 0.0;
      v_y_im = 0.0;
    } else {
      r = exp(v_y_re / 2.0);
      v_y_re = r * (r * cos(v_y_im));
      v_y_im = r * (r * sin(v_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    w_y_re = rt_powd_snf(t.re, 3.0);
    w_y_im = 0.0;
  } else if (t.re == 0.0) {
    w_y_re = 0.0;
    w_y_im = -rt_powd_snf(t.im, 3.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    w_y_re = 3.0 * x_re;
    w_y_im = 3.0 * x_im;
    if (w_y_im == 0.0) {
      w_y_re = exp(w_y_re);
      w_y_im = 0.0;
    } else if (rtIsInf(w_y_im) && rtIsInf(w_y_re) && (w_y_re < 0.0)) {
      w_y_re = 0.0;
      w_y_im = 0.0;
    } else {
      r = exp(w_y_re / 2.0);
      w_y_re = r * (r * cos(w_y_im));
      w_y_im = r * (r * sin(w_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    x_y_re = rt_powd_snf(t_s.re, 4.0);
    x_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    x_y_re = rt_powd_snf(t_s.im, 4.0);
    x_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    x_y_re = 4.0 * x_re;
    x_y_im = 4.0 * x_im;
    if (x_y_im == 0.0) {
      x_y_re = exp(x_y_re);
      x_y_im = 0.0;
    } else if (rtIsInf(x_y_im) && rtIsInf(x_y_re) && (x_y_re < 0.0)) {
      x_y_re = 0.0;
      x_y_im = 0.0;
    } else {
      r = exp(x_y_re / 2.0);
      x_y_re = r * (r * cos(x_y_im));
      x_y_im = r * (r * sin(x_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    y_y_re = rt_powd_snf(t_s.re, 3.0);
    y_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    y_y_re = 0.0;
    y_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    y_y_re = 3.0 * x_re;
    y_y_im = 3.0 * x_im;
    if (y_y_im == 0.0) {
      y_y_re = exp(y_y_re);
      y_y_im = 0.0;
    } else if (rtIsInf(y_y_im) && rtIsInf(y_y_re) && (y_y_re < 0.0)) {
      y_y_re = 0.0;
      y_y_im = 0.0;
    } else {
      r = exp(y_y_re / 2.0);
      y_y_re = r * (r * cos(y_y_im));
      y_y_im = r * (r * sin(y_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    ab_y_re = rt_powd_snf(t.re, 3.0);
    ab_y_im = 0.0;
  } else if (t.re == 0.0) {
    ab_y_re = 0.0;
    ab_y_im = -rt_powd_snf(t.im, 3.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    ab_y_re = 3.0 * x_re;
    ab_y_im = 3.0 * x_im;
    if (ab_y_im == 0.0) {
      ab_y_re = exp(ab_y_re);
      ab_y_im = 0.0;
    } else if (rtIsInf(ab_y_im) && rtIsInf(ab_y_re) && (ab_y_re < 0.0)) {
      ab_y_re = 0.0;
      ab_y_im = 0.0;
    } else {
      r = exp(ab_y_re / 2.0);
      ab_y_re = r * (r * cos(ab_y_im));
      ab_y_im = r * (r * sin(ab_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    bb_y_re = rt_powd_snf(t_s.re, 3.0);
    bb_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    bb_y_re = 0.0;
    bb_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    bb_y_re = 3.0 * x_re;
    bb_y_im = 3.0 * x_im;
    if (bb_y_im == 0.0) {
      bb_y_re = exp(bb_y_re);
      bb_y_im = 0.0;
    } else if (rtIsInf(bb_y_im) && rtIsInf(bb_y_re) && (bb_y_re < 0.0)) {
      bb_y_re = 0.0;
      bb_y_im = 0.0;
    } else {
      r = exp(bb_y_re / 2.0);
      bb_y_re = r * (r * cos(bb_y_im));
      bb_y_im = r * (r * sin(bb_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    cb_y_re = rt_powd_snf(t.re, 3.0);
    cb_y_im = 0.0;
  } else if (t.re == 0.0) {
    cb_y_re = 0.0;
    cb_y_im = -rt_powd_snf(t.im, 3.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    cb_y_re = 3.0 * x_re;
    cb_y_im = 3.0 * x_im;
    if (cb_y_im == 0.0) {
      cb_y_re = exp(cb_y_re);
      cb_y_im = 0.0;
    } else if (rtIsInf(cb_y_im) && rtIsInf(cb_y_re) && (cb_y_re < 0.0)) {
      cb_y_re = 0.0;
      cb_y_im = 0.0;
    } else {
      r = exp(cb_y_re / 2.0);
      cb_y_re = r * (r * cos(cb_y_im));
      cb_y_im = r * (r * sin(cb_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    db_y_re = rt_powd_snf(t_s.re, 4.0);
    db_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    db_y_re = rt_powd_snf(t_s.im, 4.0);
    db_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    db_y_re = 4.0 * x_re;
    db_y_im = 4.0 * x_im;
    if (db_y_im == 0.0) {
      db_y_re = exp(db_y_re);
      db_y_im = 0.0;
    } else if (rtIsInf(db_y_im) && rtIsInf(db_y_re) && (db_y_re < 0.0)) {
      db_y_re = 0.0;
      db_y_im = 0.0;
    } else {
      r = exp(db_y_re / 2.0);
      db_y_re = r * (r * cos(db_y_im));
      db_y_im = r * (r * sin(db_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    eb_y_re = rt_powd_snf(t_s.re, 3.0);
    eb_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    eb_y_re = 0.0;
    eb_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    eb_y_re = 3.0 * x_re;
    eb_y_im = 3.0 * x_im;
    if (eb_y_im == 0.0) {
      eb_y_re = exp(eb_y_re);
      eb_y_im = 0.0;
    } else if (rtIsInf(eb_y_im) && rtIsInf(eb_y_re) && (eb_y_re < 0.0)) {
      eb_y_re = 0.0;
      eb_y_im = 0.0;
    } else {
      r = exp(eb_y_re / 2.0);
      eb_y_re = r * (r * cos(eb_y_im));
      eb_y_im = r * (r * sin(eb_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    fb_y_re = rt_powd_snf(t.re, 4.0);
    fb_y_im = 0.0;
  } else if (t.re == 0.0) {
    fb_y_re = rt_powd_snf(t.im, 4.0);
    fb_y_im = 0.0;
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    fb_y_re = 4.0 * x_re;
    fb_y_im = 4.0 * x_im;
    if (fb_y_im == 0.0) {
      fb_y_re = exp(fb_y_re);
      fb_y_im = 0.0;
    } else if (rtIsInf(fb_y_im) && rtIsInf(fb_y_re) && (fb_y_re < 0.0)) {
      fb_y_re = 0.0;
      fb_y_im = 0.0;
    } else {
      r = exp(fb_y_re / 2.0);
      fb_y_re = r * (r * cos(fb_y_im));
      fb_y_im = r * (r * sin(fb_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    gb_y_re = rt_powd_snf(t_s.re, 5.0);
    gb_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    gb_y_re = 0.0;
    gb_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    gb_y_re = 5.0 * x_re;
    gb_y_im = 5.0 * x_im;
    if (gb_y_im == 0.0) {
      gb_y_re = exp(gb_y_re);
      gb_y_im = 0.0;
    } else if (rtIsInf(gb_y_im) && rtIsInf(gb_y_re) && (gb_y_re < 0.0)) {
      gb_y_re = 0.0;
      gb_y_im = 0.0;
    } else {
      r = exp(gb_y_re / 2.0);
      gb_y_re = r * (r * cos(gb_y_im));
      gb_y_im = r * (r * sin(gb_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    hb_y_re = rt_powd_snf(t_s.re, 4.0);
    hb_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    hb_y_re = rt_powd_snf(t_s.im, 4.0);
    hb_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    hb_y_re = 4.0 * x_re;
    hb_y_im = 4.0 * x_im;
    if (hb_y_im == 0.0) {
      hb_y_re = exp(hb_y_re);
      hb_y_im = 0.0;
    } else if (rtIsInf(hb_y_im) && rtIsInf(hb_y_re) && (hb_y_re < 0.0)) {
      hb_y_re = 0.0;
      hb_y_im = 0.0;
    } else {
      r = exp(hb_y_re / 2.0);
      hb_y_re = r * (r * cos(hb_y_im));
      hb_y_im = r * (r * sin(hb_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    ib_y_re = rt_powd_snf(t.re, 4.0);
    ib_y_im = 0.0;
  } else if (t.re == 0.0) {
    ib_y_re = rt_powd_snf(t.im, 4.0);
    ib_y_im = 0.0;
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    ib_y_re = 4.0 * x_re;
    ib_y_im = 4.0 * x_im;
    if (ib_y_im == 0.0) {
      ib_y_re = exp(ib_y_re);
      ib_y_im = 0.0;
    } else if (rtIsInf(ib_y_im) && rtIsInf(ib_y_re) && (ib_y_re < 0.0)) {
      ib_y_re = 0.0;
      ib_y_im = 0.0;
    } else {
      r = exp(ib_y_re / 2.0);
      ib_y_re = r * (r * cos(ib_y_im));
      ib_y_im = r * (r * sin(ib_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    jb_y_re = rt_powd_snf(t_s.re, 5.0);
    jb_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    jb_y_re = 0.0;
    jb_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    jb_y_re = 5.0 * x_re;
    jb_y_im = 5.0 * x_im;
    if (jb_y_im == 0.0) {
      jb_y_re = exp(jb_y_re);
      jb_y_im = 0.0;
    } else if (rtIsInf(jb_y_im) && rtIsInf(jb_y_re) && (jb_y_re < 0.0)) {
      jb_y_re = 0.0;
      jb_y_im = 0.0;
    } else {
      r = exp(jb_y_re / 2.0);
      jb_y_re = r * (r * cos(jb_y_im));
      jb_y_im = r * (r * sin(jb_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    kb_y_re = rt_powd_snf(t_s.re, 4.0);
    kb_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    kb_y_re = rt_powd_snf(t_s.im, 4.0);
    kb_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    kb_y_re = 4.0 * x_re;
    kb_y_im = 4.0 * x_im;
    if (kb_y_im == 0.0) {
      kb_y_re = exp(kb_y_re);
      kb_y_im = 0.0;
    } else if (rtIsInf(kb_y_im) && rtIsInf(kb_y_re) && (kb_y_re < 0.0)) {
      kb_y_re = 0.0;
      kb_y_im = 0.0;
    } else {
      r = exp(kb_y_re / 2.0);
      kb_y_re = r * (r * cos(kb_y_im));
      kb_y_im = r * (r * sin(kb_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    lb_y_re = rt_powd_snf(t.re, 5.0);
    lb_y_im = 0.0;
  } else if (t.re == 0.0) {
    lb_y_re = 0.0;
    lb_y_im = rt_powd_snf(t.im, 5.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    lb_y_re = 5.0 * x_re;
    lb_y_im = 5.0 * x_im;
    if (lb_y_im == 0.0) {
      lb_y_re = exp(lb_y_re);
      lb_y_im = 0.0;
    } else if (rtIsInf(lb_y_im) && rtIsInf(lb_y_re) && (lb_y_re < 0.0)) {
      lb_y_re = 0.0;
      lb_y_im = 0.0;
    } else {
      r = exp(lb_y_re / 2.0);
      lb_y_re = r * (r * cos(lb_y_im));
      lb_y_im = r * (r * sin(lb_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    mb_y_re = rt_powd_snf(t_s.re, 6.0);
    mb_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    mb_y_re = -rt_powd_snf(t_s.im, 6.0);
    mb_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    mb_y_re = 6.0 * x_re;
    mb_y_im = 6.0 * x_im;
    if (mb_y_im == 0.0) {
      mb_y_re = exp(mb_y_re);
      mb_y_im = 0.0;
    } else if (rtIsInf(mb_y_im) && rtIsInf(mb_y_re) && (mb_y_re < 0.0)) {
      mb_y_re = 0.0;
      mb_y_im = 0.0;
    } else {
      r = exp(mb_y_re / 2.0);
      mb_y_re = r * (r * cos(mb_y_im));
      mb_y_im = r * (r * sin(mb_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    nb_y_re = rt_powd_snf(t_s.re, 5.0);
    nb_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    nb_y_re = 0.0;
    nb_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    nb_y_re = 5.0 * x_re;
    nb_y_im = 5.0 * x_im;
    if (nb_y_im == 0.0) {
      nb_y_re = exp(nb_y_re);
      nb_y_im = 0.0;
    } else if (rtIsInf(nb_y_im) && rtIsInf(nb_y_re) && (nb_y_re < 0.0)) {
      nb_y_re = 0.0;
      nb_y_im = 0.0;
    } else {
      r = exp(nb_y_re / 2.0);
      nb_y_re = r * (r * cos(nb_y_im));
      nb_y_im = r * (r * sin(nb_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    ob_y_re = rt_powd_snf(t.re, 5.0);
    ob_y_im = 0.0;
  } else if (t.re == 0.0) {
    ob_y_re = 0.0;
    ob_y_im = rt_powd_snf(t.im, 5.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    ob_y_re = 5.0 * x_re;
    ob_y_im = 5.0 * x_im;
    if (ob_y_im == 0.0) {
      ob_y_re = exp(ob_y_re);
      ob_y_im = 0.0;
    } else if (rtIsInf(ob_y_im) && rtIsInf(ob_y_re) && (ob_y_re < 0.0)) {
      ob_y_re = 0.0;
      ob_y_im = 0.0;
    } else {
      r = exp(ob_y_re / 2.0);
      ob_y_re = r * (r * cos(ob_y_im));
      ob_y_im = r * (r * sin(ob_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    pb_y_re = rt_powd_snf(t_s.re, 6.0);
    pb_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    pb_y_re = -rt_powd_snf(t_s.im, 6.0);
    pb_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    pb_y_re = 6.0 * x_re;
    pb_y_im = 6.0 * x_im;
    if (pb_y_im == 0.0) {
      pb_y_re = exp(pb_y_re);
      pb_y_im = 0.0;
    } else if (rtIsInf(pb_y_im) && rtIsInf(pb_y_re) && (pb_y_re < 0.0)) {
      pb_y_re = 0.0;
      pb_y_im = 0.0;
    } else {
      r = exp(pb_y_re / 2.0);
      pb_y_re = r * (r * cos(pb_y_im));
      pb_y_im = r * (r * sin(pb_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    qb_y_re = rt_powd_snf(t_s.re, 5.0);
    qb_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    qb_y_re = 0.0;
    qb_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    qb_y_re = 5.0 * x_re;
    qb_y_im = 5.0 * x_im;
    if (qb_y_im == 0.0) {
      qb_y_re = exp(qb_y_re);
      qb_y_im = 0.0;
    } else if (rtIsInf(qb_y_im) && rtIsInf(qb_y_re) && (qb_y_re < 0.0)) {
      qb_y_re = 0.0;
      qb_y_im = 0.0;
    } else {
      r = exp(qb_y_re / 2.0);
      qb_y_re = r * (r * cos(qb_y_im));
      qb_y_im = r * (r * sin(qb_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    rb_y_re = rt_powd_snf(t.re, 3.0);
    rb_y_im = 0.0;
  } else if (t.re == 0.0) {
    rb_y_re = 0.0;
    rb_y_im = -rt_powd_snf(t.im, 3.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    rb_y_re = 3.0 * x_re;
    rb_y_im = 3.0 * x_im;
    if (rb_y_im == 0.0) {
      rb_y_re = exp(rb_y_re);
      rb_y_im = 0.0;
    } else if (rtIsInf(rb_y_im) && rtIsInf(rb_y_re) && (rb_y_re < 0.0)) {
      rb_y_re = 0.0;
      rb_y_im = 0.0;
    } else {
      r = exp(rb_y_re / 2.0);
      rb_y_re = r * (r * cos(rb_y_im));
      rb_y_im = r * (r * sin(rb_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    sb_y_re = rt_powd_snf(t_s.re, 4.0);
    sb_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    sb_y_re = rt_powd_snf(t_s.im, 4.0);
    sb_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    sb_y_re = 4.0 * x_re;
    sb_y_im = 4.0 * x_im;
    if (sb_y_im == 0.0) {
      sb_y_re = exp(sb_y_re);
      sb_y_im = 0.0;
    } else if (rtIsInf(sb_y_im) && rtIsInf(sb_y_re) && (sb_y_re < 0.0)) {
      sb_y_re = 0.0;
      sb_y_im = 0.0;
    } else {
      r = exp(sb_y_re / 2.0);
      sb_y_re = r * (r * cos(sb_y_im));
      sb_y_im = r * (r * sin(sb_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    tb_y_re = rt_powd_snf(t_s.re, 3.0);
    tb_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    tb_y_re = 0.0;
    tb_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    tb_y_re = 3.0 * x_re;
    tb_y_im = 3.0 * x_im;
    if (tb_y_im == 0.0) {
      tb_y_re = exp(tb_y_re);
      tb_y_im = 0.0;
    } else if (rtIsInf(tb_y_im) && rtIsInf(tb_y_re) && (tb_y_re < 0.0)) {
      tb_y_re = 0.0;
      tb_y_im = 0.0;
    } else {
      r = exp(tb_y_re / 2.0);
      tb_y_re = r * (r * cos(tb_y_im));
      tb_y_im = r * (r * sin(tb_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    ub_y_re = rt_powd_snf(t.re, 4.0);
    ub_y_im = 0.0;
  } else if (t.re == 0.0) {
    ub_y_re = rt_powd_snf(t.im, 4.0);
    ub_y_im = 0.0;
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    ub_y_re = 4.0 * x_re;
    ub_y_im = 4.0 * x_im;
    if (ub_y_im == 0.0) {
      ub_y_re = exp(ub_y_re);
      ub_y_im = 0.0;
    } else if (rtIsInf(ub_y_im) && rtIsInf(ub_y_re) && (ub_y_re < 0.0)) {
      ub_y_re = 0.0;
      ub_y_im = 0.0;
    } else {
      r = exp(ub_y_re / 2.0);
      ub_y_re = r * (r * cos(ub_y_im));
      ub_y_im = r * (r * sin(ub_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    vb_y_re = rt_powd_snf(t_s.re, 5.0);
    vb_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    vb_y_re = 0.0;
    vb_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    vb_y_re = 5.0 * x_re;
    vb_y_im = 5.0 * x_im;
    if (vb_y_im == 0.0) {
      vb_y_re = exp(vb_y_re);
      vb_y_im = 0.0;
    } else if (rtIsInf(vb_y_im) && rtIsInf(vb_y_re) && (vb_y_re < 0.0)) {
      vb_y_re = 0.0;
      vb_y_im = 0.0;
    } else {
      r = exp(vb_y_re / 2.0);
      vb_y_re = r * (r * cos(vb_y_im));
      vb_y_im = r * (r * sin(vb_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    wb_y_re = rt_powd_snf(t_s.re, 4.0);
    wb_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    wb_y_re = rt_powd_snf(t_s.im, 4.0);
    wb_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    wb_y_re = 4.0 * x_re;
    wb_y_im = 4.0 * x_im;
    if (wb_y_im == 0.0) {
      wb_y_re = exp(wb_y_re);
      wb_y_im = 0.0;
    } else if (rtIsInf(wb_y_im) && rtIsInf(wb_y_re) && (wb_y_re < 0.0)) {
      wb_y_re = 0.0;
      wb_y_im = 0.0;
    } else {
      r = exp(wb_y_re / 2.0);
      wb_y_re = r * (r * cos(wb_y_im));
      wb_y_im = r * (r * sin(wb_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    xb_y_re = rt_powd_snf(t.re, 5.0);
    xb_y_im = 0.0;
  } else if (t.re == 0.0) {
    xb_y_re = 0.0;
    xb_y_im = rt_powd_snf(t.im, 5.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    xb_y_re = 5.0 * x_re;
    xb_y_im = 5.0 * x_im;
    if (xb_y_im == 0.0) {
      xb_y_re = exp(xb_y_re);
      xb_y_im = 0.0;
    } else if (rtIsInf(xb_y_im) && rtIsInf(xb_y_re) && (xb_y_re < 0.0)) {
      xb_y_re = 0.0;
      xb_y_im = 0.0;
    } else {
      r = exp(xb_y_re / 2.0);
      xb_y_re = r * (r * cos(xb_y_im));
      xb_y_im = r * (r * sin(xb_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    yb_y_re = rt_powd_snf(t_s.re, 6.0);
    yb_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    yb_y_re = -rt_powd_snf(t_s.im, 6.0);
    yb_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    yb_y_re = 6.0 * x_re;
    yb_y_im = 6.0 * x_im;
    if (yb_y_im == 0.0) {
      yb_y_re = exp(yb_y_re);
      yb_y_im = 0.0;
    } else if (rtIsInf(yb_y_im) && rtIsInf(yb_y_re) && (yb_y_re < 0.0)) {
      yb_y_re = 0.0;
      yb_y_im = 0.0;
    } else {
      r = exp(yb_y_re / 2.0);
      yb_y_re = r * (r * cos(yb_y_im));
      yb_y_im = r * (r * sin(yb_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    ac_y_re = rt_powd_snf(t_s.re, 5.0);
    ac_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    ac_y_re = 0.0;
    ac_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    ac_y_re = 5.0 * x_re;
    ac_y_im = 5.0 * x_im;
    if (ac_y_im == 0.0) {
      ac_y_re = exp(ac_y_re);
      ac_y_im = 0.0;
    } else if (rtIsInf(ac_y_im) && rtIsInf(ac_y_re) && (ac_y_re < 0.0)) {
      ac_y_re = 0.0;
      ac_y_im = 0.0;
    } else {
      r = exp(ac_y_re / 2.0);
      ac_y_re = r * (r * cos(ac_y_im));
      ac_y_im = r * (r * sin(ac_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    bc_y_re = rt_powd_snf(t.re, 3.0);
    bc_y_im = 0.0;
  } else if (t.re == 0.0) {
    bc_y_re = 0.0;
    bc_y_im = -rt_powd_snf(t.im, 3.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    bc_y_re = 3.0 * x_re;
    bc_y_im = 3.0 * x_im;
    if (bc_y_im == 0.0) {
      bc_y_re = exp(bc_y_re);
      bc_y_im = 0.0;
    } else if (rtIsInf(bc_y_im) && rtIsInf(bc_y_re) && (bc_y_re < 0.0)) {
      bc_y_re = 0.0;
      bc_y_im = 0.0;
    } else {
      r = exp(bc_y_re / 2.0);
      bc_y_re = r * (r * cos(bc_y_im));
      bc_y_im = r * (r * sin(bc_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    cc_y_re = rt_powd_snf(t_s.re, 4.0);
    cc_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    cc_y_re = rt_powd_snf(t_s.im, 4.0);
    cc_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    cc_y_re = 4.0 * x_re;
    cc_y_im = 4.0 * x_im;
    if (cc_y_im == 0.0) {
      cc_y_re = exp(cc_y_re);
      cc_y_im = 0.0;
    } else if (rtIsInf(cc_y_im) && rtIsInf(cc_y_re) && (cc_y_re < 0.0)) {
      cc_y_re = 0.0;
      cc_y_im = 0.0;
    } else {
      r = exp(cc_y_re / 2.0);
      cc_y_re = r * (r * cos(cc_y_im));
      cc_y_im = r * (r * sin(cc_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    dc_y_re = rt_powd_snf(t_s.re, 3.0);
    dc_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    dc_y_re = 0.0;
    dc_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    dc_y_re = 3.0 * x_re;
    dc_y_im = 3.0 * x_im;
    if (dc_y_im == 0.0) {
      dc_y_re = exp(dc_y_re);
      dc_y_im = 0.0;
    } else if (rtIsInf(dc_y_im) && rtIsInf(dc_y_re) && (dc_y_re < 0.0)) {
      dc_y_re = 0.0;
      dc_y_im = 0.0;
    } else {
      r = exp(dc_y_re / 2.0);
      dc_y_re = r * (r * cos(dc_y_im));
      dc_y_im = r * (r * sin(dc_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    ec_y_re = rt_powd_snf(t.re, 3.0);
    ec_y_im = 0.0;
  } else if (t.re == 0.0) {
    ec_y_re = 0.0;
    ec_y_im = -rt_powd_snf(t.im, 3.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    ec_y_re = 3.0 * x_re;
    ec_y_im = 3.0 * x_im;
    if (ec_y_im == 0.0) {
      ec_y_re = exp(ec_y_re);
      ec_y_im = 0.0;
    } else if (rtIsInf(ec_y_im) && rtIsInf(ec_y_re) && (ec_y_re < 0.0)) {
      ec_y_re = 0.0;
      ec_y_im = 0.0;
    } else {
      r = exp(ec_y_re / 2.0);
      ec_y_re = r * (r * cos(ec_y_im));
      ec_y_im = r * (r * sin(ec_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    fc_y_re = rt_powd_snf(t_s.re, 3.0);
    fc_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    fc_y_re = 0.0;
    fc_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    fc_y_re = 3.0 * x_re;
    fc_y_im = 3.0 * x_im;
    if (fc_y_im == 0.0) {
      fc_y_re = exp(fc_y_re);
      fc_y_im = 0.0;
    } else if (rtIsInf(fc_y_im) && rtIsInf(fc_y_re) && (fc_y_re < 0.0)) {
      fc_y_re = 0.0;
      fc_y_im = 0.0;
    } else {
      r = exp(fc_y_re / 2.0);
      fc_y_re = r * (r * cos(fc_y_im));
      fc_y_im = r * (r * sin(fc_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    gc_y_re = rt_powd_snf(t.re, 3.0);
    gc_y_im = 0.0;
  } else if (t.re == 0.0) {
    gc_y_re = 0.0;
    gc_y_im = -rt_powd_snf(t.im, 3.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    gc_y_re = 3.0 * x_re;
    gc_y_im = 3.0 * x_im;
    if (gc_y_im == 0.0) {
      gc_y_re = exp(gc_y_re);
      gc_y_im = 0.0;
    } else if (rtIsInf(gc_y_im) && rtIsInf(gc_y_re) && (gc_y_re < 0.0)) {
      gc_y_re = 0.0;
      gc_y_im = 0.0;
    } else {
      r = exp(gc_y_re / 2.0);
      gc_y_re = r * (r * cos(gc_y_im));
      gc_y_im = r * (r * sin(gc_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    hc_y_re = rt_powd_snf(t_s.re, 4.0);
    hc_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    hc_y_re = rt_powd_snf(t_s.im, 4.0);
    hc_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    hc_y_re = 4.0 * x_re;
    hc_y_im = 4.0 * x_im;
    if (hc_y_im == 0.0) {
      hc_y_re = exp(hc_y_re);
      hc_y_im = 0.0;
    } else if (rtIsInf(hc_y_im) && rtIsInf(hc_y_re) && (hc_y_re < 0.0)) {
      hc_y_re = 0.0;
      hc_y_im = 0.0;
    } else {
      r = exp(hc_y_re / 2.0);
      hc_y_re = r * (r * cos(hc_y_im));
      hc_y_im = r * (r * sin(hc_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    ic_y_re = rt_powd_snf(t_s.re, 3.0);
    ic_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    ic_y_re = 0.0;
    ic_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    ic_y_re = 3.0 * x_re;
    ic_y_im = 3.0 * x_im;
    if (ic_y_im == 0.0) {
      ic_y_re = exp(ic_y_re);
      ic_y_im = 0.0;
    } else if (rtIsInf(ic_y_im) && rtIsInf(ic_y_re) && (ic_y_re < 0.0)) {
      ic_y_re = 0.0;
      ic_y_im = 0.0;
    } else {
      r = exp(ic_y_re / 2.0);
      ic_y_re = r * (r * cos(ic_y_im));
      ic_y_im = r * (r * sin(ic_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    jc_y_re = rt_powd_snf(t.re, 4.0);
    jc_y_im = 0.0;
  } else if (t.re == 0.0) {
    jc_y_re = rt_powd_snf(t.im, 4.0);
    jc_y_im = 0.0;
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    jc_y_re = 4.0 * x_re;
    jc_y_im = 4.0 * x_im;
    if (jc_y_im == 0.0) {
      jc_y_re = exp(jc_y_re);
      jc_y_im = 0.0;
    } else if (rtIsInf(jc_y_im) && rtIsInf(jc_y_re) && (jc_y_re < 0.0)) {
      jc_y_re = 0.0;
      jc_y_im = 0.0;
    } else {
      r = exp(jc_y_re / 2.0);
      jc_y_re = r * (r * cos(jc_y_im));
      jc_y_im = r * (r * sin(jc_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    kc_y_re = rt_powd_snf(t_s.re, 5.0);
    kc_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    kc_y_re = 0.0;
    kc_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    kc_y_re = 5.0 * x_re;
    kc_y_im = 5.0 * x_im;
    if (kc_y_im == 0.0) {
      kc_y_re = exp(kc_y_re);
      kc_y_im = 0.0;
    } else if (rtIsInf(kc_y_im) && rtIsInf(kc_y_re) && (kc_y_re < 0.0)) {
      kc_y_re = 0.0;
      kc_y_im = 0.0;
    } else {
      r = exp(kc_y_re / 2.0);
      kc_y_re = r * (r * cos(kc_y_im));
      kc_y_im = r * (r * sin(kc_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    lc_y_re = rt_powd_snf(t_s.re, 4.0);
    lc_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    lc_y_re = rt_powd_snf(t_s.im, 4.0);
    lc_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    lc_y_re = 4.0 * x_re;
    lc_y_im = 4.0 * x_im;
    if (lc_y_im == 0.0) {
      lc_y_re = exp(lc_y_re);
      lc_y_im = 0.0;
    } else if (rtIsInf(lc_y_im) && rtIsInf(lc_y_re) && (lc_y_re < 0.0)) {
      lc_y_re = 0.0;
      lc_y_im = 0.0;
    } else {
      r = exp(lc_y_re / 2.0);
      lc_y_re = r * (r * cos(lc_y_im));
      lc_y_im = r * (r * sin(lc_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    mc_y_re = rt_powd_snf(t.re, 3.0);
    mc_y_im = 0.0;
  } else if (t.re == 0.0) {
    mc_y_re = 0.0;
    mc_y_im = -rt_powd_snf(t.im, 3.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    mc_y_re = 3.0 * x_re;
    mc_y_im = 3.0 * x_im;
    if (mc_y_im == 0.0) {
      mc_y_re = exp(mc_y_re);
      mc_y_im = 0.0;
    } else if (rtIsInf(mc_y_im) && rtIsInf(mc_y_re) && (mc_y_re < 0.0)) {
      mc_y_re = 0.0;
      mc_y_im = 0.0;
    } else {
      r = exp(mc_y_re / 2.0);
      mc_y_re = r * (r * cos(mc_y_im));
      mc_y_im = r * (r * sin(mc_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    nc_y_re = rt_powd_snf(t_s.re, 3.0);
    nc_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    nc_y_re = 0.0;
    nc_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    nc_y_re = 3.0 * x_re;
    nc_y_im = 3.0 * x_im;
    if (nc_y_im == 0.0) {
      nc_y_re = exp(nc_y_re);
      nc_y_im = 0.0;
    } else if (rtIsInf(nc_y_im) && rtIsInf(nc_y_re) && (nc_y_re < 0.0)) {
      nc_y_re = 0.0;
      nc_y_im = 0.0;
    } else {
      r = exp(nc_y_re / 2.0);
      nc_y_re = r * (r * cos(nc_y_im));
      nc_y_im = r * (r * sin(nc_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    oc_y_re = rt_powd_snf(t_s.re, 4.0);
    oc_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    oc_y_re = rt_powd_snf(t_s.im, 4.0);
    oc_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    oc_y_re = 4.0 * x_re;
    oc_y_im = 4.0 * x_im;
    if (oc_y_im == 0.0) {
      oc_y_re = exp(oc_y_re);
      oc_y_im = 0.0;
    } else if (rtIsInf(oc_y_im) && rtIsInf(oc_y_re) && (oc_y_re < 0.0)) {
      oc_y_re = 0.0;
      oc_y_im = 0.0;
    } else {
      r = exp(oc_y_re / 2.0);
      oc_y_re = r * (r * cos(oc_y_im));
      oc_y_im = r * (r * sin(oc_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    pc_y_re = rt_powd_snf(t_s.re, 3.0);
    pc_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    pc_y_re = 0.0;
    pc_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    pc_y_re = 3.0 * x_re;
    pc_y_im = 3.0 * x_im;
    if (pc_y_im == 0.0) {
      pc_y_re = exp(pc_y_re);
      pc_y_im = 0.0;
    } else if (rtIsInf(pc_y_im) && rtIsInf(pc_y_re) && (pc_y_re < 0.0)) {
      pc_y_re = 0.0;
      pc_y_im = 0.0;
    } else {
      r = exp(pc_y_re / 2.0);
      pc_y_re = r * (r * cos(pc_y_im));
      pc_y_im = r * (r * sin(pc_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    qc_y_re = rt_powd_snf(t.re, 4.0);
    qc_y_im = 0.0;
  } else if (t.re == 0.0) {
    qc_y_re = rt_powd_snf(t.im, 4.0);
    qc_y_im = 0.0;
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    qc_y_re = 4.0 * x_re;
    qc_y_im = 4.0 * x_im;
    if (qc_y_im == 0.0) {
      qc_y_re = exp(qc_y_re);
      qc_y_im = 0.0;
    } else if (rtIsInf(qc_y_im) && rtIsInf(qc_y_re) && (qc_y_re < 0.0)) {
      qc_y_re = 0.0;
      qc_y_im = 0.0;
    } else {
      r = exp(qc_y_re / 2.0);
      qc_y_re = r * (r * cos(qc_y_im));
      qc_y_im = r * (r * sin(qc_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    rc_y_re = rt_powd_snf(t_s.re, 5.0);
    rc_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    rc_y_re = 0.0;
    rc_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    rc_y_re = 5.0 * x_re;
    rc_y_im = 5.0 * x_im;
    if (rc_y_im == 0.0) {
      rc_y_re = exp(rc_y_re);
      rc_y_im = 0.0;
    } else if (rtIsInf(rc_y_im) && rtIsInf(rc_y_re) && (rc_y_re < 0.0)) {
      rc_y_re = 0.0;
      rc_y_im = 0.0;
    } else {
      r = exp(rc_y_re / 2.0);
      rc_y_re = r * (r * cos(rc_y_im));
      rc_y_im = r * (r * sin(rc_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    sc_y_re = rt_powd_snf(t_s.re, 4.0);
    sc_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    sc_y_re = rt_powd_snf(t_s.im, 4.0);
    sc_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    sc_y_re = 4.0 * x_re;
    sc_y_im = 4.0 * x_im;
    if (sc_y_im == 0.0) {
      sc_y_re = exp(sc_y_re);
      sc_y_im = 0.0;
    } else if (rtIsInf(sc_y_im) && rtIsInf(sc_y_re) && (sc_y_re < 0.0)) {
      sc_y_re = 0.0;
      sc_y_im = 0.0;
    } else {
      r = exp(sc_y_re / 2.0);
      sc_y_re = r * (r * cos(sc_y_im));
      sc_y_im = r * (r * sin(sc_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    tc_y_re = rt_powd_snf(t.re, 5.0);
    tc_y_im = 0.0;
  } else if (t.re == 0.0) {
    tc_y_re = 0.0;
    tc_y_im = rt_powd_snf(t.im, 5.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    tc_y_re = 5.0 * x_re;
    tc_y_im = 5.0 * x_im;
    if (tc_y_im == 0.0) {
      tc_y_re = exp(tc_y_re);
      tc_y_im = 0.0;
    } else if (rtIsInf(tc_y_im) && rtIsInf(tc_y_re) && (tc_y_re < 0.0)) {
      tc_y_re = 0.0;
      tc_y_im = 0.0;
    } else {
      r = exp(tc_y_re / 2.0);
      tc_y_re = r * (r * cos(tc_y_im));
      tc_y_im = r * (r * sin(tc_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    uc_y_re = rt_powd_snf(t_s.re, 6.0);
    uc_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    uc_y_re = -rt_powd_snf(t_s.im, 6.0);
    uc_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    uc_y_re = 6.0 * x_re;
    uc_y_im = 6.0 * x_im;
    if (uc_y_im == 0.0) {
      uc_y_re = exp(uc_y_re);
      uc_y_im = 0.0;
    } else if (rtIsInf(uc_y_im) && rtIsInf(uc_y_re) && (uc_y_re < 0.0)) {
      uc_y_re = 0.0;
      uc_y_im = 0.0;
    } else {
      r = exp(uc_y_re / 2.0);
      uc_y_re = r * (r * cos(uc_y_im));
      uc_y_im = r * (r * sin(uc_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    vc_y_re = rt_powd_snf(t_s.re, 5.0);
    vc_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    vc_y_re = 0.0;
    vc_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    vc_y_re = 5.0 * x_re;
    vc_y_im = 5.0 * x_im;
    if (vc_y_im == 0.0) {
      vc_y_re = exp(vc_y_re);
      vc_y_im = 0.0;
    } else if (rtIsInf(vc_y_im) && rtIsInf(vc_y_re) && (vc_y_re < 0.0)) {
      vc_y_re = 0.0;
      vc_y_im = 0.0;
    } else {
      r = exp(vc_y_re / 2.0);
      vc_y_re = r * (r * cos(vc_y_im));
      vc_y_im = r * (r * sin(vc_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    wc_y_re = rt_powd_snf(t.re, 4.0);
    wc_y_im = 0.0;
  } else if (t.re == 0.0) {
    wc_y_re = rt_powd_snf(t.im, 4.0);
    wc_y_im = 0.0;
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    wc_y_re = 4.0 * x_re;
    wc_y_im = 4.0 * x_im;
    if (wc_y_im == 0.0) {
      wc_y_re = exp(wc_y_re);
      wc_y_im = 0.0;
    } else if (rtIsInf(wc_y_im) && rtIsInf(wc_y_re) && (wc_y_re < 0.0)) {
      wc_y_re = 0.0;
      wc_y_im = 0.0;
    } else {
      r = exp(wc_y_re / 2.0);
      wc_y_re = r * (r * cos(wc_y_im));
      wc_y_im = r * (r * sin(wc_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    xc_y_re = rt_powd_snf(t_s.re, 3.0);
    xc_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    xc_y_re = 0.0;
    xc_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    xc_y_re = 3.0 * x_re;
    xc_y_im = 3.0 * x_im;
    if (xc_y_im == 0.0) {
      xc_y_re = exp(xc_y_re);
      xc_y_im = 0.0;
    } else if (rtIsInf(xc_y_im) && rtIsInf(xc_y_re) && (xc_y_re < 0.0)) {
      xc_y_re = 0.0;
      xc_y_im = 0.0;
    } else {
      r = exp(xc_y_re / 2.0);
      xc_y_re = r * (r * cos(xc_y_im));
      xc_y_im = r * (r * sin(xc_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    yc_y_re = rt_powd_snf(t_s.re, 5.0);
    yc_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    yc_y_re = 0.0;
    yc_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    yc_y_re = 5.0 * x_re;
    yc_y_im = 5.0 * x_im;
    if (yc_y_im == 0.0) {
      yc_y_re = exp(yc_y_re);
      yc_y_im = 0.0;
    } else if (rtIsInf(yc_y_im) && rtIsInf(yc_y_re) && (yc_y_re < 0.0)) {
      yc_y_re = 0.0;
      yc_y_im = 0.0;
    } else {
      r = exp(yc_y_re / 2.0);
      yc_y_re = r * (r * cos(yc_y_im));
      yc_y_im = r * (r * sin(yc_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    ad_y_re = rt_powd_snf(t_s.re, 4.0);
    ad_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    ad_y_re = rt_powd_snf(t_s.im, 4.0);
    ad_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    ad_y_re = 4.0 * x_re;
    ad_y_im = 4.0 * x_im;
    if (ad_y_im == 0.0) {
      ad_y_re = exp(ad_y_re);
      ad_y_im = 0.0;
    } else if (rtIsInf(ad_y_im) && rtIsInf(ad_y_re) && (ad_y_re < 0.0)) {
      ad_y_re = 0.0;
      ad_y_im = 0.0;
    } else {
      r = exp(ad_y_re / 2.0);
      ad_y_re = r * (r * cos(ad_y_im));
      ad_y_im = r * (r * sin(ad_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    bd_y_re = rt_powd_snf(t.re, 5.0);
    bd_y_im = 0.0;
  } else if (t.re == 0.0) {
    bd_y_re = 0.0;
    bd_y_im = rt_powd_snf(t.im, 5.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    bd_y_re = 5.0 * x_re;
    bd_y_im = 5.0 * x_im;
    if (bd_y_im == 0.0) {
      bd_y_re = exp(bd_y_re);
      bd_y_im = 0.0;
    } else if (rtIsInf(bd_y_im) && rtIsInf(bd_y_re) && (bd_y_re < 0.0)) {
      bd_y_re = 0.0;
      bd_y_im = 0.0;
    } else {
      r = exp(bd_y_re / 2.0);
      bd_y_re = r * (r * cos(bd_y_im));
      bd_y_im = r * (r * sin(bd_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    cd_y_re = rt_powd_snf(t_s.re, 6.0);
    cd_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    cd_y_re = -rt_powd_snf(t_s.im, 6.0);
    cd_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    cd_y_re = 6.0 * x_re;
    cd_y_im = 6.0 * x_im;
    if (cd_y_im == 0.0) {
      cd_y_re = exp(cd_y_re);
      cd_y_im = 0.0;
    } else if (rtIsInf(cd_y_im) && rtIsInf(cd_y_re) && (cd_y_re < 0.0)) {
      cd_y_re = 0.0;
      cd_y_im = 0.0;
    } else {
      r = exp(cd_y_re / 2.0);
      cd_y_re = r * (r * cos(cd_y_im));
      cd_y_im = r * (r * sin(cd_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    dd_y_re = rt_powd_snf(t_s.re, 5.0);
    dd_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    dd_y_re = 0.0;
    dd_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    dd_y_re = 5.0 * x_re;
    dd_y_im = 5.0 * x_im;
    if (dd_y_im == 0.0) {
      dd_y_re = exp(dd_y_re);
      dd_y_im = 0.0;
    } else if (rtIsInf(dd_y_im) && rtIsInf(dd_y_re) && (dd_y_re < 0.0)) {
      dd_y_re = 0.0;
      dd_y_im = 0.0;
    } else {
      r = exp(dd_y_re / 2.0);
      dd_y_re = r * (r * cos(dd_y_im));
      dd_y_im = r * (r * sin(dd_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    ed_y_re = rt_powd_snf(t.re, 5.0);
    ed_y_im = 0.0;
  } else if (t.re == 0.0) {
    ed_y_re = 0.0;
    ed_y_im = rt_powd_snf(t.im, 5.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    ed_y_re = 5.0 * x_re;
    ed_y_im = 5.0 * x_im;
    if (ed_y_im == 0.0) {
      ed_y_re = exp(ed_y_re);
      ed_y_im = 0.0;
    } else if (rtIsInf(ed_y_im) && rtIsInf(ed_y_re) && (ed_y_re < 0.0)) {
      ed_y_re = 0.0;
      ed_y_im = 0.0;
    } else {
      r = exp(ed_y_re / 2.0);
      ed_y_re = r * (r * cos(ed_y_im));
      ed_y_im = r * (r * sin(ed_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    fd_y_re = rt_powd_snf(t_s.re, 3.0);
    fd_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    fd_y_re = 0.0;
    fd_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    fd_y_re = 3.0 * x_re;
    fd_y_im = 3.0 * x_im;
    if (fd_y_im == 0.0) {
      fd_y_re = exp(fd_y_re);
      fd_y_im = 0.0;
    } else if (rtIsInf(fd_y_im) && rtIsInf(fd_y_re) && (fd_y_re < 0.0)) {
      fd_y_re = 0.0;
      fd_y_im = 0.0;
    } else {
      r = exp(fd_y_re / 2.0);
      fd_y_re = r * (r * cos(fd_y_im));
      fd_y_im = r * (r * sin(fd_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    gd_y_re = rt_powd_snf(t_s.re, 6.0);
    gd_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    gd_y_re = -rt_powd_snf(t_s.im, 6.0);
    gd_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    gd_y_re = 6.0 * x_re;
    gd_y_im = 6.0 * x_im;
    if (gd_y_im == 0.0) {
      gd_y_re = exp(gd_y_re);
      gd_y_im = 0.0;
    } else if (rtIsInf(gd_y_im) && rtIsInf(gd_y_re) && (gd_y_re < 0.0)) {
      gd_y_re = 0.0;
      gd_y_im = 0.0;
    } else {
      r = exp(gd_y_re / 2.0);
      gd_y_re = r * (r * cos(gd_y_im));
      gd_y_im = r * (r * sin(gd_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    hd_y_re = rt_powd_snf(t_s.re, 5.0);
    hd_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    hd_y_re = 0.0;
    hd_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    hd_y_re = 5.0 * x_re;
    hd_y_im = 5.0 * x_im;
    if (hd_y_im == 0.0) {
      hd_y_re = exp(hd_y_re);
      hd_y_im = 0.0;
    } else if (rtIsInf(hd_y_im) && rtIsInf(hd_y_re) && (hd_y_re < 0.0)) {
      hd_y_re = 0.0;
      hd_y_im = 0.0;
    } else {
      r = exp(hd_y_re / 2.0);
      hd_y_re = r * (r * cos(hd_y_im));
      hd_y_im = r * (r * sin(hd_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    id_y_re = rt_powd_snf(t.re, 3.0);
    id_y_im = 0.0;
  } else if (t.re == 0.0) {
    id_y_re = 0.0;
    id_y_im = -rt_powd_snf(t.im, 3.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    id_y_re = 3.0 * x_re;
    id_y_im = 3.0 * x_im;
    if (id_y_im == 0.0) {
      id_y_re = exp(id_y_re);
      id_y_im = 0.0;
    } else if (rtIsInf(id_y_im) && rtIsInf(id_y_re) && (id_y_re < 0.0)) {
      id_y_re = 0.0;
      id_y_im = 0.0;
    } else {
      r = exp(id_y_re / 2.0);
      id_y_re = r * (r * cos(id_y_im));
      id_y_im = r * (r * sin(id_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    jd_y_re = rt_powd_snf(t_s.re, 3.0);
    jd_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    jd_y_re = 0.0;
    jd_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    jd_y_re = 3.0 * x_re;
    jd_y_im = 3.0 * x_im;
    if (jd_y_im == 0.0) {
      jd_y_re = exp(jd_y_re);
      jd_y_im = 0.0;
    } else if (rtIsInf(jd_y_im) && rtIsInf(jd_y_re) && (jd_y_re < 0.0)) {
      jd_y_re = 0.0;
      jd_y_im = 0.0;
    } else {
      r = exp(jd_y_re / 2.0);
      jd_y_re = r * (r * cos(jd_y_im));
      jd_y_im = r * (r * sin(jd_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    kd_y_re = rt_powd_snf(t.re, 3.0);
    kd_y_im = 0.0;
  } else if (t.re == 0.0) {
    kd_y_re = 0.0;
    kd_y_im = -rt_powd_snf(t.im, 3.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    kd_y_re = 3.0 * x_re;
    kd_y_im = 3.0 * x_im;
    if (kd_y_im == 0.0) {
      kd_y_re = exp(kd_y_re);
      kd_y_im = 0.0;
    } else if (rtIsInf(kd_y_im) && rtIsInf(kd_y_re) && (kd_y_re < 0.0)) {
      kd_y_re = 0.0;
      kd_y_im = 0.0;
    } else {
      r = exp(kd_y_re / 2.0);
      kd_y_re = r * (r * cos(kd_y_im));
      kd_y_im = r * (r * sin(kd_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    ld_y_re = rt_powd_snf(t_s.re, 4.0);
    ld_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    ld_y_re = rt_powd_snf(t_s.im, 4.0);
    ld_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    ld_y_re = 4.0 * x_re;
    ld_y_im = 4.0 * x_im;
    if (ld_y_im == 0.0) {
      ld_y_re = exp(ld_y_re);
      ld_y_im = 0.0;
    } else if (rtIsInf(ld_y_im) && rtIsInf(ld_y_re) && (ld_y_re < 0.0)) {
      ld_y_re = 0.0;
      ld_y_im = 0.0;
    } else {
      r = exp(ld_y_re / 2.0);
      ld_y_re = r * (r * cos(ld_y_im));
      ld_y_im = r * (r * sin(ld_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    md_y_re = rt_powd_snf(t_s.re, 3.0);
    md_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    md_y_re = 0.0;
    md_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    md_y_re = 3.0 * x_re;
    md_y_im = 3.0 * x_im;
    if (md_y_im == 0.0) {
      md_y_re = exp(md_y_re);
      md_y_im = 0.0;
    } else if (rtIsInf(md_y_im) && rtIsInf(md_y_re) && (md_y_re < 0.0)) {
      md_y_re = 0.0;
      md_y_im = 0.0;
    } else {
      r = exp(md_y_re / 2.0);
      md_y_re = r * (r * cos(md_y_im));
      md_y_im = r * (r * sin(md_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    nd_y_re = rt_powd_snf(t.re, 4.0);
    nd_y_im = 0.0;
  } else if (t.re == 0.0) {
    nd_y_re = rt_powd_snf(t.im, 4.0);
    nd_y_im = 0.0;
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    nd_y_re = 4.0 * x_re;
    nd_y_im = 4.0 * x_im;
    if (nd_y_im == 0.0) {
      nd_y_re = exp(nd_y_re);
      nd_y_im = 0.0;
    } else if (rtIsInf(nd_y_im) && rtIsInf(nd_y_re) && (nd_y_re < 0.0)) {
      nd_y_re = 0.0;
      nd_y_im = 0.0;
    } else {
      r = exp(nd_y_re / 2.0);
      nd_y_re = r * (r * cos(nd_y_im));
      nd_y_im = r * (r * sin(nd_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    od_y_re = rt_powd_snf(t_s.re, 5.0);
    od_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    od_y_re = 0.0;
    od_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    od_y_re = 5.0 * x_re;
    od_y_im = 5.0 * x_im;
    if (od_y_im == 0.0) {
      od_y_re = exp(od_y_re);
      od_y_im = 0.0;
    } else if (rtIsInf(od_y_im) && rtIsInf(od_y_re) && (od_y_re < 0.0)) {
      od_y_re = 0.0;
      od_y_im = 0.0;
    } else {
      r = exp(od_y_re / 2.0);
      od_y_re = r * (r * cos(od_y_im));
      od_y_im = r * (r * sin(od_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    pd_y_re = rt_powd_snf(t_s.re, 4.0);
    pd_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    pd_y_re = rt_powd_snf(t_s.im, 4.0);
    pd_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    pd_y_re = 4.0 * x_re;
    pd_y_im = 4.0 * x_im;
    if (pd_y_im == 0.0) {
      pd_y_re = exp(pd_y_re);
      pd_y_im = 0.0;
    } else if (rtIsInf(pd_y_im) && rtIsInf(pd_y_re) && (pd_y_re < 0.0)) {
      pd_y_re = 0.0;
      pd_y_im = 0.0;
    } else {
      r = exp(pd_y_re / 2.0);
      pd_y_re = r * (r * cos(pd_y_im));
      pd_y_im = r * (r * sin(pd_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    qd_y_re = rt_powd_snf(t.re, 5.0);
    qd_y_im = 0.0;
  } else if (t.re == 0.0) {
    qd_y_re = 0.0;
    qd_y_im = rt_powd_snf(t.im, 5.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    qd_y_re = 5.0 * x_re;
    qd_y_im = 5.0 * x_im;
    if (qd_y_im == 0.0) {
      qd_y_re = exp(qd_y_re);
      qd_y_im = 0.0;
    } else if (rtIsInf(qd_y_im) && rtIsInf(qd_y_re) && (qd_y_re < 0.0)) {
      qd_y_re = 0.0;
      qd_y_im = 0.0;
    } else {
      r = exp(qd_y_re / 2.0);
      qd_y_re = r * (r * cos(qd_y_im));
      qd_y_im = r * (r * sin(qd_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    rd_y_re = rt_powd_snf(t_s.re, 6.0);
    rd_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    rd_y_re = -rt_powd_snf(t_s.im, 6.0);
    rd_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    rd_y_re = 6.0 * x_re;
    rd_y_im = 6.0 * x_im;
    if (rd_y_im == 0.0) {
      rd_y_re = exp(rd_y_re);
      rd_y_im = 0.0;
    } else if (rtIsInf(rd_y_im) && rtIsInf(rd_y_re) && (rd_y_re < 0.0)) {
      rd_y_re = 0.0;
      rd_y_im = 0.0;
    } else {
      r = exp(rd_y_re / 2.0);
      rd_y_re = r * (r * cos(rd_y_im));
      rd_y_im = r * (r * sin(rd_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    sd_y_re = rt_powd_snf(t_s.re, 5.0);
    sd_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    sd_y_re = 0.0;
    sd_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    sd_y_re = 5.0 * x_re;
    sd_y_im = 5.0 * x_im;
    if (sd_y_im == 0.0) {
      sd_y_re = exp(sd_y_re);
      sd_y_im = 0.0;
    } else if (rtIsInf(sd_y_im) && rtIsInf(sd_y_re) && (sd_y_re < 0.0)) {
      sd_y_re = 0.0;
      sd_y_im = 0.0;
    } else {
      r = exp(sd_y_re / 2.0);
      sd_y_re = r * (r * cos(sd_y_im));
      sd_y_im = r * (r * sin(sd_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    td_y_re = rt_powd_snf(t.re, 3.0);
    td_y_im = 0.0;
  } else if (t.re == 0.0) {
    td_y_re = 0.0;
    td_y_im = -rt_powd_snf(t.im, 3.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    td_y_re = 3.0 * x_re;
    td_y_im = 3.0 * x_im;
    if (td_y_im == 0.0) {
      td_y_re = exp(td_y_re);
      td_y_im = 0.0;
    } else if (rtIsInf(td_y_im) && rtIsInf(td_y_re) && (td_y_re < 0.0)) {
      td_y_re = 0.0;
      td_y_im = 0.0;
    } else {
      r = exp(td_y_re / 2.0);
      td_y_re = r * (r * cos(td_y_im));
      td_y_im = r * (r * sin(td_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    ud_y_re = rt_powd_snf(t_s.re, 3.0);
    ud_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    ud_y_re = 0.0;
    ud_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    ud_y_re = 3.0 * x_re;
    ud_y_im = 3.0 * x_im;
    if (ud_y_im == 0.0) {
      ud_y_re = exp(ud_y_re);
      ud_y_im = 0.0;
    } else if (rtIsInf(ud_y_im) && rtIsInf(ud_y_re) && (ud_y_re < 0.0)) {
      ud_y_re = 0.0;
      ud_y_im = 0.0;
    } else {
      r = exp(ud_y_re / 2.0);
      ud_y_re = r * (r * cos(ud_y_im));
      ud_y_im = r * (r * sin(ud_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    vd_y_re = rt_powd_snf(t.re, 3.0);
    vd_y_im = 0.0;
  } else if (t.re == 0.0) {
    vd_y_re = 0.0;
    vd_y_im = -rt_powd_snf(t.im, 3.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    vd_y_re = 3.0 * x_re;
    vd_y_im = 3.0 * x_im;
    if (vd_y_im == 0.0) {
      vd_y_re = exp(vd_y_re);
      vd_y_im = 0.0;
    } else if (rtIsInf(vd_y_im) && rtIsInf(vd_y_re) && (vd_y_re < 0.0)) {
      vd_y_re = 0.0;
      vd_y_im = 0.0;
    } else {
      r = exp(vd_y_re / 2.0);
      vd_y_re = r * (r * cos(vd_y_im));
      vd_y_im = r * (r * sin(vd_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    wd_y_re = rt_powd_snf(t_s.re, 4.0);
    wd_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    wd_y_re = rt_powd_snf(t_s.im, 4.0);
    wd_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    wd_y_re = 4.0 * x_re;
    wd_y_im = 4.0 * x_im;
    if (wd_y_im == 0.0) {
      wd_y_re = exp(wd_y_re);
      wd_y_im = 0.0;
    } else if (rtIsInf(wd_y_im) && rtIsInf(wd_y_re) && (wd_y_re < 0.0)) {
      wd_y_re = 0.0;
      wd_y_im = 0.0;
    } else {
      r = exp(wd_y_re / 2.0);
      wd_y_re = r * (r * cos(wd_y_im));
      wd_y_im = r * (r * sin(wd_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    xd_y_re = rt_powd_snf(t_s.re, 3.0);
    xd_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    xd_y_re = 0.0;
    xd_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    xd_y_re = 3.0 * x_re;
    xd_y_im = 3.0 * x_im;
    if (xd_y_im == 0.0) {
      xd_y_re = exp(xd_y_re);
      xd_y_im = 0.0;
    } else if (rtIsInf(xd_y_im) && rtIsInf(xd_y_re) && (xd_y_re < 0.0)) {
      xd_y_re = 0.0;
      xd_y_im = 0.0;
    } else {
      r = exp(xd_y_re / 2.0);
      xd_y_re = r * (r * cos(xd_y_im));
      xd_y_im = r * (r * sin(xd_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    yd_y_re = rt_powd_snf(t.re, 4.0);
    yd_y_im = 0.0;
  } else if (t.re == 0.0) {
    yd_y_re = rt_powd_snf(t.im, 4.0);
    yd_y_im = 0.0;
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    yd_y_re = 4.0 * x_re;
    yd_y_im = 4.0 * x_im;
    if (yd_y_im == 0.0) {
      yd_y_re = exp(yd_y_re);
      yd_y_im = 0.0;
    } else if (rtIsInf(yd_y_im) && rtIsInf(yd_y_re) && (yd_y_re < 0.0)) {
      yd_y_re = 0.0;
      yd_y_im = 0.0;
    } else {
      r = exp(yd_y_re / 2.0);
      yd_y_re = r * (r * cos(yd_y_im));
      yd_y_im = r * (r * sin(yd_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    ae_y_re = rt_powd_snf(t_s.re, 5.0);
    ae_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    ae_y_re = 0.0;
    ae_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    ae_y_re = 5.0 * x_re;
    ae_y_im = 5.0 * x_im;
    if (ae_y_im == 0.0) {
      ae_y_re = exp(ae_y_re);
      ae_y_im = 0.0;
    } else if (rtIsInf(ae_y_im) && rtIsInf(ae_y_re) && (ae_y_re < 0.0)) {
      ae_y_re = 0.0;
      ae_y_im = 0.0;
    } else {
      r = exp(ae_y_re / 2.0);
      ae_y_re = r * (r * cos(ae_y_im));
      ae_y_im = r * (r * sin(ae_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    be_y_re = rt_powd_snf(t_s.re, 4.0);
    be_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    be_y_re = rt_powd_snf(t_s.im, 4.0);
    be_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    be_y_re = 4.0 * x_re;
    be_y_im = 4.0 * x_im;
    if (be_y_im == 0.0) {
      be_y_re = exp(be_y_re);
      be_y_im = 0.0;
    } else if (rtIsInf(be_y_im) && rtIsInf(be_y_re) && (be_y_re < 0.0)) {
      be_y_re = 0.0;
      be_y_im = 0.0;
    } else {
      r = exp(be_y_re / 2.0);
      be_y_re = r * (r * cos(be_y_im));
      be_y_im = r * (r * sin(be_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    ce_y_re = rt_powd_snf(t.re, 5.0);
    ce_y_im = 0.0;
  } else if (t.re == 0.0) {
    ce_y_re = 0.0;
    ce_y_im = rt_powd_snf(t.im, 5.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    ce_y_re = 5.0 * x_re;
    ce_y_im = 5.0 * x_im;
    if (ce_y_im == 0.0) {
      ce_y_re = exp(ce_y_re);
      ce_y_im = 0.0;
    } else if (rtIsInf(ce_y_im) && rtIsInf(ce_y_re) && (ce_y_re < 0.0)) {
      ce_y_re = 0.0;
      ce_y_im = 0.0;
    } else {
      r = exp(ce_y_re / 2.0);
      ce_y_re = r * (r * cos(ce_y_im));
      ce_y_im = r * (r * sin(ce_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    de_y_re = rt_powd_snf(t_s.re, 6.0);
    de_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    de_y_re = -rt_powd_snf(t_s.im, 6.0);
    de_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    de_y_re = 6.0 * x_re;
    de_y_im = 6.0 * x_im;
    if (de_y_im == 0.0) {
      de_y_re = exp(de_y_re);
      de_y_im = 0.0;
    } else if (rtIsInf(de_y_im) && rtIsInf(de_y_re) && (de_y_re < 0.0)) {
      de_y_re = 0.0;
      de_y_im = 0.0;
    } else {
      r = exp(de_y_re / 2.0);
      de_y_re = r * (r * cos(de_y_im));
      de_y_im = r * (r * sin(de_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    ee_y_re = rt_powd_snf(t_s.re, 5.0);
    ee_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    ee_y_re = 0.0;
    ee_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    ee_y_re = 5.0 * x_re;
    ee_y_im = 5.0 * x_im;
    if (ee_y_im == 0.0) {
      ee_y_re = exp(ee_y_re);
      ee_y_im = 0.0;
    } else if (rtIsInf(ee_y_im) && rtIsInf(ee_y_re) && (ee_y_re < 0.0)) {
      ee_y_re = 0.0;
      ee_y_im = 0.0;
    } else {
      r = exp(ee_y_re / 2.0);
      ee_y_re = r * (r * cos(ee_y_im));
      ee_y_im = r * (r * sin(ee_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    fe_y_re = rt_powd_snf(t.re, 3.0);
    fe_y_im = 0.0;
  } else if (t.re == 0.0) {
    fe_y_re = 0.0;
    fe_y_im = -rt_powd_snf(t.im, 3.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    fe_y_re = 3.0 * x_re;
    fe_y_im = 3.0 * x_im;
    if (fe_y_im == 0.0) {
      fe_y_re = exp(fe_y_re);
      fe_y_im = 0.0;
    } else if (rtIsInf(fe_y_im) && rtIsInf(fe_y_re) && (fe_y_re < 0.0)) {
      fe_y_re = 0.0;
      fe_y_im = 0.0;
    } else {
      r = exp(fe_y_re / 2.0);
      fe_y_re = r * (r * cos(fe_y_im));
      fe_y_im = r * (r * sin(fe_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    ge_y_re = rt_powd_snf(t_s.re, 4.0);
    ge_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    ge_y_re = rt_powd_snf(t_s.im, 4.0);
    ge_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    ge_y_re = 4.0 * x_re;
    ge_y_im = 4.0 * x_im;
    if (ge_y_im == 0.0) {
      ge_y_re = exp(ge_y_re);
      ge_y_im = 0.0;
    } else if (rtIsInf(ge_y_im) && rtIsInf(ge_y_re) && (ge_y_re < 0.0)) {
      ge_y_re = 0.0;
      ge_y_im = 0.0;
    } else {
      r = exp(ge_y_re / 2.0);
      ge_y_re = r * (r * cos(ge_y_im));
      ge_y_im = r * (r * sin(ge_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    he_y_re = rt_powd_snf(t_s.re, 3.0);
    he_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    he_y_re = 0.0;
    he_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    he_y_re = 3.0 * x_re;
    he_y_im = 3.0 * x_im;
    if (he_y_im == 0.0) {
      he_y_re = exp(he_y_re);
      he_y_im = 0.0;
    } else if (rtIsInf(he_y_im) && rtIsInf(he_y_re) && (he_y_re < 0.0)) {
      he_y_re = 0.0;
      he_y_im = 0.0;
    } else {
      r = exp(he_y_re / 2.0);
      he_y_re = r * (r * cos(he_y_im));
      he_y_im = r * (r * sin(he_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    ie_y_re = rt_powd_snf(t.re, 3.0);
    ie_y_im = 0.0;
  } else if (t.re == 0.0) {
    ie_y_re = 0.0;
    ie_y_im = -rt_powd_snf(t.im, 3.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    ie_y_re = 3.0 * x_re;
    ie_y_im = 3.0 * x_im;
    if (ie_y_im == 0.0) {
      ie_y_re = exp(ie_y_re);
      ie_y_im = 0.0;
    } else if (rtIsInf(ie_y_im) && rtIsInf(ie_y_re) && (ie_y_re < 0.0)) {
      ie_y_re = 0.0;
      ie_y_im = 0.0;
    } else {
      r = exp(ie_y_re / 2.0);
      ie_y_re = r * (r * cos(ie_y_im));
      ie_y_im = r * (r * sin(ie_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    je_y_re = rt_powd_snf(t_s.re, 3.0);
    je_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    je_y_re = 0.0;
    je_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    je_y_re = 3.0 * x_re;
    je_y_im = 3.0 * x_im;
    if (je_y_im == 0.0) {
      je_y_re = exp(je_y_re);
      je_y_im = 0.0;
    } else if (rtIsInf(je_y_im) && rtIsInf(je_y_re) && (je_y_re < 0.0)) {
      je_y_re = 0.0;
      je_y_im = 0.0;
    } else {
      r = exp(je_y_re / 2.0);
      je_y_re = r * (r * cos(je_y_im));
      je_y_im = r * (r * sin(je_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    ke_y_re = rt_powd_snf(t.re, 3.0);
    ke_y_im = 0.0;
  } else if (t.re == 0.0) {
    ke_y_re = 0.0;
    ke_y_im = -rt_powd_snf(t.im, 3.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    ke_y_re = 3.0 * x_re;
    ke_y_im = 3.0 * x_im;
    if (ke_y_im == 0.0) {
      ke_y_re = exp(ke_y_re);
      ke_y_im = 0.0;
    } else if (rtIsInf(ke_y_im) && rtIsInf(ke_y_re) && (ke_y_re < 0.0)) {
      ke_y_re = 0.0;
      ke_y_im = 0.0;
    } else {
      r = exp(ke_y_re / 2.0);
      ke_y_re = r * (r * cos(ke_y_im));
      ke_y_im = r * (r * sin(ke_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    le_y_re = rt_powd_snf(t_s.re, 4.0);
    le_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    le_y_re = rt_powd_snf(t_s.im, 4.0);
    le_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    le_y_re = 4.0 * x_re;
    le_y_im = 4.0 * x_im;
    if (le_y_im == 0.0) {
      le_y_re = exp(le_y_re);
      le_y_im = 0.0;
    } else if (rtIsInf(le_y_im) && rtIsInf(le_y_re) && (le_y_re < 0.0)) {
      le_y_re = 0.0;
      le_y_im = 0.0;
    } else {
      r = exp(le_y_re / 2.0);
      le_y_re = r * (r * cos(le_y_im));
      le_y_im = r * (r * sin(le_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    me_y_re = rt_powd_snf(t_s.re, 3.0);
    me_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    me_y_re = 0.0;
    me_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    me_y_re = 3.0 * x_re;
    me_y_im = 3.0 * x_im;
    if (me_y_im == 0.0) {
      me_y_re = exp(me_y_re);
      me_y_im = 0.0;
    } else if (rtIsInf(me_y_im) && rtIsInf(me_y_re) && (me_y_re < 0.0)) {
      me_y_re = 0.0;
      me_y_im = 0.0;
    } else {
      r = exp(me_y_re / 2.0);
      me_y_re = r * (r * cos(me_y_im));
      me_y_im = r * (r * sin(me_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    ne_y_re = rt_powd_snf(t.re, 4.0);
    ne_y_im = 0.0;
  } else if (t.re == 0.0) {
    ne_y_re = rt_powd_snf(t.im, 4.0);
    ne_y_im = 0.0;
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    ne_y_re = 4.0 * x_re;
    ne_y_im = 4.0 * x_im;
    if (ne_y_im == 0.0) {
      ne_y_re = exp(ne_y_re);
      ne_y_im = 0.0;
    } else if (rtIsInf(ne_y_im) && rtIsInf(ne_y_re) && (ne_y_re < 0.0)) {
      ne_y_re = 0.0;
      ne_y_im = 0.0;
    } else {
      r = exp(ne_y_re / 2.0);
      ne_y_re = r * (r * cos(ne_y_im));
      ne_y_im = r * (r * sin(ne_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    oe_y_re = rt_powd_snf(t_s.re, 5.0);
    oe_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    oe_y_re = 0.0;
    oe_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    oe_y_re = 5.0 * x_re;
    oe_y_im = 5.0 * x_im;
    if (oe_y_im == 0.0) {
      oe_y_re = exp(oe_y_re);
      oe_y_im = 0.0;
    } else if (rtIsInf(oe_y_im) && rtIsInf(oe_y_re) && (oe_y_re < 0.0)) {
      oe_y_re = 0.0;
      oe_y_im = 0.0;
    } else {
      r = exp(oe_y_re / 2.0);
      oe_y_re = r * (r * cos(oe_y_im));
      oe_y_im = r * (r * sin(oe_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    pe_y_re = rt_powd_snf(t_s.re, 4.0);
    pe_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    pe_y_re = rt_powd_snf(t_s.im, 4.0);
    pe_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    pe_y_re = 4.0 * x_re;
    pe_y_im = 4.0 * x_im;
    if (pe_y_im == 0.0) {
      pe_y_re = exp(pe_y_re);
      pe_y_im = 0.0;
    } else if (rtIsInf(pe_y_im) && rtIsInf(pe_y_re) && (pe_y_re < 0.0)) {
      pe_y_re = 0.0;
      pe_y_im = 0.0;
    } else {
      r = exp(pe_y_re / 2.0);
      pe_y_re = r * (r * cos(pe_y_im));
      pe_y_im = r * (r * sin(pe_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    qe_y_re = rt_powd_snf(t.re, 4.0);
    qe_y_im = 0.0;
  } else if (t.re == 0.0) {
    qe_y_re = rt_powd_snf(t.im, 4.0);
    qe_y_im = 0.0;
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    qe_y_re = 4.0 * x_re;
    qe_y_im = 4.0 * x_im;
    if (qe_y_im == 0.0) {
      qe_y_re = exp(qe_y_re);
      qe_y_im = 0.0;
    } else if (rtIsInf(qe_y_im) && rtIsInf(qe_y_re) && (qe_y_re < 0.0)) {
      qe_y_re = 0.0;
      qe_y_im = 0.0;
    } else {
      r = exp(qe_y_re / 2.0);
      qe_y_re = r * (r * cos(qe_y_im));
      qe_y_im = r * (r * sin(qe_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    re_y_re = rt_powd_snf(t_s.re, 5.0);
    re_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    re_y_re = 0.0;
    re_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    re_y_re = 5.0 * x_re;
    re_y_im = 5.0 * x_im;
    if (re_y_im == 0.0) {
      re_y_re = exp(re_y_re);
      re_y_im = 0.0;
    } else if (rtIsInf(re_y_im) && rtIsInf(re_y_re) && (re_y_re < 0.0)) {
      re_y_re = 0.0;
      re_y_im = 0.0;
    } else {
      r = exp(re_y_re / 2.0);
      re_y_re = r * (r * cos(re_y_im));
      re_y_im = r * (r * sin(re_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    se_y_re = rt_powd_snf(t_s.re, 4.0);
    se_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    se_y_re = rt_powd_snf(t_s.im, 4.0);
    se_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    se_y_re = 4.0 * x_re;
    se_y_im = 4.0 * x_im;
    if (se_y_im == 0.0) {
      se_y_re = exp(se_y_re);
      se_y_im = 0.0;
    } else if (rtIsInf(se_y_im) && rtIsInf(se_y_re) && (se_y_re < 0.0)) {
      se_y_re = 0.0;
      se_y_im = 0.0;
    } else {
      r = exp(se_y_re / 2.0);
      se_y_re = r * (r * cos(se_y_im));
      se_y_im = r * (r * sin(se_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    te_y_re = rt_powd_snf(t.re, 5.0);
    te_y_im = 0.0;
  } else if (t.re == 0.0) {
    te_y_re = 0.0;
    te_y_im = rt_powd_snf(t.im, 5.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    te_y_re = 5.0 * x_re;
    te_y_im = 5.0 * x_im;
    if (te_y_im == 0.0) {
      te_y_re = exp(te_y_re);
      te_y_im = 0.0;
    } else if (rtIsInf(te_y_im) && rtIsInf(te_y_re) && (te_y_re < 0.0)) {
      te_y_re = 0.0;
      te_y_im = 0.0;
    } else {
      r = exp(te_y_re / 2.0);
      te_y_re = r * (r * cos(te_y_im));
      te_y_im = r * (r * sin(te_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    ue_y_re = rt_powd_snf(t_s.re, 6.0);
    ue_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    ue_y_re = -rt_powd_snf(t_s.im, 6.0);
    ue_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    ue_y_re = 6.0 * x_re;
    ue_y_im = 6.0 * x_im;
    if (ue_y_im == 0.0) {
      ue_y_re = exp(ue_y_re);
      ue_y_im = 0.0;
    } else if (rtIsInf(ue_y_im) && rtIsInf(ue_y_re) && (ue_y_re < 0.0)) {
      ue_y_re = 0.0;
      ue_y_im = 0.0;
    } else {
      r = exp(ue_y_re / 2.0);
      ue_y_re = r * (r * cos(ue_y_im));
      ue_y_im = r * (r * sin(ue_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    ve_y_re = rt_powd_snf(t_s.re, 5.0);
    ve_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    ve_y_re = 0.0;
    ve_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    ve_y_re = 5.0 * x_re;
    ve_y_im = 5.0 * x_im;
    if (ve_y_im == 0.0) {
      ve_y_re = exp(ve_y_re);
      ve_y_im = 0.0;
    } else if (rtIsInf(ve_y_im) && rtIsInf(ve_y_re) && (ve_y_re < 0.0)) {
      ve_y_re = 0.0;
      ve_y_im = 0.0;
    } else {
      r = exp(ve_y_re / 2.0);
      ve_y_re = r * (r * cos(ve_y_im));
      ve_y_im = r * (r * sin(ve_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    we_y_re = rt_powd_snf(t.re, 5.0);
    we_y_im = 0.0;
  } else if (t.re == 0.0) {
    we_y_re = 0.0;
    we_y_im = rt_powd_snf(t.im, 5.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    we_y_re = 5.0 * x_re;
    we_y_im = 5.0 * x_im;
    if (we_y_im == 0.0) {
      we_y_re = exp(we_y_re);
      we_y_im = 0.0;
    } else if (rtIsInf(we_y_im) && rtIsInf(we_y_re) && (we_y_re < 0.0)) {
      we_y_re = 0.0;
      we_y_im = 0.0;
    } else {
      r = exp(we_y_re / 2.0);
      we_y_re = r * (r * cos(we_y_im));
      we_y_im = r * (r * sin(we_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    xe_y_re = rt_powd_snf(t_s.re, 6.0);
    xe_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    xe_y_re = -rt_powd_snf(t_s.im, 6.0);
    xe_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    xe_y_re = 6.0 * x_re;
    xe_y_im = 6.0 * x_im;
    if (xe_y_im == 0.0) {
      xe_y_re = exp(xe_y_re);
      xe_y_im = 0.0;
    } else if (rtIsInf(xe_y_im) && rtIsInf(xe_y_re) && (xe_y_re < 0.0)) {
      xe_y_re = 0.0;
      xe_y_im = 0.0;
    } else {
      r = exp(xe_y_re / 2.0);
      xe_y_re = r * (r * cos(xe_y_im));
      xe_y_im = r * (r * sin(xe_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    ye_y_re = rt_powd_snf(t_s.re, 5.0);
    ye_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    ye_y_re = 0.0;
    ye_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    ye_y_re = 5.0 * x_re;
    ye_y_im = 5.0 * x_im;
    if (ye_y_im == 0.0) {
      ye_y_re = exp(ye_y_re);
      ye_y_im = 0.0;
    } else if (rtIsInf(ye_y_im) && rtIsInf(ye_y_re) && (ye_y_re < 0.0)) {
      ye_y_re = 0.0;
      ye_y_im = 0.0;
    } else {
      r = exp(ye_y_re / 2.0);
      ye_y_re = r * (r * cos(ye_y_im));
      ye_y_im = r * (r * sin(ye_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    af_y_re = rt_powd_snf(t.re, 3.0);
    af_y_im = 0.0;
  } else if (t.re == 0.0) {
    af_y_re = 0.0;
    af_y_im = -rt_powd_snf(t.im, 3.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    af_y_re = 3.0 * x_re;
    af_y_im = 3.0 * x_im;
    if (af_y_im == 0.0) {
      af_y_re = exp(af_y_re);
      af_y_im = 0.0;
    } else if (rtIsInf(af_y_im) && rtIsInf(af_y_re) && (af_y_re < 0.0)) {
      af_y_re = 0.0;
      af_y_im = 0.0;
    } else {
      r = exp(af_y_re / 2.0);
      af_y_re = r * (r * cos(af_y_im));
      af_y_im = r * (r * sin(af_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    bf_y_re = rt_powd_snf(t_s.re, 4.0);
    bf_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    bf_y_re = rt_powd_snf(t_s.im, 4.0);
    bf_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    bf_y_re = 4.0 * x_re;
    bf_y_im = 4.0 * x_im;
    if (bf_y_im == 0.0) {
      bf_y_re = exp(bf_y_re);
      bf_y_im = 0.0;
    } else if (rtIsInf(bf_y_im) && rtIsInf(bf_y_re) && (bf_y_re < 0.0)) {
      bf_y_re = 0.0;
      bf_y_im = 0.0;
    } else {
      r = exp(bf_y_re / 2.0);
      bf_y_re = r * (r * cos(bf_y_im));
      bf_y_im = r * (r * sin(bf_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    cf_y_re = rt_powd_snf(t_s.re, 3.0);
    cf_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    cf_y_re = 0.0;
    cf_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    cf_y_re = 3.0 * x_re;
    cf_y_im = 3.0 * x_im;
    if (cf_y_im == 0.0) {
      cf_y_re = exp(cf_y_re);
      cf_y_im = 0.0;
    } else if (rtIsInf(cf_y_im) && rtIsInf(cf_y_re) && (cf_y_re < 0.0)) {
      cf_y_re = 0.0;
      cf_y_im = 0.0;
    } else {
      r = exp(cf_y_re / 2.0);
      cf_y_re = r * (r * cos(cf_y_im));
      cf_y_im = r * (r * sin(cf_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    df_y_re = rt_powd_snf(t.re, 4.0);
    df_y_im = 0.0;
  } else if (t.re == 0.0) {
    df_y_re = rt_powd_snf(t.im, 4.0);
    df_y_im = 0.0;
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    df_y_re = 4.0 * x_re;
    df_y_im = 4.0 * x_im;
    if (df_y_im == 0.0) {
      df_y_re = exp(df_y_re);
      df_y_im = 0.0;
    } else if (rtIsInf(df_y_im) && rtIsInf(df_y_re) && (df_y_re < 0.0)) {
      df_y_re = 0.0;
      df_y_im = 0.0;
    } else {
      r = exp(df_y_re / 2.0);
      df_y_re = r * (r * cos(df_y_im));
      df_y_im = r * (r * sin(df_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    ef_y_re = rt_powd_snf(t_s.re, 5.0);
    ef_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    ef_y_re = 0.0;
    ef_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    ef_y_re = 5.0 * x_re;
    ef_y_im = 5.0 * x_im;
    if (ef_y_im == 0.0) {
      ef_y_re = exp(ef_y_re);
      ef_y_im = 0.0;
    } else if (rtIsInf(ef_y_im) && rtIsInf(ef_y_re) && (ef_y_re < 0.0)) {
      ef_y_re = 0.0;
      ef_y_im = 0.0;
    } else {
      r = exp(ef_y_re / 2.0);
      ef_y_re = r * (r * cos(ef_y_im));
      ef_y_im = r * (r * sin(ef_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    ff_y_re = rt_powd_snf(t_s.re, 4.0);
    ff_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    ff_y_re = rt_powd_snf(t_s.im, 4.0);
    ff_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    ff_y_re = 4.0 * x_re;
    ff_y_im = 4.0 * x_im;
    if (ff_y_im == 0.0) {
      ff_y_re = exp(ff_y_re);
      ff_y_im = 0.0;
    } else if (rtIsInf(ff_y_im) && rtIsInf(ff_y_re) && (ff_y_re < 0.0)) {
      ff_y_re = 0.0;
      ff_y_im = 0.0;
    } else {
      r = exp(ff_y_re / 2.0);
      ff_y_re = r * (r * cos(ff_y_im));
      ff_y_im = r * (r * sin(ff_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    gf_y_re = rt_powd_snf(t.re, 5.0);
    gf_y_im = 0.0;
  } else if (t.re == 0.0) {
    gf_y_re = 0.0;
    gf_y_im = rt_powd_snf(t.im, 5.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    gf_y_re = 5.0 * x_re;
    gf_y_im = 5.0 * x_im;
    if (gf_y_im == 0.0) {
      gf_y_re = exp(gf_y_re);
      gf_y_im = 0.0;
    } else if (rtIsInf(gf_y_im) && rtIsInf(gf_y_re) && (gf_y_re < 0.0)) {
      gf_y_re = 0.0;
      gf_y_im = 0.0;
    } else {
      r = exp(gf_y_re / 2.0);
      gf_y_re = r * (r * cos(gf_y_im));
      gf_y_im = r * (r * sin(gf_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    hf_y_re = rt_powd_snf(t_s.re, 6.0);
    hf_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    hf_y_re = -rt_powd_snf(t_s.im, 6.0);
    hf_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    hf_y_re = 6.0 * x_re;
    hf_y_im = 6.0 * x_im;
    if (hf_y_im == 0.0) {
      hf_y_re = exp(hf_y_re);
      hf_y_im = 0.0;
    } else if (rtIsInf(hf_y_im) && rtIsInf(hf_y_re) && (hf_y_re < 0.0)) {
      hf_y_re = 0.0;
      hf_y_im = 0.0;
    } else {
      r = exp(hf_y_re / 2.0);
      hf_y_re = r * (r * cos(hf_y_im));
      hf_y_im = r * (r * sin(hf_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    if_y_re = rt_powd_snf(t_s.re, 5.0);
    if_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    if_y_re = 0.0;
    if_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    if_y_re = 5.0 * x_re;
    if_y_im = 5.0 * x_im;
    if (if_y_im == 0.0) {
      if_y_re = exp(if_y_re);
      if_y_im = 0.0;
    } else if (rtIsInf(if_y_im) && rtIsInf(if_y_re) && (if_y_re < 0.0)) {
      if_y_re = 0.0;
      if_y_im = 0.0;
    } else {
      r = exp(if_y_re / 2.0);
      if_y_re = r * (r * cos(if_y_im));
      if_y_im = r * (r * sin(if_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    jf_y_re = rt_powd_snf(t.re, 3.0);
    jf_y_im = 0.0;
  } else if (t.re == 0.0) {
    jf_y_re = 0.0;
    jf_y_im = -rt_powd_snf(t.im, 3.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    jf_y_re = 3.0 * x_re;
    jf_y_im = 3.0 * x_im;
    if (jf_y_im == 0.0) {
      jf_y_re = exp(jf_y_re);
      jf_y_im = 0.0;
    } else if (rtIsInf(jf_y_im) && rtIsInf(jf_y_re) && (jf_y_re < 0.0)) {
      jf_y_re = 0.0;
      jf_y_im = 0.0;
    } else {
      r = exp(jf_y_re / 2.0);
      jf_y_re = r * (r * cos(jf_y_im));
      jf_y_im = r * (r * sin(jf_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    kf_y_re = rt_powd_snf(t_s.re, 4.0);
    kf_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    kf_y_re = rt_powd_snf(t_s.im, 4.0);
    kf_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    kf_y_re = 4.0 * x_re;
    kf_y_im = 4.0 * x_im;
    if (kf_y_im == 0.0) {
      kf_y_re = exp(kf_y_re);
      kf_y_im = 0.0;
    } else if (rtIsInf(kf_y_im) && rtIsInf(kf_y_re) && (kf_y_re < 0.0)) {
      kf_y_re = 0.0;
      kf_y_im = 0.0;
    } else {
      r = exp(kf_y_re / 2.0);
      kf_y_re = r * (r * cos(kf_y_im));
      kf_y_im = r * (r * sin(kf_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    lf_y_re = rt_powd_snf(t_s.re, 3.0);
    lf_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    lf_y_re = 0.0;
    lf_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    lf_y_re = 3.0 * x_re;
    lf_y_im = 3.0 * x_im;
    if (lf_y_im == 0.0) {
      lf_y_re = exp(lf_y_re);
      lf_y_im = 0.0;
    } else if (rtIsInf(lf_y_im) && rtIsInf(lf_y_re) && (lf_y_re < 0.0)) {
      lf_y_re = 0.0;
      lf_y_im = 0.0;
    } else {
      r = exp(lf_y_re / 2.0);
      lf_y_re = r * (r * cos(lf_y_im));
      lf_y_im = r * (r * sin(lf_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    mf_y_re = rt_powd_snf(t.re, 3.0);
    mf_y_im = 0.0;
  } else if (t.re == 0.0) {
    mf_y_re = 0.0;
    mf_y_im = -rt_powd_snf(t.im, 3.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    mf_y_re = 3.0 * x_re;
    mf_y_im = 3.0 * x_im;
    if (mf_y_im == 0.0) {
      mf_y_re = exp(mf_y_re);
      mf_y_im = 0.0;
    } else if (rtIsInf(mf_y_im) && rtIsInf(mf_y_re) && (mf_y_re < 0.0)) {
      mf_y_re = 0.0;
      mf_y_im = 0.0;
    } else {
      r = exp(mf_y_re / 2.0);
      mf_y_re = r * (r * cos(mf_y_im));
      mf_y_im = r * (r * sin(mf_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    nf_y_re = rt_powd_snf(t_s.re, 3.0);
    nf_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    nf_y_re = 0.0;
    nf_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    nf_y_re = 3.0 * x_re;
    nf_y_im = 3.0 * x_im;
    if (nf_y_im == 0.0) {
      nf_y_re = exp(nf_y_re);
      nf_y_im = 0.0;
    } else if (rtIsInf(nf_y_im) && rtIsInf(nf_y_re) && (nf_y_re < 0.0)) {
      nf_y_re = 0.0;
      nf_y_im = 0.0;
    } else {
      r = exp(nf_y_re / 2.0);
      nf_y_re = r * (r * cos(nf_y_im));
      nf_y_im = r * (r * sin(nf_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    of_y_re = rt_powd_snf(t.re, 3.0);
    of_y_im = 0.0;
  } else if (t.re == 0.0) {
    of_y_re = 0.0;
    of_y_im = -rt_powd_snf(t.im, 3.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    of_y_re = 3.0 * x_re;
    of_y_im = 3.0 * x_im;
    if (of_y_im == 0.0) {
      of_y_re = exp(of_y_re);
      of_y_im = 0.0;
    } else if (rtIsInf(of_y_im) && rtIsInf(of_y_re) && (of_y_re < 0.0)) {
      of_y_re = 0.0;
      of_y_im = 0.0;
    } else {
      r = exp(of_y_re / 2.0);
      of_y_re = r * (r * cos(of_y_im));
      of_y_im = r * (r * sin(of_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    pf_y_re = rt_powd_snf(t_s.re, 4.0);
    pf_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    pf_y_re = rt_powd_snf(t_s.im, 4.0);
    pf_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    pf_y_re = 4.0 * x_re;
    pf_y_im = 4.0 * x_im;
    if (pf_y_im == 0.0) {
      pf_y_re = exp(pf_y_re);
      pf_y_im = 0.0;
    } else if (rtIsInf(pf_y_im) && rtIsInf(pf_y_re) && (pf_y_re < 0.0)) {
      pf_y_re = 0.0;
      pf_y_im = 0.0;
    } else {
      r = exp(pf_y_re / 2.0);
      pf_y_re = r * (r * cos(pf_y_im));
      pf_y_im = r * (r * sin(pf_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    qf_y_re = rt_powd_snf(t_s.re, 3.0);
    qf_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    qf_y_re = 0.0;
    qf_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    qf_y_re = 3.0 * x_re;
    qf_y_im = 3.0 * x_im;
    if (qf_y_im == 0.0) {
      qf_y_re = exp(qf_y_re);
      qf_y_im = 0.0;
    } else if (rtIsInf(qf_y_im) && rtIsInf(qf_y_re) && (qf_y_re < 0.0)) {
      qf_y_re = 0.0;
      qf_y_im = 0.0;
    } else {
      r = exp(qf_y_re / 2.0);
      qf_y_re = r * (r * cos(qf_y_im));
      qf_y_im = r * (r * sin(qf_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    rf_y_re = rt_powd_snf(t.re, 4.0);
    rf_y_im = 0.0;
  } else if (t.re == 0.0) {
    rf_y_re = rt_powd_snf(t.im, 4.0);
    rf_y_im = 0.0;
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    rf_y_re = 4.0 * x_re;
    rf_y_im = 4.0 * x_im;
    if (rf_y_im == 0.0) {
      rf_y_re = exp(rf_y_re);
      rf_y_im = 0.0;
    } else if (rtIsInf(rf_y_im) && rtIsInf(rf_y_re) && (rf_y_re < 0.0)) {
      rf_y_re = 0.0;
      rf_y_im = 0.0;
    } else {
      r = exp(rf_y_re / 2.0);
      rf_y_re = r * (r * cos(rf_y_im));
      rf_y_im = r * (r * sin(rf_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    sf_y_re = rt_powd_snf(t_s.re, 5.0);
    sf_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    sf_y_re = 0.0;
    sf_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    sf_y_re = 5.0 * x_re;
    sf_y_im = 5.0 * x_im;
    if (sf_y_im == 0.0) {
      sf_y_re = exp(sf_y_re);
      sf_y_im = 0.0;
    } else if (rtIsInf(sf_y_im) && rtIsInf(sf_y_re) && (sf_y_re < 0.0)) {
      sf_y_re = 0.0;
      sf_y_im = 0.0;
    } else {
      r = exp(sf_y_re / 2.0);
      sf_y_re = r * (r * cos(sf_y_im));
      sf_y_im = r * (r * sin(sf_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    tf_y_re = rt_powd_snf(t_s.re, 4.0);
    tf_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    tf_y_re = rt_powd_snf(t_s.im, 4.0);
    tf_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    tf_y_re = 4.0 * x_re;
    tf_y_im = 4.0 * x_im;
    if (tf_y_im == 0.0) {
      tf_y_re = exp(tf_y_re);
      tf_y_im = 0.0;
    } else if (rtIsInf(tf_y_im) && rtIsInf(tf_y_re) && (tf_y_re < 0.0)) {
      tf_y_re = 0.0;
      tf_y_im = 0.0;
    } else {
      r = exp(tf_y_re / 2.0);
      tf_y_re = r * (r * cos(tf_y_im));
      tf_y_im = r * (r * sin(tf_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    uf_y_re = rt_powd_snf(t.re, 3.0);
    uf_y_im = 0.0;
  } else if (t.re == 0.0) {
    uf_y_re = 0.0;
    uf_y_im = -rt_powd_snf(t.im, 3.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    uf_y_re = 3.0 * x_re;
    uf_y_im = 3.0 * x_im;
    if (uf_y_im == 0.0) {
      uf_y_re = exp(uf_y_re);
      uf_y_im = 0.0;
    } else if (rtIsInf(uf_y_im) && rtIsInf(uf_y_re) && (uf_y_re < 0.0)) {
      uf_y_re = 0.0;
      uf_y_im = 0.0;
    } else {
      r = exp(uf_y_re / 2.0);
      uf_y_re = r * (r * cos(uf_y_im));
      uf_y_im = r * (r * sin(uf_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    vf_y_re = rt_powd_snf(t_s.re, 3.0);
    vf_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    vf_y_re = 0.0;
    vf_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    vf_y_re = 3.0 * x_re;
    vf_y_im = 3.0 * x_im;
    if (vf_y_im == 0.0) {
      vf_y_re = exp(vf_y_re);
      vf_y_im = 0.0;
    } else if (rtIsInf(vf_y_im) && rtIsInf(vf_y_re) && (vf_y_re < 0.0)) {
      vf_y_re = 0.0;
      vf_y_im = 0.0;
    } else {
      r = exp(vf_y_re / 2.0);
      vf_y_re = r * (r * cos(vf_y_im));
      vf_y_im = r * (r * sin(vf_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    wf_y_re = rt_powd_snf(t_s.re, 4.0);
    wf_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    wf_y_re = rt_powd_snf(t_s.im, 4.0);
    wf_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    wf_y_re = 4.0 * x_re;
    wf_y_im = 4.0 * x_im;
    if (wf_y_im == 0.0) {
      wf_y_re = exp(wf_y_re);
      wf_y_im = 0.0;
    } else if (rtIsInf(wf_y_im) && rtIsInf(wf_y_re) && (wf_y_re < 0.0)) {
      wf_y_re = 0.0;
      wf_y_im = 0.0;
    } else {
      r = exp(wf_y_re / 2.0);
      wf_y_re = r * (r * cos(wf_y_im));
      wf_y_im = r * (r * sin(wf_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    xf_y_re = rt_powd_snf(t_s.re, 3.0);
    xf_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    xf_y_re = 0.0;
    xf_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    xf_y_re = 3.0 * x_re;
    xf_y_im = 3.0 * x_im;
    if (xf_y_im == 0.0) {
      xf_y_re = exp(xf_y_re);
      xf_y_im = 0.0;
    } else if (rtIsInf(xf_y_im) && rtIsInf(xf_y_re) && (xf_y_re < 0.0)) {
      xf_y_re = 0.0;
      xf_y_im = 0.0;
    } else {
      r = exp(xf_y_re / 2.0);
      xf_y_re = r * (r * cos(xf_y_im));
      xf_y_im = r * (r * sin(xf_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    yf_y_re = rt_powd_snf(t.re, 4.0);
    yf_y_im = 0.0;
  } else if (t.re == 0.0) {
    yf_y_re = rt_powd_snf(t.im, 4.0);
    yf_y_im = 0.0;
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    yf_y_re = 4.0 * x_re;
    yf_y_im = 4.0 * x_im;
    if (yf_y_im == 0.0) {
      yf_y_re = exp(yf_y_re);
      yf_y_im = 0.0;
    } else if (rtIsInf(yf_y_im) && rtIsInf(yf_y_re) && (yf_y_re < 0.0)) {
      yf_y_re = 0.0;
      yf_y_im = 0.0;
    } else {
      r = exp(yf_y_re / 2.0);
      yf_y_re = r * (r * cos(yf_y_im));
      yf_y_im = r * (r * sin(yf_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    ag_y_re = rt_powd_snf(t_s.re, 5.0);
    ag_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    ag_y_re = 0.0;
    ag_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    ag_y_re = 5.0 * x_re;
    ag_y_im = 5.0 * x_im;
    if (ag_y_im == 0.0) {
      ag_y_re = exp(ag_y_re);
      ag_y_im = 0.0;
    } else if (rtIsInf(ag_y_im) && rtIsInf(ag_y_re) && (ag_y_re < 0.0)) {
      ag_y_re = 0.0;
      ag_y_im = 0.0;
    } else {
      r = exp(ag_y_re / 2.0);
      ag_y_re = r * (r * cos(ag_y_im));
      ag_y_im = r * (r * sin(ag_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    bg_y_re = rt_powd_snf(t_s.re, 4.0);
    bg_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    bg_y_re = rt_powd_snf(t_s.im, 4.0);
    bg_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    bg_y_re = 4.0 * x_re;
    bg_y_im = 4.0 * x_im;
    if (bg_y_im == 0.0) {
      bg_y_re = exp(bg_y_re);
      bg_y_im = 0.0;
    } else if (rtIsInf(bg_y_im) && rtIsInf(bg_y_re) && (bg_y_re < 0.0)) {
      bg_y_re = 0.0;
      bg_y_im = 0.0;
    } else {
      r = exp(bg_y_re / 2.0);
      bg_y_re = r * (r * cos(bg_y_im));
      bg_y_im = r * (r * sin(bg_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    cg_y_re = rt_powd_snf(t.re, 5.0);
    cg_y_im = 0.0;
  } else if (t.re == 0.0) {
    cg_y_re = 0.0;
    cg_y_im = rt_powd_snf(t.im, 5.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    cg_y_re = 5.0 * x_re;
    cg_y_im = 5.0 * x_im;
    if (cg_y_im == 0.0) {
      cg_y_re = exp(cg_y_re);
      cg_y_im = 0.0;
    } else if (rtIsInf(cg_y_im) && rtIsInf(cg_y_re) && (cg_y_re < 0.0)) {
      cg_y_re = 0.0;
      cg_y_im = 0.0;
    } else {
      r = exp(cg_y_re / 2.0);
      cg_y_re = r * (r * cos(cg_y_im));
      cg_y_im = r * (r * sin(cg_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    dg_y_re = rt_powd_snf(t_s.re, 6.0);
    dg_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    dg_y_re = -rt_powd_snf(t_s.im, 6.0);
    dg_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    dg_y_re = 6.0 * x_re;
    dg_y_im = 6.0 * x_im;
    if (dg_y_im == 0.0) {
      dg_y_re = exp(dg_y_re);
      dg_y_im = 0.0;
    } else if (rtIsInf(dg_y_im) && rtIsInf(dg_y_re) && (dg_y_re < 0.0)) {
      dg_y_re = 0.0;
      dg_y_im = 0.0;
    } else {
      r = exp(dg_y_re / 2.0);
      dg_y_re = r * (r * cos(dg_y_im));
      dg_y_im = r * (r * sin(dg_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    eg_y_re = rt_powd_snf(t_s.re, 5.0);
    eg_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    eg_y_re = 0.0;
    eg_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    eg_y_re = 5.0 * x_re;
    eg_y_im = 5.0 * x_im;
    if (eg_y_im == 0.0) {
      eg_y_re = exp(eg_y_re);
      eg_y_im = 0.0;
    } else if (rtIsInf(eg_y_im) && rtIsInf(eg_y_re) && (eg_y_re < 0.0)) {
      eg_y_re = 0.0;
      eg_y_im = 0.0;
    } else {
      r = exp(eg_y_re / 2.0);
      eg_y_re = r * (r * cos(eg_y_im));
      eg_y_im = r * (r * sin(eg_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    fg_y_re = rt_powd_snf(t.re, 4.0);
    fg_y_im = 0.0;
  } else if (t.re == 0.0) {
    fg_y_re = rt_powd_snf(t.im, 4.0);
    fg_y_im = 0.0;
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    fg_y_re = 4.0 * x_re;
    fg_y_im = 4.0 * x_im;
    if (fg_y_im == 0.0) {
      fg_y_re = exp(fg_y_re);
      fg_y_im = 0.0;
    } else if (rtIsInf(fg_y_im) && rtIsInf(fg_y_re) && (fg_y_re < 0.0)) {
      fg_y_re = 0.0;
      fg_y_im = 0.0;
    } else {
      r = exp(fg_y_re / 2.0);
      fg_y_re = r * (r * cos(fg_y_im));
      fg_y_im = r * (r * sin(fg_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    gg_y_re = rt_powd_snf(t_s.re, 3.0);
    gg_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    gg_y_re = 0.0;
    gg_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    gg_y_re = 3.0 * x_re;
    gg_y_im = 3.0 * x_im;
    if (gg_y_im == 0.0) {
      gg_y_re = exp(gg_y_re);
      gg_y_im = 0.0;
    } else if (rtIsInf(gg_y_im) && rtIsInf(gg_y_re) && (gg_y_re < 0.0)) {
      gg_y_re = 0.0;
      gg_y_im = 0.0;
    } else {
      r = exp(gg_y_re / 2.0);
      gg_y_re = r * (r * cos(gg_y_im));
      gg_y_im = r * (r * sin(gg_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    hg_y_re = rt_powd_snf(t_s.re, 5.0);
    hg_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    hg_y_re = 0.0;
    hg_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    hg_y_re = 5.0 * x_re;
    hg_y_im = 5.0 * x_im;
    if (hg_y_im == 0.0) {
      hg_y_re = exp(hg_y_re);
      hg_y_im = 0.0;
    } else if (rtIsInf(hg_y_im) && rtIsInf(hg_y_re) && (hg_y_re < 0.0)) {
      hg_y_re = 0.0;
      hg_y_im = 0.0;
    } else {
      r = exp(hg_y_re / 2.0);
      hg_y_re = r * (r * cos(hg_y_im));
      hg_y_im = r * (r * sin(hg_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    ig_y_re = rt_powd_snf(t_s.re, 4.0);
    ig_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    ig_y_re = rt_powd_snf(t_s.im, 4.0);
    ig_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    ig_y_re = 4.0 * x_re;
    ig_y_im = 4.0 * x_im;
    if (ig_y_im == 0.0) {
      ig_y_re = exp(ig_y_re);
      ig_y_im = 0.0;
    } else if (rtIsInf(ig_y_im) && rtIsInf(ig_y_re) && (ig_y_re < 0.0)) {
      ig_y_re = 0.0;
      ig_y_im = 0.0;
    } else {
      r = exp(ig_y_re / 2.0);
      ig_y_re = r * (r * cos(ig_y_im));
      ig_y_im = r * (r * sin(ig_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    jg_y_re = rt_powd_snf(t.re, 5.0);
    jg_y_im = 0.0;
  } else if (t.re == 0.0) {
    jg_y_re = 0.0;
    jg_y_im = rt_powd_snf(t.im, 5.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    jg_y_re = 5.0 * x_re;
    jg_y_im = 5.0 * x_im;
    if (jg_y_im == 0.0) {
      jg_y_re = exp(jg_y_re);
      jg_y_im = 0.0;
    } else if (rtIsInf(jg_y_im) && rtIsInf(jg_y_re) && (jg_y_re < 0.0)) {
      jg_y_re = 0.0;
      jg_y_im = 0.0;
    } else {
      r = exp(jg_y_re / 2.0);
      jg_y_re = r * (r * cos(jg_y_im));
      jg_y_im = r * (r * sin(jg_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    kg_y_re = rt_powd_snf(t_s.re, 6.0);
    kg_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    kg_y_re = -rt_powd_snf(t_s.im, 6.0);
    kg_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    kg_y_re = 6.0 * x_re;
    kg_y_im = 6.0 * x_im;
    if (kg_y_im == 0.0) {
      kg_y_re = exp(kg_y_re);
      kg_y_im = 0.0;
    } else if (rtIsInf(kg_y_im) && rtIsInf(kg_y_re) && (kg_y_re < 0.0)) {
      kg_y_re = 0.0;
      kg_y_im = 0.0;
    } else {
      r = exp(kg_y_re / 2.0);
      kg_y_re = r * (r * cos(kg_y_im));
      kg_y_im = r * (r * sin(kg_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    lg_y_re = rt_powd_snf(t_s.re, 5.0);
    lg_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    lg_y_re = 0.0;
    lg_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    lg_y_re = 5.0 * x_re;
    lg_y_im = 5.0 * x_im;
    if (lg_y_im == 0.0) {
      lg_y_re = exp(lg_y_re);
      lg_y_im = 0.0;
    } else if (rtIsInf(lg_y_im) && rtIsInf(lg_y_re) && (lg_y_re < 0.0)) {
      lg_y_re = 0.0;
      lg_y_im = 0.0;
    } else {
      r = exp(lg_y_re / 2.0);
      lg_y_re = r * (r * cos(lg_y_im));
      lg_y_im = r * (r * sin(lg_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    mg_y_re = rt_powd_snf(t.re, 5.0);
    mg_y_im = 0.0;
  } else if (t.re == 0.0) {
    mg_y_re = 0.0;
    mg_y_im = rt_powd_snf(t.im, 5.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    mg_y_re = 5.0 * x_re;
    mg_y_im = 5.0 * x_im;
    if (mg_y_im == 0.0) {
      mg_y_re = exp(mg_y_re);
      mg_y_im = 0.0;
    } else if (rtIsInf(mg_y_im) && rtIsInf(mg_y_re) && (mg_y_re < 0.0)) {
      mg_y_re = 0.0;
      mg_y_im = 0.0;
    } else {
      r = exp(mg_y_re / 2.0);
      mg_y_re = r * (r * cos(mg_y_im));
      mg_y_im = r * (r * sin(mg_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    ng_y_re = rt_powd_snf(t_s.re, 3.0);
    ng_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    ng_y_re = 0.0;
    ng_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    ng_y_re = 3.0 * x_re;
    ng_y_im = 3.0 * x_im;
    if (ng_y_im == 0.0) {
      ng_y_re = exp(ng_y_re);
      ng_y_im = 0.0;
    } else if (rtIsInf(ng_y_im) && rtIsInf(ng_y_re) && (ng_y_re < 0.0)) {
      ng_y_re = 0.0;
      ng_y_im = 0.0;
    } else {
      r = exp(ng_y_re / 2.0);
      ng_y_re = r * (r * cos(ng_y_im));
      ng_y_im = r * (r * sin(ng_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    og_y_re = rt_powd_snf(t_s.re, 6.0);
    og_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    og_y_re = -rt_powd_snf(t_s.im, 6.0);
    og_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    og_y_re = 6.0 * x_re;
    og_y_im = 6.0 * x_im;
    if (og_y_im == 0.0) {
      og_y_re = exp(og_y_re);
      og_y_im = 0.0;
    } else if (rtIsInf(og_y_im) && rtIsInf(og_y_re) && (og_y_re < 0.0)) {
      og_y_re = 0.0;
      og_y_im = 0.0;
    } else {
      r = exp(og_y_re / 2.0);
      og_y_re = r * (r * cos(og_y_im));
      og_y_im = r * (r * sin(og_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    pg_y_re = rt_powd_snf(t_s.re, 5.0);
    pg_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    pg_y_re = 0.0;
    pg_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    pg_y_re = 5.0 * x_re;
    pg_y_im = 5.0 * x_im;
    if (pg_y_im == 0.0) {
      pg_y_re = exp(pg_y_re);
      pg_y_im = 0.0;
    } else if (rtIsInf(pg_y_im) && rtIsInf(pg_y_re) && (pg_y_re < 0.0)) {
      pg_y_re = 0.0;
      pg_y_im = 0.0;
    } else {
      r = exp(pg_y_re / 2.0);
      pg_y_re = r * (r * cos(pg_y_im));
      pg_y_im = r * (r * sin(pg_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    qg_y_re = rt_powd_snf(t.re, 3.0);
    qg_y_im = 0.0;
  } else if (t.re == 0.0) {
    qg_y_re = 0.0;
    qg_y_im = -rt_powd_snf(t.im, 3.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    qg_y_re = 3.0 * x_re;
    qg_y_im = 3.0 * x_im;
    if (qg_y_im == 0.0) {
      qg_y_re = exp(qg_y_re);
      qg_y_im = 0.0;
    } else if (rtIsInf(qg_y_im) && rtIsInf(qg_y_re) && (qg_y_re < 0.0)) {
      qg_y_re = 0.0;
      qg_y_im = 0.0;
    } else {
      r = exp(qg_y_re / 2.0);
      qg_y_re = r * (r * cos(qg_y_im));
      qg_y_im = r * (r * sin(qg_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    rg_y_re = rt_powd_snf(t_s.re, 3.0);
    rg_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    rg_y_re = 0.0;
    rg_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    rg_y_re = 3.0 * x_re;
    rg_y_im = 3.0 * x_im;
    if (rg_y_im == 0.0) {
      rg_y_re = exp(rg_y_re);
      rg_y_im = 0.0;
    } else if (rtIsInf(rg_y_im) && rtIsInf(rg_y_re) && (rg_y_re < 0.0)) {
      rg_y_re = 0.0;
      rg_y_im = 0.0;
    } else {
      r = exp(rg_y_re / 2.0);
      rg_y_re = r * (r * cos(rg_y_im));
      rg_y_im = r * (r * sin(rg_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    sg_y_re = rt_powd_snf(t.re, 3.0);
    sg_y_im = 0.0;
  } else if (t.re == 0.0) {
    sg_y_re = 0.0;
    sg_y_im = -rt_powd_snf(t.im, 3.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    sg_y_re = 3.0 * x_re;
    sg_y_im = 3.0 * x_im;
    if (sg_y_im == 0.0) {
      sg_y_re = exp(sg_y_re);
      sg_y_im = 0.0;
    } else if (rtIsInf(sg_y_im) && rtIsInf(sg_y_re) && (sg_y_re < 0.0)) {
      sg_y_re = 0.0;
      sg_y_im = 0.0;
    } else {
      r = exp(sg_y_re / 2.0);
      sg_y_re = r * (r * cos(sg_y_im));
      sg_y_im = r * (r * sin(sg_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    tg_y_re = rt_powd_snf(t_s.re, 4.0);
    tg_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    tg_y_re = rt_powd_snf(t_s.im, 4.0);
    tg_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    tg_y_re = 4.0 * x_re;
    tg_y_im = 4.0 * x_im;
    if (tg_y_im == 0.0) {
      tg_y_re = exp(tg_y_re);
      tg_y_im = 0.0;
    } else if (rtIsInf(tg_y_im) && rtIsInf(tg_y_re) && (tg_y_re < 0.0)) {
      tg_y_re = 0.0;
      tg_y_im = 0.0;
    } else {
      r = exp(tg_y_re / 2.0);
      tg_y_re = r * (r * cos(tg_y_im));
      tg_y_im = r * (r * sin(tg_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    ug_y_re = rt_powd_snf(t_s.re, 3.0);
    ug_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    ug_y_re = 0.0;
    ug_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    ug_y_re = 3.0 * x_re;
    ug_y_im = 3.0 * x_im;
    if (ug_y_im == 0.0) {
      ug_y_re = exp(ug_y_re);
      ug_y_im = 0.0;
    } else if (rtIsInf(ug_y_im) && rtIsInf(ug_y_re) && (ug_y_re < 0.0)) {
      ug_y_re = 0.0;
      ug_y_im = 0.0;
    } else {
      r = exp(ug_y_re / 2.0);
      ug_y_re = r * (r * cos(ug_y_im));
      ug_y_im = r * (r * sin(ug_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    vg_y_re = rt_powd_snf(t.re, 4.0);
    vg_y_im = 0.0;
  } else if (t.re == 0.0) {
    vg_y_re = rt_powd_snf(t.im, 4.0);
    vg_y_im = 0.0;
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    vg_y_re = 4.0 * x_re;
    vg_y_im = 4.0 * x_im;
    if (vg_y_im == 0.0) {
      vg_y_re = exp(vg_y_re);
      vg_y_im = 0.0;
    } else if (rtIsInf(vg_y_im) && rtIsInf(vg_y_re) && (vg_y_re < 0.0)) {
      vg_y_re = 0.0;
      vg_y_im = 0.0;
    } else {
      r = exp(vg_y_re / 2.0);
      vg_y_re = r * (r * cos(vg_y_im));
      vg_y_im = r * (r * sin(vg_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    wg_y_re = rt_powd_snf(t_s.re, 5.0);
    wg_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    wg_y_re = 0.0;
    wg_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    wg_y_re = 5.0 * x_re;
    wg_y_im = 5.0 * x_im;
    if (wg_y_im == 0.0) {
      wg_y_re = exp(wg_y_re);
      wg_y_im = 0.0;
    } else if (rtIsInf(wg_y_im) && rtIsInf(wg_y_re) && (wg_y_re < 0.0)) {
      wg_y_re = 0.0;
      wg_y_im = 0.0;
    } else {
      r = exp(wg_y_re / 2.0);
      wg_y_re = r * (r * cos(wg_y_im));
      wg_y_im = r * (r * sin(wg_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    xg_y_re = rt_powd_snf(t_s.re, 4.0);
    xg_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    xg_y_re = rt_powd_snf(t_s.im, 4.0);
    xg_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    xg_y_re = 4.0 * x_re;
    xg_y_im = 4.0 * x_im;
    if (xg_y_im == 0.0) {
      xg_y_re = exp(xg_y_re);
      xg_y_im = 0.0;
    } else if (rtIsInf(xg_y_im) && rtIsInf(xg_y_re) && (xg_y_re < 0.0)) {
      xg_y_re = 0.0;
      xg_y_im = 0.0;
    } else {
      r = exp(xg_y_re / 2.0);
      xg_y_re = r * (r * cos(xg_y_im));
      xg_y_im = r * (r * sin(xg_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    yg_y_re = rt_powd_snf(t.re, 5.0);
    yg_y_im = 0.0;
  } else if (t.re == 0.0) {
    yg_y_re = 0.0;
    yg_y_im = rt_powd_snf(t.im, 5.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    yg_y_re = 5.0 * x_re;
    yg_y_im = 5.0 * x_im;
    if (yg_y_im == 0.0) {
      yg_y_re = exp(yg_y_re);
      yg_y_im = 0.0;
    } else if (rtIsInf(yg_y_im) && rtIsInf(yg_y_re) && (yg_y_re < 0.0)) {
      yg_y_re = 0.0;
      yg_y_im = 0.0;
    } else {
      r = exp(yg_y_re / 2.0);
      yg_y_re = r * (r * cos(yg_y_im));
      yg_y_im = r * (r * sin(yg_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    ah_y_re = rt_powd_snf(t_s.re, 6.0);
    ah_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    ah_y_re = -rt_powd_snf(t_s.im, 6.0);
    ah_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    ah_y_re = 6.0 * x_re;
    ah_y_im = 6.0 * x_im;
    if (ah_y_im == 0.0) {
      ah_y_re = exp(ah_y_re);
      ah_y_im = 0.0;
    } else if (rtIsInf(ah_y_im) && rtIsInf(ah_y_re) && (ah_y_re < 0.0)) {
      ah_y_re = 0.0;
      ah_y_im = 0.0;
    } else {
      r = exp(ah_y_re / 2.0);
      ah_y_re = r * (r * cos(ah_y_im));
      ah_y_im = r * (r * sin(ah_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    bh_y_re = rt_powd_snf(t_s.re, 5.0);
    bh_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    bh_y_re = 0.0;
    bh_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    bh_y_re = 5.0 * x_re;
    bh_y_im = 5.0 * x_im;
    if (bh_y_im == 0.0) {
      bh_y_re = exp(bh_y_re);
      bh_y_im = 0.0;
    } else if (rtIsInf(bh_y_im) && rtIsInf(bh_y_re) && (bh_y_re < 0.0)) {
      bh_y_re = 0.0;
      bh_y_im = 0.0;
    } else {
      r = exp(bh_y_re / 2.0);
      bh_y_re = r * (r * cos(bh_y_im));
      bh_y_im = r * (r * sin(bh_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    ch_y_re = rt_powd_snf(t.re, 3.0);
    ch_y_im = 0.0;
  } else if (t.re == 0.0) {
    ch_y_re = 0.0;
    ch_y_im = -rt_powd_snf(t.im, 3.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    ch_y_re = 3.0 * x_re;
    ch_y_im = 3.0 * x_im;
    if (ch_y_im == 0.0) {
      ch_y_re = exp(ch_y_re);
      ch_y_im = 0.0;
    } else if (rtIsInf(ch_y_im) && rtIsInf(ch_y_re) && (ch_y_re < 0.0)) {
      ch_y_re = 0.0;
      ch_y_im = 0.0;
    } else {
      r = exp(ch_y_re / 2.0);
      ch_y_re = r * (r * cos(ch_y_im));
      ch_y_im = r * (r * sin(ch_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    dh_y_re = rt_powd_snf(t_s.re, 3.0);
    dh_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    dh_y_re = 0.0;
    dh_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    dh_y_re = 3.0 * x_re;
    dh_y_im = 3.0 * x_im;
    if (dh_y_im == 0.0) {
      dh_y_re = exp(dh_y_re);
      dh_y_im = 0.0;
    } else if (rtIsInf(dh_y_im) && rtIsInf(dh_y_re) && (dh_y_re < 0.0)) {
      dh_y_re = 0.0;
      dh_y_im = 0.0;
    } else {
      r = exp(dh_y_re / 2.0);
      dh_y_re = r * (r * cos(dh_y_im));
      dh_y_im = r * (r * sin(dh_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    eh_y_re = rt_powd_snf(t.re, 3.0);
    eh_y_im = 0.0;
  } else if (t.re == 0.0) {
    eh_y_re = 0.0;
    eh_y_im = -rt_powd_snf(t.im, 3.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    eh_y_re = 3.0 * x_re;
    eh_y_im = 3.0 * x_im;
    if (eh_y_im == 0.0) {
      eh_y_re = exp(eh_y_re);
      eh_y_im = 0.0;
    } else if (rtIsInf(eh_y_im) && rtIsInf(eh_y_re) && (eh_y_re < 0.0)) {
      eh_y_re = 0.0;
      eh_y_im = 0.0;
    } else {
      r = exp(eh_y_re / 2.0);
      eh_y_re = r * (r * cos(eh_y_im));
      eh_y_im = r * (r * sin(eh_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    fh_y_re = rt_powd_snf(t_s.re, 4.0);
    fh_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    fh_y_re = rt_powd_snf(t_s.im, 4.0);
    fh_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    fh_y_re = 4.0 * x_re;
    fh_y_im = 4.0 * x_im;
    if (fh_y_im == 0.0) {
      fh_y_re = exp(fh_y_re);
      fh_y_im = 0.0;
    } else if (rtIsInf(fh_y_im) && rtIsInf(fh_y_re) && (fh_y_re < 0.0)) {
      fh_y_re = 0.0;
      fh_y_im = 0.0;
    } else {
      r = exp(fh_y_re / 2.0);
      fh_y_re = r * (r * cos(fh_y_im));
      fh_y_im = r * (r * sin(fh_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    gh_y_re = rt_powd_snf(t_s.re, 3.0);
    gh_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    gh_y_re = 0.0;
    gh_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    gh_y_re = 3.0 * x_re;
    gh_y_im = 3.0 * x_im;
    if (gh_y_im == 0.0) {
      gh_y_re = exp(gh_y_re);
      gh_y_im = 0.0;
    } else if (rtIsInf(gh_y_im) && rtIsInf(gh_y_re) && (gh_y_re < 0.0)) {
      gh_y_re = 0.0;
      gh_y_im = 0.0;
    } else {
      r = exp(gh_y_re / 2.0);
      gh_y_re = r * (r * cos(gh_y_im));
      gh_y_im = r * (r * sin(gh_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    hh_y_re = rt_powd_snf(t.re, 4.0);
    hh_y_im = 0.0;
  } else if (t.re == 0.0) {
    hh_y_re = rt_powd_snf(t.im, 4.0);
    hh_y_im = 0.0;
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    hh_y_re = 4.0 * x_re;
    hh_y_im = 4.0 * x_im;
    if (hh_y_im == 0.0) {
      hh_y_re = exp(hh_y_re);
      hh_y_im = 0.0;
    } else if (rtIsInf(hh_y_im) && rtIsInf(hh_y_re) && (hh_y_re < 0.0)) {
      hh_y_re = 0.0;
      hh_y_im = 0.0;
    } else {
      r = exp(hh_y_re / 2.0);
      hh_y_re = r * (r * cos(hh_y_im));
      hh_y_im = r * (r * sin(hh_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    ih_y_re = rt_powd_snf(t_s.re, 5.0);
    ih_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    ih_y_re = 0.0;
    ih_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    ih_y_re = 5.0 * x_re;
    ih_y_im = 5.0 * x_im;
    if (ih_y_im == 0.0) {
      ih_y_re = exp(ih_y_re);
      ih_y_im = 0.0;
    } else if (rtIsInf(ih_y_im) && rtIsInf(ih_y_re) && (ih_y_re < 0.0)) {
      ih_y_re = 0.0;
      ih_y_im = 0.0;
    } else {
      r = exp(ih_y_re / 2.0);
      ih_y_re = r * (r * cos(ih_y_im));
      ih_y_im = r * (r * sin(ih_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    jh_y_re = rt_powd_snf(t_s.re, 4.0);
    jh_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    jh_y_re = rt_powd_snf(t_s.im, 4.0);
    jh_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    jh_y_re = 4.0 * x_re;
    jh_y_im = 4.0 * x_im;
    if (jh_y_im == 0.0) {
      jh_y_re = exp(jh_y_re);
      jh_y_im = 0.0;
    } else if (rtIsInf(jh_y_im) && rtIsInf(jh_y_re) && (jh_y_re < 0.0)) {
      jh_y_re = 0.0;
      jh_y_im = 0.0;
    } else {
      r = exp(jh_y_re / 2.0);
      jh_y_re = r * (r * cos(jh_y_im));
      jh_y_im = r * (r * sin(jh_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    kh_y_re = rt_powd_snf(t.re, 5.0);
    kh_y_im = 0.0;
  } else if (t.re == 0.0) {
    kh_y_re = 0.0;
    kh_y_im = rt_powd_snf(t.im, 5.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    kh_y_re = 5.0 * x_re;
    kh_y_im = 5.0 * x_im;
    if (kh_y_im == 0.0) {
      kh_y_re = exp(kh_y_re);
      kh_y_im = 0.0;
    } else if (rtIsInf(kh_y_im) && rtIsInf(kh_y_re) && (kh_y_re < 0.0)) {
      kh_y_re = 0.0;
      kh_y_im = 0.0;
    } else {
      r = exp(kh_y_re / 2.0);
      kh_y_re = r * (r * cos(kh_y_im));
      kh_y_im = r * (r * sin(kh_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    lh_y_re = rt_powd_snf(t_s.re, 6.0);
    lh_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    lh_y_re = -rt_powd_snf(t_s.im, 6.0);
    lh_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    lh_y_re = 6.0 * x_re;
    lh_y_im = 6.0 * x_im;
    if (lh_y_im == 0.0) {
      lh_y_re = exp(lh_y_re);
      lh_y_im = 0.0;
    } else if (rtIsInf(lh_y_im) && rtIsInf(lh_y_re) && (lh_y_re < 0.0)) {
      lh_y_re = 0.0;
      lh_y_im = 0.0;
    } else {
      r = exp(lh_y_re / 2.0);
      lh_y_re = r * (r * cos(lh_y_im));
      lh_y_im = r * (r * sin(lh_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    mh_y_re = rt_powd_snf(t_s.re, 5.0);
    mh_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    mh_y_re = 0.0;
    mh_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    mh_y_re = 5.0 * x_re;
    mh_y_im = 5.0 * x_im;
    if (mh_y_im == 0.0) {
      mh_y_re = exp(mh_y_re);
      mh_y_im = 0.0;
    } else if (rtIsInf(mh_y_im) && rtIsInf(mh_y_re) && (mh_y_re < 0.0)) {
      mh_y_re = 0.0;
      mh_y_im = 0.0;
    } else {
      r = exp(mh_y_re / 2.0);
      mh_y_re = r * (r * cos(mh_y_im));
      mh_y_im = r * (r * sin(mh_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    nh_y_re = rt_powd_snf(t.re, 3.0);
    nh_y_im = 0.0;
  } else if (t.re == 0.0) {
    nh_y_re = 0.0;
    nh_y_im = -rt_powd_snf(t.im, 3.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    nh_y_re = 3.0 * x_re;
    nh_y_im = 3.0 * x_im;
    if (nh_y_im == 0.0) {
      nh_y_re = exp(nh_y_re);
      nh_y_im = 0.0;
    } else if (rtIsInf(nh_y_im) && rtIsInf(nh_y_re) && (nh_y_re < 0.0)) {
      nh_y_re = 0.0;
      nh_y_im = 0.0;
    } else {
      r = exp(nh_y_re / 2.0);
      nh_y_re = r * (r * cos(nh_y_im));
      nh_y_im = r * (r * sin(nh_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    oh_y_re = rt_powd_snf(t_s.re, 4.0);
    oh_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    oh_y_re = rt_powd_snf(t_s.im, 4.0);
    oh_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    oh_y_re = 4.0 * x_re;
    oh_y_im = 4.0 * x_im;
    if (oh_y_im == 0.0) {
      oh_y_re = exp(oh_y_re);
      oh_y_im = 0.0;
    } else if (rtIsInf(oh_y_im) && rtIsInf(oh_y_re) && (oh_y_re < 0.0)) {
      oh_y_re = 0.0;
      oh_y_im = 0.0;
    } else {
      r = exp(oh_y_re / 2.0);
      oh_y_re = r * (r * cos(oh_y_im));
      oh_y_im = r * (r * sin(oh_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    ph_y_re = rt_powd_snf(t_s.re, 3.0);
    ph_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    ph_y_re = 0.0;
    ph_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    ph_y_re = 3.0 * x_re;
    ph_y_im = 3.0 * x_im;
    if (ph_y_im == 0.0) {
      ph_y_re = exp(ph_y_re);
      ph_y_im = 0.0;
    } else if (rtIsInf(ph_y_im) && rtIsInf(ph_y_re) && (ph_y_re < 0.0)) {
      ph_y_re = 0.0;
      ph_y_im = 0.0;
    } else {
      r = exp(ph_y_re / 2.0);
      ph_y_re = r * (r * cos(ph_y_im));
      ph_y_im = r * (r * sin(ph_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    qh_y_re = rt_powd_snf(t.re, 3.0);
    qh_y_im = 0.0;
  } else if (t.re == 0.0) {
    qh_y_re = 0.0;
    qh_y_im = -rt_powd_snf(t.im, 3.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    qh_y_re = 3.0 * x_re;
    qh_y_im = 3.0 * x_im;
    if (qh_y_im == 0.0) {
      qh_y_re = exp(qh_y_re);
      qh_y_im = 0.0;
    } else if (rtIsInf(qh_y_im) && rtIsInf(qh_y_re) && (qh_y_re < 0.0)) {
      qh_y_re = 0.0;
      qh_y_im = 0.0;
    } else {
      r = exp(qh_y_re / 2.0);
      qh_y_re = r * (r * cos(qh_y_im));
      qh_y_im = r * (r * sin(qh_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    rh_y_re = rt_powd_snf(t_s.re, 3.0);
    rh_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    rh_y_re = 0.0;
    rh_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    rh_y_re = 3.0 * x_re;
    rh_y_im = 3.0 * x_im;
    if (rh_y_im == 0.0) {
      rh_y_re = exp(rh_y_re);
      rh_y_im = 0.0;
    } else if (rtIsInf(rh_y_im) && rtIsInf(rh_y_re) && (rh_y_re < 0.0)) {
      rh_y_re = 0.0;
      rh_y_im = 0.0;
    } else {
      r = exp(rh_y_re / 2.0);
      rh_y_re = r * (r * cos(rh_y_im));
      rh_y_im = r * (r * sin(rh_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    sh_y_re = rt_powd_snf(t.re, 3.0);
    sh_y_im = 0.0;
  } else if (t.re == 0.0) {
    sh_y_re = 0.0;
    sh_y_im = -rt_powd_snf(t.im, 3.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    sh_y_re = 3.0 * x_re;
    sh_y_im = 3.0 * x_im;
    if (sh_y_im == 0.0) {
      sh_y_re = exp(sh_y_re);
      sh_y_im = 0.0;
    } else if (rtIsInf(sh_y_im) && rtIsInf(sh_y_re) && (sh_y_re < 0.0)) {
      sh_y_re = 0.0;
      sh_y_im = 0.0;
    } else {
      r = exp(sh_y_re / 2.0);
      sh_y_re = r * (r * cos(sh_y_im));
      sh_y_im = r * (r * sin(sh_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    th_y_re = rt_powd_snf(t_s.re, 4.0);
    th_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    th_y_re = rt_powd_snf(t_s.im, 4.0);
    th_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    th_y_re = 4.0 * x_re;
    th_y_im = 4.0 * x_im;
    if (th_y_im == 0.0) {
      th_y_re = exp(th_y_re);
      th_y_im = 0.0;
    } else if (rtIsInf(th_y_im) && rtIsInf(th_y_re) && (th_y_re < 0.0)) {
      th_y_re = 0.0;
      th_y_im = 0.0;
    } else {
      r = exp(th_y_re / 2.0);
      th_y_re = r * (r * cos(th_y_im));
      th_y_im = r * (r * sin(th_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    uh_y_re = rt_powd_snf(t_s.re, 3.0);
    uh_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    uh_y_re = 0.0;
    uh_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    uh_y_re = 3.0 * x_re;
    uh_y_im = 3.0 * x_im;
    if (uh_y_im == 0.0) {
      uh_y_re = exp(uh_y_re);
      uh_y_im = 0.0;
    } else if (rtIsInf(uh_y_im) && rtIsInf(uh_y_re) && (uh_y_re < 0.0)) {
      uh_y_re = 0.0;
      uh_y_im = 0.0;
    } else {
      r = exp(uh_y_re / 2.0);
      uh_y_re = r * (r * cos(uh_y_im));
      uh_y_im = r * (r * sin(uh_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    vh_y_re = rt_powd_snf(t.re, 4.0);
    vh_y_im = 0.0;
  } else if (t.re == 0.0) {
    vh_y_re = rt_powd_snf(t.im, 4.0);
    vh_y_im = 0.0;
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    vh_y_re = 4.0 * x_re;
    vh_y_im = 4.0 * x_im;
    if (vh_y_im == 0.0) {
      vh_y_re = exp(vh_y_re);
      vh_y_im = 0.0;
    } else if (rtIsInf(vh_y_im) && rtIsInf(vh_y_re) && (vh_y_re < 0.0)) {
      vh_y_re = 0.0;
      vh_y_im = 0.0;
    } else {
      r = exp(vh_y_re / 2.0);
      vh_y_re = r * (r * cos(vh_y_im));
      vh_y_im = r * (r * sin(vh_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    wh_y_re = rt_powd_snf(t_s.re, 5.0);
    wh_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    wh_y_re = 0.0;
    wh_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    wh_y_re = 5.0 * x_re;
    wh_y_im = 5.0 * x_im;
    if (wh_y_im == 0.0) {
      wh_y_re = exp(wh_y_re);
      wh_y_im = 0.0;
    } else if (rtIsInf(wh_y_im) && rtIsInf(wh_y_re) && (wh_y_re < 0.0)) {
      wh_y_re = 0.0;
      wh_y_im = 0.0;
    } else {
      r = exp(wh_y_re / 2.0);
      wh_y_re = r * (r * cos(wh_y_im));
      wh_y_im = r * (r * sin(wh_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    xh_y_re = rt_powd_snf(t_s.re, 4.0);
    xh_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    xh_y_re = rt_powd_snf(t_s.im, 4.0);
    xh_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    xh_y_re = 4.0 * x_re;
    xh_y_im = 4.0 * x_im;
    if (xh_y_im == 0.0) {
      xh_y_re = exp(xh_y_re);
      xh_y_im = 0.0;
    } else if (rtIsInf(xh_y_im) && rtIsInf(xh_y_re) && (xh_y_re < 0.0)) {
      xh_y_re = 0.0;
      xh_y_im = 0.0;
    } else {
      r = exp(xh_y_re / 2.0);
      xh_y_re = r * (r * cos(xh_y_im));
      xh_y_im = r * (r * sin(xh_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    yh_y_re = rt_powd_snf(t.re, 4.0);
    yh_y_im = 0.0;
  } else if (t.re == 0.0) {
    yh_y_re = rt_powd_snf(t.im, 4.0);
    yh_y_im = 0.0;
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    yh_y_re = 4.0 * x_re;
    yh_y_im = 4.0 * x_im;
    if (yh_y_im == 0.0) {
      yh_y_re = exp(yh_y_re);
      yh_y_im = 0.0;
    } else if (rtIsInf(yh_y_im) && rtIsInf(yh_y_re) && (yh_y_re < 0.0)) {
      yh_y_re = 0.0;
      yh_y_im = 0.0;
    } else {
      r = exp(yh_y_re / 2.0);
      yh_y_re = r * (r * cos(yh_y_im));
      yh_y_im = r * (r * sin(yh_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    ai_y_re = rt_powd_snf(t_s.re, 5.0);
    ai_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    ai_y_re = 0.0;
    ai_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    ai_y_re = 5.0 * x_re;
    ai_y_im = 5.0 * x_im;
    if (ai_y_im == 0.0) {
      ai_y_re = exp(ai_y_re);
      ai_y_im = 0.0;
    } else if (rtIsInf(ai_y_im) && rtIsInf(ai_y_re) && (ai_y_re < 0.0)) {
      ai_y_re = 0.0;
      ai_y_im = 0.0;
    } else {
      r = exp(ai_y_re / 2.0);
      ai_y_re = r * (r * cos(ai_y_im));
      ai_y_im = r * (r * sin(ai_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    bi_y_re = rt_powd_snf(t_s.re, 4.0);
    bi_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    bi_y_re = rt_powd_snf(t_s.im, 4.0);
    bi_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    bi_y_re = 4.0 * x_re;
    bi_y_im = 4.0 * x_im;
    if (bi_y_im == 0.0) {
      bi_y_re = exp(bi_y_re);
      bi_y_im = 0.0;
    } else if (rtIsInf(bi_y_im) && rtIsInf(bi_y_re) && (bi_y_re < 0.0)) {
      bi_y_re = 0.0;
      bi_y_im = 0.0;
    } else {
      r = exp(bi_y_re / 2.0);
      bi_y_re = r * (r * cos(bi_y_im));
      bi_y_im = r * (r * sin(bi_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    ci_y_re = rt_powd_snf(t.re, 5.0);
    ci_y_im = 0.0;
  } else if (t.re == 0.0) {
    ci_y_re = 0.0;
    ci_y_im = rt_powd_snf(t.im, 5.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    ci_y_re = 5.0 * x_re;
    ci_y_im = 5.0 * x_im;
    if (ci_y_im == 0.0) {
      ci_y_re = exp(ci_y_re);
      ci_y_im = 0.0;
    } else if (rtIsInf(ci_y_im) && rtIsInf(ci_y_re) && (ci_y_re < 0.0)) {
      ci_y_re = 0.0;
      ci_y_im = 0.0;
    } else {
      r = exp(ci_y_re / 2.0);
      ci_y_re = r * (r * cos(ci_y_im));
      ci_y_im = r * (r * sin(ci_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    di_y_re = rt_powd_snf(t_s.re, 6.0);
    di_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    di_y_re = -rt_powd_snf(t_s.im, 6.0);
    di_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    di_y_re = 6.0 * x_re;
    di_y_im = 6.0 * x_im;
    if (di_y_im == 0.0) {
      di_y_re = exp(di_y_re);
      di_y_im = 0.0;
    } else if (rtIsInf(di_y_im) && rtIsInf(di_y_re) && (di_y_re < 0.0)) {
      di_y_re = 0.0;
      di_y_im = 0.0;
    } else {
      r = exp(di_y_re / 2.0);
      di_y_re = r * (r * cos(di_y_im));
      di_y_im = r * (r * sin(di_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    ei_y_re = rt_powd_snf(t_s.re, 5.0);
    ei_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    ei_y_re = 0.0;
    ei_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    ei_y_re = 5.0 * x_re;
    ei_y_im = 5.0 * x_im;
    if (ei_y_im == 0.0) {
      ei_y_re = exp(ei_y_re);
      ei_y_im = 0.0;
    } else if (rtIsInf(ei_y_im) && rtIsInf(ei_y_re) && (ei_y_re < 0.0)) {
      ei_y_re = 0.0;
      ei_y_im = 0.0;
    } else {
      r = exp(ei_y_re / 2.0);
      ei_y_re = r * (r * cos(ei_y_im));
      ei_y_im = r * (r * sin(ei_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    fi_y_re = rt_powd_snf(t.re, 5.0);
    fi_y_im = 0.0;
  } else if (t.re == 0.0) {
    fi_y_re = 0.0;
    fi_y_im = rt_powd_snf(t.im, 5.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    fi_y_re = 5.0 * x_re;
    fi_y_im = 5.0 * x_im;
    if (fi_y_im == 0.0) {
      fi_y_re = exp(fi_y_re);
      fi_y_im = 0.0;
    } else if (rtIsInf(fi_y_im) && rtIsInf(fi_y_re) && (fi_y_re < 0.0)) {
      fi_y_re = 0.0;
      fi_y_im = 0.0;
    } else {
      r = exp(fi_y_re / 2.0);
      fi_y_re = r * (r * cos(fi_y_im));
      fi_y_im = r * (r * sin(fi_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    gi_y_re = rt_powd_snf(t_s.re, 6.0);
    gi_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    gi_y_re = -rt_powd_snf(t_s.im, 6.0);
    gi_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    gi_y_re = 6.0 * x_re;
    gi_y_im = 6.0 * x_im;
    if (gi_y_im == 0.0) {
      gi_y_re = exp(gi_y_re);
      gi_y_im = 0.0;
    } else if (rtIsInf(gi_y_im) && rtIsInf(gi_y_re) && (gi_y_re < 0.0)) {
      gi_y_re = 0.0;
      gi_y_im = 0.0;
    } else {
      r = exp(gi_y_re / 2.0);
      gi_y_re = r * (r * cos(gi_y_im));
      gi_y_im = r * (r * sin(gi_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    hi_y_re = rt_powd_snf(t_s.re, 5.0);
    hi_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    hi_y_re = 0.0;
    hi_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    hi_y_re = 5.0 * x_re;
    hi_y_im = 5.0 * x_im;
    if (hi_y_im == 0.0) {
      hi_y_re = exp(hi_y_re);
      hi_y_im = 0.0;
    } else if (rtIsInf(hi_y_im) && rtIsInf(hi_y_re) && (hi_y_re < 0.0)) {
      hi_y_re = 0.0;
      hi_y_im = 0.0;
    } else {
      r = exp(hi_y_re / 2.0);
      hi_y_re = r * (r * cos(hi_y_im));
      hi_y_im = r * (r * sin(hi_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    ii_y_re = rt_powd_snf(t.re, 3.0);
    ii_y_im = 0.0;
  } else if (t.re == 0.0) {
    ii_y_re = 0.0;
    ii_y_im = -rt_powd_snf(t.im, 3.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    ii_y_re = 3.0 * x_re;
    ii_y_im = 3.0 * x_im;
    if (ii_y_im == 0.0) {
      ii_y_re = exp(ii_y_re);
      ii_y_im = 0.0;
    } else if (rtIsInf(ii_y_im) && rtIsInf(ii_y_re) && (ii_y_re < 0.0)) {
      ii_y_re = 0.0;
      ii_y_im = 0.0;
    } else {
      r = exp(ii_y_re / 2.0);
      ii_y_re = r * (r * cos(ii_y_im));
      ii_y_im = r * (r * sin(ii_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    ji_y_re = rt_powd_snf(t_s.re, 4.0);
    ji_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    ji_y_re = rt_powd_snf(t_s.im, 4.0);
    ji_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    ji_y_re = 4.0 * x_re;
    ji_y_im = 4.0 * x_im;
    if (ji_y_im == 0.0) {
      ji_y_re = exp(ji_y_re);
      ji_y_im = 0.0;
    } else if (rtIsInf(ji_y_im) && rtIsInf(ji_y_re) && (ji_y_re < 0.0)) {
      ji_y_re = 0.0;
      ji_y_im = 0.0;
    } else {
      r = exp(ji_y_re / 2.0);
      ji_y_re = r * (r * cos(ji_y_im));
      ji_y_im = r * (r * sin(ji_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    ki_y_re = rt_powd_snf(t_s.re, 3.0);
    ki_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    ki_y_re = 0.0;
    ki_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    ki_y_re = 3.0 * x_re;
    ki_y_im = 3.0 * x_im;
    if (ki_y_im == 0.0) {
      ki_y_re = exp(ki_y_re);
      ki_y_im = 0.0;
    } else if (rtIsInf(ki_y_im) && rtIsInf(ki_y_re) && (ki_y_re < 0.0)) {
      ki_y_re = 0.0;
      ki_y_im = 0.0;
    } else {
      r = exp(ki_y_re / 2.0);
      ki_y_re = r * (r * cos(ki_y_im));
      ki_y_im = r * (r * sin(ki_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    li_y_re = rt_powd_snf(t.re, 4.0);
    li_y_im = 0.0;
  } else if (t.re == 0.0) {
    li_y_re = rt_powd_snf(t.im, 4.0);
    li_y_im = 0.0;
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    li_y_re = 4.0 * x_re;
    li_y_im = 4.0 * x_im;
    if (li_y_im == 0.0) {
      li_y_re = exp(li_y_re);
      li_y_im = 0.0;
    } else if (rtIsInf(li_y_im) && rtIsInf(li_y_re) && (li_y_re < 0.0)) {
      li_y_re = 0.0;
      li_y_im = 0.0;
    } else {
      r = exp(li_y_re / 2.0);
      li_y_re = r * (r * cos(li_y_im));
      li_y_im = r * (r * sin(li_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    mi_y_re = rt_powd_snf(t_s.re, 5.0);
    mi_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    mi_y_re = 0.0;
    mi_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    mi_y_re = 5.0 * x_re;
    mi_y_im = 5.0 * x_im;
    if (mi_y_im == 0.0) {
      mi_y_re = exp(mi_y_re);
      mi_y_im = 0.0;
    } else if (rtIsInf(mi_y_im) && rtIsInf(mi_y_re) && (mi_y_re < 0.0)) {
      mi_y_re = 0.0;
      mi_y_im = 0.0;
    } else {
      r = exp(mi_y_re / 2.0);
      mi_y_re = r * (r * cos(mi_y_im));
      mi_y_im = r * (r * sin(mi_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    ni_y_re = rt_powd_snf(t_s.re, 4.0);
    ni_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    ni_y_re = rt_powd_snf(t_s.im, 4.0);
    ni_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    ni_y_re = 4.0 * x_re;
    ni_y_im = 4.0 * x_im;
    if (ni_y_im == 0.0) {
      ni_y_re = exp(ni_y_re);
      ni_y_im = 0.0;
    } else if (rtIsInf(ni_y_im) && rtIsInf(ni_y_re) && (ni_y_re < 0.0)) {
      ni_y_re = 0.0;
      ni_y_im = 0.0;
    } else {
      r = exp(ni_y_re / 2.0);
      ni_y_re = r * (r * cos(ni_y_im));
      ni_y_im = r * (r * sin(ni_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    oi_y_re = rt_powd_snf(t.re, 5.0);
    oi_y_im = 0.0;
  } else if (t.re == 0.0) {
    oi_y_re = 0.0;
    oi_y_im = rt_powd_snf(t.im, 5.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    oi_y_re = 5.0 * x_re;
    oi_y_im = 5.0 * x_im;
    if (oi_y_im == 0.0) {
      oi_y_re = exp(oi_y_re);
      oi_y_im = 0.0;
    } else if (rtIsInf(oi_y_im) && rtIsInf(oi_y_re) && (oi_y_re < 0.0)) {
      oi_y_re = 0.0;
      oi_y_im = 0.0;
    } else {
      r = exp(oi_y_re / 2.0);
      oi_y_re = r * (r * cos(oi_y_im));
      oi_y_im = r * (r * sin(oi_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    pi_y_re = rt_powd_snf(t_s.re, 6.0);
    pi_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    pi_y_re = -rt_powd_snf(t_s.im, 6.0);
    pi_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    pi_y_re = 6.0 * x_re;
    pi_y_im = 6.0 * x_im;
    if (pi_y_im == 0.0) {
      pi_y_re = exp(pi_y_re);
      pi_y_im = 0.0;
    } else if (rtIsInf(pi_y_im) && rtIsInf(pi_y_re) && (pi_y_re < 0.0)) {
      pi_y_re = 0.0;
      pi_y_im = 0.0;
    } else {
      r = exp(pi_y_re / 2.0);
      pi_y_re = r * (r * cos(pi_y_im));
      pi_y_im = r * (r * sin(pi_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    qi_y_re = rt_powd_snf(t_s.re, 5.0);
    qi_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    qi_y_re = 0.0;
    qi_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    qi_y_re = 5.0 * x_re;
    qi_y_im = 5.0 * x_im;
    if (qi_y_im == 0.0) {
      qi_y_re = exp(qi_y_re);
      qi_y_im = 0.0;
    } else if (rtIsInf(qi_y_im) && rtIsInf(qi_y_re) && (qi_y_re < 0.0)) {
      qi_y_re = 0.0;
      qi_y_im = 0.0;
    } else {
      r = exp(qi_y_re / 2.0);
      qi_y_re = r * (r * cos(qi_y_im));
      qi_y_im = r * (r * sin(qi_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    ri_y_re = rt_powd_snf(t.re, 3.0);
    ri_y_im = 0.0;
  } else if (t.re == 0.0) {
    ri_y_re = 0.0;
    ri_y_im = -rt_powd_snf(t.im, 3.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    ri_y_re = 3.0 * x_re;
    ri_y_im = 3.0 * x_im;
    if (ri_y_im == 0.0) {
      ri_y_re = exp(ri_y_re);
      ri_y_im = 0.0;
    } else if (rtIsInf(ri_y_im) && rtIsInf(ri_y_re) && (ri_y_re < 0.0)) {
      ri_y_re = 0.0;
      ri_y_im = 0.0;
    } else {
      r = exp(ri_y_re / 2.0);
      ri_y_re = r * (r * cos(ri_y_im));
      ri_y_im = r * (r * sin(ri_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    si_y_re = rt_powd_snf(t_s.re, 4.0);
    si_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    si_y_re = rt_powd_snf(t_s.im, 4.0);
    si_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    si_y_re = 4.0 * x_re;
    si_y_im = 4.0 * x_im;
    if (si_y_im == 0.0) {
      si_y_re = exp(si_y_re);
      si_y_im = 0.0;
    } else if (rtIsInf(si_y_im) && rtIsInf(si_y_re) && (si_y_re < 0.0)) {
      si_y_re = 0.0;
      si_y_im = 0.0;
    } else {
      r = exp(si_y_re / 2.0);
      si_y_re = r * (r * cos(si_y_im));
      si_y_im = r * (r * sin(si_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    ti_y_re = rt_powd_snf(t_s.re, 3.0);
    ti_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    ti_y_re = 0.0;
    ti_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    ti_y_re = 3.0 * x_re;
    ti_y_im = 3.0 * x_im;
    if (ti_y_im == 0.0) {
      ti_y_re = exp(ti_y_re);
      ti_y_im = 0.0;
    } else if (rtIsInf(ti_y_im) && rtIsInf(ti_y_re) && (ti_y_re < 0.0)) {
      ti_y_re = 0.0;
      ti_y_im = 0.0;
    } else {
      r = exp(ti_y_re / 2.0);
      ti_y_re = r * (r * cos(ti_y_im));
      ti_y_im = r * (r * sin(ti_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    ui_y_re = rt_powd_snf(t.re, 3.0);
    ui_y_im = 0.0;
  } else if (t.re == 0.0) {
    ui_y_re = 0.0;
    ui_y_im = -rt_powd_snf(t.im, 3.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    ui_y_re = 3.0 * x_re;
    ui_y_im = 3.0 * x_im;
    if (ui_y_im == 0.0) {
      ui_y_re = exp(ui_y_re);
      ui_y_im = 0.0;
    } else if (rtIsInf(ui_y_im) && rtIsInf(ui_y_re) && (ui_y_re < 0.0)) {
      ui_y_re = 0.0;
      ui_y_im = 0.0;
    } else {
      r = exp(ui_y_re / 2.0);
      ui_y_re = r * (r * cos(ui_y_im));
      ui_y_im = r * (r * sin(ui_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    vi_y_re = rt_powd_snf(t_s.re, 3.0);
    vi_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    vi_y_re = 0.0;
    vi_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    vi_y_re = 3.0 * x_re;
    vi_y_im = 3.0 * x_im;
    if (vi_y_im == 0.0) {
      vi_y_re = exp(vi_y_re);
      vi_y_im = 0.0;
    } else if (rtIsInf(vi_y_im) && rtIsInf(vi_y_re) && (vi_y_re < 0.0)) {
      vi_y_re = 0.0;
      vi_y_im = 0.0;
    } else {
      r = exp(vi_y_re / 2.0);
      vi_y_re = r * (r * cos(vi_y_im));
      vi_y_im = r * (r * sin(vi_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    wi_y_re = rt_powd_snf(t.re, 3.0);
    wi_y_im = 0.0;
  } else if (t.re == 0.0) {
    wi_y_re = 0.0;
    wi_y_im = -rt_powd_snf(t.im, 3.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    wi_y_re = 3.0 * x_re;
    wi_y_im = 3.0 * x_im;
    if (wi_y_im == 0.0) {
      wi_y_re = exp(wi_y_re);
      wi_y_im = 0.0;
    } else if (rtIsInf(wi_y_im) && rtIsInf(wi_y_re) && (wi_y_re < 0.0)) {
      wi_y_re = 0.0;
      wi_y_im = 0.0;
    } else {
      r = exp(wi_y_re / 2.0);
      wi_y_re = r * (r * cos(wi_y_im));
      wi_y_im = r * (r * sin(wi_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    xi_y_re = rt_powd_snf(t_s.re, 4.0);
    xi_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    xi_y_re = rt_powd_snf(t_s.im, 4.0);
    xi_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    xi_y_re = 4.0 * x_re;
    xi_y_im = 4.0 * x_im;
    if (xi_y_im == 0.0) {
      xi_y_re = exp(xi_y_re);
      xi_y_im = 0.0;
    } else if (rtIsInf(xi_y_im) && rtIsInf(xi_y_re) && (xi_y_re < 0.0)) {
      xi_y_re = 0.0;
      xi_y_im = 0.0;
    } else {
      r = exp(xi_y_re / 2.0);
      xi_y_re = r * (r * cos(xi_y_im));
      xi_y_im = r * (r * sin(xi_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    yi_y_re = rt_powd_snf(t_s.re, 3.0);
    yi_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    yi_y_re = 0.0;
    yi_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    yi_y_re = 3.0 * x_re;
    yi_y_im = 3.0 * x_im;
    if (yi_y_im == 0.0) {
      yi_y_re = exp(yi_y_re);
      yi_y_im = 0.0;
    } else if (rtIsInf(yi_y_im) && rtIsInf(yi_y_re) && (yi_y_re < 0.0)) {
      yi_y_re = 0.0;
      yi_y_im = 0.0;
    } else {
      r = exp(yi_y_re / 2.0);
      yi_y_re = r * (r * cos(yi_y_im));
      yi_y_im = r * (r * sin(yi_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    aj_y_re = rt_powd_snf(t.re, 4.0);
    aj_y_im = 0.0;
  } else if (t.re == 0.0) {
    aj_y_re = rt_powd_snf(t.im, 4.0);
    aj_y_im = 0.0;
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    aj_y_re = 4.0 * x_re;
    aj_y_im = 4.0 * x_im;
    if (aj_y_im == 0.0) {
      aj_y_re = exp(aj_y_re);
      aj_y_im = 0.0;
    } else if (rtIsInf(aj_y_im) && rtIsInf(aj_y_re) && (aj_y_re < 0.0)) {
      aj_y_re = 0.0;
      aj_y_im = 0.0;
    } else {
      r = exp(aj_y_re / 2.0);
      aj_y_re = r * (r * cos(aj_y_im));
      aj_y_im = r * (r * sin(aj_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    bj_y_re = rt_powd_snf(t_s.re, 5.0);
    bj_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    bj_y_re = 0.0;
    bj_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    bj_y_re = 5.0 * x_re;
    bj_y_im = 5.0 * x_im;
    if (bj_y_im == 0.0) {
      bj_y_re = exp(bj_y_re);
      bj_y_im = 0.0;
    } else if (rtIsInf(bj_y_im) && rtIsInf(bj_y_re) && (bj_y_re < 0.0)) {
      bj_y_re = 0.0;
      bj_y_im = 0.0;
    } else {
      r = exp(bj_y_re / 2.0);
      bj_y_re = r * (r * cos(bj_y_im));
      bj_y_im = r * (r * sin(bj_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    cj_y_re = rt_powd_snf(t_s.re, 4.0);
    cj_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    cj_y_re = rt_powd_snf(t_s.im, 4.0);
    cj_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    cj_y_re = 4.0 * x_re;
    cj_y_im = 4.0 * x_im;
    if (cj_y_im == 0.0) {
      cj_y_re = exp(cj_y_re);
      cj_y_im = 0.0;
    } else if (rtIsInf(cj_y_im) && rtIsInf(cj_y_re) && (cj_y_re < 0.0)) {
      cj_y_re = 0.0;
      cj_y_im = 0.0;
    } else {
      r = exp(cj_y_re / 2.0);
      cj_y_re = r * (r * cos(cj_y_im));
      cj_y_im = r * (r * sin(cj_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    dj_y_re = rt_powd_snf(t.re, 3.0);
    dj_y_im = 0.0;
  } else if (t.re == 0.0) {
    dj_y_re = 0.0;
    dj_y_im = -rt_powd_snf(t.im, 3.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    dj_y_re = 3.0 * x_re;
    dj_y_im = 3.0 * x_im;
    if (dj_y_im == 0.0) {
      dj_y_re = exp(dj_y_re);
      dj_y_im = 0.0;
    } else if (rtIsInf(dj_y_im) && rtIsInf(dj_y_re) && (dj_y_re < 0.0)) {
      dj_y_re = 0.0;
      dj_y_im = 0.0;
    } else {
      r = exp(dj_y_re / 2.0);
      dj_y_re = r * (r * cos(dj_y_im));
      dj_y_im = r * (r * sin(dj_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    ej_y_re = rt_powd_snf(t_s.re, 3.0);
    ej_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    ej_y_re = 0.0;
    ej_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    ej_y_re = 3.0 * x_re;
    ej_y_im = 3.0 * x_im;
    if (ej_y_im == 0.0) {
      ej_y_re = exp(ej_y_re);
      ej_y_im = 0.0;
    } else if (rtIsInf(ej_y_im) && rtIsInf(ej_y_re) && (ej_y_re < 0.0)) {
      ej_y_re = 0.0;
      ej_y_im = 0.0;
    } else {
      r = exp(ej_y_re / 2.0);
      ej_y_re = r * (r * cos(ej_y_im));
      ej_y_im = r * (r * sin(ej_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    fj_y_re = rt_powd_snf(t_s.re, 4.0);
    fj_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    fj_y_re = rt_powd_snf(t_s.im, 4.0);
    fj_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    fj_y_re = 4.0 * x_re;
    fj_y_im = 4.0 * x_im;
    if (fj_y_im == 0.0) {
      fj_y_re = exp(fj_y_re);
      fj_y_im = 0.0;
    } else if (rtIsInf(fj_y_im) && rtIsInf(fj_y_re) && (fj_y_re < 0.0)) {
      fj_y_re = 0.0;
      fj_y_im = 0.0;
    } else {
      r = exp(fj_y_re / 2.0);
      fj_y_re = r * (r * cos(fj_y_im));
      fj_y_im = r * (r * sin(fj_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    gj_y_re = rt_powd_snf(t_s.re, 3.0);
    gj_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    gj_y_re = 0.0;
    gj_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    gj_y_re = 3.0 * x_re;
    gj_y_im = 3.0 * x_im;
    if (gj_y_im == 0.0) {
      gj_y_re = exp(gj_y_re);
      gj_y_im = 0.0;
    } else if (rtIsInf(gj_y_im) && rtIsInf(gj_y_re) && (gj_y_re < 0.0)) {
      gj_y_re = 0.0;
      gj_y_im = 0.0;
    } else {
      r = exp(gj_y_re / 2.0);
      gj_y_re = r * (r * cos(gj_y_im));
      gj_y_im = r * (r * sin(gj_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    hj_y_re = rt_powd_snf(t.re, 4.0);
    hj_y_im = 0.0;
  } else if (t.re == 0.0) {
    hj_y_re = rt_powd_snf(t.im, 4.0);
    hj_y_im = 0.0;
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    hj_y_re = 4.0 * x_re;
    hj_y_im = 4.0 * x_im;
    if (hj_y_im == 0.0) {
      hj_y_re = exp(hj_y_re);
      hj_y_im = 0.0;
    } else if (rtIsInf(hj_y_im) && rtIsInf(hj_y_re) && (hj_y_re < 0.0)) {
      hj_y_re = 0.0;
      hj_y_im = 0.0;
    } else {
      r = exp(hj_y_re / 2.0);
      hj_y_re = r * (r * cos(hj_y_im));
      hj_y_im = r * (r * sin(hj_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    ij_y_re = rt_powd_snf(t_s.re, 5.0);
    ij_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    ij_y_re = 0.0;
    ij_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    ij_y_re = 5.0 * x_re;
    ij_y_im = 5.0 * x_im;
    if (ij_y_im == 0.0) {
      ij_y_re = exp(ij_y_re);
      ij_y_im = 0.0;
    } else if (rtIsInf(ij_y_im) && rtIsInf(ij_y_re) && (ij_y_re < 0.0)) {
      ij_y_re = 0.0;
      ij_y_im = 0.0;
    } else {
      r = exp(ij_y_re / 2.0);
      ij_y_re = r * (r * cos(ij_y_im));
      ij_y_im = r * (r * sin(ij_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    jj_y_re = rt_powd_snf(t_s.re, 4.0);
    jj_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    jj_y_re = rt_powd_snf(t_s.im, 4.0);
    jj_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    jj_y_re = 4.0 * x_re;
    jj_y_im = 4.0 * x_im;
    if (jj_y_im == 0.0) {
      jj_y_re = exp(jj_y_re);
      jj_y_im = 0.0;
    } else if (rtIsInf(jj_y_im) && rtIsInf(jj_y_re) && (jj_y_re < 0.0)) {
      jj_y_re = 0.0;
      jj_y_im = 0.0;
    } else {
      r = exp(jj_y_re / 2.0);
      jj_y_re = r * (r * cos(jj_y_im));
      jj_y_im = r * (r * sin(jj_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    kj_y_re = rt_powd_snf(t.re, 5.0);
    kj_y_im = 0.0;
  } else if (t.re == 0.0) {
    kj_y_re = 0.0;
    kj_y_im = rt_powd_snf(t.im, 5.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    kj_y_re = 5.0 * x_re;
    kj_y_im = 5.0 * x_im;
    if (kj_y_im == 0.0) {
      kj_y_re = exp(kj_y_re);
      kj_y_im = 0.0;
    } else if (rtIsInf(kj_y_im) && rtIsInf(kj_y_re) && (kj_y_re < 0.0)) {
      kj_y_re = 0.0;
      kj_y_im = 0.0;
    } else {
      r = exp(kj_y_re / 2.0);
      kj_y_re = r * (r * cos(kj_y_im));
      kj_y_im = r * (r * sin(kj_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    lj_y_re = rt_powd_snf(t_s.re, 6.0);
    lj_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    lj_y_re = -rt_powd_snf(t_s.im, 6.0);
    lj_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    lj_y_re = 6.0 * x_re;
    lj_y_im = 6.0 * x_im;
    if (lj_y_im == 0.0) {
      lj_y_re = exp(lj_y_re);
      lj_y_im = 0.0;
    } else if (rtIsInf(lj_y_im) && rtIsInf(lj_y_re) && (lj_y_re < 0.0)) {
      lj_y_re = 0.0;
      lj_y_im = 0.0;
    } else {
      r = exp(lj_y_re / 2.0);
      lj_y_re = r * (r * cos(lj_y_im));
      lj_y_im = r * (r * sin(lj_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    mj_y_re = rt_powd_snf(t_s.re, 5.0);
    mj_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    mj_y_re = 0.0;
    mj_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    mj_y_re = 5.0 * x_re;
    mj_y_im = 5.0 * x_im;
    if (mj_y_im == 0.0) {
      mj_y_re = exp(mj_y_re);
      mj_y_im = 0.0;
    } else if (rtIsInf(mj_y_im) && rtIsInf(mj_y_re) && (mj_y_re < 0.0)) {
      mj_y_re = 0.0;
      mj_y_im = 0.0;
    } else {
      r = exp(mj_y_re / 2.0);
      mj_y_re = r * (r * cos(mj_y_im));
      mj_y_im = r * (r * sin(mj_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    nj_y_re = rt_powd_snf(t.re, 4.0);
    nj_y_im = 0.0;
  } else if (t.re == 0.0) {
    nj_y_re = rt_powd_snf(t.im, 4.0);
    nj_y_im = 0.0;
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    nj_y_re = 4.0 * x_re;
    nj_y_im = 4.0 * x_im;
    if (nj_y_im == 0.0) {
      nj_y_re = exp(nj_y_re);
      nj_y_im = 0.0;
    } else if (rtIsInf(nj_y_im) && rtIsInf(nj_y_re) && (nj_y_re < 0.0)) {
      nj_y_re = 0.0;
      nj_y_im = 0.0;
    } else {
      r = exp(nj_y_re / 2.0);
      nj_y_re = r * (r * cos(nj_y_im));
      nj_y_im = r * (r * sin(nj_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    oj_y_re = rt_powd_snf(t_s.re, 3.0);
    oj_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    oj_y_re = 0.0;
    oj_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    oj_y_re = 3.0 * x_re;
    oj_y_im = 3.0 * x_im;
    if (oj_y_im == 0.0) {
      oj_y_re = exp(oj_y_re);
      oj_y_im = 0.0;
    } else if (rtIsInf(oj_y_im) && rtIsInf(oj_y_re) && (oj_y_re < 0.0)) {
      oj_y_re = 0.0;
      oj_y_im = 0.0;
    } else {
      r = exp(oj_y_re / 2.0);
      oj_y_re = r * (r * cos(oj_y_im));
      oj_y_im = r * (r * sin(oj_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    pj_y_re = rt_powd_snf(t_s.re, 5.0);
    pj_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    pj_y_re = 0.0;
    pj_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    pj_y_re = 5.0 * x_re;
    pj_y_im = 5.0 * x_im;
    if (pj_y_im == 0.0) {
      pj_y_re = exp(pj_y_re);
      pj_y_im = 0.0;
    } else if (rtIsInf(pj_y_im) && rtIsInf(pj_y_re) && (pj_y_re < 0.0)) {
      pj_y_re = 0.0;
      pj_y_im = 0.0;
    } else {
      r = exp(pj_y_re / 2.0);
      pj_y_re = r * (r * cos(pj_y_im));
      pj_y_im = r * (r * sin(pj_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    qj_y_re = rt_powd_snf(t_s.re, 4.0);
    qj_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    qj_y_re = rt_powd_snf(t_s.im, 4.0);
    qj_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    qj_y_re = 4.0 * x_re;
    qj_y_im = 4.0 * x_im;
    if (qj_y_im == 0.0) {
      qj_y_re = exp(qj_y_re);
      qj_y_im = 0.0;
    } else if (rtIsInf(qj_y_im) && rtIsInf(qj_y_re) && (qj_y_re < 0.0)) {
      qj_y_re = 0.0;
      qj_y_im = 0.0;
    } else {
      r = exp(qj_y_re / 2.0);
      qj_y_re = r * (r * cos(qj_y_im));
      qj_y_im = r * (r * sin(qj_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    rj_y_re = rt_powd_snf(t.re, 5.0);
    rj_y_im = 0.0;
  } else if (t.re == 0.0) {
    rj_y_re = 0.0;
    rj_y_im = rt_powd_snf(t.im, 5.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    rj_y_re = 5.0 * x_re;
    rj_y_im = 5.0 * x_im;
    if (rj_y_im == 0.0) {
      rj_y_re = exp(rj_y_re);
      rj_y_im = 0.0;
    } else if (rtIsInf(rj_y_im) && rtIsInf(rj_y_re) && (rj_y_re < 0.0)) {
      rj_y_re = 0.0;
      rj_y_im = 0.0;
    } else {
      r = exp(rj_y_re / 2.0);
      rj_y_re = r * (r * cos(rj_y_im));
      rj_y_im = r * (r * sin(rj_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    sj_y_re = rt_powd_snf(t_s.re, 6.0);
    sj_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    sj_y_re = -rt_powd_snf(t_s.im, 6.0);
    sj_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    sj_y_re = 6.0 * x_re;
    sj_y_im = 6.0 * x_im;
    if (sj_y_im == 0.0) {
      sj_y_re = exp(sj_y_re);
      sj_y_im = 0.0;
    } else if (rtIsInf(sj_y_im) && rtIsInf(sj_y_re) && (sj_y_re < 0.0)) {
      sj_y_re = 0.0;
      sj_y_im = 0.0;
    } else {
      r = exp(sj_y_re / 2.0);
      sj_y_re = r * (r * cos(sj_y_im));
      sj_y_im = r * (r * sin(sj_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    tj_y_re = rt_powd_snf(t_s.re, 5.0);
    tj_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    tj_y_re = 0.0;
    tj_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    tj_y_re = 5.0 * x_re;
    tj_y_im = 5.0 * x_im;
    if (tj_y_im == 0.0) {
      tj_y_re = exp(tj_y_re);
      tj_y_im = 0.0;
    } else if (rtIsInf(tj_y_im) && rtIsInf(tj_y_re) && (tj_y_re < 0.0)) {
      tj_y_re = 0.0;
      tj_y_im = 0.0;
    } else {
      r = exp(tj_y_re / 2.0);
      tj_y_re = r * (r * cos(tj_y_im));
      tj_y_im = r * (r * sin(tj_y_im));
    }
  }

  if ((t.im == 0.0) && (t.re >= 0.0)) {
    uj_y_re = rt_powd_snf(t.re, 5.0);
    uj_y_im = 0.0;
  } else if (t.re == 0.0) {
    uj_y_re = 0.0;
    uj_y_im = rt_powd_snf(t.im, 5.0);
  } else {
    if (t.im == 0.0) {
      if (t.re < 0.0) {
        x_re = log(fabs(t.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t.re);
        x_im = 0.0;
      }
    } else if ((fabs(t.re) > 8.9884656743115785E+307) || (fabs(t.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t.re / 2.0, t.im / 2.0)) + 0.69314718055994529;
      x_im = rt_atan2d_snf(t.im, t.re);
    } else {
      x_re = log(rt_hypotd_snf(t.re, t.im));
      x_im = rt_atan2d_snf(t.im, t.re);
    }

    uj_y_re = 5.0 * x_re;
    uj_y_im = 5.0 * x_im;
    if (uj_y_im == 0.0) {
      uj_y_re = exp(uj_y_re);
      uj_y_im = 0.0;
    } else if (rtIsInf(uj_y_im) && rtIsInf(uj_y_re) && (uj_y_re < 0.0)) {
      uj_y_re = 0.0;
      uj_y_im = 0.0;
    } else {
      r = exp(uj_y_re / 2.0);
      uj_y_re = r * (r * cos(uj_y_im));
      uj_y_im = r * (r * sin(uj_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    vj_y_re = rt_powd_snf(t_s.re, 3.0);
    vj_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    vj_y_re = 0.0;
    vj_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    vj_y_re = 3.0 * x_re;
    vj_y_im = 3.0 * x_im;
    if (vj_y_im == 0.0) {
      vj_y_re = exp(vj_y_re);
      vj_y_im = 0.0;
    } else if (rtIsInf(vj_y_im) && rtIsInf(vj_y_re) && (vj_y_re < 0.0)) {
      vj_y_re = 0.0;
      vj_y_im = 0.0;
    } else {
      r = exp(vj_y_re / 2.0);
      vj_y_re = r * (r * cos(vj_y_im));
      vj_y_im = r * (r * sin(vj_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    wj_y_re = rt_powd_snf(t_s.re, 6.0);
    wj_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    wj_y_re = -rt_powd_snf(t_s.im, 6.0);
    wj_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    wj_y_re = 6.0 * x_re;
    wj_y_im = 6.0 * x_im;
    if (wj_y_im == 0.0) {
      wj_y_re = exp(wj_y_re);
      wj_y_im = 0.0;
    } else if (rtIsInf(wj_y_im) && rtIsInf(wj_y_re) && (wj_y_re < 0.0)) {
      wj_y_re = 0.0;
      wj_y_im = 0.0;
    } else {
      r = exp(wj_y_re / 2.0);
      wj_y_re = r * (r * cos(wj_y_im));
      wj_y_im = r * (r * sin(wj_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    xj_y_re = rt_powd_snf(t_s.re, 5.0);
    xj_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    xj_y_re = 0.0;
    xj_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    xj_y_re = 5.0 * x_re;
    xj_y_im = 5.0 * x_im;
    if (xj_y_im == 0.0) {
      xj_y_re = exp(xj_y_re);
      xj_y_im = 0.0;
    } else if (rtIsInf(xj_y_im) && rtIsInf(xj_y_re) && (xj_y_re < 0.0)) {
      xj_y_re = 0.0;
      xj_y_im = 0.0;
    } else {
      r = exp(xj_y_re / 2.0);
      xj_y_re = r * (r * cos(xj_y_im));
      xj_y_im = r * (r * sin(xj_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    yj_y_re = rt_powd_snf(t_s.re, 3.0);
    yj_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    yj_y_re = 0.0;
    yj_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    yj_y_re = 3.0 * x_re;
    yj_y_im = 3.0 * x_im;
    if (yj_y_im == 0.0) {
      yj_y_re = exp(yj_y_re);
      yj_y_im = 0.0;
    } else if (rtIsInf(yj_y_im) && rtIsInf(yj_y_re) && (yj_y_re < 0.0)) {
      yj_y_re = 0.0;
      yj_y_im = 0.0;
    } else {
      r = exp(yj_y_re / 2.0);
      yj_y_re = r * (r * cos(yj_y_im));
      yj_y_im = r * (r * sin(yj_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    ak_y_re = rt_powd_snf(t_s.re, 4.0);
    ak_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    ak_y_re = rt_powd_snf(t_s.im, 4.0);
    ak_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    ak_y_re = 4.0 * x_re;
    ak_y_im = 4.0 * x_im;
    if (ak_y_im == 0.0) {
      ak_y_re = exp(ak_y_re);
      ak_y_im = 0.0;
    } else if (rtIsInf(ak_y_im) && rtIsInf(ak_y_re) && (ak_y_re < 0.0)) {
      ak_y_re = 0.0;
      ak_y_im = 0.0;
    } else {
      r = exp(ak_y_re / 2.0);
      ak_y_re = r * (r * cos(ak_y_im));
      ak_y_im = r * (r * sin(ak_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    bk_y_re = rt_powd_snf(t_s.re, 5.0);
    bk_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    bk_y_re = 0.0;
    bk_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    bk_y_re = 5.0 * x_re;
    bk_y_im = 5.0 * x_im;
    if (bk_y_im == 0.0) {
      bk_y_re = exp(bk_y_re);
      bk_y_im = 0.0;
    } else if (rtIsInf(bk_y_im) && rtIsInf(bk_y_re) && (bk_y_re < 0.0)) {
      bk_y_re = 0.0;
      bk_y_im = 0.0;
    } else {
      r = exp(bk_y_re / 2.0);
      bk_y_re = r * (r * cos(bk_y_im));
      bk_y_im = r * (r * sin(bk_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    ck_y_re = rt_powd_snf(t_s.re, 5.0);
    ck_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    ck_y_re = 0.0;
    ck_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    ck_y_re = 5.0 * x_re;
    ck_y_im = 5.0 * x_im;
    if (ck_y_im == 0.0) {
      ck_y_re = exp(ck_y_re);
      ck_y_im = 0.0;
    } else if (rtIsInf(ck_y_im) && rtIsInf(ck_y_re) && (ck_y_re < 0.0)) {
      ck_y_re = 0.0;
      ck_y_im = 0.0;
    } else {
      r = exp(ck_y_re / 2.0);
      ck_y_re = r * (r * cos(ck_y_im));
      ck_y_im = r * (r * sin(ck_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    dk_y_re = rt_powd_snf(t_s.re, 3.0);
    dk_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    dk_y_re = 0.0;
    dk_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    dk_y_re = 3.0 * x_re;
    dk_y_im = 3.0 * x_im;
    if (dk_y_im == 0.0) {
      dk_y_re = exp(dk_y_re);
      dk_y_im = 0.0;
    } else if (rtIsInf(dk_y_im) && rtIsInf(dk_y_re) && (dk_y_re < 0.0)) {
      dk_y_re = 0.0;
      dk_y_im = 0.0;
    } else {
      r = exp(dk_y_re / 2.0);
      dk_y_re = r * (r * cos(dk_y_im));
      dk_y_im = r * (r * sin(dk_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    ek_y_re = rt_powd_snf(t_s.re, 4.0);
    ek_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    ek_y_re = rt_powd_snf(t_s.im, 4.0);
    ek_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    ek_y_re = 4.0 * x_re;
    ek_y_im = 4.0 * x_im;
    if (ek_y_im == 0.0) {
      ek_y_re = exp(ek_y_re);
      ek_y_im = 0.0;
    } else if (rtIsInf(ek_y_im) && rtIsInf(ek_y_re) && (ek_y_re < 0.0)) {
      ek_y_re = 0.0;
      ek_y_im = 0.0;
    } else {
      r = exp(ek_y_re / 2.0);
      ek_y_re = r * (r * cos(ek_y_im));
      ek_y_im = r * (r * sin(ek_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    fk_y_re = rt_powd_snf(t_s.re, 5.0);
    fk_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    fk_y_re = 0.0;
    fk_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    fk_y_re = 5.0 * x_re;
    fk_y_im = 5.0 * x_im;
    if (fk_y_im == 0.0) {
      fk_y_re = exp(fk_y_re);
      fk_y_im = 0.0;
    } else if (rtIsInf(fk_y_im) && rtIsInf(fk_y_re) && (fk_y_re < 0.0)) {
      fk_y_re = 0.0;
      fk_y_im = 0.0;
    } else {
      r = exp(fk_y_re / 2.0);
      fk_y_re = r * (r * cos(fk_y_im));
      fk_y_im = r * (r * sin(fk_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    gk_y_re = rt_powd_snf(t_s.re, 5.0);
    gk_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    gk_y_re = 0.0;
    gk_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    gk_y_re = 5.0 * x_re;
    gk_y_im = 5.0 * x_im;
    if (gk_y_im == 0.0) {
      gk_y_re = exp(gk_y_re);
      gk_y_im = 0.0;
    } else if (rtIsInf(gk_y_im) && rtIsInf(gk_y_re) && (gk_y_re < 0.0)) {
      gk_y_re = 0.0;
      gk_y_im = 0.0;
    } else {
      r = exp(gk_y_re / 2.0);
      gk_y_re = r * (r * cos(gk_y_im));
      gk_y_im = r * (r * sin(gk_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    hk_y_re = rt_powd_snf(t_s.re, 3.0);
    hk_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    hk_y_re = 0.0;
    hk_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    hk_y_re = 3.0 * x_re;
    hk_y_im = 3.0 * x_im;
    if (hk_y_im == 0.0) {
      hk_y_re = exp(hk_y_re);
      hk_y_im = 0.0;
    } else if (rtIsInf(hk_y_im) && rtIsInf(hk_y_re) && (hk_y_re < 0.0)) {
      hk_y_re = 0.0;
      hk_y_im = 0.0;
    } else {
      r = exp(hk_y_re / 2.0);
      hk_y_re = r * (r * cos(hk_y_im));
      hk_y_im = r * (r * sin(hk_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    ik_y_re = rt_powd_snf(t_s.re, 4.0);
    ik_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    ik_y_re = rt_powd_snf(t_s.im, 4.0);
    ik_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    ik_y_re = 4.0 * x_re;
    ik_y_im = 4.0 * x_im;
    if (ik_y_im == 0.0) {
      ik_y_re = exp(ik_y_re);
      ik_y_im = 0.0;
    } else if (rtIsInf(ik_y_im) && rtIsInf(ik_y_re) && (ik_y_re < 0.0)) {
      ik_y_re = 0.0;
      ik_y_im = 0.0;
    } else {
      r = exp(ik_y_re / 2.0);
      ik_y_re = r * (r * cos(ik_y_im));
      ik_y_im = r * (r * sin(ik_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    jk_y_re = rt_powd_snf(t_s.re, 5.0);
    jk_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    jk_y_re = 0.0;
    jk_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    jk_y_re = 5.0 * x_re;
    jk_y_im = 5.0 * x_im;
    if (jk_y_im == 0.0) {
      jk_y_re = exp(jk_y_re);
      jk_y_im = 0.0;
    } else if (rtIsInf(jk_y_im) && rtIsInf(jk_y_re) && (jk_y_re < 0.0)) {
      jk_y_re = 0.0;
      jk_y_im = 0.0;
    } else {
      r = exp(jk_y_re / 2.0);
      jk_y_re = r * (r * cos(jk_y_im));
      jk_y_im = r * (r * sin(jk_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    kk_y_re = rt_powd_snf(t_s.re, 5.0);
    kk_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    kk_y_re = 0.0;
    kk_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    kk_y_re = 5.0 * x_re;
    kk_y_im = 5.0 * x_im;
    if (kk_y_im == 0.0) {
      kk_y_re = exp(kk_y_re);
      kk_y_im = 0.0;
    } else if (rtIsInf(kk_y_im) && rtIsInf(kk_y_re) && (kk_y_re < 0.0)) {
      kk_y_re = 0.0;
      kk_y_im = 0.0;
    } else {
      r = exp(kk_y_re / 2.0);
      kk_y_re = r * (r * cos(kk_y_im));
      kk_y_im = r * (r * sin(kk_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    lk_y_re = rt_powd_snf(t_s.re, 5.0);
    lk_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    lk_y_re = 0.0;
    lk_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    lk_y_re = 5.0 * x_re;
    lk_y_im = 5.0 * x_im;
    if (lk_y_im == 0.0) {
      lk_y_re = exp(lk_y_re);
      lk_y_im = 0.0;
    } else if (rtIsInf(lk_y_im) && rtIsInf(lk_y_re) && (lk_y_re < 0.0)) {
      lk_y_re = 0.0;
      lk_y_im = 0.0;
    } else {
      r = exp(lk_y_re / 2.0);
      lk_y_re = r * (r * cos(lk_y_im));
      lk_y_im = r * (r * sin(lk_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    mk_y_re = rt_powd_snf(t_s.re, 4.0);
    mk_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    mk_y_re = rt_powd_snf(t_s.im, 4.0);
    mk_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    mk_y_re = 4.0 * x_re;
    mk_y_im = 4.0 * x_im;
    if (mk_y_im == 0.0) {
      mk_y_re = exp(mk_y_re);
      mk_y_im = 0.0;
    } else if (rtIsInf(mk_y_im) && rtIsInf(mk_y_re) && (mk_y_re < 0.0)) {
      mk_y_re = 0.0;
      mk_y_im = 0.0;
    } else {
      r = exp(mk_y_re / 2.0);
      mk_y_re = r * (r * cos(mk_y_im));
      mk_y_im = r * (r * sin(mk_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    nk_y_re = rt_powd_snf(t_s.re, 3.0);
    nk_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    nk_y_re = 0.0;
    nk_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    nk_y_re = 3.0 * x_re;
    nk_y_im = 3.0 * x_im;
    if (nk_y_im == 0.0) {
      nk_y_re = exp(nk_y_re);
      nk_y_im = 0.0;
    } else if (rtIsInf(nk_y_im) && rtIsInf(nk_y_re) && (nk_y_re < 0.0)) {
      nk_y_re = 0.0;
      nk_y_im = 0.0;
    } else {
      r = exp(nk_y_re / 2.0);
      nk_y_re = r * (r * cos(nk_y_im));
      nk_y_im = r * (r * sin(nk_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    ok_y_re = rt_powd_snf(t_s.re, 5.0);
    ok_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    ok_y_re = 0.0;
    ok_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    ok_y_re = 5.0 * x_re;
    ok_y_im = 5.0 * x_im;
    if (ok_y_im == 0.0) {
      ok_y_re = exp(ok_y_re);
      ok_y_im = 0.0;
    } else if (rtIsInf(ok_y_im) && rtIsInf(ok_y_re) && (ok_y_re < 0.0)) {
      ok_y_re = 0.0;
      ok_y_im = 0.0;
    } else {
      r = exp(ok_y_re / 2.0);
      ok_y_re = r * (r * cos(ok_y_im));
      ok_y_im = r * (r * sin(ok_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    pk_y_re = rt_powd_snf(t_s.re, 5.0);
    pk_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    pk_y_re = 0.0;
    pk_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    pk_y_re = 5.0 * x_re;
    pk_y_im = 5.0 * x_im;
    if (pk_y_im == 0.0) {
      pk_y_re = exp(pk_y_re);
      pk_y_im = 0.0;
    } else if (rtIsInf(pk_y_im) && rtIsInf(pk_y_re) && (pk_y_re < 0.0)) {
      pk_y_re = 0.0;
      pk_y_im = 0.0;
    } else {
      r = exp(pk_y_re / 2.0);
      pk_y_re = r * (r * cos(pk_y_im));
      pk_y_im = r * (r * sin(pk_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    qk_y_re = rt_powd_snf(t_s.re, 4.0);
    qk_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    qk_y_re = rt_powd_snf(t_s.im, 4.0);
    qk_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    qk_y_re = 4.0 * x_re;
    qk_y_im = 4.0 * x_im;
    if (qk_y_im == 0.0) {
      qk_y_re = exp(qk_y_re);
      qk_y_im = 0.0;
    } else if (rtIsInf(qk_y_im) && rtIsInf(qk_y_re) && (qk_y_re < 0.0)) {
      qk_y_re = 0.0;
      qk_y_im = 0.0;
    } else {
      r = exp(qk_y_re / 2.0);
      qk_y_re = r * (r * cos(qk_y_im));
      qk_y_im = r * (r * sin(qk_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    rk_y_re = rt_powd_snf(t_s.re, 3.0);
    rk_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    rk_y_re = 0.0;
    rk_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    rk_y_re = 3.0 * x_re;
    rk_y_im = 3.0 * x_im;
    if (rk_y_im == 0.0) {
      rk_y_re = exp(rk_y_re);
      rk_y_im = 0.0;
    } else if (rtIsInf(rk_y_im) && rtIsInf(rk_y_re) && (rk_y_re < 0.0)) {
      rk_y_re = 0.0;
      rk_y_im = 0.0;
    } else {
      r = exp(rk_y_re / 2.0);
      rk_y_re = r * (r * cos(rk_y_im));
      rk_y_im = r * (r * sin(rk_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    sk_y_re = rt_powd_snf(t_s.re, 5.0);
    sk_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    sk_y_re = 0.0;
    sk_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    sk_y_re = 5.0 * x_re;
    sk_y_im = 5.0 * x_im;
    if (sk_y_im == 0.0) {
      sk_y_re = exp(sk_y_re);
      sk_y_im = 0.0;
    } else if (rtIsInf(sk_y_im) && rtIsInf(sk_y_re) && (sk_y_re < 0.0)) {
      sk_y_re = 0.0;
      sk_y_im = 0.0;
    } else {
      r = exp(sk_y_re / 2.0);
      sk_y_re = r * (r * cos(sk_y_im));
      sk_y_im = r * (r * sin(sk_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    tk_y_re = rt_powd_snf(t_s.re, 5.0);
    tk_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    tk_y_re = 0.0;
    tk_y_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    tk_y_re = 5.0 * x_re;
    tk_y_im = 5.0 * x_im;
    if (tk_y_im == 0.0) {
      tk_y_re = exp(tk_y_re);
      tk_y_im = 0.0;
    } else if (rtIsInf(tk_y_im) && rtIsInf(tk_y_re) && (tk_y_re < 0.0)) {
      tk_y_re = 0.0;
      tk_y_im = 0.0;
    } else {
      r = exp(tk_y_re / 2.0);
      tk_y_re = r * (r * cos(tk_y_im));
      tk_y_im = r * (r * sin(tk_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    uk_y_re = rt_powd_snf(t_s.re, 4.0);
    uk_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    uk_y_re = rt_powd_snf(t_s.im, 4.0);
    uk_y_im = 0.0;
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    uk_y_re = 4.0 * x_re;
    uk_y_im = 4.0 * x_im;
    if (uk_y_im == 0.0) {
      uk_y_re = exp(uk_y_re);
      uk_y_im = 0.0;
    } else if (rtIsInf(uk_y_im) && rtIsInf(uk_y_re) && (uk_y_re < 0.0)) {
      uk_y_re = 0.0;
      uk_y_im = 0.0;
    } else {
      r = exp(uk_y_re / 2.0);
      uk_y_re = r * (r * cos(uk_y_im));
      uk_y_im = r * (r * sin(uk_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    vk_y_re = rt_powd_snf(t_s.re, 3.0);
    vk_y_im = 0.0;
  } else if (t_s.re == 0.0) {
    vk_y_re = 0.0;
    vk_y_im = -rt_powd_snf(t_s.im, 3.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    vk_y_re = 3.0 * x_re;
    vk_y_im = 3.0 * x_im;
    if (vk_y_im == 0.0) {
      vk_y_re = exp(vk_y_re);
      vk_y_im = 0.0;
    } else if (rtIsInf(vk_y_im) && rtIsInf(vk_y_re) && (vk_y_re < 0.0)) {
      vk_y_re = 0.0;
      vk_y_im = 0.0;
    } else {
      r = exp(vk_y_re / 2.0);
      vk_y_re = r * (r * cos(vk_y_im));
      vk_y_im = r * (r * sin(vk_y_im));
    }
  }

  if ((t_s.im == 0.0) && (t_s.re >= 0.0)) {
    x_re = rt_powd_snf(t_s.re, 5.0);
    x_im = 0.0;
  } else if (t_s.re == 0.0) {
    x_re = 0.0;
    x_im = rt_powd_snf(t_s.im, 5.0);
  } else {
    if (t_s.im == 0.0) {
      if (t_s.re < 0.0) {
        x_re = log(fabs(t_s.re));
        x_im = 3.1415926535897931;
      } else {
        x_re = log(t_s.re);
        x_im = 0.0;
      }
    } else if ((fabs(t_s.re) > 8.9884656743115785E+307) || (fabs(t_s.im) >
                8.9884656743115785E+307)) {
      x_re = log(rt_hypotd_snf(t_s.re / 2.0, t_s.im / 2.0)) +
        0.69314718055994529;
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    } else {
      x_re = log(rt_hypotd_snf(t_s.re, t_s.im));
      x_im = rt_atan2d_snf(t_s.im, t_s.re);
    }

    x_re *= 5.0;
    x_im *= 5.0;
    if (x_im == 0.0) {
      x_re = exp(x_re);
      x_im = 0.0;
    } else if (rtIsInf(x_im) && rtIsInf(x_re) && (x_re < 0.0)) {
      x_re = 0.0;
      x_im = 0.0;
    } else {
      r = exp(x_re / 2.0);
      x_re = r * (r * cos(x_im));
      x_im = r * (r * sin(x_im));
    }
  }

  t_re_tmp = t.re * t.re - t.im * t.im;
  t_im_tmp = t.re * t.im + t.im * t.re;
  ar = t_re_tmp * x07;
  if (t_im_tmp * x07 == 0.0) {
    t_re = ar / 2.0;
  } else if (ar == 0.0) {
    t_re = 0.0;
  } else {
    t_re = ar / 2.0;
  }

  t_s_re_tmp = t_s.re * t_s.re - t_s.im * t_s.im;
  t_s_im_tmp = t_s.re * t_s.im + t_s.im * t_s.re;
  ar = 625.0 * y_re * x01;
  ai = 625.0 * y_im * x01;
  y_re = 500.0 * b_y_re + 9.0 * t_s_re_tmp;
  r = 500.0 * b_y_im + 9.0 * t_s_im_tmp;
  if (r == 0.0) {
    if (ai == 0.0) {
      re = ar / y_re;
    } else if (ar == 0.0) {
      re = 0.0;
    } else {
      re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      re = ai / r;
    } else if (ai == 0.0) {
      re = 0.0;
    } else {
      re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        wk_y_re = 0.5;
      } else {
        wk_y_re = -0.5;
      }

      if (r > 0.0) {
        b_r = 0.5;
      } else {
        b_r = -0.5;
      }

      re = (ar * wk_y_re + ai * b_r) / brm;
    } else {
      brm = y_re / r;
      re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  ar = 15.0 * c_y_re * x01;
  ai = 15.0 * c_y_im * x01;
  y_re = 500.0 * d_y_re + 9.0 * e_y_re;
  r = 500.0 * d_y_im + 9.0 * e_y_im;
  if (r == 0.0) {
    if (ai == 0.0) {
      b_re = ar / y_re;
    } else if (ar == 0.0) {
      b_re = 0.0;
    } else {
      b_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      b_re = ai / r;
    } else if (ai == 0.0) {
      b_re = 0.0;
    } else {
      b_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      b_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        xk_y_re = 0.5;
      } else {
        xk_y_re = -0.5;
      }

      if (r > 0.0) {
        c_r = 0.5;
      } else {
        c_r = -0.5;
      }

      b_re = (ar * xk_y_re + ai * c_r) / brm;
    } else {
      brm = y_re / r;
      b_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  ar = 15.0 * f_y_re * x01;
  ai = 15.0 * f_y_im * x01;
  y_re = 2.0 * (500.0 * g_y_re + 9.0 * h_y_re);
  r = 2.0 * (500.0 * g_y_im + 9.0 * h_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      c_re = ar / y_re;
    } else if (ar == 0.0) {
      c_re = 0.0;
    } else {
      c_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      c_re = ai / r;
    } else if (ai == 0.0) {
      c_re = 0.0;
    } else {
      c_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      c_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        yk_y_re = 0.5;
      } else {
        yk_y_re = -0.5;
      }

      if (r > 0.0) {
        d_r = 0.5;
      } else {
        d_r = -0.5;
      }

      c_re = (ar * yk_y_re + ai * d_r) / brm;
    } else {
      brm = y_re / r;
      c_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  ar = 3.0 * i_y_re * x01;
  ai = 3.0 * i_y_im * x01;
  y_re = 2.0 * (500.0 * j_y_re + 9.0 * k_y_re);
  r = 2.0 * (500.0 * j_y_im + 9.0 * k_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      d_re = ar / y_re;
    } else if (ar == 0.0) {
      d_re = 0.0;
    } else {
      d_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      d_re = ai / r;
    } else if (ai == 0.0) {
      d_re = 0.0;
    } else {
      d_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      d_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        al_y_re = 0.5;
      } else {
        al_y_re = -0.5;
      }

      if (r > 0.0) {
        e_r = 0.5;
      } else {
        e_r = -0.5;
      }

      d_re = (ar * al_y_re + ai * e_r) / brm;
    } else {
      brm = y_re / r;
      d_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  ar = 625.0 * l_y_re * x11;
  ai = 625.0 * l_y_im * x11;
  y_re = 500.0 * m_y_re + 9.0 * t_s_re_tmp;
  r = 500.0 * m_y_im + 9.0 * t_s_im_tmp;
  if (r == 0.0) {
    if (ai == 0.0) {
      m_y_im = ar / y_re;
    } else if (ar == 0.0) {
      m_y_im = 0.0;
    } else {
      m_y_im = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      m_y_im = ai / r;
    } else if (ai == 0.0) {
      m_y_im = 0.0;
    } else {
      m_y_im = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      m_y_im = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        bl_y_re = 0.5;
      } else {
        bl_y_re = -0.5;
      }

      if (r > 0.0) {
        f_r = 0.5;
      } else {
        f_r = -0.5;
      }

      m_y_im = (ar * bl_y_re + ai * f_r) / brm;
    } else {
      brm = y_re / r;
      m_y_im = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  ar = 15.0 * n_y_re * x11;
  ai = 15.0 * n_y_im * x11;
  y_re = 500.0 * o_y_re + 9.0 * p_y_re;
  r = 500.0 * o_y_im + 9.0 * p_y_im;
  if (r == 0.0) {
    if (ai == 0.0) {
      m_y_re = ar / y_re;
    } else if (ar == 0.0) {
      m_y_re = 0.0;
    } else {
      m_y_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      m_y_re = ai / r;
    } else if (ai == 0.0) {
      m_y_re = 0.0;
    } else {
      m_y_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      m_y_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        cl_y_re = 0.5;
      } else {
        cl_y_re = -0.5;
      }

      if (r > 0.0) {
        g_r = 0.5;
      } else {
        g_r = -0.5;
      }

      m_y_re = (ar * cl_y_re + ai * g_r) / brm;
    } else {
      brm = y_re / r;
      m_y_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  ar = 15.0 * q_y_re * x11;
  ai = 15.0 * q_y_im * x11;
  y_re = 2.0 * (500.0 * r_y_re + 9.0 * s_y_re);
  r = 2.0 * (500.0 * r_y_im + 9.0 * s_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      l_y_im = ar / y_re;
    } else if (ar == 0.0) {
      l_y_im = 0.0;
    } else {
      l_y_im = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      l_y_im = ai / r;
    } else if (ai == 0.0) {
      l_y_im = 0.0;
    } else {
      l_y_im = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      l_y_im = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        dl_y_re = 0.5;
      } else {
        dl_y_re = -0.5;
      }

      if (r > 0.0) {
        h_r = 0.5;
      } else {
        h_r = -0.5;
      }

      l_y_im = (ar * dl_y_re + ai * h_r) / brm;
    } else {
      brm = y_re / r;
      l_y_im = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  ar = 3.0 * t_y_re * x11;
  ai = 3.0 * t_y_im * x11;
  y_re = 2.0 * (500.0 * u_y_re + 9.0 * v_y_re);
  r = 2.0 * (500.0 * u_y_im + 9.0 * v_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      l_y_re = ar / y_re;
    } else if (ar == 0.0) {
      l_y_re = 0.0;
    } else {
      l_y_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      l_y_re = ai / r;
    } else if (ai == 0.0) {
      l_y_re = 0.0;
    } else {
      l_y_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      l_y_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        el_y_re = 0.5;
      } else {
        el_y_re = -0.5;
      }

      if (r > 0.0) {
        i_r = 0.5;
      } else {
        i_r = -0.5;
      }

      l_y_re = (ar * el_y_re + ai * i_r) / brm;
    } else {
      brm = y_re / r;
      l_y_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  k_y_re = 1875.0 * w_y_re;
  k_y_im = 1875.0 * w_y_im;
  ar = (k_y_re * t_s.re - k_y_im * t_s.im) * x01;
  ai = (k_y_re * t_s.im + k_y_im * t_s.re) * x01;
  y_re = 500.0 * x_y_re + 9.0 * y_y_re;
  r = 500.0 * x_y_im + 9.0 * y_y_im;
  if (r == 0.0) {
    if (ai == 0.0) {
      k_y_re = ar / y_re;
    } else if (ar == 0.0) {
      k_y_re = 0.0;
    } else {
      k_y_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      k_y_re = ai / r;
    } else if (ai == 0.0) {
      k_y_re = 0.0;
    } else {
      k_y_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      k_y_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        fl_y_re = 0.5;
      } else {
        fl_y_re = -0.5;
      }

      if (r > 0.0) {
        j_r = 0.5;
      } else {
        j_r = -0.5;
      }

      k_y_re = (ar * fl_y_re + ai * j_r) / brm;
    } else {
      brm = y_re / r;
      k_y_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  j_y_im = 625.0 * ab_y_re;
  k_y_im = 625.0 * ab_y_im;
  ar = (j_y_im * t_s.re - k_y_im * t_s.im) * x04;
  ai = (j_y_im * t_s.im + k_y_im * t_s.re) * x04;
  y_re = 500.0 * bb_y_re + 9.0 * t_s_re_tmp;
  r = 500.0 * bb_y_im + 9.0 * t_s_im_tmp;
  if (r == 0.0) {
    if (ai == 0.0) {
      j_y_im = ar / y_re;
    } else if (ar == 0.0) {
      j_y_im = 0.0;
    } else {
      j_y_im = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      j_y_im = ai / r;
    } else if (ai == 0.0) {
      j_y_im = 0.0;
    } else {
      j_y_im = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      j_y_im = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        gl_y_re = 0.5;
      } else {
        gl_y_re = -0.5;
      }

      if (r > 0.0) {
        k_r = 0.5;
      } else {
        k_r = -0.5;
      }

      j_y_im = (ar * gl_y_re + ai * k_r) / brm;
    } else {
      brm = y_re / r;
      j_y_im = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  j_y_re = 15.0 * cb_y_re;
  k_y_im = 15.0 * cb_y_im;
  ar = (j_y_re * t_s.re - k_y_im * t_s.im) * x04;
  ai = (j_y_re * t_s.im + k_y_im * t_s.re) * x04;
  y_re = 500.0 * db_y_re + 9.0 * eb_y_re;
  r = 500.0 * db_y_im + 9.0 * eb_y_im;
  if (r == 0.0) {
    if (ai == 0.0) {
      j_y_re = ar / y_re;
    } else if (ar == 0.0) {
      j_y_re = 0.0;
    } else {
      j_y_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      j_y_re = ai / r;
    } else if (ai == 0.0) {
      j_y_re = 0.0;
    } else {
      j_y_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      j_y_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        hl_y_re = 0.5;
      } else {
        hl_y_re = -0.5;
      }

      if (r > 0.0) {
        l_r = 0.5;
      } else {
        l_r = -0.5;
      }

      j_y_re = (ar * hl_y_re + ai * l_r) / brm;
    } else {
      brm = y_re / r;
      j_y_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  i_y_im = 1875.0 * fb_y_re;
  k_y_im = 1875.0 * fb_y_im;
  ar = (i_y_im * t_s.re - k_y_im * t_s.im) * x01;
  ai = (i_y_im * t_s.im + k_y_im * t_s.re) * x01;
  y_re = 2.0 * (500.0 * gb_y_re + 9.0 * hb_y_re);
  r = 2.0 * (500.0 * gb_y_im + 9.0 * hb_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      i_y_im = ar / y_re;
    } else if (ar == 0.0) {
      i_y_im = 0.0;
    } else {
      i_y_im = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      i_y_im = ai / r;
    } else if (ai == 0.0) {
      i_y_im = 0.0;
    } else {
      i_y_im = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      i_y_im = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        il_y_re = 0.5;
      } else {
        il_y_re = -0.5;
      }

      if (r > 0.0) {
        m_r = 0.5;
      } else {
        m_r = -0.5;
      }

      i_y_im = (ar * il_y_re + ai * m_r) / brm;
    } else {
      brm = y_re / r;
      i_y_im = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  i_y_re = 15.0 * ib_y_re;
  k_y_im = 15.0 * ib_y_im;
  ar = (i_y_re * t_s.re - k_y_im * t_s.im) * x04;
  ai = (i_y_re * t_s.im + k_y_im * t_s.re) * x04;
  y_re = 2.0 * (500.0 * jb_y_re + 9.0 * kb_y_re);
  r = 2.0 * (500.0 * jb_y_im + 9.0 * kb_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      i_y_re = ar / y_re;
    } else if (ar == 0.0) {
      i_y_re = 0.0;
    } else {
      i_y_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      i_y_re = ai / r;
    } else if (ai == 0.0) {
      i_y_re = 0.0;
    } else {
      i_y_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      i_y_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        jl_y_re = 0.5;
      } else {
        jl_y_re = -0.5;
      }

      if (r > 0.0) {
        n_r = 0.5;
      } else {
        n_r = -0.5;
      }

      i_y_re = (ar * jl_y_re + ai * n_r) / brm;
    } else {
      brm = y_re / r;
      i_y_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  h_y_im = 375.0 * lb_y_re;
  k_y_im = 375.0 * lb_y_im;
  ar = (h_y_im * t_s.re - k_y_im * t_s.im) * x01;
  ai = (h_y_im * t_s.im + k_y_im * t_s.re) * x01;
  y_re = 2.0 * (500.0 * mb_y_re + 9.0 * nb_y_re);
  r = 2.0 * (500.0 * mb_y_im + 9.0 * nb_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      h_y_im = ar / y_re;
    } else if (ar == 0.0) {
      h_y_im = 0.0;
    } else {
      h_y_im = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      h_y_im = ai / r;
    } else if (ai == 0.0) {
      h_y_im = 0.0;
    } else {
      h_y_im = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      h_y_im = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        kl_y_re = 0.5;
      } else {
        kl_y_re = -0.5;
      }

      if (r > 0.0) {
        o_r = 0.5;
      } else {
        o_r = -0.5;
      }

      h_y_im = (ar * kl_y_re + ai * o_r) / brm;
    } else {
      brm = y_re / r;
      h_y_im = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  h_y_re = 3.0 * ob_y_re;
  k_y_im = 3.0 * ob_y_im;
  ar = (h_y_re * t_s.re - k_y_im * t_s.im) * x04;
  ai = (h_y_re * t_s.im + k_y_im * t_s.re) * x04;
  y_re = 2.0 * (500.0 * pb_y_re + 9.0 * qb_y_re);
  r = 2.0 * (500.0 * pb_y_im + 9.0 * qb_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      h_y_re = ar / y_re;
    } else if (ar == 0.0) {
      h_y_re = 0.0;
    } else {
      h_y_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      h_y_re = ai / r;
    } else if (ai == 0.0) {
      h_y_re = 0.0;
    } else {
      h_y_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      h_y_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        ll_y_re = 0.5;
      } else {
        ll_y_re = -0.5;
      }

      if (r > 0.0) {
        p_r = 0.5;
      } else {
        p_r = -0.5;
      }

      h_y_re = (ar * ll_y_re + ai * p_r) / brm;
    } else {
      brm = y_re / r;
      h_y_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  g_y_im = 1875.0 * rb_y_re;
  k_y_im = 1875.0 * rb_y_im;
  ar = (g_y_im * t_s.re - k_y_im * t_s.im) * x11;
  ai = (g_y_im * t_s.im + k_y_im * t_s.re) * x11;
  y_re = 500.0 * sb_y_re + 9.0 * tb_y_re;
  r = 500.0 * sb_y_im + 9.0 * tb_y_im;
  if (r == 0.0) {
    if (ai == 0.0) {
      g_y_im = ar / y_re;
    } else if (ar == 0.0) {
      g_y_im = 0.0;
    } else {
      g_y_im = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      g_y_im = ai / r;
    } else if (ai == 0.0) {
      g_y_im = 0.0;
    } else {
      g_y_im = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      g_y_im = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        ml_y_re = 0.5;
      } else {
        ml_y_re = -0.5;
      }

      if (r > 0.0) {
        q_r = 0.5;
      } else {
        q_r = -0.5;
      }

      g_y_im = (ar * ml_y_re + ai * q_r) / brm;
    } else {
      brm = y_re / r;
      g_y_im = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  g_y_re = 1875.0 * ub_y_re;
  k_y_im = 1875.0 * ub_y_im;
  ar = (g_y_re * t_s.re - k_y_im * t_s.im) * x11;
  ai = (g_y_re * t_s.im + k_y_im * t_s.re) * x11;
  y_re = 2.0 * (500.0 * vb_y_re + 9.0 * wb_y_re);
  r = 2.0 * (500.0 * vb_y_im + 9.0 * wb_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      g_y_re = ar / y_re;
    } else if (ar == 0.0) {
      g_y_re = 0.0;
    } else {
      g_y_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      g_y_re = ai / r;
    } else if (ai == 0.0) {
      g_y_re = 0.0;
    } else {
      g_y_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      g_y_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        nl_y_re = 0.5;
      } else {
        nl_y_re = -0.5;
      }

      if (r > 0.0) {
        r_r = 0.5;
      } else {
        r_r = -0.5;
      }

      g_y_re = (ar * nl_y_re + ai * r_r) / brm;
    } else {
      brm = y_re / r;
      g_y_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  f_y_im = 375.0 * xb_y_re;
  k_y_im = 375.0 * xb_y_im;
  ar = (f_y_im * t_s.re - k_y_im * t_s.im) * x11;
  ai = (f_y_im * t_s.im + k_y_im * t_s.re) * x11;
  y_re = 2.0 * (500.0 * yb_y_re + 9.0 * ac_y_re);
  r = 2.0 * (500.0 * yb_y_im + 9.0 * ac_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      f_y_im = ar / y_re;
    } else if (ar == 0.0) {
      f_y_im = 0.0;
    } else {
      f_y_im = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      f_y_im = ai / r;
    } else if (ai == 0.0) {
      f_y_im = 0.0;
    } else {
      f_y_im = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      f_y_im = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        ol_y_re = 0.5;
      } else {
        ol_y_re = -0.5;
      }

      if (r > 0.0) {
        s_r = 0.5;
      } else {
        s_r = -0.5;
      }

      f_y_im = (ar * ol_y_re + ai * s_r) / brm;
    } else {
      brm = y_re / r;
      f_y_im = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  f_y_re = 1875.0 * bc_y_re;
  k_y_im = 1875.0 * bc_y_im;
  ar = (f_y_re * t_s_re_tmp - k_y_im * t_s_im_tmp) * x04;
  ai = (f_y_re * t_s_im_tmp + k_y_im * t_s_re_tmp) * x04;
  y_re = 500.0 * cc_y_re + 9.0 * dc_y_re;
  r = 500.0 * cc_y_im + 9.0 * dc_y_im;
  if (r == 0.0) {
    if (ai == 0.0) {
      f_y_re = ar / y_re;
    } else if (ar == 0.0) {
      f_y_re = 0.0;
    } else {
      f_y_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      f_y_re = ai / r;
    } else if (ai == 0.0) {
      f_y_re = 0.0;
    } else {
      f_y_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      f_y_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        pl_y_re = 0.5;
      } else {
        pl_y_re = -0.5;
      }

      if (r > 0.0) {
        t_r = 0.5;
      } else {
        t_r = -0.5;
      }

      f_y_re = (ar * pl_y_re + ai * t_r) / brm;
    } else {
      brm = y_re / r;
      f_y_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  e_y_im = 125.0 * ec_y_re;
  k_y_im = 125.0 * ec_y_im;
  ar = (e_y_im * t_s_re_tmp - k_y_im * t_s_im_tmp) * x07;
  ai = (e_y_im * t_s_im_tmp + k_y_im * t_s_re_tmp) * x07;
  y_re = 500.0 * fc_y_re + 9.0 * t_s_re_tmp;
  r = 500.0 * fc_y_im + 9.0 * t_s_im_tmp;
  if (r == 0.0) {
    if (ai == 0.0) {
      e_y_im = ar / y_re;
    } else if (ar == 0.0) {
      e_y_im = 0.0;
    } else {
      e_y_im = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      e_y_im = ai / r;
    } else if (ai == 0.0) {
      e_y_im = 0.0;
    } else {
      e_y_im = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      e_y_im = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        ql_y_re = 0.5;
      } else {
        ql_y_re = -0.5;
      }

      if (r > 0.0) {
        u_r = 0.5;
      } else {
        u_r = -0.5;
      }

      e_y_im = (ar * ql_y_re + ai * u_r) / brm;
    } else {
      brm = y_re / r;
      e_y_im = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  e_y_re = 15.0 * gc_y_re;
  k_y_im = 15.0 * gc_y_im;
  ar = (e_y_re * t_s_re_tmp - k_y_im * t_s_im_tmp) * x07;
  ai = (e_y_re * t_s_im_tmp + k_y_im * t_s_re_tmp) * x07;
  y_re = 2.0 * (500.0 * hc_y_re + 9.0 * ic_y_re);
  r = 2.0 * (500.0 * hc_y_im + 9.0 * ic_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      e_y_re = ar / y_re;
    } else if (ar == 0.0) {
      e_y_re = 0.0;
    } else {
      e_y_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      e_y_re = ai / r;
    } else if (ai == 0.0) {
      e_y_re = 0.0;
    } else {
      e_y_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      e_y_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        rl_y_re = 0.5;
      } else {
        rl_y_re = -0.5;
      }

      if (r > 0.0) {
        v_r = 0.5;
      } else {
        v_r = -0.5;
      }

      e_y_re = (ar * rl_y_re + ai * v_r) / brm;
    } else {
      brm = y_re / r;
      e_y_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  d_y_im = 1875.0 * jc_y_re;
  k_y_im = 1875.0 * jc_y_im;
  ar = (d_y_im * t_s_re_tmp - k_y_im * t_s_im_tmp) * x04;
  ai = (d_y_im * t_s_im_tmp + k_y_im * t_s_re_tmp) * x04;
  y_re = 2.0 * (500.0 * kc_y_re + 9.0 * lc_y_re);
  r = 2.0 * (500.0 * kc_y_im + 9.0 * lc_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      d_y_im = ar / y_re;
    } else if (ar == 0.0) {
      d_y_im = 0.0;
    } else {
      d_y_im = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      d_y_im = ai / r;
    } else if (ai == 0.0) {
      d_y_im = 0.0;
    } else {
      d_y_im = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      d_y_im = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        sl_y_re = 0.5;
      } else {
        sl_y_re = -0.5;
      }

      if (r > 0.0) {
        w_r = 0.5;
      } else {
        w_r = -0.5;
      }

      d_y_im = (ar * sl_y_re + ai * w_r) / brm;
    } else {
      brm = y_re / r;
      d_y_im = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  d_y_re = 625.0 * mc_y_re;
  k_y_im = 625.0 * mc_y_im;
  ar = (d_y_re * nc_y_re - k_y_im * nc_y_im) * x07;
  ai = (d_y_re * nc_y_im + k_y_im * nc_y_re) * x07;
  y_re = 500.0 * oc_y_re + 9.0 * pc_y_re;
  r = 500.0 * oc_y_im + 9.0 * pc_y_im;
  if (r == 0.0) {
    if (ai == 0.0) {
      d_y_re = ar / y_re;
    } else if (ar == 0.0) {
      d_y_re = 0.0;
    } else {
      d_y_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      d_y_re = ai / r;
    } else if (ai == 0.0) {
      d_y_re = 0.0;
    } else {
      d_y_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      d_y_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        tl_y_re = 0.5;
      } else {
        tl_y_re = -0.5;
      }

      if (r > 0.0) {
        x_r = 0.5;
      } else {
        x_r = -0.5;
      }

      d_y_re = (ar * tl_y_re + ai * x_r) / brm;
    } else {
      brm = y_re / r;
      d_y_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  c_y_im = 15.0 * qc_y_re;
  k_y_im = 15.0 * qc_y_im;
  ar = (c_y_im * t_s_re_tmp - k_y_im * t_s_im_tmp) * x07;
  ai = (c_y_im * t_s_im_tmp + k_y_im * t_s_re_tmp) * x07;
  y_re = 4.0 * (500.0 * rc_y_re + 9.0 * sc_y_re);
  r = 4.0 * (500.0 * rc_y_im + 9.0 * sc_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      c_y_im = ar / y_re;
    } else if (ar == 0.0) {
      c_y_im = 0.0;
    } else {
      c_y_im = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      c_y_im = ai / r;
    } else if (ai == 0.0) {
      c_y_im = 0.0;
    } else {
      c_y_im = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      c_y_im = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        ul_y_re = 0.5;
      } else {
        ul_y_re = -0.5;
      }

      if (r > 0.0) {
        y_r = 0.5;
      } else {
        y_r = -0.5;
      }

      c_y_im = (ar * ul_y_re + ai * y_r) / brm;
    } else {
      brm = y_re / r;
      c_y_im = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  c_y_re = 375.0 * tc_y_re;
  k_y_im = 375.0 * tc_y_im;
  ar = (c_y_re * t_s_re_tmp - k_y_im * t_s_im_tmp) * x04;
  ai = (c_y_re * t_s_im_tmp + k_y_im * t_s_re_tmp) * x04;
  y_re = 2.0 * (500.0 * uc_y_re + 9.0 * vc_y_re);
  r = 2.0 * (500.0 * uc_y_im + 9.0 * vc_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      c_y_re = ar / y_re;
    } else if (ar == 0.0) {
      c_y_re = 0.0;
    } else {
      c_y_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      c_y_re = ai / r;
    } else if (ai == 0.0) {
      c_y_re = 0.0;
    } else {
      c_y_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      c_y_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        vl_y_re = 0.5;
      } else {
        vl_y_re = -0.5;
      }

      if (r > 0.0) {
        ab_r = 0.5;
      } else {
        ab_r = -0.5;
      }

      c_y_re = (ar * vl_y_re + ai * ab_r) / brm;
    } else {
      brm = y_re / r;
      c_y_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  b_y_im = 625.0 * wc_y_re;
  k_y_im = 625.0 * wc_y_im;
  ar = (b_y_im * xc_y_re - k_y_im * xc_y_im) * x07;
  ai = (b_y_im * xc_y_im + k_y_im * xc_y_re) * x07;
  y_re = 2.0 * (500.0 * yc_y_re + 9.0 * ad_y_re);
  r = 2.0 * (500.0 * yc_y_im + 9.0 * ad_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      b_y_im = ar / y_re;
    } else if (ar == 0.0) {
      b_y_im = 0.0;
    } else {
      b_y_im = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      b_y_im = ai / r;
    } else if (ai == 0.0) {
      b_y_im = 0.0;
    } else {
      b_y_im = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      b_y_im = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        wl_y_re = 0.5;
      } else {
        wl_y_re = -0.5;
      }

      if (r > 0.0) {
        bb_r = 0.5;
      } else {
        bb_r = -0.5;
      }

      b_y_im = (ar * wl_y_re + ai * bb_r) / brm;
    } else {
      brm = y_re / r;
      b_y_im = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  b_y_re = 3.0 * bd_y_re;
  k_y_im = 3.0 * bd_y_im;
  ar = (b_y_re * t_s_re_tmp - k_y_im * t_s_im_tmp) * x07;
  ai = (b_y_re * t_s_im_tmp + k_y_im * t_s_re_tmp) * x07;
  y_re = 4.0 * (500.0 * cd_y_re + 9.0 * dd_y_re);
  r = 4.0 * (500.0 * cd_y_im + 9.0 * dd_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      b_y_re = ar / y_re;
    } else if (ar == 0.0) {
      b_y_re = 0.0;
    } else {
      b_y_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      b_y_re = ai / r;
    } else if (ai == 0.0) {
      b_y_re = 0.0;
    } else {
      b_y_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      b_y_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        xl_y_re = 0.5;
      } else {
        xl_y_re = -0.5;
      }

      if (r > 0.0) {
        cb_r = 0.5;
      } else {
        cb_r = -0.5;
      }

      b_y_re = (ar * xl_y_re + ai * cb_r) / brm;
    } else {
      brm = y_re / r;
      b_y_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  r = 125.0 * ed_y_re;
  k_y_im = 125.0 * ed_y_im;
  ar = (r * fd_y_re - k_y_im * fd_y_im) * x07;
  ai = (r * fd_y_im + k_y_im * fd_y_re) * x07;
  y_re = 2.0 * (500.0 * gd_y_re + 9.0 * hd_y_re);
  r = 2.0 * (500.0 * gd_y_im + 9.0 * hd_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      r = ar / y_re;
    } else if (ar == 0.0) {
      r = 0.0;
    } else {
      r = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      r = ai / r;
    } else if (ai == 0.0) {
      r = 0.0;
    } else {
      r = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      r = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        yl_y_re = 0.5;
      } else {
        yl_y_re = -0.5;
      }

      if (r > 0.0) {
        db_r = 0.5;
      } else {
        db_r = -0.5;
      }

      r = (ar * yl_y_re + ai * db_r) / brm;
    } else {
      brm = y_re / r;
      r = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  states[0] = (((((((((((((((((((((((((((((x01 + t.re * x04) + t_re) + re) -
    b_re) + c_re) - d_re) - m_y_im) + m_y_re) - l_y_im) + l_y_re) - k_y_re) +
    j_y_im) - j_y_re) + i_y_im) + i_y_re) - h_y_im) - h_y_re) + g_y_im) - g_y_re)
                        + f_y_im) - f_y_re) + e_y_im) - e_y_re) + d_y_im) -
                   d_y_re) + c_y_im) - c_y_re) + b_y_im) - b_y_re) - r;
  ar = t_re_tmp * x08;
  if (t_im_tmp * x08 == 0.0) {
    t_re = ar / 2.0;
  } else if (ar == 0.0) {
    t_re = 0.0;
  } else {
    t_re = ar / 2.0;
  }

  ar = 625.0 * id_y_re * x02;
  ai = 625.0 * id_y_im * x02;
  y_re = 500.0 * jd_y_re + 9.0 * t_s_re_tmp;
  r = 500.0 * jd_y_im + 9.0 * t_s_im_tmp;
  if (r == 0.0) {
    if (ai == 0.0) {
      re = ar / y_re;
    } else if (ar == 0.0) {
      re = 0.0;
    } else {
      re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      re = ai / r;
    } else if (ai == 0.0) {
      re = 0.0;
    } else {
      re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        am_y_re = 0.5;
      } else {
        am_y_re = -0.5;
      }

      if (r > 0.0) {
        eb_r = 0.5;
      } else {
        eb_r = -0.5;
      }

      re = (ar * am_y_re + ai * eb_r) / brm;
    } else {
      brm = y_re / r;
      re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  ar = 15.0 * kd_y_re * x02;
  ai = 15.0 * kd_y_im * x02;
  y_re = 500.0 * ld_y_re + 9.0 * md_y_re;
  r = 500.0 * ld_y_im + 9.0 * md_y_im;
  if (r == 0.0) {
    if (ai == 0.0) {
      b_re = ar / y_re;
    } else if (ar == 0.0) {
      b_re = 0.0;
    } else {
      b_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      b_re = ai / r;
    } else if (ai == 0.0) {
      b_re = 0.0;
    } else {
      b_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      b_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        bm_y_re = 0.5;
      } else {
        bm_y_re = -0.5;
      }

      if (r > 0.0) {
        fb_r = 0.5;
      } else {
        fb_r = -0.5;
      }

      b_re = (ar * bm_y_re + ai * fb_r) / brm;
    } else {
      brm = y_re / r;
      b_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  ar = 15.0 * nd_y_re * x02;
  ai = 15.0 * nd_y_im * x02;
  y_re = 2.0 * (500.0 * od_y_re + 9.0 * pd_y_re);
  r = 2.0 * (500.0 * od_y_im + 9.0 * pd_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      c_re = ar / y_re;
    } else if (ar == 0.0) {
      c_re = 0.0;
    } else {
      c_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      c_re = ai / r;
    } else if (ai == 0.0) {
      c_re = 0.0;
    } else {
      c_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      c_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        cm_y_re = 0.5;
      } else {
        cm_y_re = -0.5;
      }

      if (r > 0.0) {
        gb_r = 0.5;
      } else {
        gb_r = -0.5;
      }

      c_re = (ar * cm_y_re + ai * gb_r) / brm;
    } else {
      brm = y_re / r;
      c_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  ar = 3.0 * qd_y_re * x02;
  ai = 3.0 * qd_y_im * x02;
  y_re = 2.0 * (500.0 * rd_y_re + 9.0 * sd_y_re);
  r = 2.0 * (500.0 * rd_y_im + 9.0 * sd_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      d_re = ar / y_re;
    } else if (ar == 0.0) {
      d_re = 0.0;
    } else {
      d_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      d_re = ai / r;
    } else if (ai == 0.0) {
      d_re = 0.0;
    } else {
      d_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      d_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        dm_y_re = 0.5;
      } else {
        dm_y_re = -0.5;
      }

      if (r > 0.0) {
        hb_r = 0.5;
      } else {
        hb_r = -0.5;
      }

      d_re = (ar * dm_y_re + ai * hb_r) / brm;
    } else {
      brm = y_re / r;
      d_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  ar = 625.0 * td_y_re * x12;
  ai = 625.0 * td_y_im * x12;
  y_re = 500.0 * ud_y_re + 9.0 * t_s_re_tmp;
  r = 500.0 * ud_y_im + 9.0 * t_s_im_tmp;
  if (r == 0.0) {
    if (ai == 0.0) {
      m_y_im = ar / y_re;
    } else if (ar == 0.0) {
      m_y_im = 0.0;
    } else {
      m_y_im = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      m_y_im = ai / r;
    } else if (ai == 0.0) {
      m_y_im = 0.0;
    } else {
      m_y_im = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      m_y_im = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        em_y_re = 0.5;
      } else {
        em_y_re = -0.5;
      }

      if (r > 0.0) {
        ib_r = 0.5;
      } else {
        ib_r = -0.5;
      }

      m_y_im = (ar * em_y_re + ai * ib_r) / brm;
    } else {
      brm = y_re / r;
      m_y_im = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  ar = 15.0 * vd_y_re * x12;
  ai = 15.0 * vd_y_im * x12;
  y_re = 500.0 * wd_y_re + 9.0 * xd_y_re;
  r = 500.0 * wd_y_im + 9.0 * xd_y_im;
  if (r == 0.0) {
    if (ai == 0.0) {
      m_y_re = ar / y_re;
    } else if (ar == 0.0) {
      m_y_re = 0.0;
    } else {
      m_y_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      m_y_re = ai / r;
    } else if (ai == 0.0) {
      m_y_re = 0.0;
    } else {
      m_y_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      m_y_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        fm_y_re = 0.5;
      } else {
        fm_y_re = -0.5;
      }

      if (r > 0.0) {
        jb_r = 0.5;
      } else {
        jb_r = -0.5;
      }

      m_y_re = (ar * fm_y_re + ai * jb_r) / brm;
    } else {
      brm = y_re / r;
      m_y_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  ar = 15.0 * yd_y_re * x12;
  ai = 15.0 * yd_y_im * x12;
  y_re = 2.0 * (500.0 * ae_y_re + 9.0 * be_y_re);
  r = 2.0 * (500.0 * ae_y_im + 9.0 * be_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      l_y_im = ar / y_re;
    } else if (ar == 0.0) {
      l_y_im = 0.0;
    } else {
      l_y_im = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      l_y_im = ai / r;
    } else if (ai == 0.0) {
      l_y_im = 0.0;
    } else {
      l_y_im = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      l_y_im = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        gm_y_re = 0.5;
      } else {
        gm_y_re = -0.5;
      }

      if (r > 0.0) {
        kb_r = 0.5;
      } else {
        kb_r = -0.5;
      }

      l_y_im = (ar * gm_y_re + ai * kb_r) / brm;
    } else {
      brm = y_re / r;
      l_y_im = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  ar = 3.0 * ce_y_re * x12;
  ai = 3.0 * ce_y_im * x12;
  y_re = 2.0 * (500.0 * de_y_re + 9.0 * ee_y_re);
  r = 2.0 * (500.0 * de_y_im + 9.0 * ee_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      l_y_re = ar / y_re;
    } else if (ar == 0.0) {
      l_y_re = 0.0;
    } else {
      l_y_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      l_y_re = ai / r;
    } else if (ai == 0.0) {
      l_y_re = 0.0;
    } else {
      l_y_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      l_y_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        hm_y_re = 0.5;
      } else {
        hm_y_re = -0.5;
      }

      if (r > 0.0) {
        lb_r = 0.5;
      } else {
        lb_r = -0.5;
      }

      l_y_re = (ar * hm_y_re + ai * lb_r) / brm;
    } else {
      brm = y_re / r;
      l_y_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  k_y_re = 1875.0 * fe_y_re;
  k_y_im = 1875.0 * fe_y_im;
  ar = (k_y_re * t_s.re - k_y_im * t_s.im) * x02;
  ai = (k_y_re * t_s.im + k_y_im * t_s.re) * x02;
  y_re = 500.0 * ge_y_re + 9.0 * he_y_re;
  r = 500.0 * ge_y_im + 9.0 * he_y_im;
  if (r == 0.0) {
    if (ai == 0.0) {
      k_y_re = ar / y_re;
    } else if (ar == 0.0) {
      k_y_re = 0.0;
    } else {
      k_y_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      k_y_re = ai / r;
    } else if (ai == 0.0) {
      k_y_re = 0.0;
    } else {
      k_y_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      k_y_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        im_y_re = 0.5;
      } else {
        im_y_re = -0.5;
      }

      if (r > 0.0) {
        mb_r = 0.5;
      } else {
        mb_r = -0.5;
      }

      k_y_re = (ar * im_y_re + ai * mb_r) / brm;
    } else {
      brm = y_re / r;
      k_y_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  j_y_im = 625.0 * ie_y_re;
  k_y_im = 625.0 * ie_y_im;
  ar = (j_y_im * t_s.re - k_y_im * t_s.im) * x05;
  ai = (j_y_im * t_s.im + k_y_im * t_s.re) * x05;
  y_re = 500.0 * je_y_re + 9.0 * t_s_re_tmp;
  r = 500.0 * je_y_im + 9.0 * t_s_im_tmp;
  if (r == 0.0) {
    if (ai == 0.0) {
      j_y_im = ar / y_re;
    } else if (ar == 0.0) {
      j_y_im = 0.0;
    } else {
      j_y_im = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      j_y_im = ai / r;
    } else if (ai == 0.0) {
      j_y_im = 0.0;
    } else {
      j_y_im = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      j_y_im = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        jm_y_re = 0.5;
      } else {
        jm_y_re = -0.5;
      }

      if (r > 0.0) {
        nb_r = 0.5;
      } else {
        nb_r = -0.5;
      }

      j_y_im = (ar * jm_y_re + ai * nb_r) / brm;
    } else {
      brm = y_re / r;
      j_y_im = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  j_y_re = 15.0 * ke_y_re;
  k_y_im = 15.0 * ke_y_im;
  ar = (j_y_re * t_s.re - k_y_im * t_s.im) * x05;
  ai = (j_y_re * t_s.im + k_y_im * t_s.re) * x05;
  y_re = 500.0 * le_y_re + 9.0 * me_y_re;
  r = 500.0 * le_y_im + 9.0 * me_y_im;
  if (r == 0.0) {
    if (ai == 0.0) {
      j_y_re = ar / y_re;
    } else if (ar == 0.0) {
      j_y_re = 0.0;
    } else {
      j_y_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      j_y_re = ai / r;
    } else if (ai == 0.0) {
      j_y_re = 0.0;
    } else {
      j_y_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      j_y_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        km_y_re = 0.5;
      } else {
        km_y_re = -0.5;
      }

      if (r > 0.0) {
        ob_r = 0.5;
      } else {
        ob_r = -0.5;
      }

      j_y_re = (ar * km_y_re + ai * ob_r) / brm;
    } else {
      brm = y_re / r;
      j_y_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  i_y_im = 1875.0 * ne_y_re;
  k_y_im = 1875.0 * ne_y_im;
  ar = (i_y_im * t_s.re - k_y_im * t_s.im) * x02;
  ai = (i_y_im * t_s.im + k_y_im * t_s.re) * x02;
  y_re = 2.0 * (500.0 * oe_y_re + 9.0 * pe_y_re);
  r = 2.0 * (500.0 * oe_y_im + 9.0 * pe_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      i_y_im = ar / y_re;
    } else if (ar == 0.0) {
      i_y_im = 0.0;
    } else {
      i_y_im = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      i_y_im = ai / r;
    } else if (ai == 0.0) {
      i_y_im = 0.0;
    } else {
      i_y_im = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      i_y_im = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        lm_y_re = 0.5;
      } else {
        lm_y_re = -0.5;
      }

      if (r > 0.0) {
        pb_r = 0.5;
      } else {
        pb_r = -0.5;
      }

      i_y_im = (ar * lm_y_re + ai * pb_r) / brm;
    } else {
      brm = y_re / r;
      i_y_im = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  i_y_re = 15.0 * qe_y_re;
  k_y_im = 15.0 * qe_y_im;
  ar = (i_y_re * t_s.re - k_y_im * t_s.im) * x05;
  ai = (i_y_re * t_s.im + k_y_im * t_s.re) * x05;
  y_re = 2.0 * (500.0 * re_y_re + 9.0 * se_y_re);
  r = 2.0 * (500.0 * re_y_im + 9.0 * se_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      i_y_re = ar / y_re;
    } else if (ar == 0.0) {
      i_y_re = 0.0;
    } else {
      i_y_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      i_y_re = ai / r;
    } else if (ai == 0.0) {
      i_y_re = 0.0;
    } else {
      i_y_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      i_y_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        mm_y_re = 0.5;
      } else {
        mm_y_re = -0.5;
      }

      if (r > 0.0) {
        qb_r = 0.5;
      } else {
        qb_r = -0.5;
      }

      i_y_re = (ar * mm_y_re + ai * qb_r) / brm;
    } else {
      brm = y_re / r;
      i_y_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  h_y_im = 375.0 * te_y_re;
  k_y_im = 375.0 * te_y_im;
  ar = (h_y_im * t_s.re - k_y_im * t_s.im) * x02;
  ai = (h_y_im * t_s.im + k_y_im * t_s.re) * x02;
  y_re = 2.0 * (500.0 * ue_y_re + 9.0 * ve_y_re);
  r = 2.0 * (500.0 * ue_y_im + 9.0 * ve_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      h_y_im = ar / y_re;
    } else if (ar == 0.0) {
      h_y_im = 0.0;
    } else {
      h_y_im = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      h_y_im = ai / r;
    } else if (ai == 0.0) {
      h_y_im = 0.0;
    } else {
      h_y_im = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      h_y_im = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        nm_y_re = 0.5;
      } else {
        nm_y_re = -0.5;
      }

      if (r > 0.0) {
        rb_r = 0.5;
      } else {
        rb_r = -0.5;
      }

      h_y_im = (ar * nm_y_re + ai * rb_r) / brm;
    } else {
      brm = y_re / r;
      h_y_im = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  h_y_re = 3.0 * we_y_re;
  k_y_im = 3.0 * we_y_im;
  ar = (h_y_re * t_s.re - k_y_im * t_s.im) * x05;
  ai = (h_y_re * t_s.im + k_y_im * t_s.re) * x05;
  y_re = 2.0 * (500.0 * xe_y_re + 9.0 * ye_y_re);
  r = 2.0 * (500.0 * xe_y_im + 9.0 * ye_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      h_y_re = ar / y_re;
    } else if (ar == 0.0) {
      h_y_re = 0.0;
    } else {
      h_y_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      h_y_re = ai / r;
    } else if (ai == 0.0) {
      h_y_re = 0.0;
    } else {
      h_y_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      h_y_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        om_y_re = 0.5;
      } else {
        om_y_re = -0.5;
      }

      if (r > 0.0) {
        sb_r = 0.5;
      } else {
        sb_r = -0.5;
      }

      h_y_re = (ar * om_y_re + ai * sb_r) / brm;
    } else {
      brm = y_re / r;
      h_y_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  g_y_im = 1875.0 * af_y_re;
  k_y_im = 1875.0 * af_y_im;
  ar = (g_y_im * t_s.re - k_y_im * t_s.im) * x12;
  ai = (g_y_im * t_s.im + k_y_im * t_s.re) * x12;
  y_re = 500.0 * bf_y_re + 9.0 * cf_y_re;
  r = 500.0 * bf_y_im + 9.0 * cf_y_im;
  if (r == 0.0) {
    if (ai == 0.0) {
      g_y_im = ar / y_re;
    } else if (ar == 0.0) {
      g_y_im = 0.0;
    } else {
      g_y_im = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      g_y_im = ai / r;
    } else if (ai == 0.0) {
      g_y_im = 0.0;
    } else {
      g_y_im = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      g_y_im = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        pm_y_re = 0.5;
      } else {
        pm_y_re = -0.5;
      }

      if (r > 0.0) {
        tb_r = 0.5;
      } else {
        tb_r = -0.5;
      }

      g_y_im = (ar * pm_y_re + ai * tb_r) / brm;
    } else {
      brm = y_re / r;
      g_y_im = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  g_y_re = 1875.0 * df_y_re;
  k_y_im = 1875.0 * df_y_im;
  ar = (g_y_re * t_s.re - k_y_im * t_s.im) * x12;
  ai = (g_y_re * t_s.im + k_y_im * t_s.re) * x12;
  y_re = 2.0 * (500.0 * ef_y_re + 9.0 * ff_y_re);
  r = 2.0 * (500.0 * ef_y_im + 9.0 * ff_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      g_y_re = ar / y_re;
    } else if (ar == 0.0) {
      g_y_re = 0.0;
    } else {
      g_y_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      g_y_re = ai / r;
    } else if (ai == 0.0) {
      g_y_re = 0.0;
    } else {
      g_y_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      g_y_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        qm_y_re = 0.5;
      } else {
        qm_y_re = -0.5;
      }

      if (r > 0.0) {
        ub_r = 0.5;
      } else {
        ub_r = -0.5;
      }

      g_y_re = (ar * qm_y_re + ai * ub_r) / brm;
    } else {
      brm = y_re / r;
      g_y_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  f_y_im = 375.0 * gf_y_re;
  k_y_im = 375.0 * gf_y_im;
  ar = (f_y_im * t_s.re - k_y_im * t_s.im) * x12;
  ai = (f_y_im * t_s.im + k_y_im * t_s.re) * x12;
  y_re = 2.0 * (500.0 * hf_y_re + 9.0 * if_y_re);
  r = 2.0 * (500.0 * hf_y_im + 9.0 * if_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      f_y_im = ar / y_re;
    } else if (ar == 0.0) {
      f_y_im = 0.0;
    } else {
      f_y_im = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      f_y_im = ai / r;
    } else if (ai == 0.0) {
      f_y_im = 0.0;
    } else {
      f_y_im = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      f_y_im = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        rm_y_re = 0.5;
      } else {
        rm_y_re = -0.5;
      }

      if (r > 0.0) {
        vb_r = 0.5;
      } else {
        vb_r = -0.5;
      }

      f_y_im = (ar * rm_y_re + ai * vb_r) / brm;
    } else {
      brm = y_re / r;
      f_y_im = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  f_y_re = 1875.0 * jf_y_re;
  k_y_im = 1875.0 * jf_y_im;
  ar = (f_y_re * t_s_re_tmp - k_y_im * t_s_im_tmp) * x05;
  ai = (f_y_re * t_s_im_tmp + k_y_im * t_s_re_tmp) * x05;
  y_re = 500.0 * kf_y_re + 9.0 * lf_y_re;
  r = 500.0 * kf_y_im + 9.0 * lf_y_im;
  if (r == 0.0) {
    if (ai == 0.0) {
      f_y_re = ar / y_re;
    } else if (ar == 0.0) {
      f_y_re = 0.0;
    } else {
      f_y_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      f_y_re = ai / r;
    } else if (ai == 0.0) {
      f_y_re = 0.0;
    } else {
      f_y_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      f_y_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        sm_y_re = 0.5;
      } else {
        sm_y_re = -0.5;
      }

      if (r > 0.0) {
        wb_r = 0.5;
      } else {
        wb_r = -0.5;
      }

      f_y_re = (ar * sm_y_re + ai * wb_r) / brm;
    } else {
      brm = y_re / r;
      f_y_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  e_y_im = 125.0 * mf_y_re;
  k_y_im = 125.0 * mf_y_im;
  ar = (e_y_im * t_s_re_tmp - k_y_im * t_s_im_tmp) * x08;
  ai = (e_y_im * t_s_im_tmp + k_y_im * t_s_re_tmp) * x08;
  y_re = 500.0 * nf_y_re + 9.0 * t_s_re_tmp;
  r = 500.0 * nf_y_im + 9.0 * t_s_im_tmp;
  if (r == 0.0) {
    if (ai == 0.0) {
      e_y_im = ar / y_re;
    } else if (ar == 0.0) {
      e_y_im = 0.0;
    } else {
      e_y_im = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      e_y_im = ai / r;
    } else if (ai == 0.0) {
      e_y_im = 0.0;
    } else {
      e_y_im = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      e_y_im = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        tm_y_re = 0.5;
      } else {
        tm_y_re = -0.5;
      }

      if (r > 0.0) {
        xb_r = 0.5;
      } else {
        xb_r = -0.5;
      }

      e_y_im = (ar * tm_y_re + ai * xb_r) / brm;
    } else {
      brm = y_re / r;
      e_y_im = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  e_y_re = 15.0 * of_y_re;
  k_y_im = 15.0 * of_y_im;
  ar = (e_y_re * t_s_re_tmp - k_y_im * t_s_im_tmp) * x08;
  ai = (e_y_re * t_s_im_tmp + k_y_im * t_s_re_tmp) * x08;
  y_re = 2.0 * (500.0 * pf_y_re + 9.0 * qf_y_re);
  r = 2.0 * (500.0 * pf_y_im + 9.0 * qf_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      e_y_re = ar / y_re;
    } else if (ar == 0.0) {
      e_y_re = 0.0;
    } else {
      e_y_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      e_y_re = ai / r;
    } else if (ai == 0.0) {
      e_y_re = 0.0;
    } else {
      e_y_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      e_y_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        um_y_re = 0.5;
      } else {
        um_y_re = -0.5;
      }

      if (r > 0.0) {
        yb_r = 0.5;
      } else {
        yb_r = -0.5;
      }

      e_y_re = (ar * um_y_re + ai * yb_r) / brm;
    } else {
      brm = y_re / r;
      e_y_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  d_y_im = 1875.0 * rf_y_re;
  k_y_im = 1875.0 * rf_y_im;
  ar = (d_y_im * t_s_re_tmp - k_y_im * t_s_im_tmp) * x05;
  ai = (d_y_im * t_s_im_tmp + k_y_im * t_s_re_tmp) * x05;
  y_re = 2.0 * (500.0 * sf_y_re + 9.0 * tf_y_re);
  r = 2.0 * (500.0 * sf_y_im + 9.0 * tf_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      d_y_im = ar / y_re;
    } else if (ar == 0.0) {
      d_y_im = 0.0;
    } else {
      d_y_im = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      d_y_im = ai / r;
    } else if (ai == 0.0) {
      d_y_im = 0.0;
    } else {
      d_y_im = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      d_y_im = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        vm_y_re = 0.5;
      } else {
        vm_y_re = -0.5;
      }

      if (r > 0.0) {
        ac_r = 0.5;
      } else {
        ac_r = -0.5;
      }

      d_y_im = (ar * vm_y_re + ai * ac_r) / brm;
    } else {
      brm = y_re / r;
      d_y_im = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  d_y_re = 625.0 * uf_y_re;
  k_y_im = 625.0 * uf_y_im;
  ar = (d_y_re * vf_y_re - k_y_im * vf_y_im) * x08;
  ai = (d_y_re * vf_y_im + k_y_im * vf_y_re) * x08;
  y_re = 500.0 * wf_y_re + 9.0 * xf_y_re;
  r = 500.0 * wf_y_im + 9.0 * xf_y_im;
  if (r == 0.0) {
    if (ai == 0.0) {
      d_y_re = ar / y_re;
    } else if (ar == 0.0) {
      d_y_re = 0.0;
    } else {
      d_y_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      d_y_re = ai / r;
    } else if (ai == 0.0) {
      d_y_re = 0.0;
    } else {
      d_y_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      d_y_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        wm_y_re = 0.5;
      } else {
        wm_y_re = -0.5;
      }

      if (r > 0.0) {
        bc_r = 0.5;
      } else {
        bc_r = -0.5;
      }

      d_y_re = (ar * wm_y_re + ai * bc_r) / brm;
    } else {
      brm = y_re / r;
      d_y_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  c_y_im = 15.0 * yf_y_re;
  k_y_im = 15.0 * yf_y_im;
  ar = (c_y_im * t_s_re_tmp - k_y_im * t_s_im_tmp) * x08;
  ai = (c_y_im * t_s_im_tmp + k_y_im * t_s_re_tmp) * x08;
  y_re = 4.0 * (500.0 * ag_y_re + 9.0 * bg_y_re);
  r = 4.0 * (500.0 * ag_y_im + 9.0 * bg_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      c_y_im = ar / y_re;
    } else if (ar == 0.0) {
      c_y_im = 0.0;
    } else {
      c_y_im = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      c_y_im = ai / r;
    } else if (ai == 0.0) {
      c_y_im = 0.0;
    } else {
      c_y_im = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      c_y_im = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        xm_y_re = 0.5;
      } else {
        xm_y_re = -0.5;
      }

      if (r > 0.0) {
        cc_r = 0.5;
      } else {
        cc_r = -0.5;
      }

      c_y_im = (ar * xm_y_re + ai * cc_r) / brm;
    } else {
      brm = y_re / r;
      c_y_im = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  c_y_re = 375.0 * cg_y_re;
  k_y_im = 375.0 * cg_y_im;
  ar = (c_y_re * t_s_re_tmp - k_y_im * t_s_im_tmp) * x05;
  ai = (c_y_re * t_s_im_tmp + k_y_im * t_s_re_tmp) * x05;
  y_re = 2.0 * (500.0 * dg_y_re + 9.0 * eg_y_re);
  r = 2.0 * (500.0 * dg_y_im + 9.0 * eg_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      c_y_re = ar / y_re;
    } else if (ar == 0.0) {
      c_y_re = 0.0;
    } else {
      c_y_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      c_y_re = ai / r;
    } else if (ai == 0.0) {
      c_y_re = 0.0;
    } else {
      c_y_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      c_y_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        ym_y_re = 0.5;
      } else {
        ym_y_re = -0.5;
      }

      if (r > 0.0) {
        dc_r = 0.5;
      } else {
        dc_r = -0.5;
      }

      c_y_re = (ar * ym_y_re + ai * dc_r) / brm;
    } else {
      brm = y_re / r;
      c_y_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  b_y_im = 625.0 * fg_y_re;
  k_y_im = 625.0 * fg_y_im;
  ar = (b_y_im * gg_y_re - k_y_im * gg_y_im) * x08;
  ai = (b_y_im * gg_y_im + k_y_im * gg_y_re) * x08;
  y_re = 2.0 * (500.0 * hg_y_re + 9.0 * ig_y_re);
  r = 2.0 * (500.0 * hg_y_im + 9.0 * ig_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      b_y_im = ar / y_re;
    } else if (ar == 0.0) {
      b_y_im = 0.0;
    } else {
      b_y_im = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      b_y_im = ai / r;
    } else if (ai == 0.0) {
      b_y_im = 0.0;
    } else {
      b_y_im = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      b_y_im = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        an_y_re = 0.5;
      } else {
        an_y_re = -0.5;
      }

      if (r > 0.0) {
        ec_r = 0.5;
      } else {
        ec_r = -0.5;
      }

      b_y_im = (ar * an_y_re + ai * ec_r) / brm;
    } else {
      brm = y_re / r;
      b_y_im = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  b_y_re = 3.0 * jg_y_re;
  k_y_im = 3.0 * jg_y_im;
  ar = (b_y_re * t_s_re_tmp - k_y_im * t_s_im_tmp) * x08;
  ai = (b_y_re * t_s_im_tmp + k_y_im * t_s_re_tmp) * x08;
  y_re = 4.0 * (500.0 * kg_y_re + 9.0 * lg_y_re);
  r = 4.0 * (500.0 * kg_y_im + 9.0 * lg_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      b_y_re = ar / y_re;
    } else if (ar == 0.0) {
      b_y_re = 0.0;
    } else {
      b_y_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      b_y_re = ai / r;
    } else if (ai == 0.0) {
      b_y_re = 0.0;
    } else {
      b_y_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      b_y_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        bn_y_re = 0.5;
      } else {
        bn_y_re = -0.5;
      }

      if (r > 0.0) {
        fc_r = 0.5;
      } else {
        fc_r = -0.5;
      }

      b_y_re = (ar * bn_y_re + ai * fc_r) / brm;
    } else {
      brm = y_re / r;
      b_y_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  r = 125.0 * mg_y_re;
  k_y_im = 125.0 * mg_y_im;
  ar = (r * ng_y_re - k_y_im * ng_y_im) * x08;
  ai = (r * ng_y_im + k_y_im * ng_y_re) * x08;
  y_re = 2.0 * (500.0 * og_y_re + 9.0 * pg_y_re);
  r = 2.0 * (500.0 * og_y_im + 9.0 * pg_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      r = ar / y_re;
    } else if (ar == 0.0) {
      r = 0.0;
    } else {
      r = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      r = ai / r;
    } else if (ai == 0.0) {
      r = 0.0;
    } else {
      r = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      r = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        cn_y_re = 0.5;
      } else {
        cn_y_re = -0.5;
      }

      if (r > 0.0) {
        gc_r = 0.5;
      } else {
        gc_r = -0.5;
      }

      r = (ar * cn_y_re + ai * gc_r) / brm;
    } else {
      brm = y_re / r;
      r = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  states[1] = (((((((((((((((((((((((((((((x02 + t.re * x05) + t_re) + re) -
    b_re) + c_re) - d_re) - m_y_im) + m_y_re) - l_y_im) + l_y_re) - k_y_re) +
    j_y_im) - j_y_re) + i_y_im) + i_y_re) - h_y_im) - h_y_re) + g_y_im) - g_y_re)
                        + f_y_im) - f_y_re) + e_y_im) - e_y_re) + d_y_im) -
                   d_y_re) + c_y_im) - c_y_re) + b_y_im) - b_y_re) - r;
  ar = t_re_tmp * x09;
  if (t_im_tmp * x09 == 0.0) {
    t_re = ar / 2.0;
  } else if (ar == 0.0) {
    t_re = 0.0;
  } else {
    t_re = ar / 2.0;
  }

  ar = 625.0 * qg_y_re * x03;
  ai = 625.0 * qg_y_im * x03;
  y_re = 500.0 * rg_y_re + 9.0 * t_s_re_tmp;
  r = 500.0 * rg_y_im + 9.0 * t_s_im_tmp;
  if (r == 0.0) {
    if (ai == 0.0) {
      re = ar / y_re;
    } else if (ar == 0.0) {
      re = 0.0;
    } else {
      re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      re = ai / r;
    } else if (ai == 0.0) {
      re = 0.0;
    } else {
      re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        dn_y_re = 0.5;
      } else {
        dn_y_re = -0.5;
      }

      if (r > 0.0) {
        hc_r = 0.5;
      } else {
        hc_r = -0.5;
      }

      re = (ar * dn_y_re + ai * hc_r) / brm;
    } else {
      brm = y_re / r;
      re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  ar = 15.0 * sg_y_re * x03;
  ai = 15.0 * sg_y_im * x03;
  y_re = 500.0 * tg_y_re + 9.0 * ug_y_re;
  r = 500.0 * tg_y_im + 9.0 * ug_y_im;
  if (r == 0.0) {
    if (ai == 0.0) {
      b_re = ar / y_re;
    } else if (ar == 0.0) {
      b_re = 0.0;
    } else {
      b_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      b_re = ai / r;
    } else if (ai == 0.0) {
      b_re = 0.0;
    } else {
      b_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      b_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        en_y_re = 0.5;
      } else {
        en_y_re = -0.5;
      }

      if (r > 0.0) {
        ic_r = 0.5;
      } else {
        ic_r = -0.5;
      }

      b_re = (ar * en_y_re + ai * ic_r) / brm;
    } else {
      brm = y_re / r;
      b_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  ar = 15.0 * vg_y_re * x03;
  ai = 15.0 * vg_y_im * x03;
  y_re = 2.0 * (500.0 * wg_y_re + 9.0 * xg_y_re);
  r = 2.0 * (500.0 * wg_y_im + 9.0 * xg_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      c_re = ar / y_re;
    } else if (ar == 0.0) {
      c_re = 0.0;
    } else {
      c_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      c_re = ai / r;
    } else if (ai == 0.0) {
      c_re = 0.0;
    } else {
      c_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      c_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        fn_y_re = 0.5;
      } else {
        fn_y_re = -0.5;
      }

      if (r > 0.0) {
        jc_r = 0.5;
      } else {
        jc_r = -0.5;
      }

      c_re = (ar * fn_y_re + ai * jc_r) / brm;
    } else {
      brm = y_re / r;
      c_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  ar = 3.0 * yg_y_re * x03;
  ai = 3.0 * yg_y_im * x03;
  y_re = 2.0 * (500.0 * ah_y_re + 9.0 * bh_y_re);
  r = 2.0 * (500.0 * ah_y_im + 9.0 * bh_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      d_re = ar / y_re;
    } else if (ar == 0.0) {
      d_re = 0.0;
    } else {
      d_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      d_re = ai / r;
    } else if (ai == 0.0) {
      d_re = 0.0;
    } else {
      d_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      d_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        gn_y_re = 0.5;
      } else {
        gn_y_re = -0.5;
      }

      if (r > 0.0) {
        kc_r = 0.5;
      } else {
        kc_r = -0.5;
      }

      d_re = (ar * gn_y_re + ai * kc_r) / brm;
    } else {
      brm = y_re / r;
      d_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  ar = 625.0 * ch_y_re * x13;
  ai = 625.0 * ch_y_im * x13;
  y_re = 500.0 * dh_y_re + 9.0 * t_s_re_tmp;
  r = 500.0 * dh_y_im + 9.0 * t_s_im_tmp;
  if (r == 0.0) {
    if (ai == 0.0) {
      m_y_im = ar / y_re;
    } else if (ar == 0.0) {
      m_y_im = 0.0;
    } else {
      m_y_im = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      m_y_im = ai / r;
    } else if (ai == 0.0) {
      m_y_im = 0.0;
    } else {
      m_y_im = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      m_y_im = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        hn_y_re = 0.5;
      } else {
        hn_y_re = -0.5;
      }

      if (r > 0.0) {
        lc_r = 0.5;
      } else {
        lc_r = -0.5;
      }

      m_y_im = (ar * hn_y_re + ai * lc_r) / brm;
    } else {
      brm = y_re / r;
      m_y_im = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  ar = 15.0 * eh_y_re * x13;
  ai = 15.0 * eh_y_im * x13;
  y_re = 500.0 * fh_y_re + 9.0 * gh_y_re;
  r = 500.0 * fh_y_im + 9.0 * gh_y_im;
  if (r == 0.0) {
    if (ai == 0.0) {
      m_y_re = ar / y_re;
    } else if (ar == 0.0) {
      m_y_re = 0.0;
    } else {
      m_y_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      m_y_re = ai / r;
    } else if (ai == 0.0) {
      m_y_re = 0.0;
    } else {
      m_y_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      m_y_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        in_y_re = 0.5;
      } else {
        in_y_re = -0.5;
      }

      if (r > 0.0) {
        mc_r = 0.5;
      } else {
        mc_r = -0.5;
      }

      m_y_re = (ar * in_y_re + ai * mc_r) / brm;
    } else {
      brm = y_re / r;
      m_y_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  ar = 15.0 * hh_y_re * x13;
  ai = 15.0 * hh_y_im * x13;
  y_re = 2.0 * (500.0 * ih_y_re + 9.0 * jh_y_re);
  r = 2.0 * (500.0 * ih_y_im + 9.0 * jh_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      l_y_im = ar / y_re;
    } else if (ar == 0.0) {
      l_y_im = 0.0;
    } else {
      l_y_im = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      l_y_im = ai / r;
    } else if (ai == 0.0) {
      l_y_im = 0.0;
    } else {
      l_y_im = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      l_y_im = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        jn_y_re = 0.5;
      } else {
        jn_y_re = -0.5;
      }

      if (r > 0.0) {
        nc_r = 0.5;
      } else {
        nc_r = -0.5;
      }

      l_y_im = (ar * jn_y_re + ai * nc_r) / brm;
    } else {
      brm = y_re / r;
      l_y_im = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  ar = 3.0 * kh_y_re * x13;
  ai = 3.0 * kh_y_im * x13;
  y_re = 2.0 * (500.0 * lh_y_re + 9.0 * mh_y_re);
  r = 2.0 * (500.0 * lh_y_im + 9.0 * mh_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      l_y_re = ar / y_re;
    } else if (ar == 0.0) {
      l_y_re = 0.0;
    } else {
      l_y_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      l_y_re = ai / r;
    } else if (ai == 0.0) {
      l_y_re = 0.0;
    } else {
      l_y_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      l_y_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        kn_y_re = 0.5;
      } else {
        kn_y_re = -0.5;
      }

      if (r > 0.0) {
        oc_r = 0.5;
      } else {
        oc_r = -0.5;
      }

      l_y_re = (ar * kn_y_re + ai * oc_r) / brm;
    } else {
      brm = y_re / r;
      l_y_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  k_y_re = 1875.0 * nh_y_re;
  k_y_im = 1875.0 * nh_y_im;
  ar = (k_y_re * t_s.re - k_y_im * t_s.im) * x03;
  ai = (k_y_re * t_s.im + k_y_im * t_s.re) * x03;
  y_re = 500.0 * oh_y_re + 9.0 * ph_y_re;
  r = 500.0 * oh_y_im + 9.0 * ph_y_im;
  if (r == 0.0) {
    if (ai == 0.0) {
      k_y_re = ar / y_re;
    } else if (ar == 0.0) {
      k_y_re = 0.0;
    } else {
      k_y_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      k_y_re = ai / r;
    } else if (ai == 0.0) {
      k_y_re = 0.0;
    } else {
      k_y_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      k_y_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        ln_y_re = 0.5;
      } else {
        ln_y_re = -0.5;
      }

      if (r > 0.0) {
        pc_r = 0.5;
      } else {
        pc_r = -0.5;
      }

      k_y_re = (ar * ln_y_re + ai * pc_r) / brm;
    } else {
      brm = y_re / r;
      k_y_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  j_y_im = 625.0 * qh_y_re;
  k_y_im = 625.0 * qh_y_im;
  ar = (j_y_im * t_s.re - k_y_im * t_s.im) * x06;
  ai = (j_y_im * t_s.im + k_y_im * t_s.re) * x06;
  y_re = 500.0 * rh_y_re + 9.0 * t_s_re_tmp;
  r = 500.0 * rh_y_im + 9.0 * t_s_im_tmp;
  if (r == 0.0) {
    if (ai == 0.0) {
      j_y_im = ar / y_re;
    } else if (ar == 0.0) {
      j_y_im = 0.0;
    } else {
      j_y_im = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      j_y_im = ai / r;
    } else if (ai == 0.0) {
      j_y_im = 0.0;
    } else {
      j_y_im = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      j_y_im = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        mn_y_re = 0.5;
      } else {
        mn_y_re = -0.5;
      }

      if (r > 0.0) {
        qc_r = 0.5;
      } else {
        qc_r = -0.5;
      }

      j_y_im = (ar * mn_y_re + ai * qc_r) / brm;
    } else {
      brm = y_re / r;
      j_y_im = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  j_y_re = 15.0 * sh_y_re;
  k_y_im = 15.0 * sh_y_im;
  ar = (j_y_re * t_s.re - k_y_im * t_s.im) * x06;
  ai = (j_y_re * t_s.im + k_y_im * t_s.re) * x06;
  y_re = 500.0 * th_y_re + 9.0 * uh_y_re;
  r = 500.0 * th_y_im + 9.0 * uh_y_im;
  if (r == 0.0) {
    if (ai == 0.0) {
      j_y_re = ar / y_re;
    } else if (ar == 0.0) {
      j_y_re = 0.0;
    } else {
      j_y_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      j_y_re = ai / r;
    } else if (ai == 0.0) {
      j_y_re = 0.0;
    } else {
      j_y_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      j_y_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        nn_y_re = 0.5;
      } else {
        nn_y_re = -0.5;
      }

      if (r > 0.0) {
        rc_r = 0.5;
      } else {
        rc_r = -0.5;
      }

      j_y_re = (ar * nn_y_re + ai * rc_r) / brm;
    } else {
      brm = y_re / r;
      j_y_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  i_y_im = 1875.0 * vh_y_re;
  k_y_im = 1875.0 * vh_y_im;
  ar = (i_y_im * t_s.re - k_y_im * t_s.im) * x03;
  ai = (i_y_im * t_s.im + k_y_im * t_s.re) * x03;
  y_re = 2.0 * (500.0 * wh_y_re + 9.0 * xh_y_re);
  r = 2.0 * (500.0 * wh_y_im + 9.0 * xh_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      i_y_im = ar / y_re;
    } else if (ar == 0.0) {
      i_y_im = 0.0;
    } else {
      i_y_im = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      i_y_im = ai / r;
    } else if (ai == 0.0) {
      i_y_im = 0.0;
    } else {
      i_y_im = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      i_y_im = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        on_y_re = 0.5;
      } else {
        on_y_re = -0.5;
      }

      if (r > 0.0) {
        sc_r = 0.5;
      } else {
        sc_r = -0.5;
      }

      i_y_im = (ar * on_y_re + ai * sc_r) / brm;
    } else {
      brm = y_re / r;
      i_y_im = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  i_y_re = 15.0 * yh_y_re;
  k_y_im = 15.0 * yh_y_im;
  ar = (i_y_re * t_s.re - k_y_im * t_s.im) * x06;
  ai = (i_y_re * t_s.im + k_y_im * t_s.re) * x06;
  y_re = 2.0 * (500.0 * ai_y_re + 9.0 * bi_y_re);
  r = 2.0 * (500.0 * ai_y_im + 9.0 * bi_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      i_y_re = ar / y_re;
    } else if (ar == 0.0) {
      i_y_re = 0.0;
    } else {
      i_y_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      i_y_re = ai / r;
    } else if (ai == 0.0) {
      i_y_re = 0.0;
    } else {
      i_y_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      i_y_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        pn_y_re = 0.5;
      } else {
        pn_y_re = -0.5;
      }

      if (r > 0.0) {
        tc_r = 0.5;
      } else {
        tc_r = -0.5;
      }

      i_y_re = (ar * pn_y_re + ai * tc_r) / brm;
    } else {
      brm = y_re / r;
      i_y_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  h_y_im = 375.0 * ci_y_re;
  k_y_im = 375.0 * ci_y_im;
  ar = (h_y_im * t_s.re - k_y_im * t_s.im) * x03;
  ai = (h_y_im * t_s.im + k_y_im * t_s.re) * x03;
  y_re = 2.0 * (500.0 * di_y_re + 9.0 * ei_y_re);
  r = 2.0 * (500.0 * di_y_im + 9.0 * ei_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      h_y_im = ar / y_re;
    } else if (ar == 0.0) {
      h_y_im = 0.0;
    } else {
      h_y_im = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      h_y_im = ai / r;
    } else if (ai == 0.0) {
      h_y_im = 0.0;
    } else {
      h_y_im = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      h_y_im = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        qn_y_re = 0.5;
      } else {
        qn_y_re = -0.5;
      }

      if (r > 0.0) {
        uc_r = 0.5;
      } else {
        uc_r = -0.5;
      }

      h_y_im = (ar * qn_y_re + ai * uc_r) / brm;
    } else {
      brm = y_re / r;
      h_y_im = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  h_y_re = 3.0 * fi_y_re;
  k_y_im = 3.0 * fi_y_im;
  ar = (h_y_re * t_s.re - k_y_im * t_s.im) * x06;
  ai = (h_y_re * t_s.im + k_y_im * t_s.re) * x06;
  y_re = 2.0 * (500.0 * gi_y_re + 9.0 * hi_y_re);
  r = 2.0 * (500.0 * gi_y_im + 9.0 * hi_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      h_y_re = ar / y_re;
    } else if (ar == 0.0) {
      h_y_re = 0.0;
    } else {
      h_y_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      h_y_re = ai / r;
    } else if (ai == 0.0) {
      h_y_re = 0.0;
    } else {
      h_y_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      h_y_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        rn_y_re = 0.5;
      } else {
        rn_y_re = -0.5;
      }

      if (r > 0.0) {
        vc_r = 0.5;
      } else {
        vc_r = -0.5;
      }

      h_y_re = (ar * rn_y_re + ai * vc_r) / brm;
    } else {
      brm = y_re / r;
      h_y_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  g_y_im = 1875.0 * ii_y_re;
  k_y_im = 1875.0 * ii_y_im;
  ar = (g_y_im * t_s.re - k_y_im * t_s.im) * x13;
  ai = (g_y_im * t_s.im + k_y_im * t_s.re) * x13;
  y_re = 500.0 * ji_y_re + 9.0 * ki_y_re;
  r = 500.0 * ji_y_im + 9.0 * ki_y_im;
  if (r == 0.0) {
    if (ai == 0.0) {
      g_y_im = ar / y_re;
    } else if (ar == 0.0) {
      g_y_im = 0.0;
    } else {
      g_y_im = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      g_y_im = ai / r;
    } else if (ai == 0.0) {
      g_y_im = 0.0;
    } else {
      g_y_im = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      g_y_im = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        sn_y_re = 0.5;
      } else {
        sn_y_re = -0.5;
      }

      if (r > 0.0) {
        wc_r = 0.5;
      } else {
        wc_r = -0.5;
      }

      g_y_im = (ar * sn_y_re + ai * wc_r) / brm;
    } else {
      brm = y_re / r;
      g_y_im = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  g_y_re = 1875.0 * li_y_re;
  k_y_im = 1875.0 * li_y_im;
  ar = (g_y_re * t_s.re - k_y_im * t_s.im) * x13;
  ai = (g_y_re * t_s.im + k_y_im * t_s.re) * x13;
  y_re = 2.0 * (500.0 * mi_y_re + 9.0 * ni_y_re);
  r = 2.0 * (500.0 * mi_y_im + 9.0 * ni_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      g_y_re = ar / y_re;
    } else if (ar == 0.0) {
      g_y_re = 0.0;
    } else {
      g_y_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      g_y_re = ai / r;
    } else if (ai == 0.0) {
      g_y_re = 0.0;
    } else {
      g_y_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      g_y_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        tn_y_re = 0.5;
      } else {
        tn_y_re = -0.5;
      }

      if (r > 0.0) {
        xc_r = 0.5;
      } else {
        xc_r = -0.5;
      }

      g_y_re = (ar * tn_y_re + ai * xc_r) / brm;
    } else {
      brm = y_re / r;
      g_y_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  f_y_im = 375.0 * oi_y_re;
  k_y_im = 375.0 * oi_y_im;
  ar = (f_y_im * t_s.re - k_y_im * t_s.im) * x13;
  ai = (f_y_im * t_s.im + k_y_im * t_s.re) * x13;
  y_re = 2.0 * (500.0 * pi_y_re + 9.0 * qi_y_re);
  r = 2.0 * (500.0 * pi_y_im + 9.0 * qi_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      f_y_im = ar / y_re;
    } else if (ar == 0.0) {
      f_y_im = 0.0;
    } else {
      f_y_im = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      f_y_im = ai / r;
    } else if (ai == 0.0) {
      f_y_im = 0.0;
    } else {
      f_y_im = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      f_y_im = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        un_y_re = 0.5;
      } else {
        un_y_re = -0.5;
      }

      if (r > 0.0) {
        yc_r = 0.5;
      } else {
        yc_r = -0.5;
      }

      f_y_im = (ar * un_y_re + ai * yc_r) / brm;
    } else {
      brm = y_re / r;
      f_y_im = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  f_y_re = 1875.0 * ri_y_re;
  k_y_im = 1875.0 * ri_y_im;
  ar = (f_y_re * t_s_re_tmp - k_y_im * t_s_im_tmp) * x06;
  ai = (f_y_re * t_s_im_tmp + k_y_im * t_s_re_tmp) * x06;
  y_re = 500.0 * si_y_re + 9.0 * ti_y_re;
  r = 500.0 * si_y_im + 9.0 * ti_y_im;
  if (r == 0.0) {
    if (ai == 0.0) {
      f_y_re = ar / y_re;
    } else if (ar == 0.0) {
      f_y_re = 0.0;
    } else {
      f_y_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      f_y_re = ai / r;
    } else if (ai == 0.0) {
      f_y_re = 0.0;
    } else {
      f_y_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      f_y_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        vn_y_re = 0.5;
      } else {
        vn_y_re = -0.5;
      }

      if (r > 0.0) {
        ad_r = 0.5;
      } else {
        ad_r = -0.5;
      }

      f_y_re = (ar * vn_y_re + ai * ad_r) / brm;
    } else {
      brm = y_re / r;
      f_y_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  e_y_im = 125.0 * ui_y_re;
  k_y_im = 125.0 * ui_y_im;
  ar = (e_y_im * t_s_re_tmp - k_y_im * t_s_im_tmp) * x09;
  ai = (e_y_im * t_s_im_tmp + k_y_im * t_s_re_tmp) * x09;
  y_re = 500.0 * vi_y_re + 9.0 * t_s_re_tmp;
  r = 500.0 * vi_y_im + 9.0 * t_s_im_tmp;
  if (r == 0.0) {
    if (ai == 0.0) {
      e_y_im = ar / y_re;
    } else if (ar == 0.0) {
      e_y_im = 0.0;
    } else {
      e_y_im = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      e_y_im = ai / r;
    } else if (ai == 0.0) {
      e_y_im = 0.0;
    } else {
      e_y_im = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      e_y_im = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        wn_y_re = 0.5;
      } else {
        wn_y_re = -0.5;
      }

      if (r > 0.0) {
        bd_r = 0.5;
      } else {
        bd_r = -0.5;
      }

      e_y_im = (ar * wn_y_re + ai * bd_r) / brm;
    } else {
      brm = y_re / r;
      e_y_im = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  e_y_re = 15.0 * wi_y_re;
  k_y_im = 15.0 * wi_y_im;
  ar = (e_y_re * t_s_re_tmp - k_y_im * t_s_im_tmp) * x09;
  ai = (e_y_re * t_s_im_tmp + k_y_im * t_s_re_tmp) * x09;
  y_re = 2.0 * (500.0 * xi_y_re + 9.0 * yi_y_re);
  r = 2.0 * (500.0 * xi_y_im + 9.0 * yi_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      e_y_re = ar / y_re;
    } else if (ar == 0.0) {
      e_y_re = 0.0;
    } else {
      e_y_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      e_y_re = ai / r;
    } else if (ai == 0.0) {
      e_y_re = 0.0;
    } else {
      e_y_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      e_y_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        xn_y_re = 0.5;
      } else {
        xn_y_re = -0.5;
      }

      if (r > 0.0) {
        cd_r = 0.5;
      } else {
        cd_r = -0.5;
      }

      e_y_re = (ar * xn_y_re + ai * cd_r) / brm;
    } else {
      brm = y_re / r;
      e_y_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  d_y_im = 1875.0 * aj_y_re;
  k_y_im = 1875.0 * aj_y_im;
  ar = (d_y_im * t_s_re_tmp - k_y_im * t_s_im_tmp) * x06;
  ai = (d_y_im * t_s_im_tmp + k_y_im * t_s_re_tmp) * x06;
  y_re = 2.0 * (500.0 * bj_y_re + 9.0 * cj_y_re);
  r = 2.0 * (500.0 * bj_y_im + 9.0 * cj_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      d_y_im = ar / y_re;
    } else if (ar == 0.0) {
      d_y_im = 0.0;
    } else {
      d_y_im = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      d_y_im = ai / r;
    } else if (ai == 0.0) {
      d_y_im = 0.0;
    } else {
      d_y_im = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      d_y_im = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        yn_y_re = 0.5;
      } else {
        yn_y_re = -0.5;
      }

      if (r > 0.0) {
        dd_r = 0.5;
      } else {
        dd_r = -0.5;
      }

      d_y_im = (ar * yn_y_re + ai * dd_r) / brm;
    } else {
      brm = y_re / r;
      d_y_im = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  d_y_re = 625.0 * dj_y_re;
  k_y_im = 625.0 * dj_y_im;
  ar = (d_y_re * ej_y_re - k_y_im * ej_y_im) * x09;
  ai = (d_y_re * ej_y_im + k_y_im * ej_y_re) * x09;
  y_re = 500.0 * fj_y_re + 9.0 * gj_y_re;
  r = 500.0 * fj_y_im + 9.0 * gj_y_im;
  if (r == 0.0) {
    if (ai == 0.0) {
      d_y_re = ar / y_re;
    } else if (ar == 0.0) {
      d_y_re = 0.0;
    } else {
      d_y_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      d_y_re = ai / r;
    } else if (ai == 0.0) {
      d_y_re = 0.0;
    } else {
      d_y_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      d_y_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        ao_y_re = 0.5;
      } else {
        ao_y_re = -0.5;
      }

      if (r > 0.0) {
        ed_r = 0.5;
      } else {
        ed_r = -0.5;
      }

      d_y_re = (ar * ao_y_re + ai * ed_r) / brm;
    } else {
      brm = y_re / r;
      d_y_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  c_y_im = 15.0 * hj_y_re;
  k_y_im = 15.0 * hj_y_im;
  ar = (c_y_im * t_s_re_tmp - k_y_im * t_s_im_tmp) * x09;
  ai = (c_y_im * t_s_im_tmp + k_y_im * t_s_re_tmp) * x09;
  y_re = 4.0 * (500.0 * ij_y_re + 9.0 * jj_y_re);
  r = 4.0 * (500.0 * ij_y_im + 9.0 * jj_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      c_y_im = ar / y_re;
    } else if (ar == 0.0) {
      c_y_im = 0.0;
    } else {
      c_y_im = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      c_y_im = ai / r;
    } else if (ai == 0.0) {
      c_y_im = 0.0;
    } else {
      c_y_im = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      c_y_im = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        bo_y_re = 0.5;
      } else {
        bo_y_re = -0.5;
      }

      if (r > 0.0) {
        fd_r = 0.5;
      } else {
        fd_r = -0.5;
      }

      c_y_im = (ar * bo_y_re + ai * fd_r) / brm;
    } else {
      brm = y_re / r;
      c_y_im = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  c_y_re = 375.0 * kj_y_re;
  k_y_im = 375.0 * kj_y_im;
  ar = (c_y_re * t_s_re_tmp - k_y_im * t_s_im_tmp) * x06;
  ai = (c_y_re * t_s_im_tmp + k_y_im * t_s_re_tmp) * x06;
  y_re = 2.0 * (500.0 * lj_y_re + 9.0 * mj_y_re);
  r = 2.0 * (500.0 * lj_y_im + 9.0 * mj_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      c_y_re = ar / y_re;
    } else if (ar == 0.0) {
      c_y_re = 0.0;
    } else {
      c_y_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      c_y_re = ai / r;
    } else if (ai == 0.0) {
      c_y_re = 0.0;
    } else {
      c_y_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      c_y_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        co_y_re = 0.5;
      } else {
        co_y_re = -0.5;
      }

      if (r > 0.0) {
        gd_r = 0.5;
      } else {
        gd_r = -0.5;
      }

      c_y_re = (ar * co_y_re + ai * gd_r) / brm;
    } else {
      brm = y_re / r;
      c_y_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  b_y_im = 625.0 * nj_y_re;
  k_y_im = 625.0 * nj_y_im;
  ar = (b_y_im * oj_y_re - k_y_im * oj_y_im) * x09;
  ai = (b_y_im * oj_y_im + k_y_im * oj_y_re) * x09;
  y_re = 2.0 * (500.0 * pj_y_re + 9.0 * qj_y_re);
  r = 2.0 * (500.0 * pj_y_im + 9.0 * qj_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      b_y_im = ar / y_re;
    } else if (ar == 0.0) {
      b_y_im = 0.0;
    } else {
      b_y_im = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      b_y_im = ai / r;
    } else if (ai == 0.0) {
      b_y_im = 0.0;
    } else {
      b_y_im = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      b_y_im = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        do_y_re = 0.5;
      } else {
        do_y_re = -0.5;
      }

      if (r > 0.0) {
        hd_r = 0.5;
      } else {
        hd_r = -0.5;
      }

      b_y_im = (ar * do_y_re + ai * hd_r) / brm;
    } else {
      brm = y_re / r;
      b_y_im = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  b_y_re = 3.0 * rj_y_re;
  k_y_im = 3.0 * rj_y_im;
  ar = (b_y_re * t_s_re_tmp - k_y_im * t_s_im_tmp) * x09;
  ai = (b_y_re * t_s_im_tmp + k_y_im * t_s_re_tmp) * x09;
  y_re = 4.0 * (500.0 * sj_y_re + 9.0 * tj_y_re);
  r = 4.0 * (500.0 * sj_y_im + 9.0 * tj_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      b_y_re = ar / y_re;
    } else if (ar == 0.0) {
      b_y_re = 0.0;
    } else {
      b_y_re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      b_y_re = ai / r;
    } else if (ai == 0.0) {
      b_y_re = 0.0;
    } else {
      b_y_re = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      b_y_re = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        eo_y_re = 0.5;
      } else {
        eo_y_re = -0.5;
      }

      if (r > 0.0) {
        id_r = 0.5;
      } else {
        id_r = -0.5;
      }

      b_y_re = (ar * eo_y_re + ai * id_r) / brm;
    } else {
      brm = y_re / r;
      b_y_re = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  r = 125.0 * uj_y_re;
  k_y_im = 125.0 * uj_y_im;
  ar = (r * vj_y_re - k_y_im * vj_y_im) * x09;
  ai = (r * vj_y_im + k_y_im * vj_y_re) * x09;
  y_re = 2.0 * (500.0 * wj_y_re + 9.0 * xj_y_re);
  r = 2.0 * (500.0 * wj_y_im + 9.0 * xj_y_im);
  if (r == 0.0) {
    if (ai == 0.0) {
      r = ar / y_re;
    } else if (ar == 0.0) {
      r = 0.0;
    } else {
      r = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      r = ai / r;
    } else if (ai == 0.0) {
      r = 0.0;
    } else {
      r = ai / r;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(r);
    if (brm > bim) {
      brm = r / y_re;
      r = (ar + brm * ai) / (y_re + brm * r);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        fo_y_re = 0.5;
      } else {
        fo_y_re = -0.5;
      }

      if (r > 0.0) {
        jd_r = 0.5;
      } else {
        jd_r = -0.5;
      }

      r = (ar * fo_y_re + ai * jd_r) / brm;
    } else {
      brm = y_re / r;
      r = (brm * ar + ai) / (r + brm * y_re);
    }
  }

  states[2] = (((((((((((((((((((((((((((((x03 + t.re * x06) + t_re) + re) -
    b_re) + c_re) - d_re) - m_y_im) + m_y_re) - l_y_im) + l_y_re) - k_y_re) +
    j_y_im) - j_y_re) + i_y_im) + i_y_re) - h_y_im) - h_y_re) + g_y_im) - g_y_re)
                        + f_y_im) - f_y_re) + e_y_im) - e_y_re) + d_y_im) -
                   d_y_re) + c_y_im) - c_y_re) + b_y_im) - b_y_re) - r;
  re = 5.0 * t_re_tmp;
  k_y_im = 5.0 * t_im_tmp;
  b_re = 6.0 * t_re_tmp * x01 - 6.0 * t_re_tmp * x11;
  f_y_im = 6.0 * t_im_tmp * x01 - 6.0 * t_im_tmp * x11;
  c_re = re * b_re - k_y_im * f_y_im;
  k_y_im = re * f_y_im + k_y_im * b_re;
  if (k_y_im == 0.0) {
    re = c_re / 4.0;
    k_y_im = 0.0;
  } else if (c_re == 0.0) {
    re = 0.0;
    k_y_im /= 4.0;
  } else {
    re = c_re / 4.0;
    k_y_im /= 4.0;
  }

  b_re = 5.0 * t_re_tmp;
  f_y_im = 5.0 * t_im_tmp;
  c_re = b_re * yj_y_re - f_y_im * yj_y_im;
  f_y_im = b_re * yj_y_im + f_y_im * yj_y_re;
  h_y_re = 3000.0 * t.re;
  h_y_im = 3000.0 * t.im;
  i_y_re = 12.0 * t.re;
  i_y_im = 12.0 * t.im;
  b_re = ((((3000.0 * x01 + 36.0 * x04) - 3000.0 * x11) - h_y_re * x04) - i_y_re
          * x07) + 250.0 * t_re_tmp * x07;
  d_y_im = ((0.0 - h_y_im * x04) - i_y_im * x07) + 250.0 * t_im_tmp * x07;
  d_re = c_re * b_re - f_y_im * d_y_im;
  f_y_im = c_re * d_y_im + f_y_im * b_re;
  if (f_y_im == 0.0) {
    b_re = d_re / 4.0;
    f_y_im = 0.0;
  } else if (d_re == 0.0) {
    b_re = 0.0;
    f_y_im /= 4.0;
  } else {
    b_re = d_re / 4.0;
    f_y_im /= 4.0;
  }

  c_re = 5.0 * t_re_tmp;
  d_y_im = 5.0 * t_im_tmp;
  d_re = c_re * t_s_re_tmp - d_y_im * t_s_im_tmp;
  d_y_im = c_re * t_s_im_tmp + d_y_im * t_s_re_tmp;
  j_y_re = 24.0 * t.re;
  j_y_im = 24.0 * t.im;
  c_re = (((((36.0 * x01 - 36.0 * x11) - h_y_re * x01) - j_y_re * x04) + h_y_re *
           x11) + 750.0 * t_re_tmp * x04) + 3.0 * t_re_tmp * x07;
  e_y_re = ((((0.0 - h_y_im * x01) - j_y_im * x04) + h_y_im * x11) + 750.0 *
            t_im_tmp * x04) + 3.0 * t_im_tmp * x07;
  m_y_im = d_re * c_re - d_y_im * e_y_re;
  d_y_im = d_re * e_y_re + d_y_im * c_re;
  if (d_y_im == 0.0) {
    c_re = m_y_im / 4.0;
    d_y_im = 0.0;
  } else if (m_y_im == 0.0) {
    c_re = 0.0;
    d_y_im /= 4.0;
  } else {
    c_re = m_y_im / 4.0;
    d_y_im /= 4.0;
  }

  d_re = 5.0 * t_re_tmp;
  e_y_re = 5.0 * t_im_tmp;
  m_y_im = d_re * ak_y_re - e_y_re * ak_y_im;
  e_y_re = d_re * ak_y_im + e_y_re * ak_y_re;
  e_y_im = 1000.0 * t.re;
  f_y_re = 1000.0 * t.im;
  d_re = (3000.0 * x04 + 18.0 * x07) - e_y_im * x07;
  c_y_im = 0.0 - f_y_re * x07;
  m_y_re = m_y_im * d_re - e_y_re * c_y_im;
  e_y_re = m_y_im * c_y_im + e_y_re * d_re;
  if (e_y_re == 0.0) {
    d_re = m_y_re / 4.0;
    e_y_re = 0.0;
  } else if (m_y_re == 0.0) {
    d_re = 0.0;
    e_y_re /= 4.0;
  } else {
    d_re = m_y_re / 4.0;
    e_y_re /= 4.0;
  }

  m_y_im = 1500.0 * t_re_tmp;
  c_y_im = 1500.0 * t_im_tmp;
  m_y_re = 5.0 * t_re_tmp;
  d_y_re = 5.0 * t_im_tmp;
  l_y_im = m_y_re * t_s.re - d_y_re * t_s.im;
  d_y_re = m_y_re * t_s.im + d_y_re * t_s.re;
  m_y_re = (((j_y_re * x11 - j_y_re * x01) + 750.0 * t_re_tmp * x01) + 6.0 *
            t_re_tmp * x04) - 750.0 * t_re_tmp * x11;
  c_y_re = (((j_y_im * x11 - j_y_im * x01) + 750.0 * t_im_tmp * x01) + 6.0 *
            t_im_tmp * x04) - 750.0 * t_im_tmp * x11;
  l_y_re = l_y_im * m_y_re - d_y_re * c_y_re;
  d_y_re = l_y_im * c_y_re + d_y_re * m_y_re;
  if (d_y_re == 0.0) {
    m_y_re = l_y_re / 4.0;
    d_y_re = 0.0;
  } else if (l_y_re == 0.0) {
    m_y_re = 0.0;
    d_y_re /= 4.0;
  } else {
    m_y_re = l_y_re / 4.0;
    d_y_re /= 4.0;
  }

  g_y_re = 500.0 * t_s.re + 9.0;
  g_y_im = 500.0 * t_s.im;
  y_re = ck_y_re * g_y_re - ck_y_im * g_y_im;
  y_im = ck_y_re * g_y_im + ck_y_im * g_y_re;
  ar = ((((re + b_re) + c_re) + d_re) + (m_y_im * bk_y_re - c_y_im * bk_y_im) *
        x07) + m_y_re;
  ai = ((((k_y_im + f_y_im) + d_y_im) + e_y_re) + (m_y_im * bk_y_im + c_y_im *
         bk_y_re) * x07) + d_y_re;
  if (y_im == 0.0) {
    if (ai == 0.0) {
      re = ar / y_re;
    } else if (ar == 0.0) {
      re = 0.0;
    } else {
      re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      re = ai / y_im;
    } else if (ai == 0.0) {
      re = 0.0;
    } else {
      re = ai / y_im;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(y_im);
    if (brm > bim) {
      brm = y_im / y_re;
      re = (ar + brm * ai) / (y_re + brm * y_im);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        go_y_re = 0.5;
      } else {
        go_y_re = -0.5;
      }

      if (y_im > 0.0) {
        wk_y_im = 0.5;
      } else {
        wk_y_im = -0.5;
      }

      re = (ar * go_y_re + ai * wk_y_im) / brm;
    } else {
      brm = y_re / y_im;
      re = (brm * ar + ai) / (y_im + brm * y_re);
    }
  }

  states[3] = (x04 + t.re * x07) - re;
  re = 5.0 * t_re_tmp;
  k_y_im = 5.0 * t_im_tmp;
  b_re = 6.0 * t_re_tmp * x02 - 6.0 * t_re_tmp * x12;
  f_y_im = 6.0 * t_im_tmp * x02 - 6.0 * t_im_tmp * x12;
  c_re = re * b_re - k_y_im * f_y_im;
  k_y_im = re * f_y_im + k_y_im * b_re;
  if (k_y_im == 0.0) {
    re = c_re / 4.0;
    k_y_im = 0.0;
  } else if (c_re == 0.0) {
    re = 0.0;
    k_y_im /= 4.0;
  } else {
    re = c_re / 4.0;
    k_y_im /= 4.0;
  }

  b_re = 5.0 * t_re_tmp;
  f_y_im = 5.0 * t_im_tmp;
  c_re = b_re * dk_y_re - f_y_im * dk_y_im;
  f_y_im = b_re * dk_y_im + f_y_im * dk_y_re;
  b_re = ((((3000.0 * x02 + 36.0 * x05) - 3000.0 * x12) - h_y_re * x05) - i_y_re
          * x08) + 250.0 * t_re_tmp * x08;
  d_y_im = ((0.0 - h_y_im * x05) - i_y_im * x08) + 250.0 * t_im_tmp * x08;
  d_re = c_re * b_re - f_y_im * d_y_im;
  f_y_im = c_re * d_y_im + f_y_im * b_re;
  if (f_y_im == 0.0) {
    b_re = d_re / 4.0;
    f_y_im = 0.0;
  } else if (d_re == 0.0) {
    b_re = 0.0;
    f_y_im /= 4.0;
  } else {
    b_re = d_re / 4.0;
    f_y_im /= 4.0;
  }

  c_re = 5.0 * t_re_tmp;
  d_y_im = 5.0 * t_im_tmp;
  d_re = c_re * t_s_re_tmp - d_y_im * t_s_im_tmp;
  d_y_im = c_re * t_s_im_tmp + d_y_im * t_s_re_tmp;
  c_re = (((((36.0 * x02 - 36.0 * x12) - h_y_re * x02) - j_y_re * x05) + h_y_re *
           x12) + 750.0 * t_re_tmp * x05) + 3.0 * t_re_tmp * x08;
  e_y_re = ((((0.0 - h_y_im * x02) - j_y_im * x05) + h_y_im * x12) + 750.0 *
            t_im_tmp * x05) + 3.0 * t_im_tmp * x08;
  m_y_im = d_re * c_re - d_y_im * e_y_re;
  d_y_im = d_re * e_y_re + d_y_im * c_re;
  if (d_y_im == 0.0) {
    c_re = m_y_im / 4.0;
    d_y_im = 0.0;
  } else if (m_y_im == 0.0) {
    c_re = 0.0;
    d_y_im /= 4.0;
  } else {
    c_re = m_y_im / 4.0;
    d_y_im /= 4.0;
  }

  d_re = 5.0 * t_re_tmp;
  e_y_re = 5.0 * t_im_tmp;
  m_y_im = d_re * ek_y_re - e_y_re * ek_y_im;
  e_y_re = d_re * ek_y_im + e_y_re * ek_y_re;
  d_re = (3000.0 * x05 + 18.0 * x08) - e_y_im * x08;
  c_y_im = 0.0 - f_y_re * x08;
  m_y_re = m_y_im * d_re - e_y_re * c_y_im;
  e_y_re = m_y_im * c_y_im + e_y_re * d_re;
  if (e_y_re == 0.0) {
    d_re = m_y_re / 4.0;
    e_y_re = 0.0;
  } else if (m_y_re == 0.0) {
    d_re = 0.0;
    e_y_re /= 4.0;
  } else {
    d_re = m_y_re / 4.0;
    e_y_re /= 4.0;
  }

  m_y_im = 1500.0 * t_re_tmp;
  c_y_im = 1500.0 * t_im_tmp;
  m_y_re = 5.0 * t_re_tmp;
  d_y_re = 5.0 * t_im_tmp;
  l_y_im = m_y_re * t_s.re - d_y_re * t_s.im;
  d_y_re = m_y_re * t_s.im + d_y_re * t_s.re;
  m_y_re = (((j_y_re * x12 - j_y_re * x02) + 750.0 * t_re_tmp * x02) + 6.0 *
            t_re_tmp * x05) - 750.0 * t_re_tmp * x12;
  c_y_re = (((j_y_im * x12 - j_y_im * x02) + 750.0 * t_im_tmp * x02) + 6.0 *
            t_im_tmp * x05) - 750.0 * t_im_tmp * x12;
  l_y_re = l_y_im * m_y_re - d_y_re * c_y_re;
  d_y_re = l_y_im * c_y_re + d_y_re * m_y_re;
  if (d_y_re == 0.0) {
    m_y_re = l_y_re / 4.0;
    d_y_re = 0.0;
  } else if (l_y_re == 0.0) {
    m_y_re = 0.0;
    d_y_re /= 4.0;
  } else {
    m_y_re = l_y_re / 4.0;
    d_y_re /= 4.0;
  }

  y_re = gk_y_re * g_y_re - gk_y_im * g_y_im;
  y_im = gk_y_re * g_y_im + gk_y_im * g_y_re;
  ar = ((((re + b_re) + c_re) + d_re) + (m_y_im * fk_y_re - c_y_im * fk_y_im) *
        x08) + m_y_re;
  ai = ((((k_y_im + f_y_im) + d_y_im) + e_y_re) + (m_y_im * fk_y_im + c_y_im *
         fk_y_re) * x08) + d_y_re;
  if (y_im == 0.0) {
    if (ai == 0.0) {
      re = ar / y_re;
    } else if (ar == 0.0) {
      re = 0.0;
    } else {
      re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      re = ai / y_im;
    } else if (ai == 0.0) {
      re = 0.0;
    } else {
      re = ai / y_im;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(y_im);
    if (brm > bim) {
      brm = y_im / y_re;
      re = (ar + brm * ai) / (y_re + brm * y_im);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        ho_y_re = 0.5;
      } else {
        ho_y_re = -0.5;
      }

      if (y_im > 0.0) {
        xk_y_im = 0.5;
      } else {
        xk_y_im = -0.5;
      }

      re = (ar * ho_y_re + ai * xk_y_im) / brm;
    } else {
      brm = y_re / y_im;
      re = (brm * ar + ai) / (y_im + brm * y_re);
    }
  }

  states[4] = (x05 + t.re * x08) - re;
  re = 5.0 * t_re_tmp;
  k_y_im = 5.0 * t_im_tmp;
  b_re = 6.0 * t_re_tmp * x03 - 6.0 * t_re_tmp * x13;
  f_y_im = 6.0 * t_im_tmp * x03 - 6.0 * t_im_tmp * x13;
  c_re = re * b_re - k_y_im * f_y_im;
  k_y_im = re * f_y_im + k_y_im * b_re;
  if (k_y_im == 0.0) {
    re = c_re / 4.0;
    k_y_im = 0.0;
  } else if (c_re == 0.0) {
    re = 0.0;
    k_y_im /= 4.0;
  } else {
    re = c_re / 4.0;
    k_y_im /= 4.0;
  }

  b_re = 5.0 * t_re_tmp;
  f_y_im = 5.0 * t_im_tmp;
  c_re = b_re * hk_y_re - f_y_im * hk_y_im;
  f_y_im = b_re * hk_y_im + f_y_im * hk_y_re;
  b_re = ((((3000.0 * x03 + 36.0 * x06) - 3000.0 * x13) - h_y_re * x06) - i_y_re
          * x09) + 250.0 * t_re_tmp * x09;
  d_y_im = ((0.0 - h_y_im * x06) - i_y_im * x09) + 250.0 * t_im_tmp * x09;
  d_re = c_re * b_re - f_y_im * d_y_im;
  f_y_im = c_re * d_y_im + f_y_im * b_re;
  if (f_y_im == 0.0) {
    b_re = d_re / 4.0;
    f_y_im = 0.0;
  } else if (d_re == 0.0) {
    b_re = 0.0;
    f_y_im /= 4.0;
  } else {
    b_re = d_re / 4.0;
    f_y_im /= 4.0;
  }

  c_re = 5.0 * t_re_tmp;
  d_y_im = 5.0 * t_im_tmp;
  d_re = c_re * t_s_re_tmp - d_y_im * t_s_im_tmp;
  d_y_im = c_re * t_s_im_tmp + d_y_im * t_s_re_tmp;
  c_re = (((((36.0 * x03 - 36.0 * x13) - h_y_re * x03) - j_y_re * x06) + h_y_re *
           x13) + 750.0 * t_re_tmp * x06) + 3.0 * t_re_tmp * x09;
  e_y_re = ((((0.0 - h_y_im * x03) - j_y_im * x06) + h_y_im * x13) + 750.0 *
            t_im_tmp * x06) + 3.0 * t_im_tmp * x09;
  m_y_im = d_re * c_re - d_y_im * e_y_re;
  d_y_im = d_re * e_y_re + d_y_im * c_re;
  if (d_y_im == 0.0) {
    c_re = m_y_im / 4.0;
    d_y_im = 0.0;
  } else if (m_y_im == 0.0) {
    c_re = 0.0;
    d_y_im /= 4.0;
  } else {
    c_re = m_y_im / 4.0;
    d_y_im /= 4.0;
  }

  d_re = 5.0 * t_re_tmp;
  e_y_re = 5.0 * t_im_tmp;
  m_y_im = d_re * ik_y_re - e_y_re * ik_y_im;
  e_y_re = d_re * ik_y_im + e_y_re * ik_y_re;
  d_re = (3000.0 * x06 + 18.0 * x09) - e_y_im * x09;
  c_y_im = 0.0 - f_y_re * x09;
  m_y_re = m_y_im * d_re - e_y_re * c_y_im;
  e_y_re = m_y_im * c_y_im + e_y_re * d_re;
  if (e_y_re == 0.0) {
    d_re = m_y_re / 4.0;
    e_y_re = 0.0;
  } else if (m_y_re == 0.0) {
    d_re = 0.0;
    e_y_re /= 4.0;
  } else {
    d_re = m_y_re / 4.0;
    e_y_re /= 4.0;
  }

  m_y_im = 1500.0 * t_re_tmp;
  c_y_im = 1500.0 * t_im_tmp;
  m_y_re = 5.0 * t_re_tmp;
  d_y_re = 5.0 * t_im_tmp;
  l_y_im = m_y_re * t_s.re - d_y_re * t_s.im;
  d_y_re = m_y_re * t_s.im + d_y_re * t_s.re;
  m_y_re = (((j_y_re * x13 - j_y_re * x03) + 750.0 * t_re_tmp * x03) + 6.0 *
            t_re_tmp * x06) - 750.0 * t_re_tmp * x13;
  c_y_re = (((j_y_im * x13 - j_y_im * x03) + 750.0 * t_im_tmp * x03) + 6.0 *
            t_im_tmp * x06) - 750.0 * t_im_tmp * x13;
  l_y_re = l_y_im * m_y_re - d_y_re * c_y_re;
  d_y_re = l_y_im * c_y_re + d_y_re * m_y_re;
  if (d_y_re == 0.0) {
    m_y_re = l_y_re / 4.0;
    d_y_re = 0.0;
  } else if (l_y_re == 0.0) {
    m_y_re = 0.0;
    d_y_re /= 4.0;
  } else {
    m_y_re = l_y_re / 4.0;
    d_y_re /= 4.0;
  }

  y_re = kk_y_re * g_y_re - kk_y_im * g_y_im;
  y_im = kk_y_re * g_y_im + kk_y_im * g_y_re;
  ar = ((((re + b_re) + c_re) + d_re) + (m_y_im * jk_y_re - c_y_im * jk_y_im) *
        x09) + m_y_re;
  ai = ((((k_y_im + f_y_im) + d_y_im) + e_y_re) + (m_y_im * jk_y_im + c_y_im *
         jk_y_re) * x09) + d_y_re;
  if (y_im == 0.0) {
    if (ai == 0.0) {
      re = ar / y_re;
    } else if (ar == 0.0) {
      re = 0.0;
    } else {
      re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      re = ai / y_im;
    } else if (ai == 0.0) {
      re = 0.0;
    } else {
      re = ai / y_im;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(y_im);
    if (brm > bim) {
      brm = y_im / y_re;
      re = (ar + brm * ai) / (y_re + brm * y_im);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        io_y_re = 0.5;
      } else {
        io_y_re = -0.5;
      }

      if (y_im > 0.0) {
        yk_y_im = 0.5;
      } else {
        yk_y_im = -0.5;
      }

      re = (ar * io_y_re + ai * yk_y_im) / brm;
    } else {
      brm = y_re / y_im;
      re = (brm * ar + ai) / (y_im + brm * y_re);
    }
  }

  states[5] = (x06 + t.re * x09) - re;
  re = x07 * h_y_re;
  k_y_im = x07 * h_y_im;
  i_y_re = 750.0 * t.re;
  i_y_im = 750.0 * t.im;
  j_y_re = 5.0 * t.re;
  j_y_im = 5.0 * t.im;
  b_re = (1500.0 * x04 + 9.0 * x07) - i_y_re * x07;
  f_y_im = 0.0 - i_y_im * x07;
  c_re = j_y_re * b_re - j_y_im * f_y_im;
  f_y_im = j_y_re * f_y_im + j_y_im * b_re;
  e_y_im = 2250.0 * t.re;
  f_y_re = 2250.0 * t.im;
  c_y_re = 9.0 * t.re;
  r = 9.0 * t.im;
  b_re = ((((1500.0 * x01 + 18.0 * x04) - 1500.0 * x11) - e_y_im * x04) - c_y_re
          * x07) + 250.0 * t_re_tmp * x07;
  d_y_im = ((0.0 - f_y_re * x04) - r * x07) + 250.0 * t_im_tmp * x07;
  d_re = j_y_re * b_re - j_y_im * d_y_im;
  d_y_im = j_y_re * d_y_im + j_y_im * b_re;
  b_y_re = 18.0 * t.re;
  b_y_im = 18.0 * t.im;
  b_re = (((((18.0 * x01 - 18.0 * x11) - e_y_im * x01) - b_y_re * x04) + e_y_im *
           x11) + 750.0 * t_re_tmp * x04) + 3.0 * t_re_tmp * x07;
  e_y_re = ((((0.0 - f_y_re * x01) - b_y_im * x04) + f_y_re * x11) + 750.0 *
            t_im_tmp * x04) + 3.0 * t_im_tmp * x07;
  m_y_im = j_y_re * b_re - j_y_im * e_y_re;
  e_y_re = j_y_re * e_y_re + j_y_im * b_re;
  b_re = (((b_y_re * x11 - b_y_re * x01) + 750.0 * t_re_tmp * x01) + 6.0 *
          t_re_tmp * x04) - 750.0 * t_re_tmp * x11;
  c_y_im = (((b_y_im * x11 - b_y_im * x01) + 750.0 * t_im_tmp * x01) + 6.0 *
            t_im_tmp * x04) - 750.0 * t_im_tmp * x11;
  m_y_re = j_y_re * b_re - j_y_im * c_y_im;
  c_y_im = j_y_re * c_y_im + j_y_im * b_re;
  b_re = 6.0 * t_re_tmp * x01 - 6.0 * t_re_tmp * x11;
  d_y_re = 6.0 * t_im_tmp * x01 - 6.0 * t_im_tmp * x11;
  y_re = ok_y_re * g_y_re - ok_y_im * g_y_im;
  y_im = ok_y_re * g_y_im + ok_y_im * g_y_re;
  ar = (((((re * lk_y_re - k_y_im * lk_y_im) + (c_re * mk_y_re - f_y_im *
            mk_y_im)) + (d_re * nk_y_re - d_y_im * nk_y_im)) + (m_y_im *
          t_s_re_tmp - e_y_re * t_s_im_tmp)) + (m_y_re * t_s.re - c_y_im *
         t_s.im)) + (j_y_re * b_re - j_y_im * d_y_re);
  ai = (((((re * lk_y_im + k_y_im * lk_y_re) + (c_re * mk_y_im + f_y_im *
            mk_y_re)) + (d_re * nk_y_im + d_y_im * nk_y_re)) + (m_y_im *
          t_s_im_tmp + e_y_re * t_s_re_tmp)) + (m_y_re * t_s.im + c_y_im *
         t_s.re)) + (j_y_re * d_y_re + j_y_im * b_re);
  if (y_im == 0.0) {
    if (ai == 0.0) {
      re = ar / y_re;
    } else if (ar == 0.0) {
      re = 0.0;
    } else {
      re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      re = ai / y_im;
    } else if (ai == 0.0) {
      re = 0.0;
    } else {
      re = ai / y_im;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(y_im);
    if (brm > bim) {
      brm = y_im / y_re;
      re = (ar + brm * ai) / (y_re + brm * y_im);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        jo_y_re = 0.5;
      } else {
        jo_y_re = -0.5;
      }

      if (y_im > 0.0) {
        al_y_im = 0.5;
      } else {
        al_y_im = -0.5;
      }

      re = (ar * jo_y_re + ai * al_y_im) / brm;
    } else {
      brm = y_re / y_im;
      re = (brm * ar + ai) / (y_im + brm * y_re);
    }
  }

  states[6] = x07 - re;
  re = x08 * h_y_re;
  k_y_im = x08 * h_y_im;
  b_re = (1500.0 * x05 + 9.0 * x08) - i_y_re * x08;
  f_y_im = 0.0 - i_y_im * x08;
  c_re = j_y_re * b_re - j_y_im * f_y_im;
  f_y_im = j_y_re * f_y_im + j_y_im * b_re;
  b_re = ((((1500.0 * x02 + 18.0 * x05) - 1500.0 * x12) - e_y_im * x05) - c_y_re
          * x08) + 250.0 * t_re_tmp * x08;
  d_y_im = ((0.0 - f_y_re * x05) - r * x08) + 250.0 * t_im_tmp * x08;
  d_re = j_y_re * b_re - j_y_im * d_y_im;
  d_y_im = j_y_re * d_y_im + j_y_im * b_re;
  b_re = (((((18.0 * x02 - 18.0 * x12) - e_y_im * x02) - b_y_re * x05) + e_y_im *
           x12) + 750.0 * t_re_tmp * x05) + 3.0 * t_re_tmp * x08;
  e_y_re = ((((0.0 - f_y_re * x02) - b_y_im * x05) + f_y_re * x12) + 750.0 *
            t_im_tmp * x05) + 3.0 * t_im_tmp * x08;
  m_y_im = j_y_re * b_re - j_y_im * e_y_re;
  e_y_re = j_y_re * e_y_re + j_y_im * b_re;
  b_re = (((b_y_re * x12 - b_y_re * x02) + 750.0 * t_re_tmp * x02) + 6.0 *
          t_re_tmp * x05) - 750.0 * t_re_tmp * x12;
  c_y_im = (((b_y_im * x12 - b_y_im * x02) + 750.0 * t_im_tmp * x02) + 6.0 *
            t_im_tmp * x05) - 750.0 * t_im_tmp * x12;
  m_y_re = j_y_re * b_re - j_y_im * c_y_im;
  c_y_im = j_y_re * c_y_im + j_y_im * b_re;
  b_re = 6.0 * t_re_tmp * x02 - 6.0 * t_re_tmp * x12;
  d_y_re = 6.0 * t_im_tmp * x02 - 6.0 * t_im_tmp * x12;
  y_re = sk_y_re * g_y_re - sk_y_im * g_y_im;
  y_im = sk_y_re * g_y_im + sk_y_im * g_y_re;
  ar = (((((re * pk_y_re - k_y_im * pk_y_im) + (c_re * qk_y_re - f_y_im *
            qk_y_im)) + (d_re * rk_y_re - d_y_im * rk_y_im)) + (m_y_im *
          t_s_re_tmp - e_y_re * t_s_im_tmp)) + (m_y_re * t_s.re - c_y_im *
         t_s.im)) + (j_y_re * b_re - j_y_im * d_y_re);
  ai = (((((re * pk_y_im + k_y_im * pk_y_re) + (c_re * qk_y_im + f_y_im *
            qk_y_re)) + (d_re * rk_y_im + d_y_im * rk_y_re)) + (m_y_im *
          t_s_im_tmp + e_y_re * t_s_re_tmp)) + (m_y_re * t_s.im + c_y_im *
         t_s.re)) + (j_y_re * d_y_re + j_y_im * b_re);
  if (y_im == 0.0) {
    if (ai == 0.0) {
      re = ar / y_re;
    } else if (ar == 0.0) {
      re = 0.0;
    } else {
      re = ar / y_re;
    }
  } else if (y_re == 0.0) {
    if (ar == 0.0) {
      re = ai / y_im;
    } else if (ai == 0.0) {
      re = 0.0;
    } else {
      re = ai / y_im;
    }
  } else {
    brm = fabs(y_re);
    bim = fabs(y_im);
    if (brm > bim) {
      brm = y_im / y_re;
      re = (ar + brm * ai) / (y_re + brm * y_im);
    } else if (bim == brm) {
      if (y_re > 0.0) {
        ko_y_re = 0.5;
      } else {
        ko_y_re = -0.5;
      }

      if (y_im > 0.0) {
        bl_y_im = 0.5;
      } else {
        bl_y_im = -0.5;
      }

      re = (ar * ko_y_re + ai * bl_y_im) / brm;
    } else {
      brm = y_re / y_im;
      re = (brm * ar + ai) / (y_im + brm * y_re);
    }
  }

  states[7] = x08 - re;
  re = x09 * h_y_re;
  k_y_im = x09 * h_y_im;
  b_re = (1500.0 * x06 + 9.0 * x09) - i_y_re * x09;
  f_y_im = 0.0 - i_y_im * x09;
  c_re = j_y_re * b_re - j_y_im * f_y_im;
  f_y_im = j_y_re * f_y_im + j_y_im * b_re;
  b_re = ((((1500.0 * x03 + 18.0 * x06) - 1500.0 * x13) - e_y_im * x06) - c_y_re
          * x09) + 250.0 * t_re_tmp * x09;
  d_y_im = ((0.0 - f_y_re * x06) - r * x09) + 250.0 * t_im_tmp * x09;
  d_re = j_y_re * b_re - j_y_im * d_y_im;
  d_y_im = j_y_re * d_y_im + j_y_im * b_re;
  b_re = (((((18.0 * x03 - 18.0 * x13) - e_y_im * x03) - b_y_re * x06) + e_y_im *
           x13) + 750.0 * t_re_tmp * x06) + 3.0 * t_re_tmp * x09;
  e_y_re = ((((0.0 - f_y_re * x03) - b_y_im * x06) + f_y_re * x13) + 750.0 *
            t_im_tmp * x06) + 3.0 * t_im_tmp * x09;
  m_y_im = j_y_re * b_re - j_y_im * e_y_re;
  e_y_re = j_y_re * e_y_re + j_y_im * b_re;
  b_re = (((b_y_re * x13 - b_y_re * x03) + 750.0 * t_re_tmp * x03) + 6.0 *
          t_re_tmp * x06) - 750.0 * t_re_tmp * x13;
  c_y_im = (((b_y_im * x13 - b_y_im * x03) + 750.0 * t_im_tmp * x03) + 6.0 *
            t_im_tmp * x06) - 750.0 * t_im_tmp * x13;
  m_y_re = j_y_re * b_re - j_y_im * c_y_im;
  c_y_im = j_y_re * c_y_im + j_y_im * b_re;
  b_re = 6.0 * t_re_tmp * x03 - 6.0 * t_re_tmp * x13;
  d_y_re = 6.0 * t_im_tmp * x03 - 6.0 * t_im_tmp * x13;
  c_y_re = x_re * g_y_re - x_im * g_y_im;
  x_im = x_re * g_y_im + x_im * g_y_re;
  ar = (((((re * tk_y_re - k_y_im * tk_y_im) + (c_re * uk_y_re - f_y_im *
            uk_y_im)) + (d_re * vk_y_re - d_y_im * vk_y_im)) + (m_y_im *
          t_s_re_tmp - e_y_re * t_s_im_tmp)) + (m_y_re * t_s.re - c_y_im *
         t_s.im)) + (j_y_re * b_re - j_y_im * d_y_re);
  ai = (((((re * tk_y_im + k_y_im * tk_y_re) + (c_re * uk_y_im + f_y_im *
            uk_y_re)) + (d_re * vk_y_im + d_y_im * vk_y_re)) + (m_y_im *
          t_s_im_tmp + e_y_re * t_s_re_tmp)) + (m_y_re * t_s.im + c_y_im *
         t_s.re)) + (j_y_re * d_y_re + j_y_im * b_re);
  if (x_im == 0.0) {
    if (ai == 0.0) {
      re = ar / c_y_re;
    } else if (ar == 0.0) {
      re = 0.0;
    } else {
      re = ar / c_y_re;
    }
  } else if (c_y_re == 0.0) {
    if (ar == 0.0) {
      re = ai / x_im;
    } else if (ai == 0.0) {
      re = 0.0;
    } else {
      re = ai / x_im;
    }
  } else {
    brm = fabs(c_y_re);
    bim = fabs(x_im);
    if (brm > bim) {
      brm = x_im / c_y_re;
      re = (ar + brm * ai) / (c_y_re + brm * x_im);
    } else if (bim == brm) {
      if (c_y_re > 0.0) {
        lo_y_re = 0.5;
      } else {
        lo_y_re = -0.5;
      }

      if (x_im > 0.0) {
        b_x_im = 0.5;
      } else {
        b_x_im = -0.5;
      }

      re = (ar * lo_y_re + ai * b_x_im) / brm;
    } else {
      brm = c_y_re / x_im;
      re = (brm * ar + ai) / (x_im + brm * c_y_re);
    }
  }

  states[8] = x09 - re;
}

/* End of code generation (quadf_statePFFull.c) */
