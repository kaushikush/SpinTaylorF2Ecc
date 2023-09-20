/* The equations used below have obtained from
 * References:                                                       *
 * T. Damour, A. Gopakumar, and B. Iyer                              *
 * (PRD 70 064028)                                                   *
 * Equations 71a, 71b                                                *
 *                                                                   *
 * K. G. Arun, Luc Blanchet, Bala R. Iyer and Siddharta Sinha        *
 * e-print: arXiv:0908.3854 [gr-qc]                                  *
 *                                                                   *
 * I. Hinder, F. Herrmann, P. Laguna and D. Shoemaker                *
 * e-Print: arXiv:0806.1037 [gr-qc]                                  *
 * Equations A35, A36                                                *
 *                                                                   *
 * and derived in the paper where this model is introduced:          *
 *                                                                   *
 * Huerta et al                                                      *
 * e-Print: arXiv:                                                   *
 * Equations                                                         */

static REAL8 x_dot_0pn(REAL8 e, REAL8 eta) /* Eq. (A26) */
{
  REAL8 x_0_pn;
  REAL8 e_pow_2 = e * e;
  REAL8 e_pow_4 = e_pow_2 * e_pow_2;
  REAL8 e_fact = 1.0 - e_pow_2;
  REAL8 num = 2. * eta * (37. * e_pow_4 + 292. * e_pow_2 + 96.);
  REAL8 den = 15. * e_fact * e_fact * e_fact * sqrt(e_fact);
  REAL8 e_0_lim = 2. * eta * 96. / 15.;

  if (e) {
    x_0_pn = num / den;

  } else {
    x_0_pn = e_0_lim;
    ;
  }

  return (x_0_pn);
}

static REAL8 x_dot_1pn(REAL8 e, REAL8 eta) /* Eq. (A27) */
{
  REAL8 x_1_pn;
  REAL8 e_pow_2 = e * e;
  REAL8 e_pow_4 = e_pow_2 * e_pow_2;
  REAL8 e_pow_6 = e_pow_4 * e_pow_2;
  REAL8 e_fact = 1. - e_pow_2;
  REAL8 e_0_term = 11888. + 14784. * eta;
  REAL8 e_2_term = e_pow_2 * (-87720. + 159600. * eta);
  REAL8 e_4_term = e_pow_4 * (-171038. + 141708. * eta);
  REAL8 e_6_term = e_pow_6 * (-11717. + 8288. * eta);
  REAL8 den = 420. * e_fact * e_fact * e_fact * e_fact * sqrt(e_fact);
  REAL8 num = -eta * (e_0_term + e_2_term + e_4_term + e_6_term);
  REAL8 e_0_lim = -eta * (2972. / 105 + 176. * eta / 5.);

  if (e) {
    x_1_pn = num / den;
  } else {
    x_1_pn = e_0_lim;
  }

  return (x_1_pn);
}

static REAL8 x_dot_1_5_pn(REAL8 e, REAL8 eta, REAL8 m1, REAL8 m2, REAL8 S1z,
                          REAL8 S2z) /*computed using inputs from klein et al.
                                        arXiv:1801.08542. See equation C1a.*/
{
  REAL8 x_1_5_pn;
  REAL8 pre_factor = 64. * eta / 5;
  REAL8 e_pow_2 = e * e;
  REAL8 e_pow_4 = e_pow_2 * e_pow_2;
  REAL8 e_pow_6 = e_pow_4 * e_pow_2;
  REAL8 e_fact = 1.0 - e_pow_2;

  if (e) {
    x_1_5_pn =
        -((m1 * m2 *
           ((5424 + 27608 * e_pow_2 + 16694 * e_pow_4 + 585 * e_pow_6) * m1 *
                m1 * S1z +
            (5424 + 27608 * e_pow_2 + 16694 * e_pow_4 + 585 * e_pow_6) * m2 *
                m2 * S2z +
            3 * (1200 + 6976 * e_pow_2 + 4886 * e_pow_4 + 207 * e_pow_6) * m1 *
                m2 * (S1z + S2z))) /
          (45. * e_fact * e_fact * e_fact * e_fact * e_fact * (m1 + m2) *
           (m1 + m2) * (m1 + m2) * (m1 + m2)));
  } else {

    x_1_5_pn = pre_factor * (-0.08333333333333333 *
                             (113 * pow(m1, 2) * S1z + 113 * pow(m2, 2) * S2z +
                              75 * m1 * m2 * (S1z + S2z)) /
                             pow(m1 + m2, 2));
  }
  return (x_1_5_pn);
}

static REAL8 x_dot_hereditary_1_5(REAL8 e, REAL8 eta, REAL8 x) /* Eq. (A28) */
{
  REAL8 x_dot_her_1_5;
  REAL8 x_term = 256. * M_PI * phi_e(e) / 5.;
  REAL8 x_term_e_0 = 256. * M_PI / 5.;
  REAL8 pre_factor = eta * x * x * x * x * x * x * sqrt(x);

  if (e) {
    x_dot_her_1_5 = pre_factor * x_term;
  } else {
    x_dot_her_1_5 = pre_factor * x_term_e_0;
  }

  return (x_dot_her_1_5);
}

static REAL8 x_dot_2pn_SS(REAL8 e, REAL8 eta, REAL8 m1, REAL8 m2, REAL8 S1z,
                          REAL8 S2z) /*computed using inputs from klein et al.
                                        arXiv:1801.08542. See equation C1a.*/
{
  REAL8 kappa1 = 1.0; /*for black holes kappa_{1,2} is 1*/
  REAL8 kappa2 = 1.0;
  REAL8 x_2pn_SS;
  REAL8 pre_factor = 64. * eta / 5;
  REAL8 e_2 = e * e;
  REAL8 e_4 = e_2 * e_2;
  REAL8 e_6 = e_4 * e_2;
  REAL8 e_fact = (1 - e_2);
  REAL8 e_fact_2 = e_fact * e_fact;
  REAL8 e_fact_sqrt = sqrt(e_fact);
  REAL8 e_fact_55 = e_fact_2 * e_fact_2 * e_fact * e_fact_sqrt;
  REAL8 M_fact_2 = (m1 + m2) * (m1 + m2);
  REAL8 M_fact_4 = M_fact_2 * M_fact_2;

  if (e) {
    /* Quentin Henry et al terms, arXiv:2308.13606 */
    x_2pn_SS =
        ((m1 * m2 *
          ((48 * (1 + 80 * kappa1) + 3 * e_6 * (9 + 236 * kappa1) +
            8 * e_2 * (57 + 2692 * kappa1) + 2 * e_4 * (207 + 7472 * kappa1)) *
               m1 * m1 * S1z * S1z +
           2 * (3792 + 21080 * e_2 + 14530 * e_4 + 681 * e_6) * m1 * m2 * S1z *
               S2z +
           (48 * (1 + 80 * kappa2) + 3 * e_6 * (9 + 236 * kappa2) +
            8 * e_2 * (57 + 2692 * kappa2) + 2 * e_4 * (207 + 7472 * kappa2)) *
               m2 * m2 * S2z * S2z)) /
         (60. * e_fact_55 * M_fact_4));

    /* Klein et al. terms */
    /*x_2pn_SS =
        (-0.016666666666666666 *
         (m1 * m2 *
          ((-48 * (1 + 80 * kappa1) + 3 * pow(e, 6) * (-9 + 319 * kappa1) -
            8 * pow(e, 2) * (57 + 2692 * kappa1) -
            2 * pow(e, 4) * (207 + 7472 * kappa1)) *
               pow(m1, 2) * pow(S1z, 2) -
           2 *
               (3792 + 21080 * pow(e, 2) + 14530 * pow(e, 4) +
                681 * pow(e, 6)) *
               m1 * m2 * S1z * S2z +
           (-48 * (1 + 80 * kappa2) + 3 * pow(e, 6) * (-9 + 319 * kappa2) -
            8 * pow(e, 2) * (57 + 2692 * kappa2) -
            2 * pow(e, 4) * (207 + 7472 * kappa2)) *
               pow(m2, 2) * pow(S2z, 2))) /
         (pow(1 - pow(e, 2), 5.5) * pow(m1 + m2, 4)));*/

  } else {

    x_2pn_SS = pre_factor * (((1 + 80 * kappa1) * pow(m1, 2) * pow(S1z, 2) +
                              158 * m1 * m2 * S1z * S2z +
                              (1 + 80 * kappa2) * pow(m2, 2) * pow(S2z, 2)) /
                             (16. * pow(m1 + m2, 2)));
  }

  return (x_2pn_SS);
}

static REAL8 x_dot_hereditary_2_5(REAL8 e, REAL8 eta,
                                  REAL8 x) /* See Huerta et al article */
{
  REAL8 x_2_5_pn;
  REAL8 e_fact = 1. - e * e;
  REAL8 e_fact_pow_2 = e_fact * e_fact;
  REAL8 e_pow_2 = e * e;

  REAL8 bit_1 = 256. * M_PI * phi_e(e) / e_fact;
  REAL8 bit_2 = -17599. * M_PI * psi_n(e) / 35. -
                2268. * eta * M_PI * zed_n(e) / 5. -
                788. * e_pow_2 * M_PI * phi_e_rad(e) / e_fact_pow_2;
  REAL8 x_2_5_term = bit_1 + 2 * bit_2 / 3.;

  REAL8 bit_1_e_0 = 256. * M_PI;
  REAL8 bit_2_e_0 = -17599. * M_PI / 35. - 2268. * eta * M_PI / 5.;
  REAL8 x_2_5_term_e_0 = bit_1_e_0 + 2. * bit_2_e_0 / 3.;

  REAL8 pre_factor = eta * x * x * x * x * x * x * x * sqrt(x);

  if (e) {
    x_2_5_pn = pre_factor * x_2_5_term;
  } else {
    x_2_5_pn = pre_factor * x_2_5_term_e_0;
  }

  return (x_2_5_pn);
}

static REAL8 x_dot_hereditary_3(REAL8 e, REAL8 eta,
                                REAL8 x) /* See Huerta et al article */
{
  REAL8 x_3_pn_her;
  REAL8 pi_pow_2 = M_PI * M_PI;
  REAL8 euler = LAL_GAMMA;

  const REAL8 log4 = 1.38629436111989061883446424292; // ln(4)

  REAL8 x_3_term =
      64. *
      (-116761. * kappa_e(e) +
       (19600. * pi_pow_2 - 59920. * euler - 59920. * log4 - 89880. * log(x)) *
           f_e(e)) /
      18375.;

  REAL8 x_3_term_e_0 = 64. *
                       (-59920. * euler - 116761. + 19600. * pi_pow_2 -
                        59920. * log4 - 89880. * log(x)) /
                       18375.;

  REAL8 pre_factor = eta * x * x * x * x * x * x * x * x;

  if (e) {
    x_3_pn_her = pre_factor * x_3_term;
  } else {
    x_3_pn_her = pre_factor * x_3_term_e_0;
  }

  return (x_3_pn_her);
}

static REAL8 x_dot_2_5pn_SO(REAL8 e, REAL8 eta, REAL8 m1, REAL8 m2, REAL8 S1z,
                            REAL8 S2z) /* Computed using inputs from energy and
  flux expressions. See Blanchet liv. rev.
  + Bohe et al. arXiv: 1501.01529 + Marsat et al. arXiv: 1411.4118*/
{
  REAL8 x_2_5pn;
  REAL8 pre_factor = 64. * eta / 5;
  REAL8 e_2 = e * e;
  REAL8 e_4 = e_2 * e_2;
  REAL8 e_6 = e_4 * e_2;
  REAL8 e_8 = e_6 * e_2;
  REAL8 e_fact = (1 - e_2);
  REAL8 e_fact_sqrt = sqrt(e_fact);
  REAL8 M_fact_2 = (m1 + m2) * (m1 + m2);
  REAL8 M_fact_4 = M_fact_2 * M_fact_2;
  REAL8 M_fact_6 = M_fact_4 * M_fact_2;

  if (e) {
    /* Quentin Henry et al terms, arXiv:2308.13606 */
    x_2_5pn = (-0.0000992063492063492 *
               (m1 * m2 *
                ((4008832 + 808515 * e_8 +
                  896 * e_4 * (126373 + 26748 * e_fact_sqrt) +
                  16 * e_6 * (2581907 + 30576 * e_fact_sqrt) +
                  384 * e_2 * (111403 + 61264 * e_fact_sqrt)) *
                     m1 * m1 * m1 * m1 * S1z +
                 (4008832 + 808515 * e_8 +
                  896 * e_4 * (126373 + 26748 * e_fact_sqrt) +
                  16 * e_6 * (2581907 + 30576 * e_fact_sqrt) +
                  384 * e_2 * (111403 + 61264 * e_fact_sqrt)) *
                     m2 * m2 * m2 * m2 * S2z +
                 (1962112 + 1803645 * e_8 +
                  32 * e_6 * (2542879 + 26754 * e_fact_sqrt) +
                  896 * e_4 * (202075 + 46809 * e_fact_sqrt) +
                  768 * e_2 * (51189 + 53606 * e_fact_sqrt)) *
                     m1 * m1 * m2 * m2 * (S1z + S2z) +
                 m1 * m1 * m1 * m2 *
                     ((3029504 + 1807155 * e_8 +
                       64 * e_6 * (1314169 + 16562 * e_fact_sqrt) +
                       128 * e_2 * (340325 + 398216 * e_fact_sqrt) +
                       112 * e_4 * (1732273 + 463632 * e_fact_sqrt)) *
                          S1z +
                      3 *
                          (414208 + 281775 * e_8 +
                           112 * e_6 * (103639 + 728 * e_fact_sqrt) +
                           336 * e_4 * (77531 + 11888 * e_fact_sqrt) +
                           128 * e_2 * (55999 + 30632 * e_fact_sqrt)) *
                          S2z) +
                 m1 * m2 * m2 * m2 *
                     (3 *
                          (414208 + 281775 * e_8 +
                           112 * e_6 * (103639 + 728 * e_fact_sqrt) +
                           336 * e_4 * (77531 + 11888 * e_fact_sqrt) +
                           128 * e_2 * (55999 + 30632 * e_fact_sqrt)) *
                          S1z +
                      (3029504 + 1807155 * e_8 +
                       64 * e_6 * (1314169 + 16562 * e_fact_sqrt) +
                       128 * e_2 * (340325 + 398216 * e_fact_sqrt) +
                       112 * e_4 * (1732273 + 463632 * e_fact_sqrt)) *
                          S2z))) /
               (pow(-1 + e_2, 6) * M_fact_6));

  } else {

    x_2_5pn =
        pre_factor *
        (-0.000992063492063492 *
         (31319 * m1 * m1 * m1 * m1 * S1z + 31319 * m2 * m2 * m2 * m2 * S2z +
          15329 * m1 * m1 * m2 * m2 * (S1z + S2z) +
          4 * m1 * m1 * m1 * m2 * (5917 * S1z + 2427 * S2z) +
          4 * m1 * m2 * m2 * m2 * (2427 * S1z + 5917 * S2z)) /
         M_fact_4);
  }

  return (x_2_5pn);
}

static REAL8 x_dot_3pn_SO(REAL8 e, REAL8 eta, REAL8 m1, REAL8 m2, REAL8 S1z,
                          REAL8 S2z) /* Computed using inputs from energy and
 flux expressions. See Blanchet liv. rev.
 + Bohe et al. arXiv: 1501.01529 + Marsat et al. arXiv: 1411.4118*/
{
  REAL8 x_3_pn_SO;
  REAL8 pre_factor = 64. * eta / 5;
  REAL8 e_2 = e * e;
  REAL8 e_4 = e_2 * e_2;
  REAL8 e_6 = e_4 * e_2;
  REAL8 e_8 = e_6 * e_2;
  REAL8 e_fact = (1 - e_2);
  REAL8 e_fact_sqrt = sqrt(e_fact);
  REAL8 e_fact_65 =
      e_fact * e_fact * e_fact * e_fact * e_fact * e_fact * e_fact_sqrt;
  REAL8 M_fact_2 = (m1 + m2) * (m1 + m2);
  REAL8 M_fact_4 = M_fact_2 * M_fact_2;
  REAL8 M_fact_6 = M_fact_4 * M_fact_2;

  if (e) {

    x_3_pn_SO = (-9.645061728395062e-6 *
                 (m1 * m2 * M_PI *
                  ((49766400 + 528887808 * e_2 + 814424832 * e_4 +
                    213166272 * e_6 + 3911917 * e_8) *
                       m1 * m1 * S1z +
                   (49766400 + 528887808 * e_2 + 814424832 * e_4 +
                    213166272 * e_6 + 3911917 * e_8) *
                       m2 * m2 * S2z +
                   3 *
                       (11132928 + 133936128 * e_2 + 232455168 * e_4 +
                        67616624 * e_6 + 1479919 * e_8) *
                       m1 * m2 * (S1z + S2z))) /
                 (e_fact_65 * M_fact_4));
  } else {

    x_3_pn_SO =
        pre_factor * (-0.16666666666666666 *
                      (M_PI * (225 * pow(m1, 2) * S1z + 225 * pow(m2, 2) * S2z +
                               151 * m1 * m2 * (S1z + S2z))) /
                      M_fact_2);
  }

  return (x_3_pn_SO);
}

static REAL8 x_dot_3pn_SS(REAL8 e, REAL8 eta, REAL8 m1, REAL8 m2, REAL8 S1z,
                          REAL8 S2z) /* Computed using inputs from energy and
 flux expressions. See Blanchet liv. rev.
 + Bohe et al. arXiv: 1501.01529 + Marsat et al. arXiv: 1411.4118*/
{
  REAL8 x_3pn_SS;
  REAL8 kappa1 = 1.0; /*for black holes kappa_{1,2} is 1*/
  REAL8 kappa2 = 1.0;
  REAL8 pre_factor = 64. * eta / 5;
  REAL8 e_2 = e * e;
  REAL8 e_4 = e_2 * e_2;
  REAL8 e_6 = e_4 * e_2;
  REAL8 e_8 = e_6 * e_2;
  REAL8 e_fact = (1 - e_2);
  REAL8 e_fact_sqrt = sqrt(e_fact);
  REAL8 e_fact_65 =
      e_fact * e_fact * e_fact * e_fact * e_fact * e_fact * e_fact_sqrt;
  REAL8 M_fact_2 = (m1 + m2) * (m1 + m2);
  REAL8 M_fact_4 = M_fact_2 * M_fact_2;
  REAL8 M_fact_6 = M_fact_4 * M_fact_2;

  if (e) {

    x_3pn_SS = ((m1 * m2 *
                 ((27 * e_8 * (15141 + 109852 * kappa1) +
                   8 * e_4 *
                       (22676605 + 3044160 * e_fact_sqrt + 32290722 * kappa1 +
                        5327280 * e_fact_sqrt * kappa1) +
                   84 * e_6 *
                       (437079 + 448 * e_fact_sqrt +
                        14 * (89253 + 56 * e_fact_sqrt) * kappa1) -
                   576 * (-37891 + 896 * e_fact_sqrt +
                          2 * (-8963 + 784 * e_fact_sqrt) * kappa1) +
                   224 * e_2 *
                       (633619 + 107616 * e_fact_sqrt +
                        42 * (10029 + 4484 * e_fact_sqrt) * kappa1)) *
                      m1 * m1 * m1 * m1 * S1z * S1z +
                  (27 * e_8 * (15141 + 109852 * kappa2) +
                   8 * e_4 *
                       (22676605 + 3044160 * e_fact_sqrt + 32290722 * kappa2 +
                        5327280 * e_fact_sqrt * kappa2) +
                   84 * e_6 *
                       (437079 + 448 * e_fact_sqrt +
                        14 * (89253 + 56 * e_fact_sqrt) * kappa2) -
                   576 * (-37891 + 896 * e_fact_sqrt +
                          2 * (-8963 + 784 * e_fact_sqrt) * kappa2) +
                   224 * e_2 *
                       (633619 + 107616 * e_fact_sqrt +
                        42 * (10029 + 4484 * e_fact_sqrt) * kappa2)) *
                      m2 * m2 * m2 * m2 * S2z * S2z +
                  6 * m1 * m1 * m1 * m2 * S1z *
                      ((3 * e_8 * (47103 + 202177 * kappa1) -
                        96 * (-34713 + 336 * e_fact_sqrt - 8440 * kappa1 +
                              2576 * e_fact_sqrt * kappa1) +
                        112 * e_2 *
                            (278677 + 13452 * e_fact_sqrt + 53216 * kappa1 +
                             103132 * e_fact_sqrt * kappa1) +
                        14 * e_6 *
                            (772539 + 168 * e_fact_sqrt +
                             4 * (369781 + 322 * e_fact_sqrt) * kappa1) +
                        4 * e_4 *
                            (7 * (1655621 + 54360 * e_fact_sqrt) +
                             460 * (22397 + 6342 * e_fact_sqrt) * kappa1)) *
                           S1z +
                       2 *
                           (521613 * e_8 + 96 * (31457 - 1680 * e_fact_sqrt) +
                            294 * e_6 * (72619 + 40 * e_fact_sqrt) +
                            112 * e_2 * (247661 + 67260 * e_fact_sqrt) +
                            e_4 * (60534556 + 7610400 * e_fact_sqrt)) *
                           S2z) +
                  6 * m1 * m2 * m2 * m2 * S2z *
                      (2 *
                           (521613 * e_8 + 96 * (31457 - 1680 * e_fact_sqrt) +
                            294 * e_6 * (72619 + 40 * e_fact_sqrt) +
                            112 * e_2 * (247661 + 67260 * e_fact_sqrt) +
                            e_4 * (60534556 + 7610400 * e_fact_sqrt)) *
                           S1z +
                       (3 * e_8 * (47103 + 202177 * kappa2) -
                        96 * (-34713 + 336 * e_fact_sqrt - 8440 * kappa2 +
                              2576 * e_fact_sqrt * kappa2) +
                        112 * e_2 *
                            (278677 + 13452 * e_fact_sqrt + 53216 * kappa2 +
                             103132 * e_fact_sqrt * kappa2) +
                        14 * e_6 *
                            (772539 + 168 * e_fact_sqrt +
                             4 * (369781 + 322 * e_fact_sqrt) * kappa2) +
                        4 * e_4 *
                            (7 * (1655621 + 54360 * e_fact_sqrt) +
                             460 * (22397 + 6342 * e_fact_sqrt) * kappa2)) *
                           S2z) +
                  m1 * m1 * m2 * m2 *
                      (9 *
                           (-86016 * e_fact_sqrt * kappa1 +
                            192 * (1729 + 112 * e_fact_sqrt + 1766 * kappa1) +
                            e_8 * (58359 + 325902 * kappa1) +
                            28 * e_6 *
                                (121449 - 56 * e_fact_sqrt +
                                 (388890 + 224 * e_fact_sqrt) * kappa1) +
                            224 * e_2 *
                                (31367 - 4484 * e_fact_sqrt +
                                 2 * (12189 + 8968 * e_fact_sqrt) * kappa1) +
                            8 * e_4 *
                                (1541617 - 126840 * e_fact_sqrt +
                                 6 * (482187 + 84560 * e_fact_sqrt) * kappa1)) *
                           S1z * S1z +
                       8 *
                           (8135280 + 1021401 * e_8 - 467712 * e_fact_sqrt +
                            147 * e_6 * (305971 + 232 * e_fact_sqrt) +
                            56 * e_2 * (1084885 + 390108 * e_fact_sqrt) +
                            2 * e_4 * (65163991 + 11035080 * e_fact_sqrt)) *
                           S1z * S2z +
                       9 *
                           (-86016 * e_fact_sqrt * kappa2 +
                            192 * (1729 + 112 * e_fact_sqrt + 1766 * kappa2) +
                            e_8 * (58359 + 325902 * kappa2) +
                            28 * e_6 *
                                (121449 - 56 * e_fact_sqrt +
                                 (388890 + 224 * e_fact_sqrt) * kappa2) +
                            224 * e_2 *
                                (31367 - 4484 * e_fact_sqrt +
                                 2 * (12189 + 8968 * e_fact_sqrt) * kappa2) +
                            8 * e_4 *
                                (1541617 - 126840 * e_fact_sqrt +
                                 6 * (482187 + 84560 * e_fact_sqrt) * kappa2)) *
                           S2z * S2z))) /
                (30240. * e_fact_65 * M_fact_6));

  } else {

    x_3pn_SS = pre_factor *
               (((36995 + 16358 * kappa1) * pow(m1, 4) * pow(S1z, 2) +
                 (36995 + 16358 * kappa2) * pow(m2, 4) * pow(S2z, 2) +
                 pow(m1, 3) * m2 * S1z *
                     ((34377 + 5864 * kappa1) * S1z + 59554 * S2z) +
                 m1 * pow(m2, 3) * S2z *
                     (59554 * S1z + (34377 + 5864 * kappa2) * S2z) +
                 3 * pow(m1, 2) * pow(m2, 2) *
                     ((1841 + 1318 * kappa1) * pow(S1z, 2) + 35498 * S1z * S2z +
                      (1841 + 1318 * kappa2) * pow(S2z, 2))) /
                (672. * M_fact_4));
  }

  return (x_3pn_SS);
}

static REAL8 x_dot_2pn(REAL8 e, REAL8 eta) /* Eq. (A29) */
{
  REAL8 x_2_pn;
  REAL8 e_pow_2 = e * e;
  REAL8 e_pow_4 = e_pow_2 * e_pow_2;
  REAL8 e_pow_6 = e_pow_4 * e_pow_2;
  REAL8 e_pow_8 = e_pow_4 * e_pow_4;
  REAL8 e_fact = 1.0 - e_pow_2;
  REAL8 eta_pow_2 = eta * eta;
  REAL8 den =
      45360. * e_fact * e_fact * e_fact * e_fact * e_fact * sqrt(e_fact);
  REAL8 e_0_term = -360224. + 4514976. * eta + 1903104. * eta_pow_2;
  REAL8 e_2_term =
      e_pow_2 * (-92846560. + 15464736. * eta + 61282032. * eta_pow_2);
  REAL8 e_4_term =
      e_pow_4 * (783768. - 207204264. * eta + 166506060. * eta_pow_2);
  REAL8 e_6_term =
      e_pow_6 * (83424402. - 123108426. * eta + 64828848. * eta_pow_2);
  REAL8 e_8_term = e_pow_8 * (3523113. - 3259980. * eta + 1964256. * eta_pow_2);
  REAL8 rt_e_0_term = 1451520. - 580608. * eta;
  REAL8 rt_e_2_term = e_pow_2 * (64532160. - 25812864. * eta);
  REAL8 rt_e_4_term = e_pow_4 * (66316320. - 26526528. * eta);
  REAL8 rt_e_6_term = e_pow_6 * (2646000. - 1058400. * eta);

  REAL8 num =
      eta *
      (e_0_term + e_2_term + e_4_term + e_6_term + e_8_term +
       sqrt(e_fact) * (rt_e_0_term + rt_e_2_term + rt_e_4_term + rt_e_6_term));

  REAL8 e_0_lim =
      eta * (68206. / 2835. + 27322. * eta / 315. + 1888. * eta_pow_2 / 45.);

  if (e) {
    x_2_pn = num / den;
  } else {
    x_2_pn = e_0_lim;
  }

  return (x_2_pn);
}

static REAL8 x_dot_3pn(REAL8 e, REAL8 eta,
                       REAL8 x) /* See Huerta et al article */
{
  const REAL8 log2 = 0.693147180559945309417232121458; // ln(2)
  REAL8 x_3_pn;
  REAL8 e_pow_2 = e * e;
  REAL8 e_pow_4 = e_pow_2 * e_pow_2;
  REAL8 e_pow_6 = e_pow_4 * e_pow_2;
  REAL8 e_pow_8 = e_pow_4 * e_pow_4;
  REAL8 e_pow_10 = e_pow_8 * e_pow_2;
  REAL8 e_fact = 1.0 - e_pow_2;
  REAL8 eta_pow_2 = eta * eta;
  REAL8 pi_pow_2 = M_PI * M_PI;
  REAL8 den =
      598752000. * e_fact * e_fact * e_fact * e_fact * e_fact * e_fact * e_fact;
  REAL8 bit_1 = -25. * e_pow_10 *
                (81. * (99269280. - 33332681. * sqrt(e_fact)) +
                 176. * eta *
                     (9. * (-5132796. * +1874543. * sqrt(e_fact)) +
                      4. * eta *
                          (3582684. - 2962791. * sqrt(e_fact) +
                           2320640. * sqrt(e_fact) * eta)));
  REAL8 bit_2 = -128. * (3950984268. - 12902173599. * sqrt(e_fact) +
                         275. * eta *
                             (-1066392. + 57265081. * sqrt(e_fact) +
                              81. * (-17696. + 16073. * sqrt(e_fact)) * eta +
                              470820. * sqrt(e_fact) * eta_pow_2 -
                              46494. * (-1. + 45 * sqrt(e_fact)) * pi_pow_2));
  REAL8 bit_3 = -32. * e_pow_2 *
                (-18. * (2603019496. + 19904214811. * sqrt(e_fact)) +
                 55. * eta *
                     (8147179440. - 5387647438. * sqrt(e_fact) +
                      270. * (-6909392. + 9657701. * sqrt(e_fact)) * eta +
                      901169500. * sqrt(e_fact) * eta_pow_2 -
                      193725. * (229. + 237. * sqrt(e_fact)) * pi_pow_2));
  REAL8 bit_4 = -8. * e_pow_4 *
                (-6. * (312191560692. + 8654689873. * sqrt(e_fact)) +
                 55. * eta *
                     (42004763280. - 88628306866. * sqrt(e_fact) -
                      1350. * (8601376. + 1306589. * sqrt(e_fact)) * eta +
                      23638717900. * sqrt(e_fact) * eta_pow_2 +
                      891135. * (-2. + 627. * sqrt(e_fact)) * pi_pow_2));
  REAL8 bit_5 = -2. * e_pow_8 *
                (4351589277552. - 1595548875627. * sqrt(e_fact) +
                 550. * eta *
                     (432. * (6368264. - 10627167. * sqrt(e_fact)) * eta +
                      2201124800. * sqrt(e_fact) * eta_pow_2 +
                      9. * (8. * (-134041982. + 65136045. * sqrt(e_fact)) +
                            861. * (14. + 891. * sqrt(e_fact)) * pi_pow_2)));
  REAL8 bit_6 =
      -12. * e_pow_6 *
      (589550775792. - 6005081022. * sqrt(e_fact) +
       55. * eta *
           (90. * (90130656. - 311841025. * sqrt(e_fact)) * eta +
            17925404000. * sqrt(e_fact) * eta_pow_2 +
            3. * (-2. * (5546517920. + 383583403. * sqrt(e_fact)) +
                  4305. * (9046. + 19113. * sqrt(e_fact)) * pi_pow_2)));
  REAL8 bit_7 = -40677120. * sqrt(e_fact) *
                (3072. + 43520. * e_pow_2 + 82736. * e_pow_4 +
                 28016. * e_pow_6 + 891. * e_pow_8) *
                (log2 - log((1. / e_fact) + (1. / sqrt(e_fact))) - log(x));
  REAL8 num = eta * (bit_1 + bit_2 + bit_3 + bit_4 + bit_5 + bit_6 + bit_7);
  REAL8 e_0_lim =
      eta * (426247111. / 222750. - 56198689. * eta / 17010. +
             541. * eta * eta / 70. - 2242. * eta * eta * eta / 81. +
             1804. * eta * pi_pow_2 / 15. + 109568. * log(x) / 525.);

  if (e) {
    x_3_pn = num / den;
  } else {
    x_3_pn = e_0_lim;
  }

  return (x_3_pn);
}

static REAL8 x_dot_3_5pnSO(REAL8 e, REAL8 eta, REAL8 m1, REAL8 m2, REAL8 S1z,
                           REAL8 S2z) /* Computed using inputs from energy and
flux expressions. See Blanchet liv. rev.
+ Bohe et al. arXiv: 1501.01529 + Marsat et al. arXiv: 1411.4118*/
{
  REAL8 x_3_5pnSO;

  REAL8 pre_factor = 64. * eta / 5;

  if (e) {

    x_3_5pnSO = 0.0;

  } else {

    x_3_5pnSO = pre_factor *
                (-0.00005511463844797178 *
                 (3127800 * pow(m1, 6) * S1z + 3127800 * pow(m2, 6) * S2z +
                  4914306 * pow(m1, 3) * pow(m2, 3) * (S1z + S2z) +
                  pow(m1, 5) * m2 * (6542338 * S1z + 1195759 * S2z) +
                  pow(m1, 4) * pow(m2, 2) * (6694579 * S1z + 3284422 * S2z) +
                  m1 * pow(m2, 5) * (1195759 * S1z + 6542338 * S2z) +
                  pow(m1, 2) * pow(m2, 4) * (3284422 * S1z + 6694579 * S2z)) /
                 pow(m1 + m2, 6));
  }

  return (x_3_5pnSO);
}

static REAL8 x_dot_3_5pn_SS(REAL8 e, REAL8 eta, REAL8 m1, REAL8 m2, REAL8 S1z,
                            REAL8 S2z) /* Computed using inputs from energy and
flux expressions. See Blanchet liv. rev.
+ Bohe et al. arXiv: 1501.01529 + Marsat et al. arXiv: 1411.4118*/
{
  REAL8 x_3_5pn_SS;
  REAL8 kappa1 = 1.0; /*for black holes kappa_{1,2} is 1*/
  REAL8 kappa2 = 1.0;
  REAL8 pre_factor = 64. * eta / 5;

  if (e) {

    return (0.0);
  } else {
    /* previous */ /* x_3_5pn_SS = pre_factor *
                      ((12*M_PI*(kappa1*pow(m1,2)*pow(S1z,2) + m2*S2z*(2*m1*S1z
                      + kappa2*m2*S2z)))/pow(m1 + m2,2)); */

    x_3_5pn_SS =
        pre_factor * ((M_PI * ((1 + 160 * kappa1) * pow(m1, 2) * pow(S1z, 2) +
                               318 * m1 * m2 * S1z * S2z +
                               (1 + 160 * kappa2) * pow(m2, 2) * pow(S2z, 2))) /
                      (8. * pow(m1 + m2, 2)));
  }
  return (x_3_5pn_SS);
}

static REAL8 x_dot_3_5pn_cubicSpin(REAL8 e, REAL8 eta, REAL8 m1, REAL8 m2,
                                   REAL8 S1z, REAL8 S2z) /* Computed using
inputs from energy and flux expressions. See Blanchet liv. rev.
+ Bohe et al. arXiv: 1501.01529 + Marsat et al. arXiv: 1411.4118*/
{
  REAL8 x_3_5pn_cubicSpin;
  REAL8 kappa1 = 1.0; /*for black holes kappa_{1,2} is 1*/
  REAL8 kappa2 = 1.0;
  REAL8 lambda1 = 1.0; /*for black holes lambda_{1,2} is 1*/
  REAL8 lambda2 = 1.0;
  REAL8 pre_factor = 64. * eta / 5;

  if (e) {

    x_3_5pn_cubicSpin = 0.0;

  } else {

    x_3_5pn_cubicSpin =
        pre_factor *
        (-0.020833333333333332 *
         (2 * (15 + 1016 * kappa1 + 528 * lambda1) * pow(m1, 4) * pow(S1z, 3) +
          2 * (15 + 1016 * kappa2 + 528 * lambda2) * pow(m2, 4) * pow(S2z, 3) +
          pow(m1, 3) * m2 * pow(S1z, 2) *
              ((21 + 536 * kappa1 + 1056 * lambda1) * S1z +
               (4033 + 3712 * kappa1) * S2z) +
          12 * pow(m1, 2) * pow(m2, 2) * S1z * S2z *
              ((89 + 434 * kappa1) * S1z + (89 + 434 * kappa2) * S2z) +
          m1 * pow(m2, 3) * pow(S2z, 2) *
              ((4033 + 3712 * kappa2) * S1z +
               (21 + 536 * kappa2 + 1056 * lambda2) * S2z)) /
         pow(m1 + m2, 4));
  }

  return (x_3_5pn_cubicSpin);
}

static REAL8 x_dot_3_5_pn(REAL8 e, REAL8 eta) /* See Huerta et al article */
{
  REAL8 x_3_5_pn;
  REAL8 e_0_lim =
      64. * eta * M_PI *
      (-4415. / 4032. + 358675. * eta / 6048. + 91495. * eta * eta / 1512.) /
      5.;

  if (e) {
    x_3_5_pn = e_0_lim;
  } else {
    x_3_5_pn = e_0_lim;
  }

  return (x_3_5_pn);
}

/* The 4PN and 4.5PN pieces of x_dot are computed using inputs from Blanchet et
al, arXiv:2304.11185. */

static REAL8 x_dot_4pn(REAL8 e, REAL8 eta, REAL8 x) {
  REAL8 x_4_pn;
  REAL8 pre_factor = 64. * eta / 5;
  double EulerGamma = 0.5772156649015329;

  if (e) {
    x_4_pn = 0;
  } else {
    x_4_pn =
        pre_factor *
        ((3959271176713 - 20643291551545 * eta) / 2.54270016e10 +
         (pow(eta, 2) * (2016887396 + 21 * eta * (-1909807 + 49518 * eta))) /
             1.306368e6 -
         (896 * eta * EulerGamma) / 3. +
         ((124741 + 734620 * eta) * EulerGamma) / 4410. -
         (361 * pow(M_PI, 2)) / 126. -
         (eta * (1472377 + 928158 * eta) * pow(M_PI, 2)) / 16128. +
         (127751 * log(2)) / 1470. - (47385 * log(3)) / 1568. +
         eta * ((-850042 * log(2)) / 2205. + (47385 * log(3)) / 392.) +
         ((124741 - 582500 * eta) * log(x)) / 8820.);
  }
  return (x_4_pn);
}

static REAL8 x_dot_4pnSO(REAL8 e, REAL8 eta, REAL8 m1, REAL8 m2, REAL8 S1z,
                         REAL8 S2z) {
  REAL8 x_4pn_SO;
  REAL8 pre_factor = 64. * eta / 5;

  if (e) {
    x_4pn_SO = 0;
  } else {
    x_4pn_SO = pre_factor *
               (-0.000496031746031746 *
                (M_PI * (307708 * pow(m1, 4) * S1z + 307708 * pow(m2, 4) * S2z +
                         93121 * pow(m1, 2) * pow(m2, 2) * (S1z + S2z) +
                         m1 * pow(m2, 3) * (119880 * S1z + 75131 * S2z) +
                         pow(m1, 3) * m2 * (75131 * S1z + 119880 * S2z))) /
                pow(m1 + m2, 4));
  }
  return (x_4pn_SO);
}

static REAL8 x_dot_4pnSS(REAL8 e, REAL8 eta, REAL8 m1, REAL8 m2, REAL8 S1z,
                         REAL8 S2z) {
  REAL8 x_4pn_SS;
  REAL8 pre_factor = 64. * eta / 5;
  REAL8 kappa1 = 1.0; /*for black holes kappa_{1,2} is 1*/
  REAL8 kappa2 = 1.0;

  if (e) {
    x_4pn_SS = 0;
  } else {
    x_4pn_SS =
        pre_factor *
        (((41400957 + 10676336 * kappa1) * pow(m1, 6) * pow(S1z, 2) +
          (41400957 + 10676336 * kappa2) * pow(m2, 6) * pow(S2z, 2) +
          pow(m1, 5) * m2 * S1z *
              (5 * (16862889 + 4890568 * kappa1) * S1z + 53648974 * S2z) +
          m1 * pow(m2, 5) * S2z *
              (53648974 * S1z + 5 * (16862889 + 4890568 * kappa2) * S2z) +
          pow(m1, 4) * pow(m2, 2) *
              ((82037757 + 30833184 * kappa1) * pow(S1z, 2) +
               168293278 * S1z * S2z +
               (8950581 + 4778168 * kappa2) * pow(S2z, 2)) +
          pow(m1, 2) * pow(m2, 4) *
              ((8950581 + 4778168 * kappa1) * pow(S1z, 2) +
               168293278 * S1z * S2z +
               3 * (27345919 + 10277728 * kappa2) * pow(S2z, 2)) +
          pow(m1, 3) * pow(m2, 3) *
              ((43557309 + 18675776 * kappa1) * pow(S1z, 2) +
               226571670 * S1z * S2z +
               (43557309 + 18675776 * kappa2) * pow(S2z, 2))) /
         (72576. * pow(m1 + m2, 6)));
  }
  return (x_4pn_SS);
}

static REAL8 x_dot_4_5_pn(REAL8 e, REAL8 eta, REAL8 x) {
  REAL8 x_4_5_pn;
  REAL8 pre_factor = 64. * eta / 5;
  double EulerGamma = 0.5772156649015329;

  if (e) {
    x_4_5_pn = 0;
  } else {
    x_4_5_pn =
        pre_factor *
        ((451 * eta * pow(M_PI, 3)) / 12. -
         (M_PI *
          (700 * eta * (3098001198 + eta * (525268513 + 289286988 * eta)) +
           145786798080 * EulerGamma +
           3 * (-343801320119 + 97191198720 * log(2)))) /
             2.2353408e9 -
         (3424 * M_PI * log(x)) / 105.);
  }
  return (x_4_5_pn);
}

static REAL8 e_dot_0pn(REAL8 e, REAL8 eta) /* Eq. (A31) */
{
  REAL8 e_0_pn;
  REAL8 e_pow_2 = e * e;
  REAL8 e_factor = 1. - e_pow_2;
  REAL8 num = -e * eta * (121. * e_pow_2 + 304.);
  REAL8 den = 15. * e_factor * e_factor * sqrt(e_factor);

  if (e) {
    e_0_pn = num / den;
  } else {
    e_0_pn = 0.0;
  }

  return (e_0_pn);
}

static REAL8 e_dot_1pn(REAL8 e, REAL8 eta) /* Eq. (A32) */
{
  REAL8 e_1_pn;
  REAL8 e_pow_2 = e * e;
  REAL8 e_pow_4 = e_pow_2 * e_pow_2;
  REAL8 e_factor = 1.0 - e_pow_2;
  REAL8 e_0_term = 8. * (28588. * eta + 8451.);
  REAL8 e_2_term = 12. * (54271. * eta - 59834.) * e_pow_2;
  REAL8 e_4_term = (93184. * eta - 125361.) * e_pow_4;
  REAL8 pre_factor =
      e * eta / (2520. * e_factor * e_factor * e_factor * sqrt(e_factor));

  if (e) {
    e_1_pn = pre_factor * (e_0_term + e_2_term + e_4_term);
  } else {
    e_1_pn = 0.0;
  }

  return (e_1_pn);
}

static REAL8 e_dot_1_5pn_SO(REAL8 e, REAL8 m1, REAL8 m2, REAL8 S1z,
                            REAL8 S2z) /*computed using inputs from klein et al.
                                          arXiv:1801.08542. See equation C1b.*/
{
  REAL8 e_1_5pn_SO;
  REAL8 e_pow_2 = e * e;
  REAL8 e_pow_4 = e_pow_2 * e_pow_2;
  REAL8 e_factor = 1.0 - e_pow_2;
  REAL8 e_prefactor = e / (e_factor * e_factor * e_factor * e_factor);
  REAL8 pre_factor =
      e_prefactor *
      ((m1 * m2) / (90 * (m1 + m2) * (m1 + m2) * (m1 + m2) * (m1 + m2)));

  if (e) {

    e_1_5pn_SO =
        pre_factor *
        ((19688 + 28256 * e_pow_2 + 2367 * e_pow_4) * m1 * m1 * S1z +
         (19688 + 28256 * e_pow_2 + 2367 * e_pow_4) * m2 * m2 * S2z +
         3 * (4344 + 8090 * e_pow_2 + 835 * e_pow_4) * m1 * m2 * (S1z + S2z));
  } else {
    e_1_5pn_SO = 0.0;
  }

  return (e_1_5pn_SO);
}

static REAL8 e_dot_2pn_SS(REAL8 e, REAL8 m1, REAL8 m2, REAL8 S1z,
                          REAL8 S2z) /*computed using inputs from klein et al.
                                        arXiv:1801.08542. See equation C1b.*/
{
  REAL8 e_2pn_SS;
  REAL8 kappa1 = 1.0;
  REAL8 kappa2 = 1.0;
  REAL8 e_2 = e * e;
  REAL8 e_4 = e_2 * e_2;
  REAL8 e_fact = (1 - e_2);
  REAL8 e_fact_45 = e_fact * e_fact * e_fact * e_fact * sqrt(e_fact);
  REAL8 M_fact_2 = (m1 + m2) * (m1 + m2);
  REAL8 M_fact_4 = M_fact_2 * M_fact_2;

  if (e) {
    /* Quentin Henry et al terms, arXiv:2308.13606v1 */
    e_2pn_SS = -0.008333333333333333 *
               (e * m1 * m2 *
                ((45 * (8 + 12 * e_2 + e_4) +
                  4 * (3752 + 5950 * e_2 + 555 * e_4) * kappa1) *
                     m1 * m1 * S1z * S1z +
                 2 * (14648 + 23260 * e_2 + 2175 * e_4) * m1 * m2 * S1z * S2z +
                 (45 * (8 + 12 * e_2 + e_4) +
                  4 * (3752 + 5950 * e_2 + 555 * e_4) * kappa2) *
                     m2 * m2 * S2z * S2z)) /
               (e_fact_45 * M_fact_4);
    /* Klein et al expressions */ /* -0.008333333333333333 *
     (e * m1 * m2 *
      ((8 * (45 + 1876 * kappa1) +
        5 * pow(e, 2) *
            (108 + 4760 * kappa1 + 3 * pow(e, 2) * (3 + 37 * kappa1))) *
           pow(m1, 2) * pow(S1z, 2) +
       2 * (14648 + 23260 * pow(e, 2) + 2175 * pow(e, 4)) * m1 * m2 * S1z *
           S2z +
       (8 * (45 + 1876 * kappa2) +
        5 * pow(e, 2) *
            (108 + 4760 * kappa2 + 3 * pow(e, 2) * (3 + 37 * kappa2))) *
           pow(m2, 2) * pow(S2z, 2))) /
     (pow(1 - pow(e, 2), 4.5) * pow(m1 + m2, 4)); */
  } else {
    e_2pn_SS = 0.0;
  }

  return (e_2pn_SS);
}

static REAL8 e_rad_hereditary_1_5(REAL8 e, REAL8 eta, REAL8 x) {
  REAL8 e_her_l_o;
  REAL8 x_pow_3 = x * x * x;
  REAL8 x_pow_4 = x_pow_3 * x;
  REAL8 x_pow_3_2 = x * sqrt(x);

  REAL8 pre_fact = 32. * eta * e * x_pow_4 * x_pow_3_2 / 5.;
  REAL8 a_1 = -985. * M_PI * phi_e_rad(e) / 48.;

  if (e) {
    e_her_l_o = pre_fact * a_1;
  } else {
    e_her_l_o = 0.0;
  }

  return (e_her_l_o);
}

static REAL8 e_rad_hereditary_2_5(REAL8 e, REAL8 eta, REAL8 x) {
  REAL8 e_her_2_5;
  REAL8 x_pow_3 = x * x * x;
  REAL8 x_pow_4 = x_pow_3 * x;
  REAL8 x_pow_3_2 = x * sqrt(x);
  REAL8 x_pow_5_2 = x_pow_3_2 * x;

  REAL8 pre_fact = 32. * eta * e * x_pow_4 * x_pow_5_2 / 5.;
  REAL8 a_2 =
      (55691. * psi_e_rad(e) / 1344. + 19067. * eta * zed_e_rad(e) / 126.) *
      M_PI;

  if (e) {
    e_her_2_5 = pre_fact * a_2;
    ;
  } else {
    e_her_2_5 = 0.0;
  }

  return (e_her_2_5);
}

static REAL8 e_rad_hereditary_3(REAL8 e, REAL8 eta, REAL8 x) {
  const REAL8 log2 = 0.693147180559945309417232121458; // ln(2)
  const REAL8 log3 = 1.09861228866810969139524523692;  // ln(3)
  REAL8 e_3_pn_her;
  REAL8 euler = LAL_GAMMA;
  REAL8 pi_pow_2 = M_PI * M_PI;
  REAL8 x_pow_3_2 = x * sqrt(x);
  REAL8 x_pow_3 = x * x * x;
  REAL8 x_pow_4 = x_pow_3 * x;
  REAL8 x_pow_7 = x_pow_3 * x_pow_4;
  REAL8 pre_fact = 32. * eta * e * x_pow_7 / 5.;
  REAL8 a_3 =
      (89789209. / 352800. - 87419. * log2 / 630. + 78003. * log3 / 560.) *
      kappa_e_rad(e);
  REAL8 a_4 = (-769. / 96.) *
              (16. * pi_pow_2 / 3. - 1712. * euler / 105. -
               1712. * log(4. * x_pow_3_2) / 105.) *
              capital_f_e(e);

  if (e) {
    e_3_pn_her = pre_fact * (a_3 + a_4);
  } else {
    e_3_pn_her = 0.0;
  }

  return (e_3_pn_her);
}

static REAL8 e_dot_2pn(REAL8 e, REAL8 eta) /* Eq. (A34) */
{
  REAL8 e_2_pn;
  REAL8 eta_pow_2 = eta * eta;
  REAL8 e_pow_2 = e * e;
  REAL8 e_pow_4 = e_pow_2 * e_pow_2;
  REAL8 e_pow_6 = e_pow_4 * e_pow_2;
  REAL8 e_fact = 1. - e_pow_2;
  REAL8 zero_term =
      -952397. / 1890. + 5937. * eta / 14. + 752. * eta_pow_2 / 5.;
  REAL8 e_2_term = e_pow_2 * (-3113989. / 2520. - 388419. * eta / 280. +
                              64433. * eta_pow_2 / 40.);
  REAL8 e_4_term = e_pow_4 * (4656611. / 3024. - 13057267. * eta / 5040. +
                              127411. * eta_pow_2 / 90.);
  REAL8 e_6_term = e_pow_6 * (420727. / 3360. - 362071. * eta / 2520. +
                              821. * eta_pow_2 / 9.);
  REAL8 zero_rt = 1336. / 3. - 2672. * eta / 15.;
  REAL8 e_2_rt = e_pow_2 * (2321. / 2. - 2321. * eta / 5.);
  REAL8 e_4_rt = e_pow_4 * (565. / 6. - 113. * eta / 3.);
  REAL8 pre_factor =
      -e * eta / (e_fact * e_fact * e_fact * e_fact * sqrt(e_fact));

  if (e) {
    e_2_pn = pre_factor * (zero_term + e_2_term + e_4_term + e_6_term +
                           sqrt(e_fact) * (zero_rt + e_2_rt + e_4_rt));
  } else {
    e_2_pn = 0.0;
  }

  return (e_2_pn);
}

static REAL8
e_dot_2_5pn_SO(REAL8 e, REAL8 m1, REAL8 m2, REAL8 S1z,
               REAL8 S2z) /* Quentin Henry et al terms, arXiv:2308.13606v1 */
{
  REAL8 e_2_5pn_SO;
  REAL8 e_2 = e * e;
  REAL8 e_4 = e_2 * e_2;
  REAL8 e_6 = e_4 * e_2;
  REAL8 e_fact = (1 - e_2);
  REAL8 e_fact_2 = e_fact * e_fact;
  REAL8 e_fact_sqrt = sqrt(e_fact);
  REAL8 e_fact_55 = e_fact_2 * e_fact_2 * e_fact * e_fact_sqrt;
  REAL8 M_fact_2 = (m1 + m2) * (m1 + m2);
  REAL8 M_fact_4 = M_fact_2 * M_fact_2;
  REAL8 M_fact_6 = M_fact_4 * M_fact_2;

  if (e) {
    e_2_5pn_SO =
        ((e * m1 * m2 *
          (3 *
               (192 * (82880 + 36083 * e_fact_sqrt) +
                168 * e_4 * (-211648 + 487675 * e_fact_sqrt) +
                3 * e_6 * (-560896 + 1037433 * e_fact_sqrt) +
                16 * e_2 * (1332912 + 6885449 * e_fact_sqrt)) *
               m1 * m1 * m1 * m1 * S1z +
           3 *
               (192 * (82880 + 36083 * e_fact_sqrt) +
                168 * e_4 * (-211648 + 487675 * e_fact_sqrt) +
                3 * e_6 * (-560896 + 1037433 * e_fact_sqrt) +
                16 * e_2 * (1332912 + 6885449 * e_fact_sqrt)) *
               m2 * m2 * m2 * m2 * S2z +
           3 *
               (27847680 - 9476416 * e_fact_sqrt +
                7 * e_6 * (-420672 + 929363 * e_fact_sqrt) +
                24 * e_4 * (-2592688 + 6508535 * e_fact_sqrt) +
                16 * e_2 * (2332596 + 9517267 * e_fact_sqrt)) *
               m1 * m1 * m2 * m2 * (S1z + S2z) +
           m1 * m1 * m1 * m2 *
               ((128 * (808080 - 259453 * e_fact_sqrt) +
                 3 * e_6 * (-3645824 + 6571731 * e_fact_sqrt) +
                 48 * e_2 * (2887976 + 10533923 * e_fact_sqrt) +
                 32 * e_4 * (-7222488 + 15232045 * e_fact_sqrt)) *
                    S1z +
                9 *
                    (206464 * e_fact_sqrt + 21830032 * e_2 * e_fact_sqrt +
                     22413824 * e_4 * e_fact_sqrt +
                     1071519 * e_6 * e_fact_sqrt -
                     896 * (-1 + e) * (1 + e) *
                         (2960 + 6927 * e_2 + 313 * e_4)) *
                    S2z) +
           m1 * m2 * m2 * m2 *
               (9 *
                    (128 * (20720 + 1613 * e_fact_sqrt) +
                     256 * e_4 * (-23149 + 87554 * e_fact_sqrt) +
                     112 * e_2 * (31736 + 194911 * e_fact_sqrt) +
                     e_6 * (-280448 + 1071519 * e_fact_sqrt)) *
                    S1z +
                (128 * (808080 - 259453 * e_fact_sqrt) +
                 3 * e_6 * (-3645824 + 6571731 * e_fact_sqrt) +
                 48 * e_2 * (2887976 + 10533923 * e_fact_sqrt) +
                 32 * e_4 * (-7222488 + 15232045 * e_fact_sqrt)) *
                    S2z))) /
         (60480. * e_fact_55 * M_fact_6));
  } else {
    e_2_5pn_SO = 0.0;
  }
  return (e_2_5pn_SO);
}

static REAL8 e_dot_3pn(REAL8 e, REAL8 eta, REAL8 x) /* Eq. (C10) */
{
  REAL8 e_3_pn;
  REAL8 eta_pow_2 = eta * eta;
  REAL8 eta_pow_3 = eta_pow_2 * eta;
  REAL8 e_pow_2 = e * e;
  REAL8 e_pow_4 = e_pow_2 * e_pow_2;
  REAL8 e_pow_6 = e_pow_4 * e_pow_2;
  REAL8 e_pow_8 = e_pow_6 * e_pow_2;
  REAL8 e_fact = 1. - e_pow_2;
  REAL8 pi_pow_2 = M_PI * M_PI;

  REAL8 zero_term = (7742634967. / 891000.) +
                    ((43386337. / 113400.) + (1017. * pi_pow_2 / 10.)) * eta -
                    (4148897. * eta_pow_2 / 2520.) -
                    (61001. * eta_pow_3 / 486.);
  REAL8 e_2_term =
      e_pow_2 *
      ((6556829759. / 891000.) +
       ((770214901. / 25200.) - (15727. * pi_pow_2 / 192.)) * eta -
       (80915371. * eta_pow_2 / 15120.) - (86910509. * eta_pow_3 / 19440.));
  REAL8 e_4_term =
      e_pow_4 *
      (-(17072216761. / 2376000.) +
       ((8799500893. / 907200.) - (295559. * pi_pow_2 / 1920.)) * eta +
       (351962207. * eta_pow_2 / 20160.) - (2223241. * eta_pow_3 / 180.));
  REAL8 e_6_term =
      e_pow_6 *
      ((17657772379. / 3696000.) +
       (-(91818931. / 10080.) - (6519. * pi_pow_2 / 640.)) * eta +
       (2495471. * eta_pow_2 / 252.) - (11792069. * eta_pow_3 / 2430.));
  REAL8 e_8_term =
      e_pow_8 * ((302322169. / 1774080.) - (1921387. * eta / 10080.) +
                 (41179. * eta_pow_2 / 216.) - (193396. * eta_pow_3 / 1215.));
  REAL8 zero_rt = -(22713049. / 15750.) +
                  (-(5526991. / 945.) + (8323. * pi_pow_2 / 180.)) * eta +
                  (54332. * eta_pow_2 / 45.);
  REAL8 e_2_rt =
      e_pow_2 * ((89395687. / 7875.) +
                 (-(38295557. / 1260.) + (94177. * pi_pow_2 / 960.)) * eta +
                 (681989. * eta_pow_2 / 90.));
  REAL8 e_4_rt =
      e_pow_4 * ((5321445613. / 378000.) +
                 (-(26478311. / 1512.) + (2501. * pi_pow_2 / 2880.)) * eta +
                 (225106. * eta_pow_2 / 45.));
  REAL8 e_6_rt = e_pow_6 * ((186961. / 336.) - (289691. * eta / 504) +
                            (3197. * eta_pow_2 / 18.));
  REAL8 free_term = 730168. / (23625. * (1. + sqrt(e_fact)));
  REAL8 log_term =
      (304. / 15.) *
      (82283. / 1995. + (297674. * e_pow_2 / 1995.) +
       (1147147. * e_pow_4 / 15960.) + (61311. * e_pow_6 / 21280.)) *
      log(x * (1. + sqrt(e_fact)) / (2. * e_fact));
  REAL8 pre_factor =
      -e * eta / (e_fact * e_fact * e_fact * e_fact * e_fact * sqrt(e_fact));

  if (e) {
    e_3_pn =
        pre_factor * (zero_term + e_2_term + e_4_term + e_6_term + e_8_term +
                      sqrt(e_fact) * (zero_rt + e_2_rt + e_4_rt + e_6_rt) +
                      free_term + log_term);
  } else {
    e_3_pn = 0.0;
  }
  return (e_3_pn);
}

static REAL8
e_dot_3pn_SO(REAL8 e, REAL8 m1, REAL8 m2, REAL8 S1z,
             REAL8 S2z) /* Quentin Henry et al terms, arXiv:2308.13606v1 */
{
  REAL8 e_3pn_SO;
  REAL8 e_2 = e * e;
  REAL8 e_4 = e_2 * e_2;
  REAL8 e_6 = e_4 * e_2;
  REAL8 e_fact = (1 - e_2);
  REAL8 e_fact_2 = e_fact * e_fact;
  REAL8 e_fact_sqrt = sqrt(e_fact);
  REAL8 e_fact_55 = e_fact_2 * e_fact_2 * e_fact * e_fact_sqrt;
  REAL8 M_fact_2 = (m1 + m2) * (m1 + m2);
  REAL8 M_fact_4 = M_fact_2 * M_fact_2;

  if (e) {
    e_3pn_SO =
        ((e * m1 * m2 * M_PI *
          ((64622592 + 238783104 * e_2 + 96887280 * e_4 + 2313613 * e_6) * m1 *
               m1 * S1z +
           (64622592 + 238783104 * e_2 + 96887280 * e_4 + 2313613 * e_6) * m2 *
               m2 * S2z +
           24 * (1744704 + 8150400 * e_2 + 3941409 * e_4 + 122714 * e_6) * m1 *
               m2 * (S1z + S2z))) /
         (51840. * e_fact_55 * M_fact_4));
  } else {
    e_3pn_SO = 0.0;
  }
  return (e_3pn_SO);
}

static REAL8 e_dot_3pn_SS(REAL8 e, REAL8 m1, REAL8 m2, REAL8 S1z,
             REAL8 S2z) /* Quentin Henry et al terms, arXiv:2308.13606v1 */
{
  REAL8 e_3pn_SS;
  REAL8 kappa1 = 1.0;
  REAL8 kappa2 = 1.0;
  REAL8 e_2 = e * e;
  REAL8 e_4 = e_2 * e_2;
  REAL8 e_6 = e_4 * e_2;
  REAL8 e_fact = (1 - e_2);
  REAL8 e_fact_2 = e_fact * e_fact;
  REAL8 e_fact_sqrt = sqrt(e_fact);
  REAL8 e_fact_55 = e_fact_2 * e_fact_2 * e_fact * e_fact_sqrt;
  REAL8 M_fact_2 = (m1 + m2) * (m1 + m2);
  REAL8 M_fact_4 = M_fact_2 * M_fact_2;
  REAL8 M_fact_6 = M_fact_4 * M_fact_2;

  if (e) {
    e_3pn_SS =
        (-0.000016534391534391536 *
         (e * m1 * m2 *
          ((27 * e_6 * (53837 + 368940 * kappa1) +
            12 * e_4 *
                (6350925 + 27328 * e_fact_sqrt + 16634634 * kappa1 +
                 47824 * e_fact_sqrt * kappa1) +
            32 * (2226847 + 545664 * e_fact_sqrt + 538758 * kappa1 +
                  954912 * e_fact_sqrt * kappa1) +
            32 * e_2 *
                (7292761 + 1157688 * e_fact_sqrt +
                 9 * (847619 + 225106 * e_fact_sqrt) * kappa1)) *
               m1 * m1 * m1 * m1 * S1z * S1z +
           (27 * e_6 * (53837 + 368940 * kappa2) +
            12 * e_4 *
                (6350925 + 27328 * e_fact_sqrt + 16634634 * kappa2 +
                 47824 * e_fact_sqrt * kappa2) +
            32 * (2226847 + 545664 * e_fact_sqrt + 538758 * kappa2 +
                  954912 * e_fact_sqrt * kappa2) +
            32 * e_2 *
                (7292761 + 1157688 * e_fact_sqrt +
                 9 * (847619 + 225106 * e_fact_sqrt) * kappa2)) *
               m2 * m2 * m2 * m2 * S2z * S2z +
           6 * m1 * m1 * m1 * m2 * S1z *
               ((9 * e_6 * (55293 + 221219 * kappa1) +
                 2 * e_4 *
                     (11036963 + 10248 * e_fact_sqrt + 19572200 * kappa1 +
                      78568 * e_fact_sqrt * kappa1) +
                 16 * (798819 + 68208 * e_fact_sqrt - 336460 * kappa1 +
                       522928 * e_fact_sqrt * kappa1) +
                 8 * e_2 *
                     (7063231 + 289422 * e_fact_sqrt +
                      6 * (698816 + 369817 * e_fact_sqrt) * kappa1)) *
                    S1z +
                2 *
                    (1818549 * e_6 + 240 * (31249 + 22736 * e_fact_sqrt) +
                     2 * e_4 * (20863999 + 51240 * e_fact_sqrt) +
                     8 * e_2 * (7764719 + 1447110 * e_fact_sqrt)) *
                    S2z) +
           6 * m1 * m2 * m2 * m2 * S2z *
               (2 *
                    (1818549 * e_6 + 240 * (31249 + 22736 * e_fact_sqrt) +
                     2 * e_4 * (20863999 + 51240 * e_fact_sqrt) +
                     8 * e_2 * (7764719 + 1447110 * e_fact_sqrt)) *
                    S1z +
                (9 * e_6 * (55293 + 221219 * kappa2) +
                 2 * e_4 *
                     (11036963 + 10248 * e_fact_sqrt + 19572200 * kappa2 +
                      78568 * e_fact_sqrt * kappa2) +
                 16 * (798819 + 68208 * e_fact_sqrt - 336460 * kappa2 +
                       522928 * e_fact_sqrt * kappa2) +
                 8 * e_2 *
                     (7063231 + 289422 * e_fact_sqrt +
                      6 * (698816 + 369817 * e_fact_sqrt) * kappa2)) *
                    S2z) +
           m1 * m1 * m2 * m2 *
               (9 *
                    (3 * e_6 * (63301 + 361786 * kappa1) +
                     4 * e_4 *
                         (1667883 - 3416 * e_fact_sqrt + 5133138 * kappa1 +
                          13664 * e_fact_sqrt * kappa1) +
                     32 * (69895 - 22736 * e_fact_sqrt - 37046 * kappa1 +
                           90944 * e_fact_sqrt * kappa1) +
                     32 * e_2 *
                         (439124 - 48237 * e_fact_sqrt + 616472 * kappa1 +
                          192948 * e_fact_sqrt * kappa1)) *
                    S1z * S1z +
                8 *
                    (3553929 * e_6 +
                     3 * e_4 * (29583761 + 99064 * e_fact_sqrt) +
                     116 * e_2 * (1176325 + 289422 * e_fact_sqrt) +
                     8 * (2057131 + 1978032 * e_fact_sqrt)) *
                    S1z * S2z +
                9 *
                    (3 * e_6 * (63301 + 361786 * kappa2) +
                     4 * e_4 *
                         (1667883 - 3416 * e_fact_sqrt + 5133138 * kappa2 +
                          13664 * e_fact_sqrt * kappa2) +
                     32 * (69895 - 22736 * e_fact_sqrt - 37046 * kappa2 +
                           90944 * e_fact_sqrt * kappa2) +
                     32 * e_2 *
                         (439124 - 48237 * e_fact_sqrt + 616472 * kappa2 +
                          192948 * e_fact_sqrt * kappa2)) *
                    S2z * S2z))) /
         (e_fact_55 * M_fact_6));
  } else {
    e_3pn_SS = 0.0;
  }
  return (e_3pn_SS);
}

static REAL8 e_dot_3_5pn(REAL8 e, REAL8 eta) /* Eq. (C10) */
{
  (void)eta;
  REAL8 e_3_5_pn;

  if (e) {
    e_3_5_pn = 0.0;
  } else {
    e_3_5_pn = 0.0;
  }
  return (e_3_5_pn);
}

static REAL8 l_dot_1pn(REAL8 e, REAL8 eta) /* Eq. (A2) */
{
  (void)eta;
  return (3. / (e * e - 1.));
}

static REAL8
l_dot_1_5pn_SO(REAL8 e, REAL8 m1, REAL8 m2, REAL8 S1z,
               REAL8 S2z) /*computed using inputs from klein et al.
                             arXiv:1801.08542. See equation B1e and B2e. */
{
  REAL8 e_pow_2 = e * e;
  REAL8 e_factor = 1.0 - e_pow_2;
  REAL8 e_prefactor = 1 / (e_factor * sqrt(e_factor));
  REAL8 pre_factor = e_prefactor / ((m1 + m2) * (m1 + m2));

  return (pre_factor *
          (4 * m1 * m1 * S1z + 4 * m2 * m2 * S2z + 3 * m1 * m2 * (S1z + S2z)));
}

static REAL8
l_dot_2pn_SS(REAL8 e, REAL8 m1, REAL8 m2, REAL8 S1z,
             REAL8 S2z) /*computed using inputs from klein et al.
                           arXiv:1801.08542. See equation B1e and B2e. */
{
  REAL8 kappa1 = 1.0;
  REAL8 kappa2 = 1.0;

  return ((-3 * (m1 * S1z * (2 * m2 * S2z + m1 * S1z * kappa1) +
                 pow(m2, 2) * pow(S2z, 2) * kappa2)) /
          (2. * pow(-1 + pow(e, 2), 2) * pow(m1 + m2, 2)));
}

static REAL8
l_dot_2_5pn_SO(REAL8 e, REAL8 m1, REAL8 m2, REAL8 S1z,
               REAL8 S2z) /* Quentin Henry et al terms, arXiv:2308.13606v1 */
{
  REAL8 l_2_5pn_SO;
  REAL8 e_2 = e * e;
  REAL8 e_4 = e_2 * e_2;
  REAL8 e_fact = (1 - e_2);
  REAL8 e_fact_2 = e_fact * e_fact;
  REAL8 e_fact_sqrt = sqrt(e_fact);
  REAL8 e_fact_25 = e_fact_2 * e_fact_sqrt;
  REAL8 M_fact_2 = (m1 + m2) * (m1 + m2);
  REAL8 M_fact_4 = M_fact_2 * M_fact_2;

  l_2_5pn_SO =
      ((20 * m1 * m1 * m1 * m1 * (S1z + 3 * e_2 * S1z) +
        20 * (1 + 3 * e_2) * m2 * m2 * m2 * m2 * S2z +
        (5 + 129 * e_2) * m1 * m1 * m2 * m2 * (S1z + S2z) +
        m1 * m1 * m1 * m2 * ((23 + 137 * e_2) * S1z + 45 * e_2 * S2z) +
        m1 * m2 * m2 * m2 * (23 * S2z + e_2 * (45 * S1z + 137 * S2z))) /
       (2. * e_fact_25 * M_fact_4));

  return (l_2_5pn_SO);
}

static REAL8
l_dot_3pn_SS(REAL8 e, REAL8 m1, REAL8 m2, REAL8 S1z,
             REAL8 S2z) /* Quentin Henry et al terms, arXiv:2308.13606v1 */
{
  REAL8 l_3pn_SS;
  REAL8 kappa1 = 1.0;
  REAL8 kappa2 = 1.0;
  REAL8 e_2 = e * e;
  REAL8 M_fact_2 = (m1 + m2) * (m1 + m2);
  REAL8 M_fact_4 = M_fact_2 * M_fact_2;

  l_3pn_SS =
      ((m1 * m1 *
            ((-8 + 42 * kappa1 + 6 * e_2 * (4 + 13 * kappa1)) * m1 * m1 +
             (-60 + 50 * kappa1 + 11 * e_2 * (3 + 11 * kappa1)) * m1 * m2 +
             3 * (-14 + 6 * kappa1 + 3 * e_2 * (1 + 8 * kappa1)) * m2 * m2) *
            S1z * S1z +
        2 * m1 * m2 *
            ((6 + 93 * e_2) * m1 * m1 + (12 + 157 * e_2) * m1 * m2 +
             3 * (2 + 31 * e_2) * m2 * m2) *
            S1z * S2z +
        m2 * m2 *
            (3 * (-14 + 6 * kappa2 + 3 * e_2 * (1 + 8 * kappa2)) * m1 * m1 +
             (-60 + 50 * kappa2 + 11 * e_2 * (3 + 11 * kappa2)) * m1 * m2 +
             2 * (-4 + 21 * kappa2 + 3 * e_2 * (4 + 13 * kappa2)) * m2 * m2) *
            S2z * S2z) /
       (4. * pow(-1 + e_2, 3) * M_fact_4));

  return (l_3pn_SS);
}

static REAL8 l_dot_2pn(REAL8 e, REAL8 eta) /* Eq. (A3) */
{
  REAL8 e_fact = 1 - e * e;
  REAL8 e_fact_pow_2 = e_fact * e_fact;
  return (((26. * eta - 51.) * e * e + 28. * eta - 18.) / (4. * e_fact_pow_2));
}

static REAL8 l_dot_3pn(REAL8 e, REAL8 eta) /* Eq. (A4) */
{
  REAL8 e_fact = 1 - e * e;
  REAL8 e_fact_pow_3 = e_fact * e_fact * e_fact;
  REAL8 pre_factor = -1. / (128. * sqrt(e_fact) * e_fact_pow_3);
  REAL8 eta_pow_2 = eta * eta;
  REAL8 pi_pow_2 = M_PI * M_PI;
  REAL8 e_pow_2 = e * e;
  REAL8 e_pow_4 = e_pow_2 * e_pow_2;
  REAL8 e_0_term = 1920. - 768. * eta;
  REAL8 e_2_term = (1920. - 768. * eta) * e_pow_2;
  REAL8 e_4_term = (1536. * eta - 3840.) * e_pow_4;
  REAL8 e_0_rt = 896. * eta_pow_2 - 14624. * eta + 492. * pi_pow_2 * eta - 192.;
  REAL8 e_2_rt =
      (5120. * eta_pow_2 + 123. * pi_pow_2 * eta - 17856. * eta + 8544.) *
      e_pow_2;
  REAL8 e_4_rt = (1040. * eta_pow_2 - 1760. * eta + 2496.) * e_pow_4;

  return (pre_factor * (e_0_term + e_2_term + e_4_term +
                        sqrt(e_fact) * (e_0_rt + e_2_rt + e_4_rt)));
}

static REAL8 phi_dot_0pn(REAL8 e, REAL8 eta, REAL8 u) /* Eq. (A11) */
{
  (void)eta;
  return sqrt(1.0 - e * e) / (cosu_factor(e, u) * cosu_factor(e, u));
}

static REAL8 phi_dot_1pn(REAL8 e, REAL8 eta, REAL8 u) /* Eq. (A12) */
{
  REAL8 u_factor = cosu_factor(e, u);

  return -(e * (eta - 4.0) * (e - cos(u)) /
           (sqrt(1.0 - e * e) * u_factor * u_factor * u_factor));
}

static REAL8
phi_dot_1_5_pnSO_ecc(REAL8 e, REAL8 m1, REAL8 m2, REAL8 S1z, REAL8 S2z,
                     REAL8 u) /*computed using inputs from klein et al.
                                 arXiv:1801.08542. See equation B1 and B2. */
{
  if (e) {
    return ((2 * e * (m1 * S1z + m2 * S2z) * (e - cos(u))) /
            ((-1 + e * e) * (m1 + m2) *
             ((-1 + e * cos(u)) * (-1 + e * cos(u)) * (-1 + e * cos(u)))));
  } else {

    return (0.);
  }
}

static REAL8 phi_dot_2pn(REAL8 e, REAL8 eta, REAL8 u) /* Eq. (A13) */
{
  REAL8 u_factor = cosu_factor(e, u);
  REAL8 u_factor_pow_5 = u_factor * u_factor * u_factor * u_factor * u_factor;
  REAL8 e_pow_2 = e * e;
  REAL8 e_pow_3 = e_pow_2 * e;
  REAL8 e_pow_4 = e_pow_2 * e_pow_2;
  REAL8 e_pow_5 = e_pow_4 * e;
  REAL8 e_pow_6 = e_pow_4 * e_pow_2;
  REAL8 e_factor = 1.0 - e_pow_2;
  REAL8 eta_pow_2 = eta * eta;
  REAL8 cos_u_pow_2 = cos(u) * cos(u);
  REAL8 cos_u_pow_3 = cos_u_pow_2 * cos(u);

  REAL8 e_0_term = 90. - 36. * eta;
  REAL8 e_2_term = (-2. * eta_pow_2 + 50. * eta + 75.) * e_pow_2;
  REAL8 e_4_term = (20. * eta_pow_2 - 26. * eta - 60.) * e_pow_4;
  REAL8 e_6_term = (-12. * eta - 18.) * eta * e_pow_6;

  REAL8 cosu_1 = ((-eta_pow_2 + 97. * eta + 12.) * e_pow_5 +
                  (-16. * eta_pow_2 - 74. * eta - 81.) * e_pow_3 +
                  (-eta_pow_2 + 67. * eta - 246.) * e) *
                 cos(u);
  REAL8 cosu_2 = ((17. * eta_pow_2 - 17. * eta + 48.) * e_pow_6 +
                  (-4. * eta_pow_2 - 38. * eta + 153.) * e_pow_4 +
                  (5. * eta_pow_2 - 35. * eta + 114.) * e_pow_2) *
                 cos_u_pow_2;
  REAL8 cosu_3 = ((-14. * eta_pow_2 + 8. * eta - 147.) * e_pow_5 +
                  (8. * eta_pow_2 + 22. * eta + 42.) * e_pow_3) *
                 cos_u_pow_3;

  REAL8 rt_zero = (180. - 72. * eta) * e_pow_2 + 36. * eta - 90.;
  REAL8 rt_cosu_1 =
      ((144. * eta - 360.) * e_pow_3 + (90. - 36. * eta) * e) * cos(u);
  REAL8 rt_cosu_2 =
      ((180. - 72. * eta) * e_pow_4 + (90. - 36. * eta) * e_pow_2) *
      cos_u_pow_2;
  REAL8 rt_cosu_3 = e_pow_3 * (36. * eta - 90.) * cos_u_pow_3;

  REAL8 pre_factor = 1. / (12. * sqrt(e_factor) * e_factor * u_factor_pow_5);

  return (pre_factor *
          (e_0_term + e_2_term + e_4_term + e_6_term + cosu_1 + cosu_2 +
           cosu_3 +
           sqrt(e_factor) * (rt_zero + rt_cosu_1 + rt_cosu_2 + rt_cosu_3)));
}

static REAL8
phi_dot_2_pnSS_ecc(REAL8 e, REAL8 m1, REAL8 m2, REAL8 S1z, REAL8 S2z,
                   REAL8 u) /*computed using inputs from klein et al.
                               arXiv:1801.08542. See equation B1 and B2. */
{
  REAL8 kappa1 = 1.0;
  REAL8 kappa2 = 1.0;

  if (e) {

    return (
        (e *
         (kappa2 * pow(m2, 2) * pow(S2z, 2) +
          m1 * S1z * (kappa1 * m1 * S1z + 2 * m2 * S2z)) *
         (e - cos(u))) /
        (pow(1 - pow(e, 2), 1.5) * pow(m1 + m2, 2) * pow(-1 + e * cos(u), 3)));
  } else {

    return (0.);
  }
}

static REAL8
phi_dot_2_5pn_SO(REAL8 e, REAL8 m1, REAL8 m2, REAL8 S1z, REAL8 S2z,
                 REAL8 u) /* Quentin Henry et al terms, arXiv:2308.13606v1 */
{
  REAL8 phi_2_5pn_SO;
  REAL8 e_2 = e * e;
  REAL8 e_3 = e_2 * e;
  REAL8 e_4 = e_2 * e_2;
  REAL8 e_6 = e_4 * e_2;
  REAL8 e_fact = (1 - e_2);
  REAL8 e_fact_sqrt = sqrt(e_fact);
  REAL8 M_fact_2 = (m1 + m2) * (m1 + m2);
  REAL8 M_fact_4 = M_fact_2 * M_fact_2;

  phi_2_5pn_SO =
      (((m1 - m2) * (m1 + m2) *
            (24 * (m1 * m2 - 3 * M_fact_2) -
             6 * e_6 * (m1 * m2 - 2 * M_fact_2) +
             18 * e_4 * (3 * m1 * m2 - 2 * M_fact_2) -
             e_2 * (31 * m1 * m2 + 56 * M_fact_2)) *
            (S1z - S2z) +
        (e_2 * (50 * m1 * m1 * m2 * m2 - 57 * m1 * m2 * M_fact_2 -
                56 * M_fact_4) +
         6 * e_6 *
             (2 * m1 * m1 * m2 * m2 - 3 * m1 * m2 * M_fact_2 + 2 * M_fact_4) -
         6 * e_4 *
             (8 * m1 * m1 * m2 * m2 - 27 * m1 * m2 * M_fact_2 + 6 * M_fact_4) -
         12 * (m1 * m1 * m2 * m2 - 8 * m1 * m2 * M_fact_2 + 6 * M_fact_4)) *
            (S1z + S2z) -
        e *
            (-8 * (56 + 49 * e_2 + 9 * e_4) * m1 * m1 * m1 * m1 * S1z -
             8 * (56 + 49 * e_2 + 9 * e_4) * m2 * m2 * m2 * m2 * S2z +
             4 * (-217 - 200 * e_2 + 9 * e_4) * m1 * m1 * m2 * m2 *
                 (S1z + S2z) -
             2 * m1 * m1 * m1 * m2 *
                 ((530 + 496 * e_2 + 6 * e_4) * S1z +
                  9 * (14 + 13 * e_2) * S2z) -
             2 * m1 * m2 * m2 * m2 *
                 (9 * (14 + 13 * e_2) * S1z +
                  2 * (265 + 248 * e_2 + 3 * e_4) * S2z)) *
            cos(u) +
        e_2 *
            (-8 * (2 + e_2) * (29 + 9 * e_2) * m1 * m1 * m1 * m1 * S1z -
             8 * (2 + e_2) * (29 + 9 * e_2) * m2 * m2 * m2 * m2 * S2z -
             8 * (2 + e_2) * (47 + 21 * e_2) * m1 * m1 * m2 * m2 * (S1z + S2z) -
             2 * m1 * m1 * m1 * m2 *
                 ((514 + 440 * e_2 + 78 * e_4) * S1z +
                  9 * (2 + e_2) * (5 + 4 * e_2) * S2z) -
             2 * m1 * m2 * m2 * m2 *
                 (9 * (2 + e_2) * (5 + 4 * e_2) * S1z +
                  2 * (257 + 220 * e_2 + 39 * e_4) * S2z)) *
            pow(cos(u), 2) -
        e_3 *
            (-8 * (17 + 21 * e_2) * m1 * m1 * m1 * m1 * S1z -
             8 * (17 + 21 * e_2) * m2 * m2 * m2 * m2 * S2z -
             8 * (11 + 57 * e_2) * m1 * m1 * m2 * m2 * (S1z + S2z) -
             2 * m1 * m1 * m1 * m2 *
                 (4 * (29 + 57 * e_2) * S1z - 6 * S2z + 87 * e_2 * S2z) -
             2 * m1 * m2 * m2 * m2 *
                 ((-6 + 87 * e_2) * S1z + 4 * (29 + 57 * e_2) * S2z)) *
            pow(cos(u), 3) -
        12 * e_fact_sqrt *
            (-12 * m1 * m1 * m1 * m1 * S1z - 12 * m2 * m2 * m2 * m2 * S2z -
             21 * m1 * m1 * m2 * m2 * (S1z + S2z) -
             2 * m1 * m1 * m1 * m2 * (13 * S1z + 3 * S2z) -
             2 * m1 * m2 * m2 * m2 * (3 * S1z + 13 * S2z)) *
            pow(-1 + e * cos(u), 2) * (1 - 2 * e_2 + e * cos(u))) /
       (12. * pow(-1 + e_2, 2) * M_fact_4 * pow(-1 + e * cos(u), 5)));

  return (phi_2_5pn_SO);
}

static REAL8 phi_dot_3pn(REAL8 e, REAL8 eta, REAL8 u) {
  REAL8 u_factor = cosu_factor(e, u);
  REAL8 u_factor_pow_5 = u_factor * u_factor * u_factor * u_factor * u_factor;
  REAL8 u_factor_pow_7 = u_factor_pow_5 * u_factor * u_factor;
  REAL8 pi_pow_2 = M_PI * M_PI;
  REAL8 eta_pow_2 = eta * eta;
  REAL8 eta_pow_3 = eta_pow_2 * eta;
  REAL8 e_pow_2 = e * e;
  REAL8 e_fact = 1.0 - e_pow_2;
  REAL8 e_pow_3 = e_pow_2 * e;
  REAL8 e_pow_4 = e_pow_2 * e_pow_2;
  REAL8 e_pow_5 = e_pow_4 * e;
  REAL8 e_pow_6 = e_pow_4 * e_pow_2;
  REAL8 e_pow_7 = e_pow_6 * e;
  REAL8 e_pow_8 = e_pow_7 * e;
  REAL8 e_pow_9 = e_pow_8 * e;
  REAL8 e_pow_10 = e_pow_8 * e_pow_2;
  REAL8 e_factor = 1.0 - e_pow_2;
  REAL8 cos_u_pow_2 = cos(u) * cos(u);
  REAL8 cos_u_pow_3 = cos_u_pow_2 * cos(u);
  REAL8 cos_u_pow_4 = cos_u_pow_2 * cos_u_pow_2;
  REAL8 cos_u_pow_5 = cos_u_pow_4 * cos(u);
  REAL8 pre_factor =
      1. / (13440. * sqrt(e_factor) * e_factor * e_factor * u_factor_pow_7);

  REAL8 e_0_term =
      67200. * eta_pow_2 - 761600. * eta + 8610. * eta * pi_pow_2 + 201600.;
  REAL8 e_2_term = (4480. * eta_pow_3 - 412160. * eta_pow_2 -
                    30135. * pi_pow_2 * eta + 553008. * eta + 342720.) *
                   e_pow_2;
  REAL8 e_4_term = (-52640. * eta_pow_3 + 516880. * eta_pow_2 +
                    68880. * pi_pow_2 * eta - 1916048. * eta + 262080.) *
                   e_pow_4;
  REAL8 e_6_term = (84000. * eta_pow_3 - 190400. * eta_pow_2 -
                    17220. * pi_pow_2 * eta - 50048. * eta - 241920.) *
                   e_pow_6;
  REAL8 e_8_term =
      (-52640. * eta_pow_2 - 13440. * eta + 483280.) * eta * e_pow_8;
  REAL8 e_10_term =
      (10080. * eta_pow_2 + 40320. * eta - 15120.) * eta * e_pow_10;

  REAL8 cosu_1 =
      ((-2240. * eta_pow_3 - 168000. * eta_pow_2 - 424480. * eta) * e_pow_9 +
       (28560. * eta_pow_3 + 242480. * eta_pow_2 + 34440. * pi_pow_2 * eta -
        1340224. * eta + 725760.) *
           e_pow_7 +
       (-33040. * eta_pow_3 - 754880. * eta_pow_2 - 172200. * pi_pow_2 * eta +
        5458480. * eta - 221760.) *
           e_pow_5 +
       (40880. * eta_pow_3 + 738640. * eta_pow_2 + 30135. * pi_pow_2 * eta +
        1554048. * eta - 2936640.) *
           e_pow_3 +
       (-560. * eta_pow_3 - 100240. * eta_pow_2 - 43050. * pi_pow_2 * eta +
        3284816. * eta - 389760.) *
           e) *
      cos(u);

  REAL8 cosu_2 =
      ((4480. * eta_pow_3 - 20160. * eta_pow_2 + 16800. * eta) * e_pow_10 +
       (3920. * eta_pow_3 + 475440. * eta_pow_2 - 17220. * pi_pow_2 * eta +
        831952. * eta - 7257600.) *
           e_pow_8 +
       (-75600. * eta_pow_3 + 96880. * eta_pow_2 + 154980. * pi_pow_2 * eta -
        3249488. * eta - 685440.) *
           e_pow_6 +
       (5040. * eta_pow_3 - 659120. * eta_pow_2 + 25830. * pi_pow_2 * eta -
        7356624. * eta + 6948480.) *
           e_pow_4 +
       (-5040. * eta_pow_3 + 190960. * eta_pow_2 + 137760. * pi_pow_2 * eta -
        7307920. * eta + 107520.) *
           e_pow_2) *
      cos_u_pow_2;

  REAL8 cosu_3 =
      ((560. * eta_pow_3 - 137200. * eta_pow_2 + 388640. * eta + 241920.) *
           e_pow_9 +
       (30800. * eta_pow_3 - 264880. * eta_pow_2 - 68880. * pi_pow_2 * eta +
        624128. * eta + 766080.) *
           e_pow_7 +
       (66640. * eta_pow_3 + 612080. * eta_pow_2 - 8610. * pi_pow_2 * eta +
        6666080. * eta - 6652800.) *
           e_pow_5 +
       (-30800. * eta_pow_3 - 294000. * eta_pow_2 - 223860. * pi_pow_2 * eta +
        9386432. * eta) *
           e_pow_3) *
      cos_u_pow_3;

  REAL8 cosu_4 =
      ((-16240. * eta_pow_3 + 12880. * eta_pow_2 + 18480. * eta) * e_pow_10 +
       (16240. * eta_pow_3 - 91840. * eta_pow_2 + 17220. * pi_pow_2 * eta -
        652192. * eta + 100800.) *
           e_pow_8 +
       (-55440. * eta_pow_3 + 34160. * eta_pow_2 - 30135. * pi_pow_2 * eta -
        2185040. * eta + 2493120.) *
           e_pow_6 +
       (21480. * eta_pow_3 + 86800. * eta_pow_2 + 163590. * pi_pow_2 * eta -
        5713888. * eta + 228480.) *
           e_pow_4) *
      cos_u_pow_4;

  REAL8 cosu_5 =
      ((13440. * eta_pow_3 + 94640. * eta_pow_2 - 113680. * eta - 221760.) *
           e_pow_9 +
       (-11200. * eta_pow_3 - 112000. * eta_pow_2 + 12915. * pi_pow_2 * eta +
        692928. * eta - 194880.) *
           e_pow_7 +
       (4480. * eta_pow_3 + 8960. * eta_pow_2 - 43050. * pi_pow_2 * eta +
        1127280. * eta - 147840.) *
           e_pow_5) *
      cos_u_pow_5;

  REAL8 rt_zero = -67200. * eta_pow_2 + 761600. * eta +
                  e_pow_4 * (40320. * eta_pow_2 + 309120. * eta - 672000.) +
                  e_pow_2 * (208320. * eta_pow_2 + 17220. * pi_pow_2 * eta -
                             2289280. * eta + 1680000.) -
                  8610. * pi_pow_2 * eta - 201600.;

  REAL8 rt_cosu_1 =
      ((-282240. * eta_pow_2 - 450240. * eta + 1478400.) * e_pow_5 +
       (-719040. * eta_pow_2 - 68880. * pi_pow_2 * eta + 8128960. * eta -
        5040000.) *
           e_pow_3 +
       (94080. * eta_pow_2 + 25830. * pi_pow_2 * eta - 1585920. * eta -
        470400.) *
           e) *
      cos(u);

  REAL8 rt_cosu_2 = ((604800. * eta * eta - 504000. * eta - 403200.) * e_pow_6 +
                     (1034880. * eta_pow_2 + 103320. * pi_pow_2 * eta -
                      11195520. * eta + 5779200.) *
                         e_pow_4 +
                     (174720. * eta_pow_2 - 17220. * pi_pow_2 * eta -
                      486080. * eta + 2688000.) *
                         e_pow_2) *
                    cos_u_pow_2;

  REAL8 rt_cosu_3 =
      ((-524160. * eta_pow_2 + 1122240. * eta - 940800.) * e_pow_7 +
       (-873600. * eta_pow_2 - 68880. * pi_pow_2 * eta + 7705600. * eta -
        3897600.) *
           e_pow_5 +
       (-416640. * eta_pow_2 - 17220. * pi_pow_2 * eta + 3357760. * eta -
        3225600.) *
           e_pow_3) *
      cos_u_pow_3;

  REAL8 rt_cosu_4 = ((161280. * eta_pow_2 - 477120. * eta + 537600.) * e_pow_8 +
                     (477120. * eta_pow_2 + 17220. * pi_pow_2 * eta -
                      2894080. * eta + 2217600.) *
                         e_pow_6 +
                     (268800. * eta_pow_2 + 25830. * pi_pow_2 * eta -
                      2721600. * eta + 1276800.) *
                         e_pow_4) *
                    cos_u_pow_4;

  REAL8 rt_cosu_5 =
      ((-127680. * eta_pow_2 + 544320. * eta - 739200.) * e_pow_7 +
       (-53760. * eta_pow_2 - 8610. * pi_pow_2 * eta + 674240. * eta - 67200.) *
           e_pow_5) *
      cos_u_pow_5;

  return (pre_factor * (e_0_term + e_2_term + e_4_term + e_6_term + e_8_term +
                        e_10_term + cosu_1 + cosu_2 + cosu_3 + cosu_4 + cosu_5 +
                        sqrt(e_fact) * (rt_zero + rt_cosu_1 + rt_cosu_2 +
                                        rt_cosu_3 + rt_cosu_4 + rt_cosu_5)));
}

static REAL8
phi_dot_3pn_SS(REAL8 e, REAL8 m1, REAL8 m2, REAL8 S1z, REAL8 S2z,
               REAL8 u) /* Quentin Henry et al terms, arXiv:2308.13606v1 */
{
  REAL8 phi_3pn_SS;
  REAL8 kappa1 = 1.0;
  REAL8 kappa2 = 1.0;
  REAL8 e_2 = e * e;
  REAL8 e_3 = e_2 * e;
  REAL8 e_4 = e_2 * e_2;
  REAL8 e_6 = e_4 * e_2;
  REAL8 e_fact = (1 - e_2);
  REAL8 e_fact_sqrt = sqrt(e_fact);
  REAL8 M_fact_2 = (m1 + m2) * (m1 + m2);
  REAL8 M_fact_4 = M_fact_2 * M_fact_2;

  phi_3pn_SS =
      ((6 *
            (4 * e_6 * e_fact_sqrt * kappa1 -
             2 * (-1 + e_fact_sqrt) * (4 + 7 * kappa1) -
             e_2 * (24 + 20 * e_fact_sqrt + 42 * kappa1 +
                    19 * e_fact_sqrt * kappa1) -
             4 * e_4 * (-4 + (-7 + 3 * e_fact_sqrt) * kappa1)) *
            m1 * m1 * m1 * m1 * S1z * S1z +
        6 *
            (4 * e_6 * e_fact_sqrt * kappa2 -
             2 * (-1 + e_fact_sqrt) * (4 + 7 * kappa2) -
             e_2 * (24 + 20 * e_fact_sqrt + 42 * kappa2 +
                    19 * e_fact_sqrt * kappa2) -
             4 * e_4 * (-4 + (-7 + 3 * e_fact_sqrt) * kappa2)) *
            m2 * m2 * m2 * m2 * S2z * S2z +
        m1 * m1 * m1 * m2 * S1z *
            ((48 * e_6 * e_fact_sqrt * (-1 + kappa1) -
              6 * (-1 + e_fact_sqrt) * (3 + 23 * kappa1) -
              4 * e_4 *
                  (-9 - 60 * e_fact_sqrt - 69 * kappa1 +
                   32 * e_fact_sqrt * kappa1) -
              e_2 * (54 + 297 * e_fact_sqrt + 414 * kappa1 +
                     145 * e_fact_sqrt * kappa1)) *
                 S1z +
             6 *
                 (60 * e_4 + 4 * e_6 * e_fact_sqrt - 30 * (-1 + e_fact_sqrt) -
                  e_2 * (90 + 67 * e_fact_sqrt)) *
                 S2z) +
        m1 * m2 * m2 * m2 * S2z *
            (6 *
                 (60 * e_4 + 4 * e_6 * e_fact_sqrt - 30 * (-1 + e_fact_sqrt) -
                  e_2 * (90 + 67 * e_fact_sqrt)) *
                 S1z +
             (48 * e_6 * e_fact_sqrt * (-1 + kappa2) -
              6 * (-1 + e_fact_sqrt) * (3 + 23 * kappa2) -
              4 * e_4 *
                  (-9 - 60 * e_fact_sqrt - 69 * kappa2 +
                   32 * e_fact_sqrt * kappa2) -
              e_2 * (54 + 297 * e_fact_sqrt + 414 * kappa2 +
                     145 * e_fact_sqrt * kappa2)) *
                 S2z) +
        m1 * m1 * m2 * m2 *
            (3 *
                 (4 * e_6 * e_fact_sqrt * (-4 + 3 * kappa1) -
                  6 * (-1 + e_fact_sqrt) * (-1 + 4 * kappa1) -
                  e_2 * (-18 + 55 * e_fact_sqrt +
                         4 * (18 + e_fact_sqrt) * kappa1) +
                  e_4 * (-12 + 68 * e_fact_sqrt -
                         8 * (-6 + 5 * e_fact_sqrt) * kappa1)) *
                 S1z * S1z +
             2 *
                 (12 * e_6 * e_fact_sqrt - 174 * (-1 + e_fact_sqrt) +
                  4 * e_4 * (87 + 7 * e_fact_sqrt) -
                  e_2 * (522 + 409 * e_fact_sqrt)) *
                 S1z * S2z +
             3 *
                 (4 * e_6 * e_fact_sqrt * (-4 + 3 * kappa2) -
                  6 * (-1 + e_fact_sqrt) * (-1 + 4 * kappa2) -
                  e_2 * (-18 + 55 * e_fact_sqrt +
                         4 * (18 + e_fact_sqrt) * kappa2) +
                  e_4 * (-12 + 68 * e_fact_sqrt -
                         8 * (-6 + 5 * e_fact_sqrt) * kappa2)) *
                 S2z * S2z) +
        e *
            (6 *
                 (-8 + 40 * e_fact_sqrt + 2 * (-7 + 25 * e_fact_sqrt) * kappa1 +
                  4 * e_4 * (-8 + (-14 + 3 * e_fact_sqrt) * kappa1) +
                  e_2 * (40 + 44 * e_fact_sqrt +
                         (70 + 61 * e_fact_sqrt) * kappa1)) *
                 m1 * m1 * m1 * m1 * S1z * S1z +
             6 *
                 (-8 + 40 * e_fact_sqrt + 2 * (-7 + 25 * e_fact_sqrt) * kappa2 +
                  4 * e_4 * (-8 + (-14 + 3 * e_fact_sqrt) * kappa2) +
                  e_2 * (40 + 44 * e_fact_sqrt +
                         (70 + 61 * e_fact_sqrt) * kappa2)) *
                 m2 * m2 * m2 * m2 * S2z * S2z +
             m1 * m1 * m1 * m2 * S1z *
                 ((-18 + 246 * e_fact_sqrt +
                   46 * (-3 + 10 * e_fact_sqrt) * kappa1 +
                   2 * e_4 *
                       (-36 - 60 * e_fact_sqrt - 276 * kappa1 +
                        47 * e_fact_sqrt * kappa1) +
                   e_2 * (90 + 243 * e_fact_sqrt +
                          (690 + 535 * e_fact_sqrt) * kappa1)) *
                      S1z +
                  18 *
                      (-10 + 42 * e_fact_sqrt + 4 * e_4 * (-10 + e_fact_sqrt) +
                       e_2 * (50 + 47 * e_fact_sqrt)) *
                      S2z) +
             m1 * m2 * m2 * m2 * S2z *
                 (18 *
                      (-10 + 42 * e_fact_sqrt + 4 * e_4 * (-10 + e_fact_sqrt) +
                       e_2 * (50 + 47 * e_fact_sqrt)) *
                      S1z +
                  (-18 + 246 * e_fact_sqrt +
                   46 * (-3 + 10 * e_fact_sqrt) * kappa2 +
                   2 * e_4 *
                       (-36 - 60 * e_fact_sqrt - 276 * kappa2 +
                        47 * e_fact_sqrt * kappa2) +
                   e_2 * (90 + 243 * e_fact_sqrt +
                          (690 + 535 * e_fact_sqrt) * kappa2)) *
                      S2z) +
             m1 * m1 * m2 * m2 *
                 (3 *
                      (6 + 14 * e_fact_sqrt +
                       8 * (-3 + 8 * e_fact_sqrt) * kappa1 +
                       4 * e_4 *
                           (6 - 7 * e_fact_sqrt +
                            8 * (-3 + e_fact_sqrt) * kappa1) +
                       e_2 * (5 * (-6 + e_fact_sqrt) +
                              24 * (5 + 3 * e_fact_sqrt) * kappa1)) *
                      S1z * S1z +
                  2 *
                      (-174 + 760 * e_fact_sqrt +
                       e_4 * (-696 + 34 * e_fact_sqrt) +
                       e_2 * (870 + 835 * e_fact_sqrt)) *
                      S1z * S2z +
                  3 *
                      (6 + 14 * e_fact_sqrt +
                       8 * (-3 + 8 * e_fact_sqrt) * kappa2 +
                       4 * e_4 *
                           (6 - 7 * e_fact_sqrt +
                            8 * (-3 + e_fact_sqrt) * kappa2) +
                       e_2 * (5 * (-6 + e_fact_sqrt) +
                              24 * (5 + 3 * e_fact_sqrt) * kappa2)) *
                      S2z * S2z)) *
            cos(u) -
        e_2 *
            (6 *
                 (8 + 56 * e_fact_sqrt + 14 * kappa1 +
                  58 * e_fact_sqrt * kappa1 +
                  4 * e_4 * (-4 - 7 * kappa1 + 3 * e_fact_sqrt * kappa1) +
                  e_2 * (8 + 28 * e_fact_sqrt + 14 * kappa1 +
                         53 * e_fact_sqrt * kappa1)) *
                 m1 * m1 * m1 * m1 * S1z * S1z +
             6 *
                 (8 + 56 * e_fact_sqrt + 14 * kappa2 +
                  58 * e_fact_sqrt * kappa2 +
                  4 * e_4 * (-4 - 7 * kappa2 + 3 * e_fact_sqrt * kappa2) +
                  e_2 * (8 + 28 * e_fact_sqrt + 14 * kappa2 +
                         53 * e_fact_sqrt * kappa2)) *
                 m2 * m2 * m2 * m2 * S2z * S2z +
             m1 * m1 * m1 * m2 * S1z *
                 ((2 * e_4 *
                       (-18 - 12 * e_fact_sqrt - 138 * kappa1 +
                        55 * e_fact_sqrt * kappa1) +
                   2 * (9 + 111 * e_fact_sqrt + 69 * kappa1 +
                        262 * e_fact_sqrt * kappa1) +
                   e_2 * (18 + 171 * e_fact_sqrt + 138 * kappa1 +
                          455 * e_fact_sqrt * kappa1)) *
                      S1z +
                  18 *
                      (10 + 46 * e_fact_sqrt +
                       4 * e_4 * (-5 + 2 * e_fact_sqrt) +
                       e_2 * (10 + 39 * e_fact_sqrt)) *
                      S2z) +
             m1 * m2 * m2 * m2 * S2z *
                 (18 *
                      (10 + 46 * e_fact_sqrt +
                       4 * e_4 * (-5 + 2 * e_fact_sqrt) +
                       e_2 * (10 + 39 * e_fact_sqrt)) *
                      S1z +
                  (2 * e_4 *
                       (-18 - 12 * e_fact_sqrt - 138 * kappa2 +
                        55 * e_fact_sqrt * kappa2) +
                   2 * (9 + 111 * e_fact_sqrt + 69 * kappa2 +
                        262 * e_fact_sqrt * kappa2) +
                   e_2 * (18 + 171 * e_fact_sqrt + 138 * kappa2 +
                          455 * e_fact_sqrt * kappa2)) *
                      S2z) +
             m1 * m1 * m2 * m2 *
                 (3 *
                      (-6 - 14 * e_fact_sqrt + 24 * kappa1 +
                       68 * e_fact_sqrt * kappa1 +
                       4 * e_4 *
                           (3 - 2 * e_fact_sqrt - 12 * kappa1 +
                            7 * e_fact_sqrt * kappa1) +
                       e_2 * (-6 + 13 * e_fact_sqrt + 24 * kappa1 +
                              72 * e_fact_sqrt * kappa1)) *
                      S1z * S1z +
                  2 *
                      (174 + 872 * e_fact_sqrt +
                       e_4 * (-348 + 98 * e_fact_sqrt) +
                       e_2 * (174 + 659 * e_fact_sqrt)) *
                      S1z * S2z +
                  3 *
                      (-6 - 14 * e_fact_sqrt + 24 * kappa2 +
                       68 * e_fact_sqrt * kappa2 +
                       4 * e_4 *
                           (3 - 2 * e_fact_sqrt - 12 * kappa2 +
                            7 * e_fact_sqrt * kappa2) +
                       e_2 * (-6 + 13 * e_fact_sqrt + 24 * kappa2 +
                              72 * e_fact_sqrt * kappa2)) *
                      S2z * S2z)) *
            pow(cos(u), 2) +
        e_3 *
            (6 *
                 (8 + 24 * e_fact_sqrt + 2 * (7 + 9 * e_fact_sqrt) * kappa1 +
                  e_2 * (-8 + 4 * e_fact_sqrt - 14 * kappa1 +
                         23 * e_fact_sqrt * kappa1)) *
                 m1 * m1 * m1 * m1 * S1z * S1z +
             6 *
                 (8 + 24 * e_fact_sqrt + 2 * (7 + 9 * e_fact_sqrt) * kappa2 +
                  e_2 * (-8 + 4 * e_fact_sqrt - 14 * kappa2 +
                         23 * e_fact_sqrt * kappa2)) *
                 m2 * m2 * m2 * m2 * S2z * S2z +
             m1 * m1 * m1 * m2 * S1z *
                 ((18 + 42 * e_fact_sqrt +
                   2 * (69 + 77 * e_fact_sqrt) * kappa1 +
                   e_2 * (-18 + 81 * e_fact_sqrt - 138 * kappa1 +
                          209 * e_fact_sqrt * kappa1)) *
                      S1z +
                  6 *
                      (30 + 38 * e_fact_sqrt +
                       5 * e_2 * (-6 + 11 * e_fact_sqrt)) *
                      S2z) +
             m1 * m2 * m2 * m2 * S2z *
                 (6 *
                      (30 + 38 * e_fact_sqrt +
                       5 * e_2 * (-6 + 11 * e_fact_sqrt)) *
                      S1z +
                  (18 + 42 * e_fact_sqrt +
                   2 * (69 + 77 * e_fact_sqrt) * kappa2 +
                   e_2 * (-18 + 81 * e_fact_sqrt - 138 * kappa2 +
                          209 * e_fact_sqrt * kappa2)) *
                      S2z) +
             m1 * m1 * m2 * m2 *
                 (3 *
                      (-6 - 18 * e_fact_sqrt +
                       8 * (3 + 2 * e_fact_sqrt) * kappa1 +
                       e_2 * (6 + 15 * e_fact_sqrt +
                              8 * (-3 + 5 * e_fact_sqrt) * kappa1)) *
                      S1z * S1z +
                  2 * (174 + 274 * e_fact_sqrt + e_2 * (-174 + 269 * e_fact_sqrt)) *
                      S1z * S2z +
                  3 * (-6 - 18 * e_fact_sqrt + 8 * (3 + 2 * e_fact_sqrt) * kappa2 + e_2 * (6 + 15 * e_fact_sqrt + 8 * (-3 + 5 * e_fact_sqrt) * kappa2)) *
                      S2z * S2z)) *
            pow(cos(u), 3)) /
       (12. * pow(-1 + e_2, 3) * M_fact_4 * pow(-1 + e * cos(u), 5)));

  return (phi_3pn_SS);
}

// static REAL8 phi_dot_3_5pn_SO(REAL8 e, REAL8 m1, REAL8 m2, REAL8 S1z, REAL8
// S2z) /* Computed using inputs from energy and flux expressions. See Blanchet
// liv. rev.
// + Bohe et al. arXiv: 1501.01529 + Marsat et al. arXiv: 1411.4118 */
// {
//   REAL8 kappa1 = 1.0; /*for black holes kappa_{1,2} is 1*/
//   REAL8 kappa2 = 1.0;

//   if(e){

//     return 0;
//   } else{

//    return (-0.5*(M_PI*((1 + 56*kappa1)*pow(m1,2)*pow(S1z,2) +
//         110*m1*m2*S1z*S2z + (1 + 56*kappa2)*pow(m2,2)*pow(S2z,2)))/pow(m1 +
//         m2,2));
//   }
// }

// The non-spinning quasi-circular piece of phi_dot at 4PN and 4.5PN order is 0.

static REAL8 phi_dot_4pn_SS(REAL8 e, REAL8 m1, REAL8 m2, REAL8 S1z,
                            REAL8 S2z) /* Computed using inputs from energy and
flux expressions. See Blanchet liv. rev.
+ Bohe et al. arXiv: 1501.01529 + Marsat et al. arXiv: 1411.4118*/
{
  REAL8 kappa1 = 1.0; /*for black holes kappa_{1,2} is 1*/
  REAL8 kappa2 = 1.0;

  if (e) {

    return 0;
  } else {
    return 0;

    /* previous */ /* return((pow(m1,6)*S1z*(2481357312*M_PI +
         S1z*(-15106969445 + 2637038096*kappa1 -
            63504*pow(S1z + 80*kappa1*S1z,2))) +
      pow(m1,4)*pow(m2,2)*
       (6*pow(S1z,2)*(-7961856981 + 2428186544*kappa1 -
            10584*pow(S1z + 80*kappa1*S1z,2)) -
         8*S1z*(5859558809 + 5016816*(1 + 80*kappa1)*pow(S1z,2))*S2z -
         (5567940197 - 1309909328*kappa2 +
            127008*(12483 + 80*kappa2 + 80*kappa1*(1 + 80*kappa2))*
             pow(S1z,2))*pow(S2z,2) +
         8064*M_PI*(551091*S1z + 332881*S2z)) -
      2*pow(m1,5)*m2*(-4032*M_PI*(690547*S1z + 119880*S2z) +
         S1z*(2*S1z*(10158136535 - 2928982112*kappa1 +
               31752*pow(S1z + 80*kappa1*S1z,2)) +
            (6916951019 + 10033632*(1 + 80*kappa1)*pow(S1z,2))*S2z)) +
      4*pow(m1,3)*pow(m2,3)*
       ((-6683990615 + 1852811072*kappa2)*pow(S2z,2) +
         768606048*M_PI*(S1z + S2z) -
         5016816*pow(S1z,3)*(S2z + 80*kappa1*S2z) +
         3*S1z*S2z*(-5429355259 -
            1672272*(1 + 80*kappa2)*pow(S2z,2)) +
         pow(S1z,2)*(-6683990615 + 1852811072*kappa1 -
            63504*(12483 + 80*kappa1 + 80*(1 + 80*kappa1)*kappa2)*
             pow(S2z,2))) +
      pow(m2,6)*S2z*(2481357312*M_PI +
         S2z*(-15106969445 + 2637038096*kappa2 -
            63504*pow(S2z + 80*kappa2*S2z,2))) +
      pow(m1,2)*pow(m2,4)*
       (8064*M_PI*(332881*S1z + 551091*S2z) +
         8*S1z*S2z*(-5859558809 -
            5016816*(1 + 80*kappa2)*pow(S2z,2)) +
         pow(S1z,2)*(-5567940197 + 1309909328*kappa1 -
            127008*(12483 + 80*kappa1 + 80*(1 + 80*kappa1)*kappa2)*
             pow(S2z,2)) + 6*pow(S2z,2)*
          (-7961856981 + 2428186544*kappa2 -
            10584*pow(S2z + 80*kappa2*S2z,2))) +
      2*m1*pow(m2,5)*(4032*M_PI*(119880*S1z + 690547*S2z) +
         S2z*(S1z*(-6916951019 -
               10033632*(1 + 80*kappa2)*pow(S2z,2)) +
            2*S2z*(-10158136535 + 2928982112*kappa2 -
               31752*pow(S2z + 80*kappa2*S2z,2)))))/
      (1.6257024e7*pow(m1 + m2,6))); */

    /*return (-0.00390625 *
            (((1 + 32 * kappa1) * pow(m1, 2) * pow(S1z, 2) +
              62 * m1 * m2 * S1z * S2z +
              (1 + 32 * kappa2) * pow(m2, 2) * pow(S2z, 2)) *
             ((1 + 80 * kappa1) * pow(m1, 2) * pow(S1z, 2) +
              158 * m1 * m2 * S1z * S2z +
              (1 + 80 * kappa2) * pow(m2, 2) * pow(S2z, 2))) /
            pow(m1 + m2, 4));*/
  }
}

static REAL8 phi_dot_4_5_pn(REAL8 e, REAL8 eta, REAL8 x) {
  REAL8 phi_4_5pn;

  if(e){
    phi_4_5pn=0.0;
  }
  else{
    phi_4_5pn=0.0;
  }
  return(phi_4_5pn);
}


static int eccentric_x_model_odes(REAL8 t, const REAL8 y[], REAL8 dydt[],
                                  void *params) {
  (void)t;
  /* parse the paramater structure */
  struct ode_parameters *op = (struct ode_parameters *)params;
  struct kepler_vars UNUSED *orb_vars = op->orbital_vars;
  REAL8 eta = op->eta;
  REAL8 m1 = op->m1;
  REAL8 m2 = op->m2;
  REAL8 S1z = op->S1z;
  REAL8 S2z = op->S2z;
  int radiation_pn_order = op->radiation_pn_order;
  int retval = GSL_SUCCESS;

  /* input variables */
  REAL8 x = y[0];
  REAL8 e = y[1];
  REAL8 l = y[2];

  REAL8 u = pn_kepler_equation(eta, x, e, l);

  if (e) {
    dydt[0] = dx_dt(radiation_pn_order, eta, m1, m2, S1z, S2z, x, e);
    dydt[1] = de_dt(radiation_pn_order, eta, m1, m2, S1z, S2z, x, e);
    dydt[2] = dl_dt(eta, m1, m2, S1z, S2z, x, e);
    dydt[3] = dphi_dt(u, eta, m1, m2, S1z, S2z, x, e);
  } else {
    /* zero eccentricity limit *
     * arXiv:0909.0066         */
    dydt[0] = dx_dt(radiation_pn_order, eta, m1, m2, S1z, S2z, x, e);
    dydt[1] = de_dt(radiation_pn_order, eta, m1, m2, S1z, S2z, x, e);
    dydt[2] = dl_dt(eta, m1, m2, S1z, S2z, x, e);
    dydt[3] = x * sqrt(x);
  }
  if (dydt[0] == XLAL_REAL8_FAIL_NAN || dydt[1] == XLAL_REAL8_FAIL_NAN) {
    retval = GSL_FAILURE;
    XLAL_ERROR_FAIL(XLAL_EFUNC);
  }

XLAL_FAIL:
  return retval;
}

static REAL8 dx_dt(int radiation_pn_order, REAL8 eta, REAL8 m1, REAL8 m2,
                   REAL8 S1z, REAL8 S2z, REAL8 x, REAL8 e) {

  REAL8 x_pow_5 = x * x * x * x * x;

  REAL8 xdot = XLAL_REAL8_FAIL_NAN;

  if (radiation_pn_order == 0) /* 0 pN term */
  {
    xdot = x_dot_0pn(e, eta) * x_pow_5;

  } else if (radiation_pn_order == 1) /* 0.5 pN term */
  {
    xdot = x_dot_0pn(e, eta) * x_pow_5;

  } else if (radiation_pn_order == 2) /* 1 pN term */
  {
    xdot = (x_dot_0pn(e, eta) + x_dot_1pn(e, eta) * x) * x_pow_5;

  } else if (radiation_pn_order == 3) /* Hereditary terms at 1.5PN order */
  {
    xdot = (x_dot_0pn(e, eta) + x_dot_1pn(e, eta) * x +
            x_dot_1_5_pn(e, eta, m1, m2, S1z, S2z) * x * sqrt(x)) *
               x_pow_5 +
           x_dot_hereditary_1_5(e, eta, x);

  } else if (radiation_pn_order == 4) /* 2 pN term */
  {
    xdot = (x_dot_0pn(e, eta) + x_dot_1pn(e, eta) * x +
            x_dot_1_5_pn(e, eta, m1, m2, S1z, S2z) * x * sqrt(x) +
            x_dot_2pn(e, eta) * x * x +
            x_dot_2pn_SS(e, eta, m1, m2, S1z, S2z) * x * x) *
               x_pow_5 +
           x_dot_hereditary_1_5(e, eta, x);

  } else if (radiation_pn_order ==
             5) /* 2 pN term + hereditary terms up to 2.5PN */
  {
    xdot = (x_dot_0pn(e, eta) + x_dot_1pn(e, eta) * x +
            x_dot_1_5_pn(e, eta, m1, m2, S1z, S2z) * x * sqrt(x) +
            x_dot_2pn(e, eta) * x * x +
            x_dot_2pn_SS(e, eta, m1, m2, S1z, S2z) * x * x +
            x_dot_2_5pn_SO(e, eta, m1, m2, S1z, S2z) * x * x * sqrt(x)) *
               x_pow_5 +
           x_dot_hereditary_1_5(e, eta, x) + x_dot_hereditary_2_5(e, eta, x);

  } else if (radiation_pn_order ==
             6) /* 3 pN term + hereditary terms up to 3PN */
  {
    xdot = (x_dot_0pn(e, eta) + x_dot_1pn(e, eta) * x +
            x_dot_1_5_pn(e, eta, m1, m2, S1z, S2z) * x * sqrt(x) +
            x_dot_2pn(e, eta) * x * x +
            x_dot_2pn_SS(e, eta, m1, m2, S1z, S2z) * x * x +
            x_dot_2_5pn_SO(e, eta, m1, m2, S1z, S2z) * x * x * sqrt(x) +
            x_dot_3pn(e, eta, x) * x * x * x +
            x_dot_3pn_SO(e, eta, m1, m2, S1z, S2z) * x * x * x +
            x_dot_3pn_SS(e, eta, m1, m2, S1z, S2z) * x * x * x) *
               x_pow_5 +
           x_dot_hereditary_1_5(e, eta, x) + x_dot_hereditary_2_5(e, eta, x) +
           x_dot_hereditary_3(e, eta, x);

  } else if (radiation_pn_order ==
             7) /* 3.5 pN term + hereditary terms up to 3PN */
  {
    xdot = (x_dot_0pn(e, eta) + x_dot_1pn(e, eta) * x +
            x_dot_1_5_pn(e, eta, m1, m2, S1z, S2z) * x * sqrt(x) +
            x_dot_2pn(e, eta) * x * x +
            x_dot_2pn_SS(e, eta, m1, m2, S1z, S2z) * x * x +
            x_dot_2_5pn_SO(e, eta, m1, m2, S1z, S2z) * x * x * sqrt(x) +
            x_dot_3pn(e, eta, x) * x * x * x +
            x_dot_3pn_SO(e, eta, m1, m2, S1z, S2z) * x * x * x +
            x_dot_3pn_SS(e, eta, m1, m2, S1z, S2z) * x * x * x +
            x_dot_3_5pnSO(e, eta, m1, m2, S1z, S2z) * x * x * x * sqrt(x) +
            x_dot_3_5_pn(e, eta) * x * x * x * sqrt(x) +
            x_dot_3_5pn_SS(e, eta, m1, m2, S1z, S2z) * x * x * x * sqrt(x) +
            x_dot_3_5pn_cubicSpin(e, eta, m1, m2, S1z, S2z) * x * x * x *
                sqrt(x)) *
               x_pow_5 +
           x_dot_hereditary_1_5(e, eta, x) + x_dot_hereditary_2_5(e, eta, x) +
           x_dot_hereditary_3(e, eta, x);

  } else if (radiation_pn_order == 8) /* 3PN eccentric terms + hereditary terms
                                         up to 3PN + 5PN flux + SF corrections*/
  {
    xdot =
        (x_dot_0pn(e, eta) + x_dot_1pn(e, eta) * x +
         x_dot_1_5_pn(e, eta, m1, m2, S1z, S2z) * x * sqrt(x) +
         x_dot_2pn(e, eta) * x * x +
         x_dot_2pn_SS(e, eta, m1, m2, S1z, S2z) * x * x +
         x_dot_2_5pn_SO(e, eta, m1, m2, S1z, S2z) * x * x * sqrt(x) +
         x_dot_3pn(e, eta, x) * x * x * x +
         x_dot_3pn_SO(e, eta, m1, m2, S1z, S2z) * x * x * x +
         x_dot_3pn_SS(e, eta, m1, m2, S1z, S2z) * x * x * x +
         x_dot_3_5pnSO(e, eta, m1, m2, S1z, S2z) * x * x * x * sqrt(x) +
         x_dot_3_5_pn(e, eta) * x * x * x * sqrt(x) +
         x_dot_3_5pn_SS(e, eta, m1, m2, S1z, S2z) * x * x * x * sqrt(x) +
         x_dot_3_5pn_cubicSpin(e, eta, m1, m2, S1z, S2z) * x * x * x * sqrt(x) +
         (x_dot_4pn(e, eta, x) + x_dot_4pnSO(e, eta, m1, m2, S1z, S2z) +
          x_dot_4pnSS(e, eta, m1, m2, S1z, S2z)) *
             x * x * x * x) *
            x_pow_5 +
        x_dot_hereditary_1_5(e, eta, x) + x_dot_hereditary_2_5(e, eta, x) +
        x_dot_hereditary_3(e, eta, x) + dxdt_4pn(x, eta);

  } else if (radiation_pn_order == 9) /* 3PN eccentric terms + hereditary terms
                                         up to 3PN + 6PN flux + SF corrections*/
  {
    xdot =
        (x_dot_0pn(e, eta) + x_dot_1pn(e, eta) * x +
         x_dot_1_5_pn(e, eta, m1, m2, S1z, S2z) * x * sqrt(x) +
         x_dot_2pn(e, eta) * x * x +
         x_dot_2pn_SS(e, eta, m1, m2, S1z, S2z) * x * x +
         x_dot_2_5pn_SO(e, eta, m1, m2, S1z, S2z) * x * x * sqrt(x) +
         x_dot_3pn(e, eta, x) * x * x * x +
         x_dot_3pn_SO(e, eta, m1, m2, S1z, S2z) * x * x * x +
         x_dot_3pn_SS(e, eta, m1, m2, S1z, S2z) * x * x * x +
         x_dot_3_5pnSO(e, eta, m1, m2, S1z, S2z) * x * x * x * sqrt(x) +
         x_dot_3_5_pn(e, eta) * x * x * x * sqrt(x) +
         x_dot_3_5pn_SS(e, eta, m1, m2, S1z, S2z) * x * x * x * sqrt(x) +
         x_dot_3_5pn_cubicSpin(e, eta, m1, m2, S1z, S2z) * x * x * x * sqrt(x) +
         (x_dot_4pn(e, eta, x) + x_dot_4pnSO(e, eta, m1, m2, S1z, S2z) +
          x_dot_4pnSS(e, eta, m1, m2, S1z, S2z)) *
             x * x * x * x +
         x_dot_4_5_pn(e, eta, x) * x * x * x * x * sqrt(x)) *
            x_pow_5 +
        x_dot_hereditary_1_5(e, eta, x) + x_dot_hereditary_2_5(e, eta, x) +
        x_dot_hereditary_3(e, eta, x) + dxdt_4_5pn(x, eta);

  } else if (radiation_pn_order ==
             10) /* 3PN eccentric terms + hereditary terms up to 3PN + 6PN flux
                    +  SF corrections*/
  {
    xdot =
        (x_dot_0pn(e, eta) + x_dot_1pn(e, eta) * x +
         x_dot_1_5_pn(e, eta, m1, m2, S1z, S2z) * x * sqrt(x) +
         x_dot_2pn(e, eta) * x * x +
         x_dot_2pn_SS(e, eta, m1, m2, S1z, S2z) * x * x +
         x_dot_2_5pn_SO(e, eta, m1, m2, S1z, S2z) * x * x * sqrt(x) +
         x_dot_3pn(e, eta, x) * x * x * x +
         x_dot_3pn_SO(e, eta, m1, m2, S1z, S2z) * x * x * x +
         x_dot_3pn_SS(e, eta, m1, m2, S1z, S2z) * x * x * x +
         x_dot_3_5pnSO(e, eta, m1, m2, S1z, S2z) * x * x * x * sqrt(x) +
         x_dot_3_5_pn(e, eta) * x * x * x * sqrt(x) +
         x_dot_3_5pn_SS(e, eta, m1, m2, S1z, S2z) * x * x * x * sqrt(x) +
         x_dot_3_5pn_cubicSpin(e, eta, m1, m2, S1z, S2z) * x * x * x * sqrt(x) +
         (x_dot_4pn(e, eta, x) + x_dot_4pnSO(e, eta, m1, m2, S1z, S2z) +
          x_dot_4pnSS(e, eta, m1, m2, S1z, S2z)) *
             x * x * x * x +
         x_dot_4_5_pn(e, eta, x) * x * x * x * x * sqrt(x)) *
            x_pow_5 +
        x_dot_hereditary_1_5(e, eta, x) + x_dot_hereditary_2_5(e, eta, x) +
        x_dot_hereditary_3(e, eta, x) + dxdt_5pn(x, eta);

  } else if (radiation_pn_order ==
             11) /* 3PN eccentric terms + hereditary terms up to 3PN + 6PN flux
                    +  SF corrections*/
  {
    xdot =
        (x_dot_0pn(e, eta) + x_dot_1pn(e, eta) * x +
         x_dot_1_5_pn(e, eta, m1, m2, S1z, S2z) * x * sqrt(x) +
         x_dot_2pn(e, eta) * x * x +
         x_dot_2pn_SS(e, eta, m1, m2, S1z, S2z) * x * x +
         x_dot_2_5pn_SO(e, eta, m1, m2, S1z, S2z) * x * x * sqrt(x) +
         x_dot_3pn(e, eta, x) * x * x * x +
         x_dot_3pn_SO(e, eta, m1, m2, S1z, S2z) * x * x * x +
         x_dot_3pn_SS(e, eta, m1, m2, S1z, S2z) * x * x * x +
         x_dot_3_5pnSO(e, eta, m1, m2, S1z, S2z) * x * x * x * sqrt(x) +
         x_dot_3_5_pn(e, eta) * x * x * x * sqrt(x) +
         x_dot_3_5pn_SS(e, eta, m1, m2, S1z, S2z) * x * x * x * sqrt(x) +
         x_dot_3_5pn_cubicSpin(e, eta, m1, m2, S1z, S2z) * x * x * x * sqrt(x) +
         (x_dot_4pn(e, eta, x) + x_dot_4pnSO(e, eta, m1, m2, S1z, S2z) +
          x_dot_4pnSS(e, eta, m1, m2, S1z, S2z)) *
             x * x * x * x +
         x_dot_4_5_pn(e, eta, x) * x * x * x * x * sqrt(x)) *
            x_pow_5 +
        x_dot_hereditary_1_5(e, eta, x) + x_dot_hereditary_2_5(e, eta, x) +
        x_dot_hereditary_3(e, eta, x) + dxdt_5_5pn(x, eta);

  } else if (radiation_pn_order ==
             12) /* 3PN eccentric terms + hereditary terms up to 3PN + 6PN flux
                    +  SF corrections*/
  {
    xdot =
        (x_dot_0pn(e, eta) + x_dot_1pn(e, eta) * x +
         x_dot_1_5_pn(e, eta, m1, m2, S1z, S2z) * x * sqrt(x) +
         x_dot_2pn(e, eta) * x * x +
         x_dot_2pn_SS(e, eta, m1, m2, S1z, S2z) * x * x +
         x_dot_2_5pn_SO(e, eta, m1, m2, S1z, S2z) * x * x * sqrt(x) +
         x_dot_3pn(e, eta, x) * x * x * x +
         x_dot_3pn_SO(e, eta, m1, m2, S1z, S2z) * x * x * x +
         x_dot_3pn_SS(e, eta, m1, m2, S1z, S2z) * x * x * x +
         x_dot_3_5pnSO(e, eta, m1, m2, S1z, S2z) * x * x * x * sqrt(x) +
         x_dot_3_5_pn(e, eta) * x * x * x * sqrt(x) +
         x_dot_3_5pn_SS(e, eta, m1, m2, S1z, S2z) * x * x * x * sqrt(x) +
         x_dot_3_5pn_cubicSpin(e, eta, m1, m2, S1z, S2z) * x * x * x * sqrt(x) +
         (x_dot_4pn(e, eta, x) + x_dot_4pnSO(e, eta, m1, m2, S1z, S2z) +
          x_dot_4pnSS(e, eta, m1, m2, S1z, S2z)) *
             x * x * x * x +
         x_dot_4_5_pn(e, eta, x) * x * x * x * x * sqrt(x)) *
            x_pow_5 +
        x_dot_hereditary_1_5(e, eta, x) + x_dot_hereditary_2_5(e, eta, x) +
        x_dot_hereditary_3(e, eta, x) + dxdt_6pn(x, eta);

  } else {
    // throw_exception( ECC_VALUE_ERROR, "error in conservative pn order" );
    XLAL_ERROR_REAL8(XLAL_EINVAL,
                     "error in radiation pn order, order %d not supported",
                     radiation_pn_order);
    return XLAL_REAL8_FAIL_NAN;
  }

  return xdot;
}

static REAL8 de_dt(int radiation_pn_order, REAL8 eta, REAL8 m1, REAL8 m2,
                   REAL8 S1z, REAL8 S2z, REAL8 x, REAL8 e) {
  REAL8 x_pow_4 = x * x * x * x;
  REAL8 edot = XLAL_REAL8_FAIL_NAN;

  if (radiation_pn_order == 0) /* 0 pN term */
  {
    edot = e_dot_0pn(e, eta) * x_pow_4;
  } else if (radiation_pn_order == 1) /* 0.5 pN term */
  {
    edot = e_dot_0pn(e, eta) * x_pow_4;
  } else if (radiation_pn_order == 2) /* 1 pN term */
  {
    edot = (e_dot_0pn(e, eta) + e_dot_1pn(e, eta) * x) * x_pow_4;
  } else if (radiation_pn_order == 3) /* 1.5 pN term */
  {
    edot = (e_dot_0pn(e, eta) + e_dot_1pn(e, eta) * x +
            e_dot_1_5pn_SO(e, m1, m2, S1z, S2z) * x * sqrt(x)) *
               x_pow_4 +
           e_rad_hereditary_1_5(e, eta, x);
  } else if (radiation_pn_order == 4) /* 2 pN term */
  {
    edot =
        (e_dot_0pn(e, eta) + e_dot_1pn(e, eta) * x + e_dot_2pn(e, eta) * x * x +
         e_dot_1_5pn_SO(e, m1, m2, S1z, S2z) * x * sqrt(x) +
         e_dot_2pn_SS(e, m1, m2, S1z, S2z) * x * x) *
            x_pow_4 +
        e_rad_hereditary_1_5(e, eta, x);
  } else if (radiation_pn_order ==
             5) /* 2 pN term + heriditary terms up to 2.5PN */
  {
    edot =
        (e_dot_0pn(e, eta) + e_dot_1pn(e, eta) * x + e_dot_2pn(e, eta) * x * x +
         e_dot_1_5pn_SO(e, m1, m2, S1z, S2z) * x * sqrt(x) +
         e_dot_2pn_SS(e, m1, m2, S1z, S2z) * x * x +
         e_dot_2_5pn_SO(e, m1, m2, S1z, S2z) * x * x * sqrt(x)) *
            x_pow_4 +
        e_rad_hereditary_1_5(e, eta, x) + e_rad_hereditary_2_5(e, eta, x);
  } else if (radiation_pn_order ==
             6) /* 3 pN term + heriditary terms up to 3PN */
  {
    edot = (e_dot_0pn(e, eta) + e_dot_1pn(e, eta) * x +
            e_dot_2pn(e, eta) * x * x + e_dot_3pn(e, eta, x) * x * x * x +
            e_dot_1_5pn_SO(e, m1, m2, S1z, S2z) * x * sqrt(x) +
            e_dot_2pn_SS(e, m1, m2, S1z, S2z) * x * x +
            e_dot_2_5pn_SO(e, m1, m2, S1z, S2z) * x * x * sqrt(x) +
            (e_dot_3pn_SO(e, m1, m2, S1z, S2z) +
             e_dot_3pn_SS(e, m1, m2, S1z, S2z)) *
                x * x * x) *
               x_pow_4 +
           e_rad_hereditary_1_5(e, eta, x) + e_rad_hereditary_2_5(e, eta, x) +
           e_rad_hereditary_3(e, eta, x);
  } else if (radiation_pn_order ==
             7) /* 3.5 pN term + heriditary terms up to 3PN */
  {
    edot = (e_dot_0pn(e, eta) + e_dot_1pn(e, eta) * x +
            e_dot_2pn(e, eta) * x * x + e_dot_3pn(e, eta, x) * x * x * x +
            e_dot_3_5pn(e, eta) * x * x * x * sqrt(x) +
            e_dot_1_5pn_SO(e, m1, m2, S1z, S2z) * x * sqrt(x) +
            e_dot_2pn_SS(e, m1, m2, S1z, S2z) * x * x +
            e_dot_2_5pn_SO(e, m1, m2, S1z, S2z) * x * x * sqrt(x) +
            (e_dot_3pn_SO(e, m1, m2, S1z, S2z) +
             e_dot_3pn_SS(e, m1, m2, S1z, S2z)) *
                x * x * x) *
               x_pow_4 +
           e_rad_hereditary_1_5(e, eta, x) + e_rad_hereditary_2_5(e, eta, x) +
           e_rad_hereditary_3(e, eta, x);
  } else if (radiation_pn_order ==
             8) /* 3.5 pN term + heriditary terms up to 3PN */
  {
    edot = (e_dot_0pn(e, eta) + e_dot_1pn(e, eta) * x +
            e_dot_2pn(e, eta) * x * x + e_dot_3pn(e, eta, x) * x * x * x +
            e_dot_3_5pn(e, eta) * x * x * x * sqrt(x) +
            e_dot_1_5pn_SO(e, m1, m2, S1z, S2z) * x * sqrt(x) +
            e_dot_2pn_SS(e, m1, m2, S1z, S2z) * x * x +
            e_dot_2_5pn_SO(e, m1, m2, S1z, S2z) * x * x * sqrt(x) +
            (e_dot_3pn_SO(e, m1, m2, S1z, S2z) +
             e_dot_3pn_SS(e, m1, m2, S1z, S2z)) *
                x * x * x) *
               x_pow_4 +
           e_rad_hereditary_1_5(e, eta, x) + e_rad_hereditary_2_5(e, eta, x) +
           e_rad_hereditary_3(e, eta, x);
  } else if (radiation_pn_order ==
             9) /* 3.5 pN term + heriditary terms up to 3PN */
  {
    edot = (e_dot_0pn(e, eta) + e_dot_1pn(e, eta) * x +
            e_dot_2pn(e, eta) * x * x + e_dot_3pn(e, eta, x) * x * x * x +
            e_dot_3_5pn(e, eta) * x * x * x * sqrt(x) +
            e_dot_1_5pn_SO(e, m1, m2, S1z, S2z) * x * sqrt(x) +
            e_dot_2pn_SS(e, m1, m2, S1z, S2z) * x * x +
            e_dot_2_5pn_SO(e, m1, m2, S1z, S2z) * x * x * sqrt(x) +
            (e_dot_3pn_SO(e, m1, m2, S1z, S2z) +
             +e_dot_3pn_SS(e, m1, m2, S1z, S2z)) *
                x * x * x) *
               x_pow_4 +
           e_rad_hereditary_1_5(e, eta, x) + e_rad_hereditary_2_5(e, eta, x) +
           e_rad_hereditary_3(e, eta, x);
  } else if (radiation_pn_order ==
             10) /* 3.5 pN term + heriditary terms up to 3PN */
  {
    edot = (e_dot_0pn(e, eta) + e_dot_1pn(e, eta) * x +
            e_dot_2pn(e, eta) * x * x + e_dot_3pn(e, eta, x) * x * x * x +
            e_dot_3_5pn(e, eta) * x * x * x * sqrt(x) +
            e_dot_1_5pn_SO(e, m1, m2, S1z, S2z) * x * sqrt(x) +
            e_dot_2pn_SS(e, m1, m2, S1z, S2z) * x * x +
            e_dot_2_5pn_SO(e, m1, m2, S1z, S2z) * x * x * sqrt(x) +
            (e_dot_3pn_SO(e, m1, m2, S1z, S2z) +
             e_dot_3pn_SS(e, m1, m2, S1z, S2z)) *
                x * x * x) *
               x_pow_4 +
           e_rad_hereditary_1_5(e, eta, x) + e_rad_hereditary_2_5(e, eta, x) +
           e_rad_hereditary_3(e, eta, x);
  } else if (radiation_pn_order ==
             11) /* 3.5 pN term + heriditary terms up to 3PN */
  {
    edot = (e_dot_0pn(e, eta) + e_dot_1pn(e, eta) * x +
            e_dot_2pn(e, eta) * x * x + e_dot_3pn(e, eta, x) * x * x * x +
            e_dot_3_5pn(e, eta) * x * x * x * sqrt(x) +
            e_dot_1_5pn_SO(e, m1, m2, S1z, S2z) * x * sqrt(x) +
            e_dot_2pn_SS(e, m1, m2, S1z, S2z) * x * x +
            e_dot_2_5pn_SO(e, m1, m2, S1z, S2z) * x * x * sqrt(x) +
            (e_dot_3pn_SO(e, m1, m2, S1z, S2z) +
             e_dot_3pn_SS(e, m1, m2, S1z, S2z)) *
                x * x * x) *
               x_pow_4 +
           e_rad_hereditary_1_5(e, eta, x) + e_rad_hereditary_2_5(e, eta, x) +
           e_rad_hereditary_3(e, eta, x);
  } else if (radiation_pn_order ==
             12) /* 3.5 pN term + heriditary terms up to 3PN */
  {
    edot = (e_dot_0pn(e, eta) + e_dot_1pn(e, eta) * x +
            e_dot_2pn(e, eta) * x * x + e_dot_3pn(e, eta, x) * x * x * x +
            e_dot_3_5pn(e, eta) * x * x * x * sqrt(x) +
            e_dot_1_5pn_SO(e, m1, m2, S1z, S2z) * x * sqrt(x) +
            e_dot_2pn_SS(e, m1, m2, S1z, S2z) * x * x +
            e_dot_2_5pn_SO(e, m1, m2, S1z, S2z) * x * x * sqrt(x) +
            (e_dot_3pn_SO(e, m1, m2, S1z, S2z) +
             e_dot_3pn_SS(e, m1, m2, S1z, S2z)) *
                x * x * x) *
               x_pow_4 +
           e_rad_hereditary_1_5(e, eta, x) + e_rad_hereditary_2_5(e, eta, x) +
           e_rad_hereditary_3(e, eta, x);
  } else {
    /*fprintf( stderr, "Error in PN order: %d\n", radiation_pn_order );*/
    // throw_exception( ECC_VALUE_ERROR, "error in conservative pn order" );
    XLAL_ERROR_REAL8(XLAL_EINVAL, "error in radiation pn order: %d",
                     radiation_pn_order);
    return XLAL_REAL8_FAIL_NAN;
  }

  return edot;
}

static REAL8 dl_dt(REAL8 eta, REAL8 m1, REAL8 m2, REAL8 S1z, REAL8 S2z, REAL8 x,
                   REAL8 e) {
  REAL8 x_pow_3_2 = sqrt(x) * x;
  REAL8 ldot = 0;

  // 3PN accurate with spin corrections

  ldot =
      (1.0 + x * l_dot_1pn(e, eta) +
       x_pow_3_2 * l_dot_1_5pn_SO(e, m1, m2, S1z, S2z) +
       x * x * l_dot_2pn(e, eta) + x * x * l_dot_2pn_SS(e, m1, m2, S1z, S2z) +
       l_dot_2_5pn_SO(e, m1, m2, S1z, S2z) * x_pow_3_2 * x +
       x * x * x * l_dot_3pn(e, eta) +
       x * x * x * l_dot_3pn_SS(e, m1, m2, S1z, S2z)) *
      x_pow_3_2;

  return ldot;
}

static REAL8 dphi_dt(REAL8 u, REAL8 eta, REAL8 m1, REAL8 m2, REAL8 S1z,
                     REAL8 S2z, REAL8 x, REAL8 e) {
  REAL8 x_pow_3_2 = sqrt(x) * x;
  REAL8 phidot = 0;

  // 4PN accurate

  phidot =
      ((phi_dot_0pn(e, eta, u) + x * phi_dot_1pn(e, eta, u) +
        x_pow_3_2 * phi_dot_1_5_pnSO_ecc(e, m1, m2, S1z, S2z, u) +
        x * x * phi_dot_2_pnSS_ecc(e, m1, m2, S1z, S2z, u) +
        x * x * phi_dot_2pn(e, eta, u) +
        x_pow_3_2 * x * phi_dot_2_5pn_SO(e, m1, m2, S1z, S2z, u) +
        x * x * x * phi_dot_3pn(e, eta, u) +
        x_pow_3_2 * x_pow_3_2 * phi_dot_3pn_SS(e, m1, m2, S1z, S2z, u)
        /* + phi_dot_3_5pn_SO(e,m1,m2,S1z,S2z) * x_pow_3_2 * x * x */
        + phi_dot_4pn_SS(e, m1, m2, S1z, S2z) * x_pow_3_2 * x * x * sqrt(x)
        + phi_dot_4_5_pn(e, eta, x)*x_pow_3_2*x_pow_3_2*x*sqrt(x)) *
       x_pow_3_2);

  return phidot;
}

static REAL8 rel_sep_0pn(REAL8 e, REAL8 u) { return (1.0 - e * cos(u)); }

static REAL8 rel_sep_1pn(REAL8 e, REAL8 u, REAL8 eta) {
  REAL8 e_fact = 1.0 - e * e;
  REAL8 bit_1 = 2. * (1.0 - e * cos(u)) / e_fact;
  REAL8 bit_2 = (-18.0 + 2. * eta - e * (6. - 7. * eta) * cos(u)) / 6.;

  return (bit_1 + bit_2);
}

static REAL8
rel_sep_1_5pn(REAL8 e, REAL8 u, REAL8 m1, REAL8 m2, REAL8 S1z,
              REAL8 S2z) /*computed using inputs from klein et al.
                            arXiv:1801.08542. See equation B1a and B2a, B2b. */
{

  return (-0.3333333333333333 *
          (2 * pow(m1, 2) * (S1z + 3 * pow(e, 2) * S1z) +
           2 * (1 + 3 * pow(e, 2)) * pow(m2, 2) * S2z +
           3 * (1 + pow(e, 2)) * m1 * m2 * (S1z + S2z) -
           2 * e *
               (4 * pow(m1, 2) * S1z + 4 * pow(m2, 2) * S2z +
                3 * m1 * m2 * (S1z + S2z)) *
               cos(u)) /
          (pow(1 - pow(e, 2), 1.5) * pow(m1 + m2, 2)));
}

static REAL8 rel_sep_2pn(REAL8 e, REAL8 u, REAL8 eta) {
  REAL8 eta_pow_2 = eta * eta;
  REAL8 e_pow_2 = e * e;
  REAL8 e_fact = 1.0 - e * e;

  REAL8 num_1 =
      (-48. + 28. * eta + e_pow_2 * (-51. + 26. * eta)) * (-1. + e * cos(u));
  REAL8 den_1 = 6. * e_fact * e_fact;
  REAL8 num_2 = 72. * (-4. + 7. * eta) +
                36. * sqrt(e_fact) * (-5. + 2. * eta) * (2. + e * cos(u)) +
                e_fact * (72. + 30. * eta + 8. * eta_pow_2 +
                          e * (-72. + 7 * (33. - 5. * eta) * eta) * cos(u));
  REAL8 den_2 = 72. * e_fact;

  return (num_1 / den_1 + num_2 / den_2);
}

static REAL8 rel_sep_2_5pn_SO(REAL8 e, REAL8 u, REAL8 m1, REAL8 m2, REAL8 S1z,
                              REAL8 S2z) {
  REAL8 r_2_5pn_SO;
  REAL8 e_2 = e * e;
  REAL8 e_4 = e_2 * e_2;
  REAL8 e_fact = (1 - e_2);
  REAL8 e_fact_sqrt = sqrt(e_fact);
  REAL8 M_fact_2 = (m1 + m2) * (m1 + m2);
  REAL8 M_fact_4 = M_fact_2 * M_fact_2;

  r_2_5pn_SO =
      ((2 * pow(-1 + e_2, 2) *
            (-12 * m1 * m1 * m1 * m1 * S1z - 12 * m2 * m2 * m2 * m2 * S2z -
             21 * m1 * m1 * m2 * m2 * (S1z + S2z) -
             2 * m1 * m1 * m1 * m2 * (13 * S1z + 3 * S2z) -
             2 * m1 * m2 * m2 * m2 * (3 * S1z + 13 * S2z)) *
            (2 + e * cos(u)) +
        2 * e_fact_sqrt *
            (12 * (2 + 7 * e_2 + e_4) * m1 * m1 * m1 * m1 * S1z +
             12 * (2 + 7 * e_2 + e_4) * m2 * m2 * m2 * m2 * S2z +
             2 * (22 + 85 * e_2 + 10 * e_4) * m1 * m1 * m2 * m2 * (S1z + S2z) +
             m1 * m1 * m1 * m2 *
                 (2 * (28 + 96 * e_2 + 11 * e_4) * S1z +
                  3 * (4 + 19 * e_2 + 2 * e_4) * S2z) +
             m1 * m2 * m2 * m2 *
                 (3 * (4 + 19 * e_2 + 2 * e_4) * S1z +
                  2 * (28 + 96 * e_2 + 11 * e_4) * S2z) -
             e *
                 (60 * (1 + e_2) * m1 * m1 * m1 * m1 * S1z +
                  60 * (1 + e_2) * m2 * m2 * m2 * m2 * S2z +
                  3 * (35 + 43 * e_2) * m1 * m1 * m2 * m2 * (S1z + S2z) +
                  m1 * m1 * m1 * m2 *
                      (133 * S1z + 137 * e_2 * S1z + 30 * S2z +
                       45 * e_2 * S2z) +
                  m1 * m2 * m2 * m2 *
                      (30 * S1z + 45 * e_2 * S1z + 133 * S2z +
                       137 * e_2 * S2z)) *
                 cos(u))) /
       (6. * pow(-1 + e_2, 3) * M_fact_4));

  return (r_2_5pn_SO);
}

static REAL8
rel_sep_2pnSS(REAL8 e, REAL8 u, REAL8 m1, REAL8 m2, REAL8 S1z,
              REAL8 S2z) /*computed using inputs from klein et al.
                            arXiv:1801.08542. See equation B1a and B2a, B2b. */
{
  REAL8 kappa1 = 1.0;
  REAL8 kappa2 = 1.0;

  return (((m1 * S1z * (2 * m2 * S2z + m1 * S1z * kappa1) +
            pow(m2, 2) * pow(S2z, 2) * kappa2) *
           (1 + pow(e, 2) - 2 * e * cos(u))) /
          (2. * pow(-1 + pow(e, 2), 2) * pow(m1 + m2, 2)));
}

static REAL8 rel_sep_3pn(REAL8 e, REAL8 u, REAL8 eta) {
  REAL8 pi_pow_2 = M_PI * M_PI;
  REAL8 e_factor = 1.0 - e * e;
  REAL8 eta_pow_2 = eta * eta;
  REAL8 e_pow_2 = e * e;
  REAL8 e_pow_4 = e_pow_2 * e_pow_2;

  return (1. / (181440. * pow7_2(e_factor))) *
         (((-665280. * eta_pow_2 + 1753920. * eta - 1814400.) * e_pow_4 +
           ((725760. * eta_pow_2 - 77490. * pi_pow_2 + 5523840.) * eta -
            3628800.) *
               e_pow_2 +
           ((544320. * eta_pow_2 + 154980. * pi_pow_2 - 14132160.) * eta +
            7257600.)) *
              e_pow_2 -
          604800. * eta_pow_2 + 6854400. * eta +
          (((302400. * eta_pow_2 - 1254960. * eta + 453600.) * e_pow_4 +
            ((-1542240. * eta_pow_2 - 38745. * pi_pow_2 + 6980400.) * eta -
             453600.) *
                e_pow_2 +
            ((2177280. * eta_pow_2 + 77490. * pi_pow_2 - 12373200.) * eta +
             4989600.)) *
               e * e_pow_2 +
           ((-937440. * eta_pow_2 - 37845. * pi_pow_2 + 6647760.) * eta -
            4989600.) *
               e) *
              cos(u) +
          sqrt(e_factor) *
              ((((-4480. * eta_pow_2 - 25200. * eta + 22680.) * eta - 120960.) *
                    e_pow_4 +
                (13440. * eta_pow_2 * eta + 4404960. * eta * eta +
                 116235. * pi_pow_2 - 12718296. * eta + 5261760.) *
                    e_pow_2 +
                ((-13440. * eta_pow_2 + 2242800. * eta + 348705. * pi_pow_2 -
                  19225080.) *
                     eta +
                 1614160.)) *
                   e_pow_2 +
               (4480. * eta_pow_2 + 45360. * eta - 8600904.) * eta +
               ((((-6860. * eta_pow_2 - 550620. * eta - 986580.) * eta +
                  120960.) *
                     e_pow_4 +
                 ((20580. * eta_pow_2 - 2458260. * eta + 3458700.) * eta -
                  2358720.) *
                     e_pow_2 +
                 ((-20580. * eta_pow_2 - 3539340. * eta - 116235. * pi_pow_2 +
                   20173860.) *
                      eta -
                  16148160.)) *
                    e * e_pow_2 +
                ((6860. * eta_pow_2 - 1220940. * eta + 464940. * pi_pow_2 +
                  17875620.) *
                     eta -
                 417440.) *
                    e) *
                   cos(u) +
               116235. * pi_pow_2 * eta + 1814400.) -
          77490. * pi_pow_2 * eta - 1814400.);
}

static REAL8 rel_sep_3pn_SS(REAL8 e, REAL8 u, REAL8 m1, REAL8 m2, REAL8 S1z,
                            REAL8 S2z) {
  REAL8 r_3pn_SS;
  REAL8 kappa1 = 1.0;
  REAL8 kappa2 = 1.0;
  REAL8 e_2 = e * e;
  REAL8 e_4 = e_2 * e_2;
  REAL8 e_fact = (1 - e_2);
  REAL8 e_fact_sqrt = sqrt(e_fact);
  REAL8 e_fact_35 = e_fact_sqrt * e_fact * e_fact * e_fact;
  REAL8 M_fact_2 = (m1 + m2) * (m1 + m2);
  REAL8 M_fact_4 = M_fact_2 * M_fact_2;

  r_3pn_SS =
      (-0.05555555555555555 *
       (-((-48 + 40 * e_fact_sqrt - 84 * kappa1 + 99 * e_fact_sqrt * kappa1 +
           6 * e_2 *
               (16 + 28 * e_fact_sqrt + 28 * kappa1 +
                51 * e_fact_sqrt * kappa1) +
           e_4 * (-48 + (-84 + 45 * e_fact_sqrt) * kappa1)) *
          m1 * m1 * m1 * m1 * S1z * S1z) -
        (-48 + 40 * e_fact_sqrt - 84 * kappa2 + 99 * e_fact_sqrt * kappa2 +
         6 * e_2 *
             (16 + 28 * e_fact_sqrt + 28 * kappa2 + 51 * e_fact_sqrt * kappa2) +
         e_4 * (-48 + (-84 + 45 * e_fact_sqrt) * kappa2)) *
            m2 * m2 * m2 * m2 * S2z * S2z -
        3 * m1 * m1 * m1 * m2 * S1z *
            ((-6 + 4 * e_fact_sqrt - 46 * kappa1 + 48 * e_fact_sqrt * kappa1 +
              e_4 * (-6 - 6 * e_fact_sqrt - 46 * kappa1 +
                     25 * e_fact_sqrt * kappa1) +
              e_2 * (12 + 55 * e_fact_sqrt + 92 * kappa1 +
                     158 * e_fact_sqrt * kappa1)) *
                 S1z +
             2 *
                 (-30 + 29 * e_fact_sqrt +
                  5 * e_2 *
                      (12 + 25 * e_fact_sqrt + 3 * e_2 * (-2 + e_fact_sqrt))) *
                 S2z) -
        3 * m1 * m2 * m2 * m2 * S2z *
            (2 *
                 (-30 + 29 * e_fact_sqrt +
                  5 * e_2 *
                      (12 + 25 * e_fact_sqrt + 3 * e_2 * (-2 + e_fact_sqrt))) *
                 S1z +
             (-6 + 4 * e_fact_sqrt - 46 * kappa2 + 48 * e_fact_sqrt * kappa2 +
              e_4 * (-6 - 6 * e_fact_sqrt - 46 * kappa2 +
                     25 * e_fact_sqrt * kappa2) +
              e_2 * (12 + 55 * e_fact_sqrt + 92 * kappa2 +
                     158 * e_fact_sqrt * kappa2)) *
                 S2z) -
        m1 * m1 * m2 * m2 *
            (9 *
                 (2 - 2 * e_fact_sqrt - 8 * kappa1 + 6 * e_fact_sqrt * kappa1 +
                  e_4 * (2 - 2 * e_fact_sqrt - 8 * kappa1 +
                         6 * e_fact_sqrt * kappa1) +
                  e_2 * (-4 + 3 * e_fact_sqrt +
                         4 * (4 + 7 * e_fact_sqrt) * kappa1)) *
                 S1z * S1z +
             2 *
                 (-174 + 175 * e_fact_sqrt + 348 * e_2 * (1 + 2 * e_fact_sqrt) +
                  6 * e_4 * (-29 + 11 * e_fact_sqrt)) *
                 S1z * S2z +
             9 *
                 (2 - 2 * e_fact_sqrt - 8 * kappa2 + 6 * e_fact_sqrt * kappa2 +
                  e_4 * (2 - 2 * e_fact_sqrt - 8 * kappa2 +
                         6 * e_fact_sqrt * kappa2) +
                  e_2 * (-4 + 3 * e_fact_sqrt +
                         4 * (4 + 7 * e_fact_sqrt) * kappa2)) *
                 S2z * S2z) +
        e *
            (2 *
                 (12 + 68 * e_fact_sqrt + 3 * (7 + 36 * e_fact_sqrt) * kappa1 +
                  3 * e_4 * (4 + 7 * kappa1) +
                  3 * e_2 *
                      (-8 + 12 * e_fact_sqrt - 14 * kappa1 +
                       39 * e_fact_sqrt * kappa1)) *
                 m1 * m1 * m1 * m1 * S1z * S1z +
             2 *
                 (12 + 68 * e_fact_sqrt + 3 * (7 + 36 * e_fact_sqrt) * kappa2 +
                  3 * e_4 * (4 + 7 * kappa2) +
                  3 * e_2 *
                      (-8 + 12 * e_fact_sqrt - 14 * kappa2 +
                       39 * e_fact_sqrt * kappa2)) *
                 m2 * m2 * m2 * m2 * S2z * S2z +
             3 * m1 * m1 * m1 * m2 * S1z *
                 (20 * e_fact_sqrt * S1z + 33 * e_2 * e_fact_sqrt * S1z +
                  110 * e_fact_sqrt * kappa1 * S1z +
                  121 * e_2 * e_fact_sqrt * kappa1 * S1z +
                  pow(-1 + e_2, 2) * (3 + 23 * kappa1) * S1z +
                  152 * e_fact_sqrt * S2z + 186 * e_2 * e_fact_sqrt * S2z +
                  30 * pow(-1 + e_2, 2) * S2z) +
             3 * m1 * m2 * m2 * m2 * S2z *
                 (152 * e_fact_sqrt * S1z + 186 * e_2 * e_fact_sqrt * S1z +
                  30 * pow(-1 + e_2, 2) * S1z + 20 * e_fact_sqrt * S2z +
                  33 * e_2 * e_fact_sqrt * S2z +
                  110 * e_fact_sqrt * kappa2 * S2z +
                  121 * e_2 * e_fact_sqrt * kappa2 * S2z +
                  pow(-1 + e_2, 2) * (3 + 23 * kappa2) * S2z) +
             m1 * m1 * m2 * m2 *
                 (9 *
                      (e_4 * (-1 + 4 * kappa1) +
                       (1 + 4 * e_fact_sqrt) * (-1 + 4 * kappa1) +
                       e_2 * (2 + 3 * e_fact_sqrt - 8 * kappa1 +
                              24 * e_fact_sqrt * kappa1)) *
                      S1z * S1z +
                  2 *
                      (87 + 87 * e_4 + 466 * e_fact_sqrt +
                       3 * e_2 * (-58 + 157 * e_fact_sqrt)) *
                      S1z * S2z +
                  9 *
                      (e_4 * (-1 + 4 * kappa2) +
                       (1 + 4 * e_fact_sqrt) * (-1 + 4 * kappa2) +
                       e_2 * (2 + 3 * e_fact_sqrt - 8 * kappa2 +
                              24 * e_fact_sqrt * kappa2)) *
                      S2z * S2z)) *
            cos(u)) /
       (e_fact_35 * M_fact_4));

  return (r_3pn_SS);
}

static REAL8 separation(REAL8 u, REAL8 eta, REAL8 x, REAL8 e, REAL8 m1,
                        REAL8 m2, REAL8 S1z, REAL8 S2z) {
  // 3PN accurate

  return ((1.0 / x) * rel_sep_0pn(e, u) + rel_sep_1pn(e, u, eta) +
          rel_sep_1_5pn(e, u, m1, m2, S1z, S2z) * sqrt(x) +
          rel_sep_2pnSS(e, u, m1, m2, S1z, S2z) * x +
          rel_sep_2pn(e, u, eta) * x +
          rel_sep_2_5pn_SO(e, u, m1, m2, S1z, S2z) * x * sqrt(x) +
          rel_sep_3pn(e, u, eta) * x * x +
          rel_sep_3pn_SS(e, u, m1, m2, S1z, S2z) * x * x);
}

static REAL8 dxdt_4pn(REAL8 x, REAL8 eta) {
  const REAL8 log2 = 0.693147180559945309417232121458; // ln(2)
  const REAL8 log3 = 1.09861228866810969139524523692;  // ln(3)

  REAL8 x_2 = x * x;

  REAL8 x_3 = x_2 * x;

  REAL8 x_4 = x_3 * x;

  REAL8 x_5 = x_4 * x;

  REAL8 eta_pow_2 = eta * eta;
  REAL8 eta_pow_3 = eta_pow_2 * eta;
  REAL8 eta_pow_4 = eta_pow_2 * eta_pow_2;

  REAL8 pi_pow_2 = M_PI * M_PI;

  REAL8 euler = LAL_GAMMA;

  REAL8 pre_fact = 64. * x_5 * eta / 5.;

  REAL8 bit_4_pn =
      (x_4 *
       (-((-891 + 477 * eta - 11 * eta_pow_2) *
          (-4.928461199294532 + (9271 * eta) / 504. + (65 * eta_pow_2) / 18.)) /
            36. +
        (5 * (-3.7113095238095237 - (35 * eta) / 12.) *
         (19683 - 66375 * eta + 1188 * eta_pow_2 + 19 * eta_pow_3 +
          2214 * eta * pi_pow_2)) /
            648. -
        (-3 - eta / 3.) * (95.10839000836025 - (134543 * eta) / 7776. -
                           (94403 * eta_pow_2) / 3024. -
                           (775 * eta_pow_3) / 324. - (1712 * euler) / 105. +
                           (16 * pi_pow_2) / 3. + (41 * eta * pi_pow_2) / 48. -
                           (3424 * log2) / 105. - (856 * log(x)) / 105.) +
        2 * (-101.65745990813167 + (232597 * euler) / 4410. -
             (1369 * pi_pow_2) / 126. + (39931 * log2) / 294. -
             (47385 * log3) / 1568. + (232597 * log(x)) / 8820.) +
        0.071630658436214 *
            (12774.514362657092 - 39782.929590804066 * eta +
             346.61235705608044 * eta_pow_2 + 2.068222621184919 * eta_pow_3 +
             eta_pow_4 - 4169.536804308797 * eta * log(x)))) /
      2.;

  REAL8 at_4_pn;

  at_4_pn = pre_fact * (bit_4_pn);

  return (at_4_pn);
}

static REAL8 dxdt_4_5pn(REAL8 x, REAL8 eta) {
  const REAL8 log2 = 0.693147180559945309417232121458; // ln(2)
  const REAL8 log3 = 1.09861228866810969139524523692;  // ln(3)

  REAL8 x_1_2 = sqrt(x);

  REAL8 x_2 = x * x;

  REAL8 x_3 = x_2 * x;

  REAL8 x_4 = x_3 * x;
  REAL8 x_9_2 = x_4 * x_1_2;
  REAL8 x_5 = x_4 * x;

  REAL8 eta_pow_2 = eta * eta;
  REAL8 eta_pow_3 = eta_pow_2 * eta;
  REAL8 eta_pow_4 = eta_pow_2 * eta_pow_2;
  REAL8 pi_pow_2 = M_PI * M_PI;

  REAL8 euler = LAL_GAMMA;

  REAL8 pre_fact = 64. * x_5 * eta / 5.;

  REAL8 bit_4_pn =
      (x_4 *
       (-((-891 + 477 * eta - 11 * eta_pow_2) *
          (-4.928461199294532 + (9271 * eta) / 504. + (65 * eta_pow_2) / 18.)) /
            36. +
        (5 * (-3.7113095238095237 - (35 * eta) / 12.) *
         (19683 - 66375 * eta + 1188 * eta_pow_2 + 19 * eta_pow_3 +
          2214 * eta * pi_pow_2)) /
            648. -
        (-3 - eta / 3.) * (95.10839000836025 - (134543 * eta) / 7776. -
                           (94403 * eta_pow_2) / 3024. -
                           (775 * eta_pow_3) / 324. - (1712 * euler) / 105. +
                           (16 * pi_pow_2) / 3. + (41 * eta * pi_pow_2) / 48. -
                           (3424 * log2) / 105. - (856 * log(x)) / 105.) +
        2 * (-101.65745990813167 + (232597 * euler) / 4410. -
             (1369 * pi_pow_2) / 126. + (39931 * log2) / 294. -
             (47385 * log3) / 1568. + (232597 * log(x)) / 8820.) +
        0.071630658436214 *
            (12774.514362657092 - 39782.929590804066 * eta +
             346.61235705608044 * eta_pow_2 + 2.068222621184919 * eta_pow_3 +
             eta_pow_4 - 4169.536804308797 * eta * log(x)))) /
      2.;

  REAL8 bit_at_9_2_pn =
      (x_9_2 *
       (-((-891 + 477 * eta - 11 * eta_pow_2) *
          ((-8191 * M_PI) / 672. - (583 * eta * M_PI) / 24.)) /
            36. -
        (-3 - eta / 3.) *
            ((-16285 * M_PI) / 504. + (214745 * eta * M_PI) / 1728. +
             (193385 * eta_pow_2 * M_PI) / 3024.) +
        (5 * M_PI *
         (19683 - 66375 * eta + 1188 * eta_pow_2 + 19 * eta_pow_3 +
          2214 * eta * pi_pow_2)) /
            162. +
        2 * ((265978667519 * M_PI) / 7.451136e8 - (6848 * euler * M_PI) / 105. -
             (13696 * M_PI * log2) / 105. - (3424 * M_PI * log(x)) / 105.))) /
      2.;

  REAL8 at_4_5_pn;

  at_4_5_pn = pre_fact * (bit_4_pn + bit_at_9_2_pn);

  return (at_4_5_pn);
}

static REAL8 dxdt_5pn(REAL8 x, REAL8 eta) {
  const REAL8 log2 = 0.693147180559945309417232121458; // ln(2)
  const REAL8 log3 = 1.09861228866810969139524523692;  // ln(3)

  REAL8 x_1_2 = sqrt(x);

  REAL8 x_2 = x * x;

  REAL8 x_3 = x_2 * x;

  REAL8 x_4 = x_3 * x;
  REAL8 x_9_2 = x_4 * x_1_2;
  REAL8 x_5 = x_4 * x;

  REAL8 eta_pow_2 = eta * eta;
  REAL8 eta_pow_3 = eta_pow_2 * eta;
  REAL8 eta_pow_4 = eta_pow_2 * eta_pow_2;
  REAL8 eta_pow_5 = eta_pow_3 * eta_pow_2;

  REAL8 pi_pow_2 = M_PI * M_PI;

  REAL8 euler = LAL_GAMMA;

  REAL8 pre_fact = 64. * x_5 * eta / 5.;

  REAL8 bit_4_pn =
      (x_4 *
       (-((-891 + 477 * eta - 11 * eta_pow_2) *
          (-4.928461199294532 + (9271 * eta) / 504. + (65 * eta_pow_2) / 18.)) /
            36. +
        (5 * (-3.7113095238095237 - (35 * eta) / 12.) *
         (19683 - 66375 * eta + 1188 * eta_pow_2 + 19 * eta_pow_3 +
          2214 * eta * pi_pow_2)) /
            648. -
        (-3 - eta / 3.) * (95.10839000836025 - (134543 * eta) / 7776. -
                           (94403 * eta_pow_2) / 3024. -
                           (775 * eta_pow_3) / 324. - (1712 * euler) / 105. +
                           (16 * pi_pow_2) / 3. + (41 * eta * pi_pow_2) / 48. -
                           (3424 * log2) / 105. - (856 * log(x)) / 105.) +
        2 * (-101.65745990813167 + (232597 * euler) / 4410. -
             (1369 * pi_pow_2) / 126. + (39931 * log2) / 294. -
             (47385 * log3) / 1568. + (232597 * log(x)) / 8820.) +
        0.071630658436214 *
            (12774.514362657092 - 39782.929590804066 * eta +
             346.61235705608044 * eta_pow_2 + 2.068222621184919 * eta_pow_3 +
             eta_pow_4 - 4169.536804308797 * eta * log(x)))) /
      2.;

  REAL8 bit_at_9_2_pn =
      (x_9_2 *
       (-((-891 + 477 * eta - 11 * eta_pow_2) *
          ((-8191 * M_PI) / 672. - (583 * eta * M_PI) / 24.)) /
            36. -
        (-3 - eta / 3.) *
            ((-16285 * M_PI) / 504. + (214745 * eta * M_PI) / 1728. +
             (193385 * eta_pow_2 * M_PI) / 3024.) +
        (5 * M_PI *
         (19683 - 66375 * eta + 1188 * eta_pow_2 + 19 * eta_pow_3 +
          2214 * eta * pi_pow_2)) /
            162. +
        2 * ((265978667519 * M_PI) / 7.451136e8 - (6848 * euler * M_PI) / 105. -
             (13696 * M_PI * log2) / 105. - (3424 * M_PI * log(x)) / 105.))) /
      2.;

  REAL8 bit_at_5_pn =
      (x_5 *
       ((5 *
         (-4.928461199294532 + (9271 * eta) / 504. + (65 * eta_pow_2) / 18.) *
         (19683 - 66375 * eta + 1188 * eta_pow_2 + 19 * eta_pow_3 +
          2214 * eta * pi_pow_2)) /
            648. -
        ((-891 + 477 * eta - 11 * eta_pow_2) *
         (95.10839000836025 - (134543 * eta) / 7776. -
          (94403 * eta_pow_2) / 3024. - (775 * eta_pow_3) / 324. -
          (1712 * euler) / 105. + (16 * pi_pow_2) / 3. +
          (41 * eta * pi_pow_2) / 48. - (3424 * log2) / 105. -
          (856 * log(x)) / 105.)) /
            36. -
        (-3 - eta / 3.) * (-101.65745990813167 + (232597 * euler) / 4410. -
                           (1369 * pi_pow_2) / 126. + (39931 * log2) / 294. -
                           (47385 * log3) / 1568. + (232597 * log(x)) / 8820.) +
        2 * (-883.093730029416 + (916628467 * euler) / 7.85862e6 -
             (424223 * pi_pow_2) / 6804. - (83217611 * log2) / 1.12266e6 +
             (47385 * log3) / 196. + (916628467 * log(x)) / 1.571724e7) +
        0.071630658436214 * (-3.7113095238095237 - (35 * eta) / 12.) *
            (12774.514362657092 - 39782.929590804066 * eta +
             346.61235705608044 * eta_pow_2 + 2.068222621184919 * eta_pow_3 +
             eta_pow_4 - 4169.536804308797 * eta * log(x)) +
        0.03851594650205761 *
            (142693.53505843072 - 375812.922000909 * eta +
             112536.48031329196 * eta_pow_2 - 6092.1750154249185 * eta_pow_3 +
             48.00500834724541 * eta_pow_4 + eta_pow_5 +
             21138.554352492254 * eta * log(x) +
             36962.46811352254 * eta_pow_2 * log(x)))) /
      2.;

  REAL8 at_5_pn;

  at_5_pn = pre_fact * (bit_4_pn + bit_at_9_2_pn + bit_at_5_pn);

  return (at_5_pn);
}

static REAL8 dxdt_5_5pn(REAL8 x, REAL8 eta) {
  const REAL8 log2 = 0.693147180559945309417232121458; // ln(2)
  const REAL8 log3 = 1.09861228866810969139524523692;  // ln(3)

  REAL8 x_1_2 = sqrt(x);

  REAL8 x_2 = x * x;

  REAL8 x_3 = x_2 * x;

  REAL8 x_4 = x_3 * x;
  REAL8 x_9_2 = x_4 * x_1_2;
  REAL8 x_5 = x_4 * x;
  REAL8 x_11_2 = x_5 * x_1_2;

  REAL8 eta_pow_2 = eta * eta;
  REAL8 eta_pow_3 = eta_pow_2 * eta;
  REAL8 eta_pow_4 = eta_pow_2 * eta_pow_2;
  REAL8 eta_pow_5 = eta_pow_3 * eta_pow_2;

  REAL8 pi_pow_2 = M_PI * M_PI;

  REAL8 euler = LAL_GAMMA;

  REAL8 pre_fact = 64. * x_5 * eta / 5.;

  REAL8 bit_4_pn =
      (x_4 *
       (-((-891 + 477 * eta - 11 * eta_pow_2) *
          (-4.928461199294532 + (9271 * eta) / 504. + (65 * eta_pow_2) / 18.)) /
            36. +
        (5 * (-3.7113095238095237 - (35 * eta) / 12.) *
         (19683 - 66375 * eta + 1188 * eta_pow_2 + 19 * eta_pow_3 +
          2214 * eta * pi_pow_2)) /
            648. -
        (-3 - eta / 3.) * (95.10839000836025 - (134543 * eta) / 7776. -
                           (94403 * eta_pow_2) / 3024. -
                           (775 * eta_pow_3) / 324. - (1712 * euler) / 105. +
                           (16 * pi_pow_2) / 3. + (41 * eta * pi_pow_2) / 48. -
                           (3424 * log2) / 105. - (856 * log(x)) / 105.) +
        2 * (-101.65745990813167 + (232597 * euler) / 4410. -
             (1369 * pi_pow_2) / 126. + (39931 * log2) / 294. -
             (47385 * log3) / 1568. + (232597 * log(x)) / 8820.) +
        0.071630658436214 *
            (12774.514362657092 - 39782.929590804066 * eta +
             346.61235705608044 * eta_pow_2 + 2.068222621184919 * eta_pow_3 +
             eta_pow_4 - 4169.536804308797 * eta * log(x)))) /
      2.;

  REAL8 bit_at_9_2_pn =
      (x_9_2 *
       (-((-891 + 477 * eta - 11 * eta_pow_2) *
          ((-8191 * M_PI) / 672. - (583 * eta * M_PI) / 24.)) /
            36. -
        (-3 - eta / 3.) *
            ((-16285 * M_PI) / 504. + (214745 * eta * M_PI) / 1728. +
             (193385 * eta_pow_2 * M_PI) / 3024.) +
        (5 * M_PI *
         (19683 - 66375 * eta + 1188 * eta_pow_2 + 19 * eta_pow_3 +
          2214 * eta * pi_pow_2)) /
            162. +
        2 * ((265978667519 * M_PI) / 7.451136e8 - (6848 * euler * M_PI) / 105. -
             (13696 * M_PI * log2) / 105. - (3424 * M_PI * log(x)) / 105.))) /
      2.;

  REAL8 bit_at_5_pn =
      (x_5 *
       ((5 *
         (-4.928461199294532 + (9271 * eta) / 504. + (65 * eta_pow_2) / 18.) *
         (19683 - 66375 * eta + 1188 * eta_pow_2 + 19 * eta_pow_3 +
          2214 * eta * pi_pow_2)) /
            648. -
        ((-891 + 477 * eta - 11 * eta_pow_2) *
         (95.10839000836025 - (134543 * eta) / 7776. -
          (94403 * eta_pow_2) / 3024. - (775 * eta_pow_3) / 324. -
          (1712 * euler) / 105. + (16 * pi_pow_2) / 3. +
          (41 * eta * pi_pow_2) / 48. - (3424 * log2) / 105. -
          (856 * log(x)) / 105.)) /
            36. -
        (-3 - eta / 3.) * (-101.65745990813167 + (232597 * euler) / 4410. -
                           (1369 * pi_pow_2) / 126. + (39931 * log2) / 294. -
                           (47385 * log3) / 1568. + (232597 * log(x)) / 8820.) +
        2 * (-883.093730029416 + (916628467 * euler) / 7.85862e6 -
             (424223 * pi_pow_2) / 6804. - (83217611 * log2) / 1.12266e6 +
             (47385 * log3) / 196. + (916628467 * log(x)) / 1.571724e7) +
        0.071630658436214 * (-3.7113095238095237 - (35 * eta) / 12.) *
            (12774.514362657092 - 39782.929590804066 * eta +
             346.61235705608044 * eta_pow_2 + 2.068222621184919 * eta_pow_3 +
             eta_pow_4 - 4169.536804308797 * eta * log(x)) +
        0.03851594650205761 *
            (142693.53505843072 - 375812.922000909 * eta +
             112536.48031329196 * eta_pow_2 - 6092.1750154249185 * eta_pow_3 +
             48.00500834724541 * eta_pow_4 + eta_pow_5 +
             21138.554352492254 * eta * log(x) +
             36962.46811352254 * eta_pow_2 * log(x)))) /
      2.;

  REAL8 bit_at_11_2_pn =
      (x_11_2 *
       (-((-891 + 477 * eta - 11 * eta_pow_2) *
          ((-16285 * M_PI) / 504. + (214745 * eta * M_PI) / 1728. +
           (193385 * eta_pow_2 * M_PI) / 3024.)) /
            36. +
        (5 * ((-8191 * M_PI) / 672. - (583 * eta * M_PI) / 24.) *
         (19683 - 66375 * eta + 1188 * eta_pow_2 + 19 * eta_pow_3 +
          2214 * eta * pi_pow_2)) /
            648. +
        0.9001374012600385 *
            (12774.514362657092 - 39782.929590804066 * eta +
             346.61235705608044 * eta_pow_2 + 2.068222621184919 * eta_pow_3 +
             eta_pow_4 - 4169.536804308797 * eta * log(x)) -
        (-3 - eta / 3.) *
            ((265978667519 * M_PI) / 7.451136e8 - (6848 * euler * M_PI) / 105. -
             (13696 * M_PI * log2) / 105. - (3424 * M_PI * log(x)) / 105.) +
        2 * ((8399309750401 * M_PI) / 1.017080064e11 +
             (177293 * euler * M_PI) / 1176. +
             (8521283 * M_PI * log2) / 17640. - (142155 * M_PI * log3) / 784. +
             (177293 * M_PI * log(x)) / 2352.))) /
      2.;

  REAL8 at_5_5_pn;

  at_5_5_pn =
      pre_fact * (bit_4_pn + bit_at_9_2_pn + bit_at_5_pn + bit_at_11_2_pn);

  return (at_5_5_pn);
}

static REAL8 dxdt_6pn(REAL8 x, REAL8 eta) {
  const REAL8 log2 = 0.693147180559945309417232121458; // ln(2)
  const REAL8 log3 = 1.09861228866810969139524523692;  // ln(3)
  const REAL8 log5 = 1.60943791243410037460075933323;  // ln(5)

  REAL8 x_1_2 = sqrt(x);

  REAL8 x_2 = x * x;

  REAL8 x_3 = x_2 * x;

  REAL8 x_4 = x_3 * x;
  REAL8 x_9_2 = x_4 * x_1_2;
  REAL8 x_5 = x_4 * x;
  REAL8 x_11_2 = x_5 * x_1_2;
  REAL8 x_6 = x_5 * x;

  REAL8 eta_pow_2 = eta * eta;
  REAL8 eta_pow_3 = eta_pow_2 * eta;
  REAL8 eta_pow_4 = eta_pow_2 * eta_pow_2;
  REAL8 eta_pow_5 = eta_pow_3 * eta_pow_2;
  REAL8 eta_pow_6 = eta_pow_4 * eta_pow_2;
  REAL8 log_x_2 = log(x) * log(x);
  REAL8 zeta_3 = 1.2020569031595942;
  REAL8 log_2_2 = log2 * log2;

  REAL8 pi_pow_2 = M_PI * M_PI;
  REAL8 pi_pow_4 = pi_pow_2 * pi_pow_2;
  REAL8 euler = LAL_GAMMA;
  REAL8 euler_2 = euler * euler;
  REAL8 pre_fact = 64. * x_5 * eta / 5.;

  REAL8 bit_4_pn =
      (x_4 *
       (-((-891 + 477 * eta - 11 * eta_pow_2) *
          (-4.928461199294532 + (9271 * eta) / 504. + (65 * eta_pow_2) / 18.)) /
            36. +
        (5 * (-3.7113095238095237 - (35 * eta) / 12.) *
         (19683 - 66375 * eta + 1188 * eta_pow_2 + 19 * eta_pow_3 +
          2214 * eta * pi_pow_2)) /
            648. -
        (-3 - eta / 3.) * (95.10839000836025 - (134543 * eta) / 7776. -
                           (94403 * eta_pow_2) / 3024. -
                           (775 * eta_pow_3) / 324. - (1712 * euler) / 105. +
                           (16 * pi_pow_2) / 3. + (41 * eta * pi_pow_2) / 48. -
                           (3424 * log2) / 105. - (856 * log(x)) / 105.) +
        2 * (-101.65745990813167 + (232597 * euler) / 4410. -
             (1369 * pi_pow_2) / 126. + (39931 * log2) / 294. -
             (47385 * log3) / 1568. + (232597 * log(x)) / 8820.) +
        0.071630658436214 *
            (12774.514362657092 - 39782.929590804066 * eta +
             346.61235705608044 * eta_pow_2 + 2.068222621184919 * eta_pow_3 +
             eta_pow_4 - 4169.536804308797 * eta * log(x)))) /
      2.;

  REAL8 bit_at_9_2_pn =
      (x_9_2 *
       (-((-891 + 477 * eta - 11 * eta_pow_2) *
          ((-8191 * M_PI) / 672. - (583 * eta * M_PI) / 24.)) /
            36. -
        (-3 - eta / 3.) *
            ((-16285 * M_PI) / 504. + (214745 * eta * M_PI) / 1728. +
             (193385 * eta_pow_2 * M_PI) / 3024.) +
        (5 * M_PI *
         (19683 - 66375 * eta + 1188 * eta_pow_2 + 19 * eta_pow_3 +
          2214 * eta * pi_pow_2)) /
            162. +
        2 * ((265978667519 * M_PI) / 7.451136e8 - (6848 * euler * M_PI) / 105. -
             (13696 * M_PI * log2) / 105. - (3424 * M_PI * log(x)) / 105.))) /
      2.;

  REAL8 bit_at_5_pn =
      (x_5 *
       ((5 *
         (-4.928461199294532 + (9271 * eta) / 504. + (65 * eta_pow_2) / 18.) *
         (19683 - 66375 * eta + 1188 * eta_pow_2 + 19 * eta_pow_3 +
          2214 * eta * pi_pow_2)) /
            648. -
        ((-891 + 477 * eta - 11 * eta_pow_2) *
         (95.10839000836025 - (134543 * eta) / 7776. -
          (94403 * eta_pow_2) / 3024. - (775 * eta_pow_3) / 324. -
          (1712 * euler) / 105. + (16 * pi_pow_2) / 3. +
          (41 * eta * pi_pow_2) / 48. - (3424 * log2) / 105. -
          (856 * log(x)) / 105.)) /
            36. -
        (-3 - eta / 3.) * (-101.65745990813167 + (232597 * euler) / 4410. -
                           (1369 * pi_pow_2) / 126. + (39931 * log2) / 294. -
                           (47385 * log3) / 1568. + (232597 * log(x)) / 8820.) +
        2 * (-883.093730029416 + (916628467 * euler) / 7.85862e6 -
             (424223 * pi_pow_2) / 6804. - (83217611 * log2) / 1.12266e6 +
             (47385 * log3) / 196. + (916628467 * log(x)) / 1.571724e7) +
        0.071630658436214 * (-3.7113095238095237 - (35 * eta) / 12.) *
            (12774.514362657092 - 39782.929590804066 * eta +
             346.61235705608044 * eta_pow_2 + 2.068222621184919 * eta_pow_3 +
             eta_pow_4 - 4169.536804308797 * eta * log(x)) +
        0.03851594650205761 *
            (142693.53505843072 - 375812.922000909 * eta +
             112536.48031329196 * eta_pow_2 - 6092.1750154249185 * eta_pow_3 +
             48.00500834724541 * eta_pow_4 + eta_pow_5 +
             21138.554352492254 * eta * log(x) +
             36962.46811352254 * eta_pow_2 * log(x)))) /
      2.;

  REAL8 bit_at_11_2_pn =
      (x_11_2 *
       (-((-891 + 477 * eta - 11 * eta_pow_2) *
          ((-16285 * M_PI) / 504. + (214745 * eta * M_PI) / 1728. +
           (193385 * eta_pow_2 * M_PI) / 3024.)) /
            36. +
        (5 * ((-8191 * M_PI) / 672. - (583 * eta * M_PI) / 24.) *
         (19683 - 66375 * eta + 1188 * eta_pow_2 + 19 * eta_pow_3 +
          2214 * eta * pi_pow_2)) /
            648. +
        0.9001374012600385 *
            (12774.514362657092 - 39782.929590804066 * eta +
             346.61235705608044 * eta_pow_2 + 2.068222621184919 * eta_pow_3 +
             eta_pow_4 - 4169.536804308797 * eta * log(x)) -
        (-3 - eta / 3.) *
            ((265978667519 * M_PI) / 7.451136e8 - (6848 * euler * M_PI) / 105. -
             (13696 * M_PI * log2) / 105. - (3424 * M_PI * log(x)) / 105.) +
        2 * ((8399309750401 * M_PI) / 1.017080064e11 +
             (177293 * euler * M_PI) / 1176. +
             (8521283 * M_PI * log2) / 17640. - (142155 * M_PI * log3) / 784. +
             (177293 * M_PI * log(x)) / 2352.))) /
      2.;

  REAL8 bit_at_6_pn =
      (x_6 *
       ((5 *
         (19683 - 66375 * eta + 1188 * eta_pow_2 + 19 * eta_pow_3 +
          2214 * eta * pi_pow_2) *
         (95.10839000836025 - (134543 * eta) / 7776. -
          (94403 * eta_pow_2) / 3024. - (775 * eta_pow_3) / 324. -
          (1712 * euler) / 105. + (16 * pi_pow_2) / 3. +
          (41 * eta * pi_pow_2) / 48. - (3424 * log2) / 105. -
          (856 * log(x)) / 105.)) /
            648. -
        ((-891 + 477 * eta - 11 * eta_pow_2) *
         (-101.65745990813167 + (232597 * euler) / 4410. -
          (1369 * pi_pow_2) / 126. + (39931 * log2) / 294. -
          (47385 * log3) / 1568. + (232597 * log(x)) / 8820.)) /
            36. -
        (-3 - eta / 3.) *
            (-883.093730029416 + (916628467 * euler) / 7.85862e6 -
             (424223 * pi_pow_2) / 6804. - (83217611 * log2) / 1.12266e6 +
             (47385 * log3) / 196. + (916628467 * log(x)) / 1.571724e7) +
        0.071630658436214 *
            (-4.928461199294532 + (9271 * eta) / 504. +
             (65 * eta_pow_2) / 18.) *
            (12774.514362657092 - 39782.929590804066 * eta +
             346.61235705608044 * eta_pow_2 + 2.068222621184919 * eta_pow_3 +
             eta_pow_4 - 4169.536804308797 * eta * log(x)) +
        0.03851594650205761 * (-3.7113095238095237 - (35 * eta) / 12.) *
            (142693.53505843072 - 375812.922000909 * eta +
             112536.48031329196 * eta_pow_2 - 6092.1750154249185 * eta_pow_3 +
             48.00500834724541 * eta_pow_4 + eta_pow_5 +
             21138.554352492254 * eta * log(x) +
             36962.46811352254 * eta_pow_2 * log(x)) +
        0.01933239502362445 *
            (1.706246232376582e6 - 6.226251580758891e6 * eta +
             4.682175752835639e6 * eta_pow_2 - 206463.70780984825 * eta_pow_3 -
             328.2772748400016 * eta_pow_4 + 55.268054571771735 * eta_pow_5 +
             eta_pow_6 + 676714.6165659907 * eta * log(x) +
             462832.31497820036 * eta_pow_2 * log(x) +
             21113.668393335593 * eta_pow_3 * log(x)) +
        2 * (3432.3197889543108 - (246137536815857 * euler) / 1.573295724e11 +
             (1465472 * euler_2) / 11025. +
             (3803225263 * pi_pow_2) / 1.047816e7 -
             (27392 * euler * pi_pow_2) / 315. - (256 * pi_pow_4) / 45. -
             (271272899815409 * log2) / 1.573295724e11 +
             (5861888 * euler * log2) / 11025. -
             (54784 * pi_pow_2 * log2) / 315. + (5861888 * log_2_2) / 11025. -
             (437114506833 * log3) / 7.8926848e8 -
             (37744140625 * log5) / 2.60941824e8 -
             (246137536815857 * log(x)) / 3.146591448e11 +
             (1465472 * euler * log(x)) / 11025. -
             (13696 * pi_pow_2 * log(x)) / 315. +
             (2930944 * log2 * log(x)) / 11025. + (366368 * log_x_2) / 11025. -
             (27392 * zeta_3) / 105.))) /
      2.;

  REAL8 at_6_pn;

  at_6_pn = pre_fact * (bit_4_pn + bit_at_9_2_pn + bit_at_5_pn +
                        bit_at_11_2_pn + bit_at_6_pn);

  return (at_6_pn);
}

/* Hereditary terms are obtained from the article
   http://arxiv.org/pdf/0908.3854v2.pdf. The
   analytical terms presented below reproduce
   the numerical data presented therein to better
   that 0.1% */

/* Eccentricity enhancement factors */

static REAL8 phi_e(REAL8 e) {
  REAL8 e_pow_2 = e * e;
  REAL8 e_fact = 1. - e_pow_2;
  REAL8 e_pow_4 = e_pow_2 * e_pow_2;
  REAL8 e_pow_6 = e_pow_4 * e_pow_2;
  REAL8 e_pow_8 = e_pow_6 * e_pow_2;
  REAL8 e_pow_10 = e_pow_8 * e_pow_2;
  REAL8 e_pow_12 = e_pow_10 * e_pow_2;

  REAL8 num =
      1.0 + (18970894028. * e_pow_2) / 2649026657. +
      (157473274. * e_pow_4) / 30734301. + (48176523. * e_pow_6) / 177473701. +
      (9293260. * e_pow_8) / 3542508891. - (5034498. * e_pow_10) / 7491716851. +
      (428340. * e_pow_12) / 9958749469.;

  REAL8 den = e_fact * e_fact * e_fact * e_fact * e_fact;

  return (num / den);
}

static REAL8 psi_e(REAL8 e) {

  REAL8 e_pow_2 = e * e;
  REAL8 e_fact = 1. - e_pow_2;
  REAL8 e_pow_4 = e_pow_2 * e_pow_2;
  REAL8 e_pow_6 = e_pow_4 * e_pow_2;
  REAL8 num = 1. - (185. * e_pow_2) / 21. - (3733. * e_pow_4) / 99. -
              (1423. * e_pow_6) / 104.;
  REAL8 den = e_fact * e_fact * e_fact * e_fact * e_fact * e_fact;

  return (num / den);
}

static REAL8 zed_e(REAL8 e) {
  REAL8 e_pow_2 = e * e;
  REAL8 e_fact = 1. - e_pow_2;
  REAL8 e_pow_4 = e_pow_2 * e_pow_2;
  REAL8 e_pow_6 = e_pow_4 * e_pow_2;

  REAL8 num = 1. + (2095. * e_pow_2) / 143. + (1590. * e_pow_4) / 59. +
              (977. * e_pow_6) / 113.;
  REAL8 den = e_fact * e_fact * e_fact * e_fact * e_fact * e_fact;

  return (num / den);
}

static REAL8 kappa_e(REAL8 e) {

  REAL8 e_pow_2 = e * e;
  REAL8 e_fact = 1. - e_pow_2;
  REAL8 e_pow_4 = e_pow_2 * e_pow_2;
  REAL8 e_pow_6 = e_pow_4 * e_pow_2;
  REAL8 e_pow_8 = e_pow_6 * e_pow_2;
  REAL8 e_pow_10 = e_pow_8 * e_pow_2;
  REAL8 num = 1. + (1497. * e_pow_2) / 79. + (7021. * e_pow_4) / 143. +
              (997. * e_pow_6) / 98. + (463. * e_pow_8) / 51. -
              (3829. * e_pow_10) / 120.;
  REAL8 den = e_fact * e_fact * e_fact * e_fact * e_fact * e_fact * e_fact;

  return (num / den);
}

static REAL8 phi_e_tilde(REAL8 e) {
  REAL8 e_pow_2 = e * e;
  REAL8 e_fact = 1. - e_pow_2;
  REAL8 e_pow_4 = e_pow_2 * e_pow_2;
  REAL8 e_pow_6 = e_pow_4 * e_pow_2;
  REAL8 e_pow_8 = e_pow_6 * e_pow_2;
  REAL8 e_pow_10 = e_pow_8 * e_pow_2;

  REAL8 num =
      1. + (413137256. * e_pow_2) / 136292703. +
      (37570495. * e_pow_4) / 98143337. - (2640201. * e_pow_6) / 993226448. -
      (4679700. * e_pow_8) / 6316712563. - (328675. * e_pow_10) / 8674876481.;
  REAL8 den = e_fact * e_fact * e_fact * sqrt(e_fact);

  return (num / den);
}

static REAL8 psi_e_tilde(REAL8 e) {
  REAL8 e_pow_2 = e * e;
  REAL8 e_fact = 1. - e_pow_2;
  REAL8 e_pow_4 = e_pow_2 * e_pow_2;
  REAL8 e_pow_6 = e_pow_4 * e_pow_2;
  REAL8 e_pow_8 = e_pow_6 * e_pow_2;
  REAL8 e_pow_10 = e_pow_8 * e_pow_2;
  REAL8 num = 1. - (2022. * e_pow_2) / 305. - (249. * e_pow_4) / 26. -
              (193. * e_pow_6) / 239. + (23. * e_pow_8) / 43. -
              (102. * e_pow_10) / 463.;
  REAL8 den = e_fact * e_fact * e_fact * e_fact * sqrt(e_fact);

  return (num / den);
}

static REAL8 zed_e_tilde(REAL8 e) {
  REAL8 e_pow_2 = e * e;
  REAL8 e_fact = 1. - e_pow_2;
  REAL8 e_pow_4 = e_pow_2 * e_pow_2;
  REAL8 e_pow_6 = e_pow_4 * e_pow_2;
  REAL8 e_pow_8 = e_pow_6 * e_pow_2;
  REAL8 num = 1. + (1563. * e_pow_2) / 194. + (1142. * e_pow_4) / 193. +
              (123. * e_pow_6) / 281. - (27. * e_pow_8) / 328.;
  REAL8 den = e_fact * e_fact * e_fact * e_fact * sqrt(e_fact);

  return (num / den);
}

static REAL8 kappa_e_tilde(REAL8 e) {
  REAL8 e_pow_2 = e * e;
  REAL8 e_fact = 1. - e_pow_2;
  REAL8 e_pow_4 = e_pow_2 * e_pow_2;
  REAL8 e_pow_6 = e_pow_4 * e_pow_2;
  REAL8 e_pow_8 = e_pow_6 * e_pow_2;
  REAL8 e_pow_10 = e_pow_8 * e_pow_2;
  REAL8 num = 1. + (1789. * e_pow_2) / 167. + (5391. * e_pow_4) / 340. +
              (2150. * e_pow_6) / 219. - (1007. * e_pow_8) / 320. +
              (2588. * e_pow_10) / 189.;
  REAL8 den = e_fact * e_fact * e_fact * e_fact * e_fact;

  return (num / den);
}

static REAL8 phi_e_rad(REAL8 e) {
  REAL8 e_fact = 1. - e * e;
  REAL8 e_pow_2 = e * e;
  REAL8 pre_fact = 192. * sqrt(e_fact) / (985. * e_pow_2);

  return (pre_fact * (sqrt(e_fact) * phi_e(e) - phi_e_tilde(e)));
}

static REAL8 psi_e_rad(REAL8 e) {
  REAL8 e_fact = 1. - e * e;
  REAL8 e_pow_2 = e * e;
  REAL8 pre_f_1 = 18816. / (55691. * e_pow_2 * sqrt(e_fact));
  REAL8 pre_f_2 = 16382. * sqrt(e_fact) / (55691. * e_pow_2);
  REAL8 bit_1 = sqrt(e_fact) * (1. - (11. / 7.) * e_pow_2) * phi_e(e) -
                (1. - (3. / 7.) * e_pow_2) * phi_e_tilde(e);
  REAL8 bit_2 = sqrt(e_fact) * psi_e(e) - psi_e_tilde(e);

  return (pre_f_1 * bit_1 + pre_f_2 * bit_2);
}

static REAL8 zed_e_rad(REAL8 e) {
  REAL8 e_fact = 1. - e * e;
  REAL8 e_pow_2 = e * e;
  REAL8 pre_f_1 = 924. / (19067. * e_pow_2 * sqrt(e_fact));
  REAL8 pre_f_2 = 12243. * sqrt(e_fact) / (76268. * e_pow_2);
  REAL8 bit_1 = -e_fact * sqrt(e_fact) * phi_e(e) +
                (1. - (5. / 11.) * e_pow_2) * phi_e_tilde(e);
  REAL8 bit_2 = sqrt(e_fact) * zed_e(e) - zed_e_tilde(e);

  return (pre_f_1 * bit_1 + pre_f_2 * bit_2);
}

static REAL8 kappa_e_rad(REAL8 e) {
  const REAL8 log2 = 0.693147180559945309417232121458; // ln(2)
  const REAL8 log3 = 1.09861228866810969139524523692;  // ln(3)

  REAL8 e_fact = 1. - e * e;
  REAL8 e_pow_2 = e * e;
  REAL8 den =
      769. / 96. - 3059665. * log2 / 700566. + 8190315. * log3 / 1868176.;
  REAL8 pre_f = sqrt(e_fact) / e_pow_2;
  REAL8 bit_1 = sqrt(e_fact) * kappa_e(e) - kappa_e_tilde(e);

  return (pre_f * bit_1 / den);
}

static REAL8 f_e(REAL8 e) {
  REAL8 e_fact = 1. - e * e;
  REAL8 e_pow_2 = e * e;
  REAL8 e_pow_4 = e_pow_2 * e_pow_2;
  REAL8 e_pow_6 = e_pow_4 * e_pow_2;
  REAL8 e_pow_8 = e_pow_6 * e_pow_2;
  REAL8 denom =
      sqrt(e_fact) * e_fact * e_fact * e_fact * e_fact * e_fact * e_fact;
  REAL8 num = 1. + 85. * e_pow_2 / 6. + 5171. * e_pow_4 / 192. +
              1751. * e_pow_6 / 192. + 297. * e_pow_8 / 1024.;

  return (num / denom);
}

static REAL8 capital_f_e(REAL8 e) {
  REAL8 e_fact = 1. - e * e;
  REAL8 e_pow_2 = e * e;
  REAL8 e_pow_4 = e_pow_2 * e_pow_2;
  REAL8 e_pow_6 = e_pow_4 * e_pow_2;
  REAL8 denom = sqrt(e_fact) * e_fact * e_fact * e_fact * e_fact * e_fact;
  REAL8 num = 1. + 2782. * e_pow_2 / 769. + 10721. * e_pow_4 / 6152. +
              1719. * e_pow_6 / 24608.;

  return (num / denom);
}

static REAL8 psi_n(REAL8 e) {
  REAL8 e_fact = 1. - e * e;
  REAL8 e_pow_2 = e * e;
  REAL8 pre_f_1 = 1344. * (7. - 5. * e_pow_2) / (17599. * e_fact);
  REAL8 pre_f_2 = 8191. / 17599.;

  return (pre_f_1 * phi_e(e) + pre_f_2 * psi_e(e));
}

static REAL8 zed_n(REAL8 e) {
  REAL8 pre_f_1 = 583. / 567.;
  REAL8 pre_f_2 = 16. / 567.;

  return (pre_f_1 * zed_e(e) - pre_f_2 * phi_e(e));
}

// static REAL8 hPlus(REAL8 x, REAL8 x0, REAL8 m1, REAL8 m2, REAL8 i, REAL8 phi,
// UINT4 vpn) {
//   const REAL8 log2 = 0.693147180559945309417232121458;   // ln(2)
//   const REAL8 log3_2 = 0.405465108108164381978013115464; // ln(3/2)
//   /* some math:
//    * sin(2*a) = 2*sin(a)*cos(a)
//    * sin(3*a) = 4*sin(a)*cos(a)^2-sin(a)
//    * sin(4*a) = 8*sin(a)*cos(a)^3-4*sin(a)*cos(a)
//    * sin(5*a) = 16*sin(a)*cos(a)^4-12*sin(a)*cos(a)^2+sin(a)
//    * sin(6*a) = 32*sin(a)*cos(a)^5-32*sin(a)*cos(a)^3+6*sin(a)*cos(a)
//    * sin(7*a) =
//    64*sin(a)*cos(a)^6-80*sin(a)*cos(a)^4+24*sin(a)*cos(a)^2-sin(a)
//    * cos(2*a) = 2*cos(a)^2  - 1
//    * cos(3*a) = 4*cos(a)^3 - 3*cos(a)
//    * cos(4*a) = 8*cos(a)^4 - 8*cos(a)^2  + 1
//    * cos(5*a) = 16*cos(a)^5 - 20*cos(a)^3 + 5*cos(a)
//    * cos(6*a) = 32*cos(a)^6 - 48*cos(a)^4 + 18*cos(a)^2 - 1
//    * cos(7*a) = 64*cos(a)^7 - 112*cos(a)^5 + 56*cos(a)^3 - 7*cos(a)
//    */
//   const REAL8 a = phi - 2 * pow3_2(x) * log(pow3_2(x / x0));
//   const REAL8 cosa = cos(a), sina = sin(a);
//   const REAL8 cos2a = 2 * pow2(cosa) - 1;
//   const REAL8 cos3a = 4 * pow3(cosa) - 3 * cosa;
//   const REAL8 cos4a = 8 * pow4(cosa) - 8 * pow2(cosa) + 1;
//   const REAL8 cos5a = 16 * pow5(cosa) - 20 * pow3(cosa) + 5 * cosa;
//   const REAL8 cos6a = 32 * pow6(cosa) - 48 * pow4(cosa) + 18 * pow2(cosa) -
//   1; const REAL8 cos7a =
//       64 * pow7(cosa) - 112 * pow5(cosa) + 56 * pow3(cosa) - 7 * cosa;
//   const REAL8 sin2a = 2 * sina * cosa;
//   const REAL8 sin3a = 4 * sina * pow2(cosa) - sina;
//   const REAL8 sin4a = 8 * sina * pow3(cosa) - 4 * sina * cosa;
//   const REAL8 sin5a = 16 * sina * pow4(cosa) - 12 * sina * pow2(cosa) + sina;
//   double EulerGamma = 0.5772156649015329;

//   REAL8 Nu = (m1*m2)/(pow2(m1+m2));
//   REAL8 delta = (m1-m2)/(m1+m2);

//   if(vpn==1){
//   //Note : 0.5PN term is removed here
//   return (2 * x *
//           (x * (((3.1666666666666665 + (3 * pow2(cos(i))) / 2. -
//                 pow4(cos(i)) / 3. +
//                 (m1 * m2 *
//                  (-3.1666666666666665 + (11 * pow2(cos(i))) / 6. +
//                   pow4(cos(i)))) /
//                     pow2(m1 + m2)) *
//                    cos2a -
//                (4 * (1 - (3 * m1 * m2) / pow2(m1 + m2)) * (1 + pow2(cos(i)))
//                *
//                 cos4a * pow2(sin(i))) /
//                    3.) )
//                     + pow3_2(x) *
//             ((-2 * M_PI * (1 + pow2(cos(i))) * cos2a +
//                ((m1 - m2) *
//                 (0.296875 + (5 * pow2(cos(i))) / 16. - pow4(cos(i)) / 192. +
//                  (m1 * m2 *
//                   (-0.5104166666666666 + pow2(cos(i)) / 8. +
//                    pow4(cos(i)) / 96.)) /
//                      pow2(m1 + m2)) *
//                 cos(a) * sin(i)) /
//                    (m1 + m2) +
//                ((m1 - m2) *
//                 (-5.1328125 - (45 * pow2(cos(i))) / 16. +
//                  (81 * pow4(cos(i))) / 128. +
//                  (m1 * m2 *
//                   (3.515625 - (9 * pow2(cos(i))) / 8. -
//                    (81 * pow4(cos(i))) / 64.)) /
//                      pow2(m1 + m2)) *
//                 cos3a * sin(i)) /
//                    (m1 + m2) +
//                (625 * (m1 - m2) * (1 - (2 * m1 * m2) / pow2(m1 + m2)) *
//                 (1 + pow2(cos(i))) * cos5a * pow3(sin(i))) /
//                    (384. * (m1 + m2))) )
//                    + pow2(x) *
//               ((0.18333333333333332 + (33 * pow2(cos(i))) / 10. +
//                 (29 * pow4(cos(i))) / 24. - pow6(cos(i)) / 24. +
//                 (pow2(m1) * pow2(m2) *
//                  (-4.083333333333333 + (9 * pow2(cos(i))) / 2. -
//                   (7 * pow4(cos(i))) / 24. - (5 * pow6(cos(i))) / 24.)) /
//                     pow4(m1 + m2) +
//                 (m1 * m2 *
//                  (9.805555555555555 - 3 * pow2(cos(i)) -
//                   (251 * pow4(cos(i))) / 72. + (5 * pow6(cos(i))) / 24.)) /
//                     pow2(m1 + m2)) *
//                    cos2a +
//                ((m1 - m2) * M_PI * (-0.625 - pow2(cos(i)) / 8.) * cos(a) *
//                 sin(i)) /
//                    (m1 + m2) +
//                (27 * (m1 - m2) * M_PI * (1 + pow2(cos(i))) * cos3a * sin(i))
//                /
//                    (8. * (m1 + m2)) +
//                (2 *
//                 (59 + 35 * pow2(cos(i)) - 8 * pow4(cos(i)) -
//                  (5 * m1 * m2 * (131 + 59 * pow2(cos(i)) - 24 *
//                  pow4(cos(i)))) /
//                      (3. * pow2(m1 + m2)) +
//                  (5 * pow2(m1) * pow2(m2) *
//                   (21 - 3 * pow2(cos(i)) - 8 * pow4(cos(i)))) /
//                      pow4(m1 + m2)) *
//                 cos4a * pow2(sin(i))) /
//                    15. -
//                (81 *
//                 (1 + (5 * pow2(m1) * pow2(m2)) / pow4(m1 + m2) -
//                  (5 * m1 * m2) / pow2(m1 + m2)) *
//                 (1 + pow2(cos(i))) * cos6a * pow4(sin(i))) /
//                    40. +
//                ((m1 - m2) *
//                 (0.275 + pow2(cos(i)) * (0.175 + log2 / 4.) + (5 * log2)
//                 / 4.) * sin(i) * sin(a)) /
//                    (m1 + m2) +
//                ((m1 - m2) * (1 + pow2(cos(i))) * (-4.725 + (27 * log3_2)
//                / 4.) *
//                 sin(i) * sin3a) /
//                    (m1 + m2) ) +
//           pow5_2(x) *
//               (M_PI *
//                    (6.333333333333333 + 3 * pow2(cos(i)) -
//                     (2 * pow4(cos(i))) / 3. +
//                     (m1 * m2 *
//                      (-5.333333333333333 + (14 * pow2(cos(i))) / 3. +
//                       2 * pow4(cos(i)))) /
//                         pow2(m1 + m2)) *
//                    cos2a +
//                ((m1 - m2) *
//                 (0.3458984375 - (1667 * pow2(cos(i))) / 5120. +
//                  (217 * pow4(cos(i))) / 9216. - pow6(cos(i)) / 9216. +
//                  (pow2(m1) * pow2(m2) *
//                   (-0.3744574652777778 + (673 * pow2(cos(i))) / 3072. -
//                    (5 * pow4(cos(i))) / 9216. - pow6(cos(i)) / 3072.)) /
//                      pow4(m1 + m2) +
//                  (m1 * m2 *
//                   (2.66015625 + (13 * pow2(cos(i))) / 768. -
//                    (35 * pow4(cos(i))) / 768. + pow6(cos(i)) / 2304.)) /
//                      pow2(m1 + m2)) *
//                 cos(a) * sin(i)) /
//                    (m1 + m2) +
//                ((m1 - m2) *
//                 (3.4541015625 - (22977 * pow2(cos(i))) / 5120. -
//                  (15309 * pow4(cos(i))) / 5120. + (729 * pow6(cos(i))) /
//                  5120. + (m1 * m2 *
//                   (-18.61640625 + (5529 * pow2(cos(i))) / 1280. +
//                    (7749 * pow4(cos(i))) / 1280. -
//                    (729 * pow6(cos(i))) / 1280.)) /
//                      pow2(m1 + m2) +
//                  (pow2(m1) * pow2(m2) *
//                   (5.6888671875 - (27267 * pow2(cos(i))) / 5120. -
//                    (1647 * pow4(cos(i))) / 5120. +
//                    (2187 * pow6(cos(i))) / 5120.)) /
//                      pow4(m1 + m2)) *
//                 cos3a * sin(i)) /
//                    (m1 + m2) +
//                ((m1 - m2) *
//                 (-11.732313368055555 + (40625 * pow2(cos(i))) / 9216. +
//                  (83125 * pow4(cos(i))) / 9216. -
//                  (15625 * pow6(cos(i))) / 9216. +
//                  (pow2(m1) * pow2(m2) *
//                   (-12.953016493055555 + (40625 * pow2(cos(i))) / 3072. +
//                    (44375 * pow4(cos(i))) / 9216. -
//                    (15625 * pow6(cos(i))) / 3072.)) /
//                      pow4(m1 + m2) +
//                  (m1 * m2 *
//                   (31.73828125 - (40625 * pow2(cos(i))) / 2304. -
//                    (48125 * pow4(cos(i))) / 2304. +
//                    (15625 * pow6(cos(i))) / 2304.)) /
//                      pow2(m1 + m2)) *
//                 cos5a * sin(i)) /
//                    (m1 + m2) -
//                (16 * (1 - (3 * m1 * m2) / pow2(m1 + m2)) * M_PI *
//                 (1 + pow2(cos(i))) * cos4a * pow2(sin(i))) /
//                    3. +
//                (117649 * (m1 - m2) *
//                 (1 + (3 * pow2(m1) * pow2(m2)) / pow4(m1 + m2) -
//                  (4 * m1 * m2) / pow2(m1 + m2)) *
//                 (1 + pow2(cos(i))) * cos7a * pow5(sin(i))) /
//                    (46080. * (m1 + m2)) +
//                (-1.8 + (14 * pow2(cos(i))) / 5. + (7 * pow4(cos(i))) / 5. +
//                 (m1 * m2 *
//                  (32 + (56 * pow2(cos(i))) / 5. - (28 * pow4(cos(i))) / 5.))
//                  /
//                     pow2(m1 + m2)) *
//                    sin2a +
//                (1 + pow2(cos(i))) *
//                    (11.2 - (32 * log2) / 3. +
//                     (m1 * m2 * (-39.766666666666666 + 32 * log2)) /
//                         pow2(m1 + m2)) *
//                    pow2(sin(i)) * sin4a ) ));}

//          else if(vpn==2){
//     //Till 1PN terms are removed
//     return (2 * x *
//          ( pow3_2(x) *
//               ((-2 * M_PI * (1 + pow2(cos(i))) * cos2a +
//                ((m1 - m2) *
//                 (0.296875 + (5 * pow2(cos(i))) / 16. - pow4(cos(i)) / 192. +
//                  (m1 * m2 *
//                   (-0.5104166666666666 + pow2(cos(i)) / 8. +
//                    pow4(cos(i)) / 96.)) /
//                      pow2(m1 + m2)) *
//                 cos(a) * sin(i)) /
//                    (m1 + m2) +
//                ((m1 - m2) *
//                 (-5.1328125 - (45 * pow2(cos(i))) / 16. +
//                  (81 * pow4(cos(i))) / 128. +
//                  (m1 * m2 *
//                   (3.515625 - (9 * pow2(cos(i))) / 8. -
//                    (81 * pow4(cos(i))) / 64.)) /
//                      pow2(m1 + m2)) *
//                 cos3a * sin(i)) /
//                    (m1 + m2) +
//                (625 * (m1 - m2) * (1 - (2 * m1 * m2) / pow2(m1 + m2)) *
//                 (1 + pow2(cos(i))) * cos5a * pow3(sin(i))) /
//                    (384. * (m1 + m2))) ) +
//           pow2(x) *
//               ((0.18333333333333332 + (33 * pow2(cos(i))) / 10. +
//                 (29 * pow4(cos(i))) / 24. - pow6(cos(i)) / 24. +
//                 (pow2(m1) * pow2(m2) *
//                  (-4.083333333333333 + (9 * pow2(cos(i))) / 2. -
//                   (7 * pow4(cos(i))) / 24. - (5 * pow6(cos(i))) / 24.)) /
//                     pow4(m1 + m2) +
//                 (m1 * m2 *
//                  (9.805555555555555 - 3 * pow2(cos(i)) -
//                   (251 * pow4(cos(i))) / 72. + (5 * pow6(cos(i))) / 24.)) /
//                     pow2(m1 + m2)) *
//                    cos2a +
//                ((m1 - m2) * M_PI * (-0.625 - pow2(cos(i)) / 8.) * cos(a) *
//                 sin(i)) /
//                    (m1 + m2) +
//                (27 * (m1 - m2) * M_PI * (1 + pow2(cos(i))) * cos3a * sin(i))
//                /
//                    (8. * (m1 + m2)) +
//                (2 *
//                 (59 + 35 * pow2(cos(i)) - 8 * pow4(cos(i)) -
//                  (5 * m1 * m2 * (131 + 59 * pow2(cos(i)) - 24 *
//                  pow4(cos(i)))) /
//                      (3. * pow2(m1 + m2)) +
//                  (5 * pow2(m1) * pow2(m2) *
//                   (21 - 3 * pow2(cos(i)) - 8 * pow4(cos(i)))) /
//                      pow4(m1 + m2)) *
//                 cos4a * pow2(sin(i))) /
//                    15. -
//                (81 *
//                 (1 + (5 * pow2(m1) * pow2(m2)) / pow4(m1 + m2) -
//                  (5 * m1 * m2) / pow2(m1 + m2)) *
//                 (1 + pow2(cos(i))) * cos6a * pow4(sin(i))) /
//                    40. +
//                ((m1 - m2) *
//                 (0.275 + pow2(cos(i)) * (0.175 + log2 / 4.) + (5 * log2)
//                 / 4.) * sin(i) * sin(a)) /
//                    (m1 + m2) +
//                ((m1 - m2) * (1 + pow2(cos(i))) * (-4.725 + (27 * log3_2)
//                / 4.) *
//                 sin(i) * sin3a) /
//                    (m1 + m2) ) +
//           pow5_2(x) *
//               (M_PI *
//                    (6.333333333333333 + 3 * pow2(cos(i)) -
//                     (2 * pow4(cos(i))) / 3. +
//                     (m1 * m2 *
//                      (-5.333333333333333 + (14 * pow2(cos(i))) / 3. +
//                       2 * pow4(cos(i)))) /
//                         pow2(m1 + m2)) *
//                    cos2a +
//                ((m1 - m2) *
//                 (0.3458984375 - (1667 * pow2(cos(i))) / 5120. +
//                  (217 * pow4(cos(i))) / 9216. - pow6(cos(i)) / 9216. +
//                  (pow2(m1) * pow2(m2) *
//                   (-0.3744574652777778 + (673 * pow2(cos(i))) / 3072. -
//                    (5 * pow4(cos(i))) / 9216. - pow6(cos(i)) / 3072.)) /
//                      pow4(m1 + m2) +
//                  (m1 * m2 *
//                   (2.66015625 + (13 * pow2(cos(i))) / 768. -
//                    (35 * pow4(cos(i))) / 768. + pow6(cos(i)) / 2304.)) /
//                      pow2(m1 + m2)) *
//                 cos(a) * sin(i)) /
//                    (m1 + m2) +
//                ((m1 - m2) *
//                 (3.4541015625 - (22977 * pow2(cos(i))) / 5120. -
//                  (15309 * pow4(cos(i))) / 5120. + (729 * pow6(cos(i))) /
//                  5120. + (m1 * m2 *
//                   (-18.61640625 + (5529 * pow2(cos(i))) / 1280. +
//                    (7749 * pow4(cos(i))) / 1280. -
//                    (729 * pow6(cos(i))) / 1280.)) /
//                      pow2(m1 + m2) +
//                  (pow2(m1) * pow2(m2) *
//                   (5.6888671875 - (27267 * pow2(cos(i))) / 5120. -
//                    (1647 * pow4(cos(i))) / 5120. +
//                    (2187 * pow6(cos(i))) / 5120.)) /
//                      pow4(m1 + m2)) *
//                 cos3a * sin(i)) /
//                    (m1 + m2) +
//                ((m1 - m2) *
//                 (-11.732313368055555 + (40625 * pow2(cos(i))) / 9216. +
//                  (83125 * pow4(cos(i))) / 9216. -
//                  (15625 * pow6(cos(i))) / 9216. +
//                  (pow2(m1) * pow2(m2) *
//                   (-12.953016493055555 + (40625 * pow2(cos(i))) / 3072. +
//                    (44375 * pow4(cos(i))) / 9216. -
//                    (15625 * pow6(cos(i))) / 3072.)) /
//                      pow4(m1 + m2) +
//                  (m1 * m2 *
//                   (31.73828125 - (40625 * pow2(cos(i))) / 2304. -
//                    (48125 * pow4(cos(i))) / 2304. +
//                    (15625 * pow6(cos(i))) / 2304.)) /
//                      pow2(m1 + m2)) *
//                 cos5a * sin(i)) /
//                    (m1 + m2) -
//                (16 * (1 - (3 * m1 * m2) / pow2(m1 + m2)) * M_PI *
//                 (1 + pow2(cos(i))) * cos4a * pow2(sin(i))) /
//                    3. +
//                (117649 * (m1 - m2) *
//                 (1 + (3 * pow2(m1) * pow2(m2)) / pow4(m1 + m2) -
//                  (4 * m1 * m2) / pow2(m1 + m2)) *
//                 (1 + pow2(cos(i))) * cos7a * pow5(sin(i))) /
//                    (46080. * (m1 + m2)) +
//                (-1.8 + (14 * pow2(cos(i))) / 5. + (7 * pow4(cos(i))) / 5. +
//                 (m1 * m2 *
//                  (32 + (56 * pow2(cos(i))) / 5. - (28 * pow4(cos(i))) / 5.))
//                  /
//                     pow2(m1 + m2)) *
//                    sin2a +
//                (1 + pow2(cos(i))) *
//                    (11.2 - (32 * log2) / 3. +
//                     (m1 * m2 * (-39.766666666666666 + 32 * log2)) /
//                         pow2(m1 + m2)) *
//                    pow2(sin(i)) * sin4a ) ));}

//   else if(vpn==3){
//     //keeping only 1.5PN hereditary term
//     return (2 * x *
//          ( pow3_2(x) *
//               (-2 * M_PI * (1 + pow2(cos(i))) * cos2a) +
//           pow2(x) *
//               ((0.18333333333333332 + (33 * pow2(cos(i))) / 10. +
//                 (29 * pow4(cos(i))) / 24. - pow6(cos(i)) / 24. +
//                 (pow2(m1) * pow2(m2) *
//                  (-4.083333333333333 + (9 * pow2(cos(i))) / 2. -
//                   (7 * pow4(cos(i))) / 24. - (5 * pow6(cos(i))) / 24.)) /
//                     pow4(m1 + m2) +
//                 (m1 * m2 *
//                  (9.805555555555555 - 3 * pow2(cos(i)) -
//                   (251 * pow4(cos(i))) / 72. + (5 * pow6(cos(i))) / 24.)) /
//                     pow2(m1 + m2)) *
//                    cos2a +
//                ((m1 - m2) * M_PI * (-0.625 - pow2(cos(i)) / 8.) * cos(a) *
//                 sin(i)) /
//                    (m1 + m2) +
//                (27 * (m1 - m2) * M_PI * (1 + pow2(cos(i))) * cos3a * sin(i))
//                /
//                    (8. * (m1 + m2)) +
//                (2 *
//                 (59 + 35 * pow2(cos(i)) - 8 * pow4(cos(i)) -
//                  (5 * m1 * m2 * (131 + 59 * pow2(cos(i)) - 24 *
//                  pow4(cos(i)))) /
//                      (3. * pow2(m1 + m2)) +
//                  (5 * pow2(m1) * pow2(m2) *
//                   (21 - 3 * pow2(cos(i)) - 8 * pow4(cos(i)))) /
//                      pow4(m1 + m2)) *
//                 cos4a * pow2(sin(i))) /
//                    15. -
//                (81 *
//                 (1 + (5 * pow2(m1) * pow2(m2)) / pow4(m1 + m2) -
//                  (5 * m1 * m2) / pow2(m1 + m2)) *
//                 (1 + pow2(cos(i))) * cos6a * pow4(sin(i))) /
//                    40. +
//                ((m1 - m2) *
//                 (0.275 + pow2(cos(i)) * (0.175 + log2 / 4.) + (5 * log2)
//                 / 4.) * sin(i) * sin(a)) /
//                    (m1 + m2) +
//                ((m1 - m2) * (1 + pow2(cos(i))) * (-4.725 + (27 * log3_2)
//                / 4.) *
//                 sin(i) * sin3a) /
//                    (m1 + m2) ) +
//           pow5_2(x) *
//               (M_PI *
//                    (6.333333333333333 + 3 * pow2(cos(i)) -
//                     (2 * pow4(cos(i))) / 3. +
//                     (m1 * m2 *
//                      (-5.333333333333333 + (14 * pow2(cos(i))) / 3. +
//                       2 * pow4(cos(i)))) /
//                         pow2(m1 + m2)) *
//                    cos2a +
//                ((m1 - m2) *
//                 (0.3458984375 - (1667 * pow2(cos(i))) / 5120. +
//                  (217 * pow4(cos(i))) / 9216. - pow6(cos(i)) / 9216. +
//                  (pow2(m1) * pow2(m2) *
//                   (-0.3744574652777778 + (673 * pow2(cos(i))) / 3072. -
//                    (5 * pow4(cos(i))) / 9216. - pow6(cos(i)) / 3072.)) /
//                      pow4(m1 + m2) +
//                  (m1 * m2 *
//                   (2.66015625 + (13 * pow2(cos(i))) / 768. -
//                    (35 * pow4(cos(i))) / 768. + pow6(cos(i)) / 2304.)) /
//                      pow2(m1 + m2)) *
//                 cos(a) * sin(i)) /
//                    (m1 + m2) +
//                ((m1 - m2) *
//                 (3.4541015625 - (22977 * pow2(cos(i))) / 5120. -
//                  (15309 * pow4(cos(i))) / 5120. + (729 * pow6(cos(i))) /
//                  5120. + (m1 * m2 *
//                   (-18.61640625 + (5529 * pow2(cos(i))) / 1280. +
//                    (7749 * pow4(cos(i))) / 1280. -
//                    (729 * pow6(cos(i))) / 1280.)) /
//                      pow2(m1 + m2) +
//                  (pow2(m1) * pow2(m2) *
//                   (5.6888671875 - (27267 * pow2(cos(i))) / 5120. -
//                    (1647 * pow4(cos(i))) / 5120. +
//                    (2187 * pow6(cos(i))) / 5120.)) /
//                      pow4(m1 + m2)) *
//                 cos3a * sin(i)) /
//                    (m1 + m2) +
//                ((m1 - m2) *
//                 (-11.732313368055555 + (40625 * pow2(cos(i))) / 9216. +
//                  (83125 * pow4(cos(i))) / 9216. -
//                  (15625 * pow6(cos(i))) / 9216. +
//                  (pow2(m1) * pow2(m2) *
//                   (-12.953016493055555 + (40625 * pow2(cos(i))) / 3072. +
//                    (44375 * pow4(cos(i))) / 9216. -
//                    (15625 * pow6(cos(i))) / 3072.)) /
//                      pow4(m1 + m2) +
//                  (m1 * m2 *
//                   (31.73828125 - (40625 * pow2(cos(i))) / 2304. -
//                    (48125 * pow4(cos(i))) / 2304. +
//                    (15625 * pow6(cos(i))) / 2304.)) /
//                      pow2(m1 + m2)) *
//                 cos5a * sin(i)) /
//                    (m1 + m2) -
//                (16 * (1 - (3 * m1 * m2) / pow2(m1 + m2)) * M_PI *
//                 (1 + pow2(cos(i))) * cos4a * pow2(sin(i))) /
//                    3. +
//                (117649 * (m1 - m2) *
//                 (1 + (3 * pow2(m1) * pow2(m2)) / pow4(m1 + m2) -
//                  (4 * m1 * m2) / pow2(m1 + m2)) *
//                 (1 + pow2(cos(i))) * cos7a * pow5(sin(i))) /
//                    (46080. * (m1 + m2)) +
//                (-1.8 + (14 * pow2(cos(i))) / 5. + (7 * pow4(cos(i))) / 5. +
//                 (m1 * m2 *
//                  (32 + (56 * pow2(cos(i))) / 5. - (28 * pow4(cos(i))) / 5.))
//                  /
//                     pow2(m1 + m2)) *
//                    sin2a +
//                (1 + pow2(cos(i))) *
//                    (11.2 - (32 * log2) / 3. +
//                     (m1 * m2 * (-39.766666666666666 + 32 * log2)) /
//                         pow2(m1 + m2)) *
//                    pow2(sin(i)) * sin4a ) ));}

//    else if(vpn==4){
//     //keeping only 1.5PN and 2PN hereditary term
//     return (2 * x *
//          ( pow3_2(x) *
//               (-2 * M_PI * (1 + pow2(cos(i))) * cos2a) +
//           pow2(x) *
//               (((m1 - m2) * M_PI * (-0.625 - pow2(cos(i)) / 8.) * cos(a) *
//                 sin(i)) /
//                    (m1 + m2) +
//                (27 * (m1 - m2) * M_PI * (1 + pow2(cos(i))) * cos3a * sin(i))
//                /
//                    (8. * (m1 + m2)) +
//                ((m1 - m2) *
//                 (0.275 + pow2(cos(i)) * (0.175 + log2 / 4.) + (5 * log2)
//                 / 4.) * sin(i) * sin(a)) /
//                    (m1 + m2) +
//                ((m1 - m2) * (1 + pow2(cos(i))) * (-4.725 + (27 * log3_2)
//                / 4.) *
//                 sin(i) * sin3a) /
//                    (m1 + m2)) +
//           pow5_2(x) *
//               (M_PI *
//                    (6.333333333333333 + 3 * pow2(cos(i)) -
//                     (2 * pow4(cos(i))) / 3. +
//                     (m1 * m2 *
//                      (-5.333333333333333 + (14 * pow2(cos(i))) / 3. +
//                       2 * pow4(cos(i)))) /
//                         pow2(m1 + m2)) *
//                    cos2a +
//                ((m1 - m2) *
//                 (0.3458984375 - (1667 * pow2(cos(i))) / 5120. +
//                  (217 * pow4(cos(i))) / 9216. - pow6(cos(i)) / 9216. +
//                  (pow2(m1) * pow2(m2) *
//                   (-0.3744574652777778 + (673 * pow2(cos(i))) / 3072. -
//                    (5 * pow4(cos(i))) / 9216. - pow6(cos(i)) / 3072.)) /
//                      pow4(m1 + m2) +
//                  (m1 * m2 *
//                   (2.66015625 + (13 * pow2(cos(i))) / 768. -
//                    (35 * pow4(cos(i))) / 768. + pow6(cos(i)) / 2304.)) /
//                      pow2(m1 + m2)) *
//                 cos(a) * sin(i)) /
//                    (m1 + m2) +
//                ((m1 - m2) *
//                 (3.4541015625 - (22977 * pow2(cos(i))) / 5120. -
//                  (15309 * pow4(cos(i))) / 5120. + (729 * pow6(cos(i))) /
//                  5120. + (m1 * m2 *
//                   (-18.61640625 + (5529 * pow2(cos(i))) / 1280. +
//                    (7749 * pow4(cos(i))) / 1280. -
//                    (729 * pow6(cos(i))) / 1280.)) /
//                      pow2(m1 + m2) +
//                  (pow2(m1) * pow2(m2) *
//                   (5.6888671875 - (27267 * pow2(cos(i))) / 5120. -
//                    (1647 * pow4(cos(i))) / 5120. +
//                    (2187 * pow6(cos(i))) / 5120.)) /
//                      pow4(m1 + m2)) *
//                 cos3a * sin(i)) /
//                    (m1 + m2) +
//                ((m1 - m2) *
//                 (-11.732313368055555 + (40625 * pow2(cos(i))) / 9216. +
//                  (83125 * pow4(cos(i))) / 9216. -
//                  (15625 * pow6(cos(i))) / 9216. +
//                  (pow2(m1) * pow2(m2) *
//                   (-12.953016493055555 + (40625 * pow2(cos(i))) / 3072. +
//                    (44375 * pow4(cos(i))) / 9216. -
//                    (15625 * pow6(cos(i))) / 3072.)) /
//                      pow4(m1 + m2) +
//                  (m1 * m2 *
//                   (31.73828125 - (40625 * pow2(cos(i))) / 2304. -
//                    (48125 * pow4(cos(i))) / 2304. +
//                    (15625 * pow6(cos(i))) / 2304.)) /
//                      pow2(m1 + m2)) *
//                 cos5a * sin(i)) /
//                    (m1 + m2) -
//                (16 * (1 - (3 * m1 * m2) / pow2(m1 + m2)) * M_PI *
//                 (1 + pow2(cos(i))) * cos4a * pow2(sin(i))) /
//                    3. +
//                (117649 * (m1 - m2) *
//                 (1 + (3 * pow2(m1) * pow2(m2)) / pow4(m1 + m2) -
//                  (4 * m1 * m2) / pow2(m1 + m2)) *
//                 (1 + pow2(cos(i))) * cos7a * pow5(sin(i))) /
//                    (46080. * (m1 + m2)) +
//                (-1.8 + (14 * pow2(cos(i))) / 5. + (7 * pow4(cos(i))) / 5. +
//                 (m1 * m2 *
//                  (32 + (56 * pow2(cos(i))) / 5. - (28 * pow4(cos(i))) / 5.))
//                  /
//                     pow2(m1 + m2)) *
//                    sin2a +
//                (1 + pow2(cos(i))) *
//                    (11.2 - (32 * log2) / 3. +
//                     (m1 * m2 * (-39.766666666666666 + 32 * log2)) /
//                         pow2(m1 + m2)) *
//                    pow2(sin(i)) * sin4a ) ));}

//   else{
//     //keeping only 1.5PN, 2PN and 2.5PN term
//      return (2 * x *
//          ( pow3_2(x) *
//               (-2 * M_PI * (1 + pow2(cos(i))) * cos2a) +
//           pow2(x) *
//               (((m1 - m2) * M_PI * (-0.625 - pow2(cos(i)) / 8.) * cos(a) *
//                 sin(i)) /
//                    (m1 + m2) +
//                (27 * (m1 - m2) * M_PI * (1 + pow2(cos(i))) * cos3a * sin(i))
//                /
//                    (8. * (m1 + m2)) +
//                ((m1 - m2) *
//                 (0.275 + pow2(cos(i)) * (0.175 + log2 / 4.) + (5 * log2)
//                 / 4.) * sin(i) * sin(a)) /
//                    (m1 + m2) +
//                ((m1 - m2) * (1 + pow2(cos(i))) * (-4.725 + (27 * log3_2)
//                / 4.) *
//                 sin(i) * sin3a) /
//                    (m1 + m2)) +
//           pow5_2(x) *
//               (M_PI *
//                    (6.333333333333333 + 3 * pow2(cos(i)) -
//                     (2 * pow4(cos(i))) / 3. +
//                     (m1 * m2 *
//                      (-5.333333333333333 + (14 * pow2(cos(i))) / 3. +
//                       2 * pow4(cos(i)))) /
//                         pow2(m1 + m2)) *
//                    cos2a   -
//                (16 * (1 - (3 * m1 * m2) / pow2(m1 + m2)) * M_PI *
//                 (1 + pow2(cos(i))) * cos4a * pow2(sin(i)))/3. +
//                (-1.8 + (14 * pow2(cos(i))) / 5. + (7 * pow4(cos(i))) / 5. +
//                 (m1 * m2 *
//                  (32 + (56 * pow2(cos(i))) / 5. - (28 * pow4(cos(i))) / 5. -
//                  (1435 - 5*(7 + 24*sqrt(35))*cos(i) +
//      12*(35 + 16*sqrt(35))*cos(2*i) -
//      21*cos(3*i) - 72*sqrt(35)*cos(3*i) -
//      7*cos(4*i))/40.)) /
//                     pow2(m1 + m2)) *
//                    sin2a +
//                (1 + pow2(cos(i))) *
//                    (11.2 - (32 * log2) / 3. +
//                     (m1 * m2 * (-39.766666666666666 + 32 * log2 +
//                     18*2*sqrt(1.4) )) /
//                         pow2(m1 + m2)) *
//                    pow2(sin(i)) * sin4a
//                    ) + pow3(x) *
//                    (cos5a * ((3125*delta*M_PI*(-1 + 2*Nu)*(3 +
//                    cos(2*i))*pow3(sin(i)))/768.) + sin5a * ((delta*(3 +
//                    cos(2*i))*(565625 - 1129522*Nu + 437500*(-1 +
//                    2*Nu)*log(2.5))*pow3(sin(i)))/
//                     53760.) +
//                     cos3a * ((27*delta*M_PI*(717 - 186*Nu + 4*(31 +
//                     42*Nu)*cos(2*i) + 9*(-1 + 2*Nu)*cos(4*i))*sin(i))/1024.)
//                     + sin3a * (-(delta*(cos(4*i)*(Nu*(791374 + 612360*log(2)
//                     - 612360*log(3)) - 2187*(181 + 140*log(2) - 140*log(3)))
//                     - 21*(2*Nu*(69623 + 150660*log(2) - 150660*log(3)) +
//                     81*(-9607 + 1020*log(2) + 6660*log(3) +
//                     768*log(57.6650390625) - 768*log(1024))) +
//                     28*cos(2*i)*(Nu*(270178 + 204120*log(2) - 204120*log(3))
//                     - 243*(-173 + 660*log(2) - 20*log(3) +
//                     64*log(57.6650390625) - 64*log(1024))))*sin(i)) /645120.)
//                     + cosa * (-(delta*M_PI*(77*(9 - 2*Nu) + 4*(59 +
//                     38*Nu)*cos(2*i) + (-1 + 2*Nu)*cos(4*i))*sin(i))/1536.) +
//                     sina * ((delta*(84043 + 385418*Nu + 291060*log(2) -
//                     64680*Nu*log(2) +
//                             3*cos(4*i)*(-181 + 4522*Nu - 140*log(2) +
//                             280*Nu*log(2)) + 12*cos(2*i)*(5519 + 8260*log(2)
//                             + 70*Nu*(847 + 76*log(2))))*sin(i))/322560.) +
//                     cos2a * (((3 + cos(2*i))*(-116761 - 59920*EulerGamma +
//                             4900*pow2(M_PI) - 119840*log(2) -
//                             14980*log(pow2(x))))/14700.) +
//                     sin2a * ((856*(M_PI)*cos(i))/105.)
//                     )  ));
//     }

// }

// static REAL8 hCross(REAL8 x, REAL8 x0, REAL8 m1, REAL8 m2, REAL8 i, REAL8
// phi, UINT4 vpn) {
//   const REAL8 log2 = 0.693147180559945309417232121458;   // ln(2)
//   const REAL8 log3_2 = 0.405465108108164381978013115464; // ln(3/2)
//   /* some math:
//    * sin(2*a) = 2*sin(a)*cos(a)
//    * sin(3*a) = 4*sin(a)*cos(a)^2-sin(a)
//    * sin(4*a) = 8*sin(a)*cos(a)^3-4*sin(a)*cos(a)
//    * sin(5*a) = 16*sin(a)*cos(a)^4-12*sin(a)*cos(a)^2+sin(a)
//    * sin(6*a) = 32*sin(a)*cos(a)^5-32*sin(a)*cos(a)^3+6*sin(a)*cos(a)
//    * sin(7*a) =
//    64*sin(a)*cos(a)^6-80*sin(a)*cos(a)^4+24*sin(a)*cos(a)^2-sin(a)
//    * cos(2*a) = 2*cos(a)^2  - 1
//    * cos(3*a) = 4*cos(a)^3 - 3*cos(a)
//    * cos(4*a) = 8*cos(a)^4 - 8*cos(a)^2  + 1
//    * cos(5*a) = 16*cos(a)^5 - 20*cos(a)^3 + 5*cos(a)
//    * cos(6*a) = 32*cos(a)^6 - 48*cos(a)^4 + 18*cos(a)^2 - 1
//    * cos(7*a) = 64*cos(a)^7 - 112*cos(a)^5 + 56*cos(a)^3 - 7*cos(a)
//    */
//   const REAL8 a = phi - 2 * pow3_2(x) * log(pow3_2(x / x0));
//   const REAL8 cosa = cos(a), sina = sin(a);
//   const REAL8 cos2a = 2 * pow2(cosa) - 1;
//   const REAL8 cos3a = 4 * pow3(cosa) - 3 * cosa;
//   const REAL8 cos4a = 8 * pow4(cosa) - 8 * pow2(cosa) + 1;
//   const REAL8 cos5a = 16 * pow5(cosa) - 20 * pow3(cosa) + 5 * cosa;
//   const REAL8 sin2a = 2 * sina * cosa;
//   const REAL8 sin3a = 4 * sina * pow2(cosa) - sina;
//   const REAL8 sin4a = 8 * sina * pow3(cosa) - 4 * sina * cosa;
//   const REAL8 sin5a = 16 * sina * pow4(cosa) - 12 * sina * pow2(cosa) + sina;
//   const REAL8 sin6a =
//       32 * sina * pow5(cosa) - 32 * sina * pow3(cosa) + 6 * sina * cosa;
//   const REAL8 sin7a = 64 * sina * pow6(cosa) - 80 * sina * pow4(cosa) +
//                       24 * sina * pow2(cosa) - sina;

//   double EulerGamma = 0.5772156649015329;

//  REAL8 Nu = (m1*m2)/(pow2(m1+m2));
//  REAL8 delta = (m1-m2)/(m1+m2);

//   if(vpn==1){
//   //Note : 0.5PN term is removed here
//   return (2 * x *
//           (x *((cos(i) *
//                    (5.666666666666667 - (4 * pow2(cos(i))) / 3. +
//                     (m1 * m2 * (-4.333333333333333 + 4 * pow2(cos(i)))) /
//                         pow2(m1 + m2)) *
//                    sin2a -
//                (8 * (1 - (3 * m1 * m2) / pow2(m1 + m2)) * cos(i) *
//                 pow2(sin(i)) * sin4a) /
//                    3.) ) +
//           pow3_2(x) *
//               ((((m1 - m2) * cos(i) *
//                 (0.65625 - (5 * pow2(cos(i))) / 96. +
//                  (m1 * m2 * (-0.4791666666666667 + (5 * pow2(cos(i))) / 48.))
//                  /
//                      pow2(m1 + m2)) *
//                 sin(i) * sina) /
//                    (m1 + m2) -
//                4 * M_PI * cos(i) * sin2a +
//                ((m1 - m2) * cos(i) *
//                 (-9.421875 + (135 * pow2(cos(i))) / 64. +
//                  (m1 * m2 * (5.34375 - (135 * pow2(cos(i))) / 32.)) /
//                      pow2(m1 + m2)) *
//                 sin(i) * sin3a) /
//                    (m1 + m2) +
//                (625 * (m1 - m2) * (1 - (2 * m1 * m2) / pow2(m1 + m2)) *
//                cos(i) *
//                 pow3(sin(i)) * sin5a) /
//                    (192. * (m1 + m2))) ) +
//           pow2(x) * (((m1 - m2) * cos(i) * cos3a * (9.45 - (27 * log3_2)
//           / 2.) *
//                       sin(i)) /
//                          (m1 + m2) +
//                      ((m1 - m2) * cos(i) * cosa * (-0.45 - (3 * log2) / 2.) *
//                       sin(i)) /
//                          (m1 + m2) -
//                      (3 * (m1 - m2) * M_PI * cos(i) * sin(i) * sina) /
//                          (4. * (m1 + m2)) +
//                      cos(i) *
//                          (1.1333333333333333 + (113 * pow2(cos(i))) / 30. -
//                           pow4(cos(i)) / 4. +
//                           (pow2(m1) * pow2(m2) *
//                            (-4.666666666666667 + (35 * pow2(cos(i))) / 6. -
//                             (5 * pow4(cos(i))) / 4.)) /
//                               pow4(m1 + m2) +
//                           (m1 * m2 *
//                            (15.88888888888889 - (245 * pow2(cos(i))) / 18. +
//                             (5 * pow4(cos(i))) / 4.)) /
//                               pow2(m1 + m2)) *
//                          sin2a +
//                      (27 * (m1 - m2) * M_PI * cos(i) * sin(i) * sin3a) /
//                          (4. * (m1 + m2)) +
//                      (4 * cos(i) *
//                       (55 - 12 * pow2(cos(i)) -
//                        (5 * m1 * m2 * (119 - 36 * pow2(cos(i)))) /
//                            (3. * pow2(m1 + m2)) +
//                        (5 * pow2(m1) * pow2(m2) * (17 - 12 * pow2(cos(i)))) /
//                            pow4(m1 + m2)) *
//                       pow2(sin(i)) * sin4a) /
//                          15. -
//                      (81 *
//                       (1 + (5 * pow2(m1) * pow2(m2)) / pow4(m1 + m2) -
//                        (5 * m1 * m2) / pow2(m1 + m2)) *
//                       cos(i) * pow4(sin(i)) * sin6a) /
//                          20. ) +
//           pow5_2(x) *
//               (cos(i) *
//                    (2 - (22 * pow2(cos(i))) / 5. +
//                     (m1 * m2 * (-56.4 + (94 * pow2(cos(i))) / 5.)) /
//                         pow2(m1 + m2)) *
//                    cos2a +
//                (6 * m1 * m2 * cos(i) * pow2(sin(i))) / (5. * pow2(m1 + m2)) +
//                cos(i) * cos4a *
//                    (-22.4 +
//                     (m1 * m2 * (79.53333333333333 - 64 * log2)) /
//                         pow2(m1 + m2) +
//                     (64 * log2) / 3.) *
//                    pow2(sin(i)) +
//                ((m1 - m2) * cos(i) *
//                 (-0.11888020833333333 + (1891 * pow2(cos(i))) / 11520. -
//                  (7 * pow4(cos(i))) / 4608. +
//                  (pow2(m1) * pow2(m2) *
//                   (-0.2823350694444444 + (301 * pow2(cos(i))) / 2304. -
//                    (7 * pow4(cos(i))) / 1536.)) /
//                      pow4(m1 + m2) +
//                  (m1 * m2 *
//                   (3.0338541666666665 - (235 * pow2(cos(i))) / 576. +
//                    (7 * pow4(cos(i))) / 1152.)) /
//                      pow2(m1 + m2)) *
//                 sin(i) * sina) /
//                    (m1 + m2) +
//                M_PI * cos(i) *
//                    (11.333333333333334 - (8 * pow2(cos(i))) / 3. +
//                     (m1 * m2 * (-6.666666666666667 + 8 * pow2(cos(i)))) /
//                         pow2(m1 + m2)) *
//                    sin2a +
//                ((m1 - m2) * cos(i) *
//                 (4.883203125 - (12069 * pow2(cos(i))) / 1280. +
//                  (1701 * pow4(cos(i))) / 2560. +
//                  (m1 * m2 *
//                   (-30.5953125 + (7821 * pow2(cos(i))) / 320. -
//                    (1701 * pow4(cos(i))) / 640.)) /
//                      pow2(m1 + m2) +
//                  (pow2(m1) * pow2(m2) *
//                   (7.383984375 - (11403 * pow2(cos(i))) / 1280. +
//                    (5103 * pow4(cos(i))) / 2560.)) /
//                      pow4(m1 + m2)) *
//                 sin(i) * sin3a) /
//                    (m1 + m2) -
//                (32 * (1 - (3 * m1 * m2) / pow2(m1 + m2)) * M_PI * cos(i) *
//                 pow2(sin(i)) * sin4a) /
//                    3. +
//                ((m1 - m2) * cos(i) *
//                 (-22.108289930555557 + (6875 * pow2(cos(i))) / 256. -
//                  (21875 * pow4(cos(i))) / 4608. +
//                  (pow2(m1) * pow2(m2) *
//                   (-21.837022569444443 + (83125 * pow2(cos(i))) / 2304. -
//                    (21875 * pow4(cos(i))) / 1536.)) /
//                      pow4(m1 + m2) +
//                  (m1 * m2 *
//                   (58.05121527777778 - (44375 * pow2(cos(i))) / 576. +
//                    (21875 * pow4(cos(i))) / 1152.)) /
//                      pow2(m1 + m2)) *
//                 sin(i) * sin5a) /
//                    (m1 + m2) +
//                (117649 * (m1 - m2) *
//                 (1 + (3 * pow2(m1) * pow2(m2)) / pow4(m1 + m2) -
//                  (4 * m1 * m2) / pow2(m1 + m2)) *
//                 cos(i) * pow5(sin(i)) * sin7a) /
//                    (23040. * (m1 + m2)) )));}

//       else if(vpn==2){
//     //1PN term removed
//     return (2 * x *
//          ( pow3_2(x) *
//               ((((m1 - m2) * cos(i) *
//                 (0.65625 - (5 * pow2(cos(i))) / 96. +
//                  (m1 * m2 * (-0.4791666666666667 + (5 * pow2(cos(i))) / 48.))
//                  /
//                      pow2(m1 + m2)) *
//                 sin(i) * sina) /
//                    (m1 + m2) -
//                4 * M_PI * cos(i) * sin2a +
//                ((m1 - m2) * cos(i) *
//                 (-9.421875 + (135 * pow2(cos(i))) / 64. +
//                  (m1 * m2 * (5.34375 - (135 * pow2(cos(i))) / 32.)) /
//                      pow2(m1 + m2)) *
//                 sin(i) * sin3a) /
//                    (m1 + m2) +
//                (625 * (m1 - m2) * (1 - (2 * m1 * m2) / pow2(m1 + m2)) *
//                cos(i) *
//                 pow3(sin(i)) * sin5a) /
//                    (192. * (m1 + m2))) ) +
//           pow2(x) * (((m1 - m2) * cos(i) * cos3a * (9.45 - (27 * log3_2)
//           / 2.) *
//                       sin(i)) /
//                          (m1 + m2) +
//                      ((m1 - m2) * cos(i) * cosa * (-0.45 - (3 * log2) / 2.) *
//                       sin(i)) /
//                          (m1 + m2) -
//                      (3 * (m1 - m2) * M_PI * cos(i) * sin(i) * sina) /
//                          (4. * (m1 + m2)) +
//                      cos(i) *
//                          (1.1333333333333333 + (113 * pow2(cos(i))) / 30. -
//                           pow4(cos(i)) / 4. +
//                           (pow2(m1) * pow2(m2) *
//                            (-4.666666666666667 + (35 * pow2(cos(i))) / 6. -
//                             (5 * pow4(cos(i))) / 4.)) /
//                               pow4(m1 + m2) +
//                           (m1 * m2 *
//                            (15.88888888888889 - (245 * pow2(cos(i))) / 18. +
//                             (5 * pow4(cos(i))) / 4.)) /
//                               pow2(m1 + m2)) *
//                          sin2a +
//                      (27 * (m1 - m2) * M_PI * cos(i) * sin(i) * sin3a) /
//                          (4. * (m1 + m2)) +
//                      (4 * cos(i) *
//                       (55 - 12 * pow2(cos(i)) -
//                        (5 * m1 * m2 * (119 - 36 * pow2(cos(i)))) /
//                            (3. * pow2(m1 + m2)) +
//                        (5 * pow2(m1) * pow2(m2) * (17 - 12 * pow2(cos(i)))) /
//                            pow4(m1 + m2)) *
//                       pow2(sin(i)) * sin4a) /
//                          15. -
//                      (81 *
//                       (1 + (5 * pow2(m1) * pow2(m2)) / pow4(m1 + m2) -
//                        (5 * m1 * m2) / pow2(m1 + m2)) *
//                       cos(i) * pow4(sin(i)) * sin6a) /
//                          20. ) +
//           pow5_2(x) *
//               (cos(i) *
//                    (2 - (22 * pow2(cos(i))) / 5. +
//                     (m1 * m2 * (-56.4 + (94 * pow2(cos(i))) / 5.)) /
//                         pow2(m1 + m2)) *
//                    cos2a +
//                (6 * m1 * m2 * cos(i) * pow2(sin(i))) / (5. * pow2(m1 + m2)) +
//                cos(i) * cos4a *
//                    (-22.4 +
//                     (m1 * m2 * (79.53333333333333 - 64 * log2)) /
//                         pow2(m1 + m2) +
//                     (64 * log2) / 3.) *
//                    pow2(sin(i)) +
//                ((m1 - m2) * cos(i) *
//                 (-0.11888020833333333 + (1891 * pow2(cos(i))) / 11520. -
//                  (7 * pow4(cos(i))) / 4608. +
//                  (pow2(m1) * pow2(m2) *
//                   (-0.2823350694444444 + (301 * pow2(cos(i))) / 2304. -
//                    (7 * pow4(cos(i))) / 1536.)) /
//                      pow4(m1 + m2) +
//                  (m1 * m2 *
//                   (3.0338541666666665 - (235 * pow2(cos(i))) / 576. +
//                    (7 * pow4(cos(i))) / 1152.)) /
//                      pow2(m1 + m2)) *
//                 sin(i) * sina) /
//                    (m1 + m2) +
//                M_PI * cos(i) *
//                    (11.333333333333334 - (8 * pow2(cos(i))) / 3. +
//                     (m1 * m2 * (-6.666666666666667 + 8 * pow2(cos(i)))) /
//                         pow2(m1 + m2)) *
//                    sin2a +
//                ((m1 - m2) * cos(i) *
//                 (4.883203125 - (12069 * pow2(cos(i))) / 1280. +
//                  (1701 * pow4(cos(i))) / 2560. +
//                  (m1 * m2 *
//                   (-30.5953125 + (7821 * pow2(cos(i))) / 320. -
//                    (1701 * pow4(cos(i))) / 640.)) /
//                      pow2(m1 + m2) +
//                  (pow2(m1) * pow2(m2) *
//                   (7.383984375 - (11403 * pow2(cos(i))) / 1280. +
//                    (5103 * pow4(cos(i))) / 2560.)) /
//                      pow4(m1 + m2)) *
//                 sin(i) * sin3a) /
//                    (m1 + m2) -
//                (32 * (1 - (3 * m1 * m2) / pow2(m1 + m2)) * M_PI * cos(i) *
//                 pow2(sin(i)) * sin4a) /
//                    3. +
//                ((m1 - m2) * cos(i) *
//                 (-22.108289930555557 + (6875 * pow2(cos(i))) / 256. -
//                  (21875 * pow4(cos(i))) / 4608. +
//                  (pow2(m1) * pow2(m2) *
//                   (-21.837022569444443 + (83125 * pow2(cos(i))) / 2304. -
//                    (21875 * pow4(cos(i))) / 1536.)) /
//                      pow4(m1 + m2) +
//                  (m1 * m2 *
//                   (58.05121527777778 - (44375 * pow2(cos(i))) / 576. +
//                    (21875 * pow4(cos(i))) / 1152.)) /
//                      pow2(m1 + m2)) *
//                 sin(i) * sin5a) /
//                    (m1 + m2) +
//                (117649 * (m1 - m2) *
//                 (1 + (3 * pow2(m1) * pow2(m2)) / pow4(m1 + m2) -
//                  (4 * m1 * m2) / pow2(m1 + m2)) *
//                 cos(i) * pow5(sin(i)) * sin7a) /
//                    (23040. * (m1 + m2)) ) ));}

//   else if(vpn==3){
//     return (2 * x *
//          (pow3_2(x) *
//               ( -4 * M_PI * cos(i) * sin2a ) +
//           pow2(x) * (((m1 - m2) * cos(i) * cos3a * (9.45 - (27 * log3_2)
//           / 2.) *
//                       sin(i)) /
//                          (m1 + m2) +
//                      ((m1 - m2) * cos(i) * cosa * (-0.45 - (3 * log2) / 2.) *
//                       sin(i)) /
//                          (m1 + m2) -
//                      (3 * (m1 - m2) * M_PI * cos(i) * sin(i) * sina) /
//                          (4. * (m1 + m2)) +
//                      cos(i) *
//                          (1.1333333333333333 + (113 * pow2(cos(i))) / 30. -
//                           pow4(cos(i)) / 4. +
//                           (pow2(m1) * pow2(m2) *
//                            (-4.666666666666667 + (35 * pow2(cos(i))) / 6. -
//                             (5 * pow4(cos(i))) / 4.)) /
//                               pow4(m1 + m2) +
//                           (m1 * m2 *
//                            (15.88888888888889 - (245 * pow2(cos(i))) / 18. +
//                             (5 * pow4(cos(i))) / 4.)) /
//                               pow2(m1 + m2)) *
//                          sin2a +
//                      (27 * (m1 - m2) * M_PI * cos(i) * sin(i) * sin3a) /
//                          (4. * (m1 + m2)) +
//                      (4 * cos(i) *
//                       (55 - 12 * pow2(cos(i)) -
//                        (5 * m1 * m2 * (119 - 36 * pow2(cos(i)))) /
//                            (3. * pow2(m1 + m2)) +
//                        (5 * pow2(m1) * pow2(m2) * (17 - 12 * pow2(cos(i)))) /
//                            pow4(m1 + m2)) *
//                       pow2(sin(i)) * sin4a) /
//                          15. -
//                      (81 *
//                       (1 + (5 * pow2(m1) * pow2(m2)) / pow4(m1 + m2) -
//                        (5 * m1 * m2) / pow2(m1 + m2)) *
//                       cos(i) * pow4(sin(i)) * sin6a) /
//                          20. ) +
//           pow5_2(x) *
//               (cos(i) *
//                    (2 - (22 * pow2(cos(i))) / 5. +
//                     (m1 * m2 * (-56.4 + (94 * pow2(cos(i))) / 5.)) /
//                         pow2(m1 + m2)) *
//                    cos2a +
//                (6 * m1 * m2 * cos(i) * pow2(sin(i))) / (5. * pow2(m1 + m2)) +
//                cos(i) * cos4a *
//                    (-22.4 +
//                     (m1 * m2 * (79.53333333333333 - 64 * log2)) /
//                         pow2(m1 + m2) +
//                     (64 * log2) / 3.) *
//                    pow2(sin(i)) +
//                ((m1 - m2) * cos(i) *
//                 (-0.11888020833333333 + (1891 * pow2(cos(i))) / 11520. -
//                  (7 * pow4(cos(i))) / 4608. +
//                  (pow2(m1) * pow2(m2) *
//                   (-0.2823350694444444 + (301 * pow2(cos(i))) / 2304. -
//                    (7 * pow4(cos(i))) / 1536.)) /
//                      pow4(m1 + m2) +
//                  (m1 * m2 *
//                   (3.0338541666666665 - (235 * pow2(cos(i))) / 576. +
//                    (7 * pow4(cos(i))) / 1152.)) /
//                      pow2(m1 + m2)) *
//                 sin(i) * sina) /
//                    (m1 + m2) +
//                M_PI * cos(i) *
//                    (11.333333333333334 - (8 * pow2(cos(i))) / 3. +
//                     (m1 * m2 * (-6.666666666666667 + 8 * pow2(cos(i)))) /
//                         pow2(m1 + m2)) *
//                    sin2a +
//                ((m1 - m2) * cos(i) *
//                 (4.883203125 - (12069 * pow2(cos(i))) / 1280. +
//                  (1701 * pow4(cos(i))) / 2560. +
//                  (m1 * m2 *
//                   (-30.5953125 + (7821 * pow2(cos(i))) / 320. -
//                    (1701 * pow4(cos(i))) / 640.)) /
//                      pow2(m1 + m2) +
//                  (pow2(m1) * pow2(m2) *
//                   (7.383984375 - (11403 * pow2(cos(i))) / 1280. +
//                    (5103 * pow4(cos(i))) / 2560.)) /
//                      pow4(m1 + m2)) *
//                 sin(i) * sin3a) /
//                    (m1 + m2) -
//                (32 * (1 - (3 * m1 * m2) / pow2(m1 + m2)) * M_PI * cos(i) *
//                 pow2(sin(i)) * sin4a) /
//                    3. +
//                ((m1 - m2) * cos(i) *
//                 (-22.108289930555557 + (6875 * pow2(cos(i))) / 256. -
//                  (21875 * pow4(cos(i))) / 4608. +
//                  (pow2(m1) * pow2(m2) *
//                   (-21.837022569444443 + (83125 * pow2(cos(i))) / 2304. -
//                    (21875 * pow4(cos(i))) / 1536.)) /
//                      pow4(m1 + m2) +
//                  (m1 * m2 *
//                   (58.05121527777778 - (44375 * pow2(cos(i))) / 576. +
//                    (21875 * pow4(cos(i))) / 1152.)) /
//                      pow2(m1 + m2)) *
//                 sin(i) * sin5a) /
//                    (m1 + m2) +
//                (117649 * (m1 - m2) *
//                 (1 + (3 * pow2(m1) * pow2(m2)) / pow4(m1 + m2) -
//                  (4 * m1 * m2) / pow2(m1 + m2)) *
//                 cos(i) * pow5(sin(i)) * sin7a) /
//                    (23040. * (m1 + m2)) ) ));}

//   else if(vpn==4){
//     //keeping only 1.5PN and 2PN hereditary term
//     return (2 * x *
//          (pow3_2(x) *
//               ( -4 * M_PI * cos(i) * sin2a ) +
//           pow2(x) * (((m1 - m2) * cos(i) * cos3a * (9.45 - (27 * log3_2)
//           / 2.) *
//                       sin(i)) /
//                          (m1 + m2) +
//                      ((m1 - m2) * cos(i) * cosa * (-0.45 - (3 * log2) / 2.) *
//                       sin(i)) /
//                          (m1 + m2) -
//                      (3 * (m1 - m2) * M_PI * cos(i) * sin(i) * sina) /
//                          (4. * (m1 + m2))  +
//                      (27 * (m1 - m2) * M_PI * cos(i) * sin(i) * sin3a) /
//                          (4. * (m1 + m2))) +
//           pow5_2(x) *
//               (cos(i) *
//                    (2 - (22 * pow2(cos(i))) / 5. +
//                     (m1 * m2 * (-56.4 + (94 * pow2(cos(i))) / 5.)) /
//                         pow2(m1 + m2)) *
//                    cos2a +
//                (6 * m1 * m2 * cos(i) * pow2(sin(i))) / (5. * pow2(m1 + m2)) +
//                cos(i) * cos4a *
//                    (-22.4 +
//                     (m1 * m2 * (79.53333333333333 - 64 * log2)) /
//                         pow2(m1 + m2) +
//                     (64 * log2) / 3.) *
//                    pow2(sin(i)) +
//                ((m1 - m2) * cos(i) *
//                 (-0.11888020833333333 + (1891 * pow2(cos(i))) / 11520. -
//                  (7 * pow4(cos(i))) / 4608. +
//                  (pow2(m1) * pow2(m2) *
//                   (-0.2823350694444444 + (301 * pow2(cos(i))) / 2304. -
//                    (7 * pow4(cos(i))) / 1536.)) /
//                      pow4(m1 + m2) +
//                  (m1 * m2 *
//                   (3.0338541666666665 - (235 * pow2(cos(i))) / 576. +
//                    (7 * pow4(cos(i))) / 1152.)) /
//                      pow2(m1 + m2)) *
//                 sin(i) * sina) /
//                    (m1 + m2) +
//                M_PI * cos(i) *
//                    (11.333333333333334 - (8 * pow2(cos(i))) / 3. +
//                     (m1 * m2 * (-6.666666666666667 + 8 * pow2(cos(i)))) /
//                         pow2(m1 + m2)) *
//                    sin2a +
//                ((m1 - m2) * cos(i) *
//                 (4.883203125 - (12069 * pow2(cos(i))) / 1280. +
//                  (1701 * pow4(cos(i))) / 2560. +
//                  (m1 * m2 *
//                   (-30.5953125 + (7821 * pow2(cos(i))) / 320. -
//                    (1701 * pow4(cos(i))) / 640.)) /
//                      pow2(m1 + m2) +
//                  (pow2(m1) * pow2(m2) *
//                   (7.383984375 - (11403 * pow2(cos(i))) / 1280. +
//                    (5103 * pow4(cos(i))) / 2560.)) /
//                      pow4(m1 + m2)) *
//                 sin(i) * sin3a) /
//                    (m1 + m2) -
//                (32 * (1 - (3 * m1 * m2) / pow2(m1 + m2)) * M_PI * cos(i) *
//                 pow2(sin(i)) * sin4a) /
//                    3. +
//                ((m1 - m2) * cos(i) *
//                 (-22.108289930555557 + (6875 * pow2(cos(i))) / 256. -
//                  (21875 * pow4(cos(i))) / 4608. +
//                  (pow2(m1) * pow2(m2) *
//                   (-21.837022569444443 + (83125 * pow2(cos(i))) / 2304. -
//                    (21875 * pow4(cos(i))) / 1536.)) /
//                      pow4(m1 + m2) +
//                  (m1 * m2 *
//                   (58.05121527777778 - (44375 * pow2(cos(i))) / 576. +
//                    (21875 * pow4(cos(i))) / 1152.)) /
//                      pow2(m1 + m2)) *
//                 sin(i) * sin5a) /
//                    (m1 + m2) +
//                (117649 * (m1 - m2) *
//                 (1 + (3 * pow2(m1) * pow2(m2)) / pow4(m1 + m2) -
//                  (4 * m1 * m2) / pow2(m1 + m2)) *
//                 cos(i) * pow5(sin(i)) * sin7a) /
//                    (23040. * (m1 + m2)) ) ));}

//   else {
//     //keeping only 1.5PN and 2PN hereditary term
//     return (2 * x *
//          (pow3_2(x) *
//               ( -4 * M_PI * cos(i) * sin2a ) +
//           pow2(x) * (((m1 - m2) * cos(i) * cos3a * (9.45 - (27 * log3_2)
//           / 2.) *
//                       sin(i)) /
//                          (m1 + m2) +
//                      ((m1 - m2) * cos(i) * cosa * (-0.45 - (3 * log2) / 2.) *
//                       sin(i)) /
//                          (m1 + m2) -
//                      (3 * (m1 - m2) * M_PI * cos(i) * sin(i) * sina) /
//                          (4. * (m1 + m2))  +
//                      (27 * (m1 - m2) * M_PI * cos(i) * sin(i) * sin3a) /
//                          (4. * (m1 + m2))) +
//           pow5_2(x) *
//               (cos(i) *
//                    (2 - (22 * pow2(cos(i))) / 5. +
//                     (m1 * m2 * (-56.4 + (94 * pow2(cos(i))) / 5.) + (((1883 +
//                     120*sqrt(35))*cos(i) -
//      8*(7 + 24*sqrt(35))*cos(2*i) +
//      (-35 + 72*sqrt(35))*cos(3*i))/40)) /
//                         pow2(m1 + m2)) *
//                    cos2a +
//                (6 * m1 * m2 * cos(i) * pow2(sin(i))) / (5. * pow2(m1 + m2)) +
//                cos(i) * cos4a *
//                    (-22.4 +
//                     (m1 * m2 * (79.53333333333333 - 64 * log2 +
//                     72*sqrt(1.4))) /
//                         pow2(m1 + m2) +
//                     (64 * log2) / 3.) *
//                    pow2(sin(i)) +
//                M_PI * cos(i) *
//                    (11.333333333333334 - (8 * pow2(cos(i))) / 3. +
//                     (m1 * m2 * (-6.666666666666667 + 8 * pow2(cos(i)))) /
//                         pow2(m1 + m2)) *
//                    sin2a  -
//                (32 * (1 - (3 * m1 * m2) / pow2(m1 + m2)) * M_PI * cos(i) *
//                 pow2(sin(i)) * sin4a) /
//                    3. )
//           + pow3(x) * ( cos5a* ((delta*cos(i)*(565625 - 1129522*Nu +
//        437500*(-1 + 2*Nu)*log(2.5))*
//      pow3(sin(i)))/13440.) +
//           sin5a* ((-3125*delta*(-1 + 2*Nu)*M_PI*cos(i)*
//      pow3(sin(i)))/192.)  +
//           cos3a * ((delta*(3*cos(3*i)*(81*
//            (2411 + 2100*log(2) - 2100*log(3)) +
//           Nu*(-390518 - 340200*log(2) +
//              340200*log(3))) +
//        7*cos(i)*(2*Nu*
//            (-13321 + 72900*log(2) - 72900*log(3))
//            + 243*(-2861 + 660*log(2) +
//              1900*log(3) +
//              256*log(57.6650390625) -
//              256*log(1024))))*sin(i))/161280.)  +
//           sin3a * ((27*delta*M_PI*(-119 + 30*Nu +
//        (15 - 30*Nu)*cos(2*i))*sin(2*i))/256.)  +
//           cosa * ((delta*(20975 + Nu*(131794 - 4200*log(2)) +
//        50820*log(2) +
//        3*cos(2*i)*(-753 - 700*log(2) +
//           14*Nu*(167 + 100*log(2))))*sin(2*i))/
//    80640.) +
//           sina * ((delta*M_PI*(121 - 10*Nu + 5*(-1 + 2*Nu)*cos(2*i))*
//      sin(2*i))/384.) +
//           cos2a * ((856*M_PI*cos(i))/105.) +
//           sin2a * ((cos(i)*(116761 + 59920*EulerGamma - 4900*pow2(M_PI) +
//        119840*log(2) + 14980*log(pow2(x))))/
//    3675.) )));
//   }
// }

// q is the mass ratio (it can be either >1 or <1, it doesn't matter). Function
// returns the symmetric mass ratio in the range [ 0 , 0.25 ]
static REAL8 SymMassRatio(REAL8 q) { return (1.0 * q) / pow2(1.0 + 1.0 * q); }

// q<1. Eta is the symmetric mass ratio in the range [ 0 , 0.25 ]
static REAL8 SmallMassRatio(REAL8 eta) {
  return (1.0 - 2.0 * eta - sqrt(1.0 - 4.0 * eta)) / (2.0 * eta);
}

static REAL8 cosu_factor(REAL8 e, REAL8 u) { return (e * cos(u) - 1.); }

static REAL8 pn_kepler_equation(REAL8 eta, REAL8 x, REAL8 e, REAL8 l) {
  // 3PN accurate
  (void)eta;
  (void)x;
  int mean_anom_negative = 0;
  REAL8 newt_err;
  REAL8 newt_thresh;
  int newton_iterations = 0;
  REAL8 u = 0.0;
  REAL8 tol = 1.0e-12;

  /* zero mean anomaly case */
  if (l == 0) {
    return u;
  }

  /* range reduction of the l */
  while (l > M_PI) {
    l -= 2 * M_PI;
  }
  while (l < -M_PI) {
    l += 2 * M_PI;
  }

  /* solve in the positive part of the orbit */
  if (l < 0.0) {
    l = -l;
    mean_anom_negative = 1;
  }

  newt_thresh = tol * fabs(1.0 - e);
  /* use mikkola for a guess */
  u = mikkola_finder(e, l);

  /* high eccentricty case */
  if ((e > 0.8) && (l < M_PI / 3.0)) {
    REAL8 trial = l / fabs(1.0 - e);
    if (trial * trial > 6.0 * fabs(1.0 - e)) {
      /* cubic term is dominant */
      if (l < M_PI)
        trial = pow1_3(6.0 * l);
    }
    u = trial;
  }

  /* iterarte using Newton's method to get solution */
  if (e < 1.0) {
    newt_err = u - e * sin(u) - l;
    while (fabs(newt_err) > newt_thresh) {
      ++newton_iterations;
      u -= newt_err / (1.0 - e * cos(u));
      newt_err = u - e * sin(u) - l;
    }
  }
  return (mean_anom_negative ? -u : u);
}

/* function to solve Keppler's equation and return the eccentric anomaly */
/* inputs are the eccentricity and the mean anomaly                      */
static REAL8 mikkola_finder(REAL8 eccentricity, REAL8 mean_anomaly) {
  /* variables needed for Mikkola's method */
  REAL8 a, b, sgn_b, z, s;
  REAL8 sgn_mean_anomaly;
  REAL8 ecc_anomaly;

  /* range reduction of mean_anomaly  */
  while (mean_anomaly > M_PI) {
    mean_anomaly -= 2 * M_PI;
  }
  while (mean_anomaly < -M_PI) {
    mean_anomaly += 2 * M_PI;
  }

  /* compute the sign of l */
  if (mean_anomaly >= 0.0)
    sgn_mean_anomaly = 1.0;
  else {
    sgn_mean_anomaly = -1.0;
  }

  mean_anomaly *= sgn_mean_anomaly;

  /* compute alpha and beta of Minkkola Eq. (9a) */
  a = (1.0 - eccentricity) / (4.0 * eccentricity + 0.5);
  b = (0.5 * mean_anomaly) / (4.0 * eccentricity + 0.5);

  /* compute the sign of beta needed in Eq. (9b) */
  if (b >= 0.0)
    sgn_b = 1.0;
  else
    sgn_b = -1.0;

  /* Mikkola Eq. (9b) */
  z = pow1_3((b + sgn_b * sqrt(b * b + a * a * a)));
  /* Mikkola Eq. (9c) */
  s = z - a / z;
  /* add the correction given in Mikkola Eq. (7) */
  s = s - 0.078 * pow5(s) / (1.0 + eccentricity);
  /* finally Mikkola Eq. (8) gives u */
  ecc_anomaly = mean_anomaly + eccentricity * (3.0 * s - 4.0 * pow3(s));
  /* correct the sign of u */
  ecc_anomaly *= sgn_mean_anomaly;

  return (ecc_anomaly);
}
// End of the file