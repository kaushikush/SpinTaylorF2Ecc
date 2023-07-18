#include <lal/SphericalHarmonics.h>
#include <math.h>

static COMPLEX16 Complex(REAL8 real, REAL8 imag) {
  COMPLEX16 z = real + imag * I;
  return z;
}
/* All the hGO_l_m fucntions contain the 3PN non-spinning general orbit terms
 from Mishra et al. arXiv:1501.07096 and 3.5PN quasi-circular spinning
 corrections as given in Henry et al, arXiv:2209.00374v2 */

/* The (2,2) mode has newly computed 4PN non-spinning quasi-circular piece from
 * arXiv:2304.11185 */

// H22

static COMPLEX16 hGO_2_m_2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z,
                           REAL8 x, struct kepler_vars params) {
  REAL8 kappa1 = 1.0; /*for black holes kappa and lambda is 1*/
  REAL8 kappa2 = 1.0;
  REAL8 lambda1 = 1.0;
  REAL8 lambda2 = 1.0;
  REAL8 delta = sqrt(1 - 4 * Nu);
  REAL8 combination_a = (PhiDOT * r + Complex(0, 1) * rDOT);
  REAL8 combination_a2 = combination_a * combination_a;
  REAL8 combination_a3 = combination_a2 * combination_a;
  REAL8 combination_a4 = combination_a3 * combination_a;
  REAL8 combination_a5 = combination_a4 * combination_a;

  REAL8 combination_b = (PhiDOT * r - Complex(0, 1) * rDOT);
  REAL8 combination_b2 = combination_b * combination_b;
  REAL8 combination_b3 = combination_b2 * combination_b;

  
  
  

  if (vpnorder == 0) {
    return (mass / r + params.PhiDOT2 * params.r2 +
            Complex(0, 2) * PhiDOT * r * rDOT - params.rDOT2);
  }

  else if (vpnorder == 2) {
    return (
        (21 * params.Mtot2 * (-10 + Nu) -
         27 * (-1 + 3 * Nu) * params.r2 * combination_b *
             combination_a3 +
         mass * r *
             ((11 + 156 * Nu) * params.PhiDOT2 * params.r2 +
              Complex(0, 10) * (5 + 27 * Nu) * PhiDOT * r * rDOT -
              3 * (15 + 32 * Nu) * params.rDOT2)) /
        (42. * params.r2));
  }

  else if (vpnorder == 3) {
    return (
        (params.Mtot2 * (Complex(0, -1) * rDOT *
                             ((3 + 3 * delta - 8 * Nu) * S1z +
                              (3 - 3 * delta - 8 * Nu) * S2z) +
                         PhiDOT * r *
                             ((-3 - 3 * delta + 5 * Nu) * S1z +
                              (-3 + 3 * delta + 5 * Nu) * S2z))) /
            (3. * params.r2) /* (<--This is the general orbit term) */
        /* (This is the quasi-circular limit of the general orbit term-->) */
        - ((-4 * ((1 + delta - Nu) * S1z + S2z - (delta + Nu) * S2z) *
            params.x2p5) /
           3.) +
        ((-4 * (S1z + delta * S1z + S2z - delta * S2z - Nu * (S1z + S2z)) *
          params.x2p5) /
         3.)); /* (<--This is Quentins quasi-circular term) */
  }

  else if (vpnorder == 4) {
    return ((6 * params.Mtot3 * (3028 + 1267 * Nu + 158 * params.eta2) +
             9 * (83 - 589 * Nu + 1111 * params.eta2) * params.r3 *
                 combination_b2 *
                 combination_a4 +
             params.Mtot2 * r *
                 ((-11891 - 36575 * Nu + 13133 * params.eta2) * params.PhiDOT2 *
                      params.r2 +
                  Complex(0, 8) * (-773 - 3767 * Nu + 2852 * params.eta2) *
                      PhiDOT * r * rDOT -
                  6 * (-619 + 2789 * Nu + 934 * params.eta2) * params.rDOT2) -
             3 * mass * params.r2 *
                 (2 * (-835 - 19 * Nu + 2995 * params.eta2) * params.PhiDOT4 *
                      params.r4 +
                  Complex(0, 6) * (-433 - 721 * Nu + 1703 * params.eta2) *
                      params.PhiDOT3 * params.r3 * rDOT +
                  6 * (-33 + 1014 * Nu + 232 * params.eta2) * params.PhiDOT2 *
                      params.r2 * params.rDOT2 +
                  Complex(0, 4) * (-863 + 1462 * Nu + 2954 * params.eta2) *
                      PhiDOT * r * params.rDOT3 -
                  3 * (-557 + 664 * Nu + 1712 * params.eta2) * params.rDOT4)) /
                (1512. * params.r3) +
            (3 * params.Mtot3 *
             (S1z * (4 * Nu * S2z + (1 + delta - 2 * Nu) * S1z * kappa1) -
              (-1 + delta + 2 * Nu) * params.S2z2 * kappa2)) /
                (4. * params.r3) -
            ((kappa1 * (1 + delta - 2 * Nu) * params.S1z2 +
              S2z * (4 * Nu * S1z - kappa2 * (-1 + delta + 2 * Nu) * S2z)) *
             params.x3) +
            ((kappa1 * (1 + delta - 2 * Nu) * params.S1z2 +
              S2z * (4 * Nu * S1z - kappa2 * (-1 + delta + 2 * Nu) * S2z)) *
             params.x3));
  }

  else if (vpnorder == 5) {
    return (
        (params.Mtot2 * Nu *
         (2 * mass * (Complex(0, -702) * PhiDOT * r + rDOT) +
          3 * r *
              (Complex(0, -316) * params.PhiDOT3 * params.r3 -
               847 * params.PhiDOT2 * params.r2 * rDOT +
               Complex(0, 184) * PhiDOT * r * params.rDOT2 -
               122 * params.rDOT3))) /
            (105. * params.r3) /* Henry et al. QC spin terms */ /* +
((2*(56*delta*Nu*(-S1z + S2z)
+ 101*Nu*(S1z + S2z) + 132*params.eta2*(S1z + S2z) - 80*(S1z + delta*S1z + S2z -
delta*S2z))*params.x3p5)/63.) */
        +
        /* Henry et al. ecc spin terms */ (
            (params.Mtot2 *
             (mass *
                  ((238 + delta * (238 - 141 * Nu) + Nu * (-181 + 474 * Nu)) *
                       PhiDOT * r * S1z +
                   Complex(0, 8) *
                       (55 + delta * (55 - 19 * Nu) +
                        2 * Nu * (-50 + 43 * Nu)) *
                       rDOT * S1z +
                   (238 + delta * (-238 + 141 * Nu) + Nu * (-181 + 474 * Nu)) *
                       PhiDOT * r * S2z +
                   Complex(0, 8) *
                       (55 + delta * (-55 + 19 * Nu) +
                        2 * Nu * (-50 + 43 * Nu)) *
                       rDOT * S2z) -
              r * (PhiDOT * r * params.rDOT2 *
                       (-((18 * (1 + delta) + 5 * (-63 + 55 * delta) * Nu +
                           188 * params.eta2) *
                          S1z) +
                        (18 * (-1 + delta) + 5 * (63 + 55 * delta) * Nu -
                         188 * params.eta2) *
                            S2z) -
                   Complex(0, 2) * params.rDOT3 *
                       ((-27 * (1 + delta) + 6 * (5 + 7 * delta) * Nu -
                         4 * params.eta2) *
                            S1z +
                        (-27 + 27 * delta + 30 * Nu - 42 * delta * Nu -
                         4 * params.eta2) *
                            S2z) +
                   Complex(0, 2) * params.PhiDOT2 * params.r2 * rDOT *
                       ((51 + 88 * Nu * (-3 + 5 * Nu) +
                         delta * (51 + 62 * Nu)) *
                            S1z +
                        (51 + 88 * Nu * (-3 + 5 * Nu) -
                         delta * (51 + 62 * Nu)) *
                            S2z) +
                   params.PhiDOT3 * params.r3 *
                       ((120 * (1 + delta) + (-483 + 83 * delta) * Nu +
                         234 * params.eta2) *
                            S1z +
                        (120 + 3 * Nu * (-161 + 78 * Nu) -
                         delta * (120 + 83 * Nu)) *
                            S2z)))) /
            (84. * params.r3)));
  }

  else if (vpnorder == 6) {
    return (
        (4 * params.Mtot4 *
             (-8203424 + 2180250 * params.eta2 + 592600 * params.eta3 +
              15 * Nu * (-5503804 + 142065 * M_PI2)) -
         2700 * (-507 + 6101 * Nu - 25050 * params.eta2 + 34525 * params.eta3) *
             params.r4 * combination_b3 *
             combination_a5 +
         params.Mtot3 * r *
             (params.PhiDOT2 *
                  (337510808 - 198882000 * params.eta2 + 56294600 * params.eta3 +
                   Nu * (183074880 - 6392925 * M_PI2)) *
                  params.r2 +
              Complex(0, 110) * PhiDOT *
                  (-5498800 - 785120 * params.eta2 + 909200 * params.eta3 +
                   3 * Nu * (-1849216 + 38745 * M_PI2)) *
                  r * rDOT +
              2 *
                  (51172744 - 94929000 * params.eta2 - 5092400 * params.eta3 +
                   45 * Nu * (2794864 + 142065 * M_PI2)) *
                  params.rDOT2) -
         20 * params.Mtot2 * params.r2 *
             ((-986439 + 1873255 * Nu - 9961400 * params.eta2 +
               6704345 * params.eta3) *
                  params.PhiDOT4 * params.r4 +
              Complex(0, 4) *
                  (-273687 - 978610 * Nu - 4599055 * params.eta2 +
                   2783005 * params.eta3) *
                  params.PhiDOT3 * params.r3 * rDOT +
              (-181719 + 19395325 * Nu + 8237980 * params.eta2 +
               2612735 * params.eta3) *
                  params.PhiDOT2 * params.r2 * params.rDOT2 +
              Complex(0, 8) *
                  (-234312 + 1541140 * Nu + 1230325 * params.eta2 +
                   1828625 * params.eta3) *
                  PhiDOT * r * params.rDOT3 -
              3 *
                  (-370268 + 1085140 * Nu + 2004715 * params.eta2 +
                   1810425 * params.eta3) *
                  params.rDOT4) +
         300 * mass * params.r3 *
             (4 *
                  (12203 - 36427 * Nu - 27334 * params.eta2 +
                   149187 * params.eta3) *
                  params.PhiDOT6 * params.r6 +
              Complex(0, 2) *
                  (44093 - 68279 * Nu - 295346 * params.eta2 +
                   541693 * params.eta3) *
                  params.PhiDOT5 * params.r5 * rDOT +
              2 *
                  (27432 - 202474 * Nu + 247505 * params.eta2 +
                   394771 * params.eta3) *
                  params.PhiDOT4 * params.r4 * params.rDOT2 +
              Complex(0, 2) *
                  (97069 - 383990 * Nu - 8741 * params.eta2 +
                   1264800 * params.eta3) *
                  params.PhiDOT3 * params.r3 * params.rDOT3 +
              (-42811 + 53992 * Nu + 309136 * params.eta2 -
               470840 * params.eta3) *
                  params.PhiDOT2 * params.r2 * params.rDOT4 +
              Complex(0, 2) *
                  (51699 - 252256 * Nu + 131150 * params.eta2 +
                   681160 * params.eta3) *
                  PhiDOT * r * params.rDOT5 -
              3 *
                  (16743 - 75104 * Nu + 26920 * params.eta2 +
                   207200 * params.eta3) *
                  params.rDOT6)) /
            (3.3264e6 * params.r4) /* Henry et al. QC spin terms */ /* + (((4*(1
+ delta)*(-7 + 9*kappa1)
- 7*(9 + 17*delta)*Nu - 9*(15 + 7*delta)*kappa1*Nu + 12*(7 -
17*kappa1)*params.eta2)*params.S1z2 + 2*S1z*(Complex(0,-42)*(1 + delta - 2*Nu) -
84*(1 + delta - Nu)*M_PI + Nu*(-271 + 288*Nu)*S2z) + S2z*(12*(7 -
17*kappa2)*params.eta2*S2z + 4*(-1 + delta)*(Complex(0,21) + 42*M_PI + 7*S2z -
9*kappa2*S2z) + Nu*(168*(Complex(0,1) + M_PI) + 7*delta*(17 + 9*kappa2)*S2z -
9*(7 + 15*kappa2)*S2z)))*params.x4)/63. */
        + /* Henry et al. ecc spin terms */ (
              -0.005952380952380952 *
              (params.Mtot3 *
               (2 * mass *
                    (S1z * (Complex(0, 14) * (1 + delta - 2 * Nu) +
                            42 * (1 + delta - 2 * Nu) * Nu * S1z +
                            kappa1 *
                                (438 + delta * (438 + 103 * Nu) +
                                 Nu * (-773 + 108 * Nu)) *
                                S1z) +
                     2 *
                         (Complex(0, -7) * (-1 + delta + 2 * Nu) +
                          (995 - 192 * Nu) * Nu * S1z) *
                         S2z -
                     (42 * Nu * (-1 + delta + 2 * Nu) +
                      kappa2 * (-438 + (773 - 108 * Nu) * Nu +
                                delta * (438 + 103 * Nu))) *
                         params.S2z2) +
                r * (params.rDOT2 *
                         (Complex(0, 56) * (1 + delta - 2 * Nu) * S1z +
                          (-56 * (1 + delta - 2 * Nu) * Nu +
                           kappa1 * (291 * (1 + delta) -
                                     2 * (445 + 154 * delta) * Nu +
                                     24 * params.eta2)) *
                              params.S1z2 +
                          4 * Nu * (-3 + 44 * Nu) * S1z * S2z +
                          S2z * (Complex(0, -56) * (-1 + delta + 2 * Nu) +
                                 56 * Nu * (-1 + delta + 2 * Nu) * S2z +
                                 kappa2 *
                                     (291 - 890 * Nu + 24 * params.eta2 +
                                      delta * (-291 + 308 * Nu)) *
                                     S2z)) +
                     params.PhiDOT2 * params.r2 *
                         (Complex(0, 196) * (1 + delta - 2 * Nu) * S1z +
                          (56 * Nu * (7 + 7 * delta + Nu) +
                           kappa1 * (-153 * (1 + delta) -
                                     2 * (62 + 215 * delta) * Nu +
                                     804 * params.eta2)) *
                              params.S1z2 +
                          8 * (60 - 187 * Nu) * Nu * S1z * S2z +
                          S2z * (Complex(0, -196) * (-1 + delta + 2 * Nu) +
                                 56 * Nu * (7 - 7 * delta + Nu) * S2z +
                                 kappa2 *
                                     (-153 + 4 * Nu * (-31 + 201 * Nu) +
                                      delta * (153 + 430 * Nu)) *
                                     S2z)) +
                     2 * PhiDOT * r * rDOT *
                         (Complex(0, -1) *
                              (117 * (1 + delta) * kappa1 -
                               434 * (1 + delta) * Nu +
                               (-23 + 211 * delta) * kappa1 * Nu +
                               2 * (14 - 195 * kappa1) * params.eta2) *
                              params.S1z2 +
                          4 * S1z *
                              (35 * (1 + delta - 2 * Nu) +
                               Complex(0, 1) * (9 - 209 * Nu) * Nu * S2z) +
                          S2z * (-140 * (-1 + delta + 2 * Nu) +
                                 Complex(0, 1) *
                                     (-14 * Nu * (-31 + 31 * delta + 2 * Nu) +
                                      kappa2 * (117 * (-1 + delta) +
                                                (23 + 211 * delta) * Nu +
                                                390 * params.eta2)) *
                                     S2z))))) /
              params.r4) + /* Henry et al. QC spinning hereditary terms */
        (((-8 * M_PI * ((1 + delta - Nu) * S1z + S2z - (delta + Nu) * S2z) *
           params.x4) /
          3.)));
  } else if (vpnorder == 7) {

    return (
        /* Henry et al QC spin terms */ /* ((3318*params.eta3*(S1z + S2z) +
    Nu*(-504*((7 + delta)*kappa1 - 3*(3 + delta)*lambda1)*params.S1z3 -
    1008*params.S1z2*(3*kappa1*M_PI - 3*(1 + delta)*S2z + 2*(1 +
    delta)*kappa1*S2z) + S1z*(17387 + 20761*delta + 1008*S2z*(6*M_PI + (-1 +
    delta)*(-3 + 2*kappa2)*S2z)) + S2z*(17387 - 20761*delta +
    504*S2z*(-6*kappa2*M_PI + (-7 + delta)*kappa2*S2z - 3*(-3 +
    delta)*lambda2*S2z))) + 2*(2809*(1 + delta)*S1z + 756*(1 +
    delta)*kappa1*M_PI*params.S1z2 + 756*(1 + delta)*(kappa1 -
    lambda1)*params.S1z3 -
    (-1 + delta)*S2z*(2809 + 756*S2z*(-(lambda2*S2z) + kappa2*(M_PI + S2z)))) -
    2*params.eta2*(708*delta*(-S1z + S2z) + (S1z + S2z)*(4427 +
    1008*(kappa1*params.S1z2 + S2z*(-2*S1z + kappa2*S2z)))))*params.x4p5)/756. */
        /* + */ /* Henry et al. ecc+spin terms */ (
            (params.Mtot2 *
             (-3 * mass * r *
                  (Complex(0, -16) * params.rDOT3 *
                       (Complex(0, 12) * Nu * (-16703 + 4427 * Nu) +
                        35 * Nu *
                            (4578 + Nu * (-4288 + 765 * Nu) +
                             delta * (3748 + 802 * Nu)) *
                            S1z +
                        35 *
                            (delta * (942 - 2 * Nu * (1874 + 401 * Nu)) +
                             Nu * (4578 + Nu * (-4288 + 765 * Nu))) *
                            S2z -
                        32970 * (S1z + delta * S1z + S2z)) +
                   4 * PhiDOT * r * params.rDOT2 *
                       (-338520 * params.eta3 * (S1z + S2z) +
                        48930 * (S1z + delta * S1z + S2z - delta * S2z) +
                        Nu * (Complex(0, 3420696) -
                              35 * (14154 + 21167 * delta) * S1z -
                              495390 * S2z + 740845 * delta * S2z) +
                        params.eta2 * (Complex(0, -612336) +
                                      245 * (3566 - 1565 * delta) * S1z +
                                      245 * (3566 + 1565 * delta) * S2z)) +
                   params.PhiDOT3 * params.r3 *
                       (2515380 * params.eta3 * (S1z + S2z) -
                        5 * Nu *
                            (Complex(0, 1859936) +
                             7 * (82329 + 37061 * delta) * S1z +
                             7 * (82329 - 37061 * delta) * S2z) -
                        128100 * (S1z + delta * S1z + S2z - delta * S2z) +
                        4 * params.eta2 *
                            (Complex(0, 381348) +
                             35 * (-18505 + 1777 * delta) * S1z -
                             35 * (18505 + 1777 * delta) * S2z)) +
                   Complex(0, 8) * params.PhiDOT2 * params.r2 * rDOT *
                       (779100 * params.eta3 * (S1z + S2z) +
                        5 * Nu *
                            (Complex(0, 828806) +
                             7 * (4839 + 5971 * delta) * S1z +
                             7 * (4839 - 5971 * delta) * S2z) -
                        62475 * (S1z + delta * S1z + S2z - delta * S2z) +
                        params.eta2 * (Complex(0, -976002) +
                                      35 * (-29599 + 3109 * delta) * S1z -
                                      35 * (29599 + 3109 * delta) * S2z))) +
              3 * params.r2 *
                  (Complex(0, 4) * params.PhiDOT2 * params.r2 * params.rDOT3 *
                       (Complex(0, 2) * (65451 - 350563 * Nu) * Nu +
                        105 * Nu *
                            (6020 + delta * (3114 + 411 * Nu) +
                             Nu * (-4513 + 9136 * Nu)) *
                            S1z +
                        105 *
                            (-3 * delta * (-408 + Nu * (1038 + 137 * Nu)) +
                             Nu * (6020 + Nu * (-4513 + 9136 * Nu))) *
                            S2z -
                        128520 * (S1z + delta * S1z + S2z)) +
                   PhiDOT * r * params.rDOT4 *
                       (Complex(0, -128) * Nu * (2487 + 18334 * Nu) -
                        105 * Nu *
                            (8689 + 8 * Nu * (-687 + 1402 * Nu) +
                             delta * (2143 + 8212 * Nu)) *
                            S1z +
                        105 *
                            (Nu * (-8689 + 8 * (687 - 1402 * Nu) * Nu) +
                             delta * (-1470 + Nu * (2143 + 8212 * Nu))) *
                            S2z +
                        154350 * (S1z + delta * S1z + S2z)) -
                   Complex(0, 6) * params.rDOT5 *
                       (-11760 * params.eta3 * (S1z + S2z) +
                        4 * params.eta2 *
                            (Complex(0, 57338) +
                             35 * (569 + 301 * delta) * S1z +
                             35 * (569 - 301 * delta) * S2z) -
                        16 * Nu *
                            (Complex(0, -2517) + 35 * (72 + 71 * delta) * S1z +
                             35 * (72 - 71 * delta) * S2z) +
                        8715 * (S1z + delta * S1z + S2z - delta * S2z)) +
                   params.PhiDOT5 * params.r5 *
                       (2263380 * params.eta3 * (S1z + S2z) +
                        3 * Nu *
                            (Complex(0, 653432) +
                             35 * (13283 + 7839 * delta) * S1z +
                             35 * (13283 - 7839 * delta) * S2z) -
                        219240 * (S1z + delta * S1z + S2z - delta * S2z) +
                        4 * params.eta2 *
                            (Complex(0, -291268) +
                             105 * (-7669 + 165 * delta) * S1z -
                             105 * (7669 + 165 * delta) * S2z)) +
                   Complex(0, 14) * params.PhiDOT4 * params.r4 * rDOT *
                       (385170 * params.eta3 * (S1z + S2z) +
                        15 * Nu *
                            (Complex(0, 16914) + (8633 + 2267 * delta) * S1z +
                             8633 * S2z - 2267 * delta * S2z) -
                        18630 * (S1z + delta * S1z + S2z - delta * S2z) +
                        2 * params.eta2 *
                            (Complex(0, -67904) +
                             15 * (-13932 + 679 * delta) * S1z -
                             15 * (13932 + 679 * delta) * S2z)) +
                   6 * params.PhiDOT3 * params.r3 * params.rDOT2 *
                       (-6720 * params.eta3 * (S1z + S2z) +
                        5 * Nu *
                            (Complex(0, -241704) +
                             7 * (4377 + 3083 * delta) * S1z +
                             7 * (4377 - 3083 * delta) * S2z) -
                        36960 * (S1z + delta * S1z + S2z - delta * S2z) +
                        params.eta2 * (Complex(0, -364384) +
                                      35 * (5059 - 4753 * delta) * S1z +
                                      35 * (5059 + 4753 * delta) * S2z))) +
              14 * params.Mtot2 *
                  (Complex(0, 30) * rDOT *
                       (15280 * params.eta3 * (S1z + S2z) -
                        4 * params.eta2 *
                            (Complex(0, -111) + 3562 * S1z + 682 * delta * S1z +
                             945 * kappa1 * params.S1z3 +
                             (3562 - 682 * delta +
                              945 * (-2 + kappa1) * params.S1z2) *
                                 S2z +
                             945 * (-2 + kappa2) * S1z * params.S2z2 +
                             945 * kappa2 * params.S2z3) +
                        Nu * (Complex(0, -27520) +
                              378 *
                                  (5 * (1 + delta) * kappa1 +
                                   3 * (3 + delta) * lambda1) *
                                  params.S1z3 +
                              (29749 - 13605 * delta) * S2z -
                              1512 * (1 + delta) * kappa1 * params.S1z2 * S2z -
                              378 *
                                  (5 * (-1 + delta) * kappa2 +
                                   3 * (-3 + delta) * lambda2) *
                                  params.S2z3 +
                              S1z * (29749 + 13605 * delta +
                                     1512 * (-1 + delta) * kappa2 *
                                         params.S2z2)) +
                        2 * (-8009 * (1 + delta) * S1z -
                             567 * (1 + delta) * lambda1 * params.S1z3 +
                             (-1 + delta) * S2z *
                                 (8009 + 567 * lambda2 * params.S2z2))) +
                   PhiDOT * r *
                       (285840 * params.eta3 * (S1z + S2z) -
                        30 * params.eta2 *
                            (Complex(0, 11504) + 1890 * kappa1 * params.S1z3 +
                             23090 * S2z - 1823 * delta * S2z +
                             1890 * (-2 + kappa1) * params.S1z2 * S2z +
                             1890 * kappa2 * params.S2z3 +
                             S1z * (23090 + 1823 * delta +
                                    1890 * (-2 + kappa2) * params.S2z2)) +
                        30 * (689 * (1 + delta) * S1z -
                              1134 * (1 + delta) * lambda1 * params.S1z3 +
                              (-1 + delta) * S2z *
                                  (-689 + 1134 * lambda2 * params.S2z2)) +
                        2 * Nu *
                            (Complex(0, 415432) - 66840 * S2z +
                             15 * (8 * (-557 + 1342 * delta) * S1z +
                                   189 *
                                       (5 * (1 + delta) * kappa1 +
                                        6 * (3 + delta) * lambda1) *
                                       params.S1z3 -
                                   10736 * delta * S2z -
                                   2457 * (1 + delta) * kappa1 * params.S1z2 *
                                       S2z +
                                   2457 * (-1 + delta) * kappa2 * S1z *
                                       params.S2z2 -
                                   189 *
                                       (5 * (-1 + delta) * kappa2 +
                                        6 * (-3 + delta) * lambda2) *
                                       params.S2z3)))))) /
            (317520. *
             params.r4)) + /* Henry et al. QC spinning hereditary terms */
        (2 * M_PI *
         (kappa1 * (1 + delta - 2 * Nu) * params.S1z2 +
          S2z * (4 * Nu * S1z - kappa2 * (-1 + delta + 2 * Nu) * S2z)) *
         params.x4p5));
  }

  else {
    return 0;
  }
}
// hQC_l_m() functions contains only the non-spinning hereditary terms at
// particular PN order.

static COMPLEX16 hQC_2_m_2(REAL8 Nu, UINT4 vpnorder, REAL8 x, struct kepler_vars params) {
  double EulerGamma = 0.5772156649015329;
  // keeping only the hereditary terms
  // if(vpnorder == 0){
  //   return(2*x);
  // }

  /* else if(vpnorder == 2){
     return((-5.095238095238095 + (55*Nu)/21.)*params.x2);
   } */

  /* else */ if (vpnorder == 3) {
    return (4 * M_PI * params.x2p5);
  }

  /* else if(vpnorder == 4){
      return((-2.874338624338624 - (1069*Nu)/108. + (2047*params.eta2)/756.)*
      params.x3);
   } */

  else if (vpnorder == 5) {
    // return((Complex(0,-48)*Nu - (214*M_PI)/21. +
    // (68*Nu*M_PI)/21.)*params.x3p5);
    return ((2 * (-107 + 34 * Nu) * M_PI * params.x3p5) / 21.);
  }

  else if (vpnorder == 6) {
    return (
        (params.x4 * (-27392 * EulerGamma +
                      M_PI * (Complex(0, 13696) + 35 * (64 + 41 * Nu) * M_PI) -
                      13696 * log(16 * x))) /
        1680.);
  }

  else if (vpnorder == 7) {
    return (((-2173 - 4990 * Nu + 1120 * params.eta2) * M_PI * params.x4p5) /
            378.);
  }

  // 4PN non-spinning quasi-circular (2,2) mode has been obtained from Blanchet
  // et al. arXiv:2304.11185

  else if (vpnorder == 8) {
    return (
        (params.x5 *
         (276756480 * EulerGamma * (11449 + 19105 * Nu) -
          12 * (846557506853 + 1008017482431 * Nu) +
          35 * (28 * params.eta2 *
                    (5385456111 + 5 * Nu * (-163158374 + 26251249 * Nu)) -
                Complex(0, 3953664) * (11449 + 109657 * Nu) * M_PI -
                135135 * (54784 + 5 * Nu * (1951 + 6560 * Nu)) * M_PI2) +
          138378240 * (11449 + 19105 * Nu) * log(16 * x))) /
        7.62810048e10);
  }

  else {
    return 0;
  }
}

static COMPLEX16 hl_2_m_2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z, REAL8 x, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_2_m_2: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_2_m_2(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x, params) +
            hQC_2_m_2(Nu, vpnorder, x, params)) *
           cpolar(1, -2 * Phi);
  }
}

static COMPLEX16 hl_2_m_min2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z, REAL8 x, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_2_m_min2: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R)  *
           conj((hGO_2_m_2(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x, params) +
                 hQC_2_m_2(Nu, vpnorder, x, params))) *
           cpolar(1, 2 * Phi);
  }
}

// H21

static COMPLEX16 hGO_2_m_1(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z,
                           REAL8 x, struct kepler_vars params) {
  REAL8 delta = sqrt(1 - 4 * Nu);
  REAL8 kappa1 = 1.0;
  REAL8 kappa2 = 1.0;
  REAL8 r0 = 1.0;

  if (vpnorder == 1) {
    return (Complex(0, 0.6666666666666666) * delta * mass * PhiDOT);
  }

  else if (vpnorder == 2) {
    return ((Complex(0, -0.5) * params.Mtot2 *
             ((1 + delta) * S1z + (-1 + delta) * S2z)) /
                params.r2 -
            (Complex(0, -0.5) * (S1z - S2z + delta * (S1z + S2z)) * params.x2) +
            (Complex(0, -0.5) * (S1z - S2z + delta * (S1z + S2z)) * params.x2));
  }

  else if (vpnorder == 3) {
    return ((Complex(0, 0.023809523809523808) * delta * mass * PhiDOT *
             (4 * mass * (-9 + 11 * Nu) +
              r * ((19 - 24 * Nu) * params.PhiDOT2 * params.r2 +
                   Complex(0, 2) * (83 + 2 * Nu) * PhiDOT * r * rDOT +
                   2 * (-33 + 10 * Nu) * params.rDOT2))) /
            r);
  }

  else if (vpnorder == 4) {
    return (
        (Complex(0, 0.011904761904761904) * params.Mtot2 *
         (2 * mass *
              ((77 + 59 * Nu + 11 * delta * (7 + Nu)) * S1z +
               (-77 - 59 * Nu + 11 * delta * (7 + Nu)) * S2z) +
          r * (Complex(0, -2) * PhiDOT * r * rDOT *
                   (147 * (1 + delta) * S1z + (-83 + 13 * delta) * Nu * S1z +
                    147 * (-1 + delta) * S2z + (83 + 13 * delta) * Nu * S2z) +
               params.rDOT2 *
                   ((105 * (1 + delta) - 4 * (13 + 15 * delta) * Nu) * S1z +
                    (-105 + 15 * delta * (7 - 4 * Nu) + 52 * Nu) * S2z) +
               4 * params.PhiDOT2 * params.r2 *
                   ((-21 - 21 * delta + 66 * Nu + 4 * delta * Nu) * S1z +
                    (21 - 21 * delta - 66 * Nu + 4 * delta * Nu) * S2z)))) /
            params.r3 -
        (Complex(0, 0.023809523809523808) *
         ((-7 + 205 * Nu + delta * (-7 + 33 * Nu)) * S1z +
          (7 - 205 * Nu + delta * (-7 + 33 * Nu)) * S2z) *
         params.x3) +
        (Complex(0, 0.023809523809523808) * params.x3 *
         ((-7 + 205 * Nu + delta * (-7 + 33 * Nu)) * S1z +
          (7 - 205 * Nu + delta * (-7 + 33 * Nu)) * S2z)));
  }

  else if (vpnorder == 5) {
    return (
        (Complex(0, 0.0013227513227513227) * delta * mass * PhiDOT *
         (10 * params.Mtot2 * (31 - 205 * Nu + 111 * params.eta2) -
          2 * mass * r *
              ((-197 + 5 * Nu + 660 * params.eta2) * params.PhiDOT2 * params.r2 +
               Complex(0, 1) * (-3167 - 5278 * Nu + 201 * params.eta2) * PhiDOT *
                   r * rDOT +
               8 * (202 + 587 * Nu - 177 * params.eta2) * params.rDOT2) +
          3 * params.r2 *
              ((152 - 692 * Nu + 333 * params.eta2) * params.PhiDOT4 *
                   params.r4 +
               Complex(0, 2) * (308 - 1607 * Nu + 111 * params.eta2) *
                   params.PhiDOT3 * params.r3 * rDOT -
               3 * (75 - 560 * Nu + 68 * params.eta2) * params.PhiDOT2 *
                   params.r2 * params.rDOT2 -
               Complex(0, 2) * (-265 + 526 * Nu + 18 * params.eta2) * PhiDOT *
                   r * params.rDOT3 +
               (-241 + 550 * Nu - 264 * params.eta2) * params.rDOT4))) /
            params.r2
        /* Henry et al QC spin terms */ /* +
  (Complex(0,-0.08333333333333333)*params.x3p5*(2*(-((1 + delta)*(5 + kappa1))
  + 2*(6 + delta + (4 + 3*delta)*kappa1)*Nu)*params.S1z2 + S1z*(Complex(0,-3) -
  Complex(0,3)*delta + 6*(1 + delta)*M_PI
  - 16*delta*Nu*S2z - Complex(0,3)*(1 + delta)*log(16)) + S2z*(6*(-1 +
  delta)*M_PI - 2*(-1 + delta)*(5 + kappa2)*S2z
  + 4*(-6 + delta - 4*kappa2 + 3*delta*kappa2)*Nu*S2z - Complex(0,3)*(-1 +
  delta)*(1 + log(16))))) */
        +
        /* Henry et al. ecc spin terms */ (
            (params.Mtot3 *
             (-4 * rDOT *
                  (kappa1 * (1 + delta - 2 * Nu) * params.S1z2 +
                   kappa2 * (-1 + delta + 2 * Nu) * params.S2z2) -
              Complex(0, 1) * PhiDOT * r *
                  ((-((1 + delta) * (9 + kappa1)) +
                    2 * (9 + (4 + 3 * delta) * kappa1) * Nu) *
                       params.S1z2 -
                   12 * delta * Nu * S1z * S2z +
                   (9 + kappa2 - 2 * (9 + 4 * kappa2) * Nu +
                    delta * (-9 - kappa2 + 6 * kappa2 * Nu)) *
                       params.S2z2))) /
            (6. * params.r3)) + /* Henry et al. QC spinning hereditary terms */
        (Complex(0, -0.5) * ((1 + delta) * S1z + (-1 + delta) * S2z) *
         params.x3p5 * (M_PI - Complex(0, 2) * log(2))));
  }

  else if (vpnorder == 6) {
    return (
        (delta * params.Mtot2 * Nu * PhiDOT *
         (mass * (195 * PhiDOT * r - Complex(0, 946) * rDOT) +
          9 * r *
              (270 * params.PhiDOT3 * params.r3 -
               Complex(0, 483) * params.PhiDOT2 * params.r2 * rDOT -
               580 * PhiDOT * r * params.rDOT2 +
               Complex(0, 42) * params.rDOT3))) /
            (315. * params.r2)
        /* Henry et al. QC spin terms */
        /* +(Complex(0,-0.0006613756613756613)*(-4*params.eta2*(-1365*(S1z - S2z)
   + 179*delta*(S1z + S2z)) + 6*(208*(1 + delta)*S1z + 63*(1 +
delta)*kappa1*params.S1z3 + (-1 + delta)*S2z*(208 + 63*kappa2*params.S2z2)) -
Nu*(378*(3 + delta)*kappa1*params.S1z3 + 378*(1 + delta)*(-2 +
kappa1)*params.S1z2*S2z + S1z*(8351 + 7027*delta + 378*(-1 + delta)*(-2 +
kappa2)*params.S2z2) + S2z*(-8351 + 7027*delta + 378*(-3 +
delta)*kappa2*params.S2z2)))*params.x4) */
        +
        /* Henry et al. ecc spin terms */ (
            Complex(0, 0.00033068783068783067) * params.Mtot2 *
            (3 * params.r2 *
                 (Complex(0, 4) * PhiDOT * r * params.rDOT3 *
                      ((-315 * (1 + delta) + 2 * (251 + 463 * delta) * Nu +
                        4 * (15 + delta) * params.eta2) *
                           S1z +
                       (315 - 315 * delta - 502 * Nu + 926 * delta * Nu +
                        4 * (-15 + delta) * params.eta2) *
                           S2z) +
                  12 * params.PhiDOT2 * params.r2 * params.rDOT2 *
                      ((189 * (1 + delta) - 2 * (521 + 293 * delta) * Nu +
                        7 * (55 + 23 * delta) * params.eta2) *
                           S1z +
                       (189 * (-1 + delta) + 2 * (521 - 293 * delta) * Nu +
                        7 * (-55 + 23 * delta) * params.eta2) *
                           S2z) +
                  params.rDOT4 *
                      ((567 * (1 + delta) - 16 * (77 + 64 * delta) * Nu +
                        8 * (177 + 173 * delta) * params.eta2) *
                           S1z +
                       (567 * (-1 + delta) + 16 * (77 - 64 * delta) * Nu +
                        8 * (-177 + 173 * delta) * params.eta2) *
                           S2z) -
                  Complex(0, 4) * params.PhiDOT3 * params.r3 * rDOT *
                      ((936 * (1 + delta) - 5 * (979 + 215 * delta) * Nu +
                        2 * (1353 + 293 * delta) * params.eta2) *
                           S1z +
                       (936 * (-1 + delta) + 5 * (979 - 215 * delta) * Nu +
                        2 * (-1353 + 293 * delta) * params.eta2) *
                           S2z) +
                  4 * params.PhiDOT4 * params.r4 *
                      ((-252 * (1 + delta) + (1315 + 857 * delta) * Nu +
                        4 * (-285 + 43 * delta) * params.eta2) *
                           S1z +
                       (252 + 5 * Nu * (-263 + 228 * Nu) +
                        delta * (-252 + Nu * (857 + 172 * Nu))) *
                           S2z)) -
             2 * mass * r *
                 (Complex(0, -1) * PhiDOT * r * rDOT *
                      ((2043 * (1 + delta) + (37 + 2597 * delta) * Nu +
                        (10635 + 139 * delta) * params.eta2) *
                           S1z +
                       (2043 * (-1 + delta) + (-37 + 2597 * delta) * Nu +
                        (-10635 + 139 * delta) * params.eta2) *
                           S2z) +
                  params.PhiDOT2 * params.r2 *
                      ((-765 - Nu * (667 + 7773 * Nu) +
                        delta * (-765 + 7 * Nu * (-533 + 245 * Nu))) *
                           S1z +
                       (765 + Nu * (667 + 7773 * Nu) +
                        delta * (-765 + 7 * Nu * (-533 + 245 * Nu))) *
                           S2z) +
                  4 * params.rDOT2 *
                      ((-234 * (1 + delta) - 4 * (560 + 901 * delta) * Nu +
                        (483 + 1111 * delta) * params.eta2) *
                           S1z +
                       (234 + 7 * (320 - 69 * Nu) * Nu +
                        delta * (-234 + Nu * (-3604 + 1111 * Nu))) *
                           S2z)) +
             2 * params.Mtot2 *
                 (1134 * kappa1 * (-1 - delta + (3 + delta) * Nu) *
                      params.S1z3 +
                  1134 * (1 + delta) * (-2 + kappa1) * Nu * params.S1z2 * S2z +
                  S1z *
                      (-5661 - 5661 * delta - 17156 * Nu - 9172 * delta * Nu +
                       231 * params.eta2 + 775 * delta * params.eta2 +
                       1134 * (-1 + delta) * (-2 + kappa2) * Nu * params.S2z2) +
                  S2z * (5661 - 5661 * delta + 17156 * Nu - 9172 * delta * Nu -
                         231 * params.eta2 + 775 * delta * params.eta2 +
                         1134 * kappa2 * (1 - delta + (-3 + delta) * Nu) *
                             params.S2z2)))) /
            params.r4);
  }

  else if (vpnorder == 7) {

    return (
        /* Henry et al. QC spin terms */ /* Complex(0,0.0003968253968253968)*params.x4p5*(10*(6*(1
  + delta)*(63 + 26*kappa1) - (84*(26 + 17*delta) + (1069 +
  757*delta)*kappa1)*Nu + (2844 + 312*delta + (1409 +
  36*delta)*kappa1)*params.eta2)*params.S1z2 + S2z*(30*(14 - 431*Nu + delta*(-14 +
  87*Nu))*M_PI + 10*(6*(-1 + delta)*(63 + 26*kappa2) + (2184 - 1428*delta +
  1069*kappa2 - 757*delta*kappa2)*Nu +
     (-2844 + 312*delta - 1409*kappa2 + 36*delta*kappa2)*params.eta2)*S2z +
  Complex(0,3)*(66 - 66*delta - 2957*Nu + 4405*delta*Nu - 20*(14 - 14*delta -
  431*Nu + 87*delta*Nu)*log(2))) + 3*S1z*(Complex(0,-1)*(66 + 66*delta - 2957*Nu
  - 280*log(2)) + 5*
   (-28*(1 + delta)*M_PI + 2*(431 + 87*delta)*Nu*M_PI + delta*Nu*(Complex(0,881)
  + 16*(7 + 23*Nu)*S2z) - Complex(0,4)*(431*Nu + delta*(-14 + 87*Nu))*log(2))))
   + */
        /* Henry et al. ecc+spin terms */ (
            (Complex(0, -1.0020843354176688e-7) *
             (-23760 * params.Mtot3 *
                  (params.PhiDOT3 * params.r4 *
                       (Complex(0, -12) * (-35 + 107 * Nu) * S1z +
                        5 *
                            (126 - 462 * Nu + 80 * params.eta2 +
                             kappa1 * (-2 - 79 * Nu + 153 * params.eta2)) *
                            params.S1z2 +
                        S2z * (Complex(0, -420) + Complex(0, 1284) * Nu +
                               10 * (-63 + kappa2) * S2z +
                               5 * (462 + 79 * kappa2) * Nu * S2z -
                               5 * (80 + 153 * kappa2) * params.eta2 * S2z)) -
                   Complex(0, 1) * params.PhiDOT2 * params.r3 * rDOT *
                       (Complex(0, -24) * (-35 + 202 * Nu) * S1z +
                        5 *
                            (-14 * Nu * (3 + 58 * Nu) +
                             kappa1 * (-409 + 1096 * Nu + 6 * params.eta2)) *
                            params.S1z2 +
                        S2z * (Complex(0, -840) + 2045 * kappa2 * S2z +
                               10 * (406 - 3 * kappa2) * params.eta2 * S2z +
                               Nu * (Complex(0, 4848) + 210 * S2z -
                                     5480 * kappa2 * S2z))) -
                   Complex(0, 4) * r * params.rDOT3 *
                       (Complex(0, 3) * (26 + 57 * Nu) * S1z +
                        5 *
                            (4 * params.eta2 +
                             kappa1 * (-7 + 49 * Nu + 15 * params.eta2)) *
                            params.S1z2 -
                        S2z * (Complex(0, 78) - 35 * kappa2 * S2z +
                               5 * (4 + 15 * kappa2) * params.eta2 * S2z +
                               Nu * (Complex(0, 171) + 245 * kappa2 * S2z))) +
                   Complex(0, 2) * mass * rDOT *
                       (Complex(0, -4) * (-78 + 769 * Nu) * S1z +
                        5 *
                            ((14 - 59 * Nu) * Nu +
                             kappa1 * (-70 - 77 * Nu + 18 * params.eta2)) *
                            params.S1z2 +
                        S2z * (Complex(0, -312) + 350 * kappa2 * S2z +
                               5 * (59 - 18 * kappa2) * params.eta2 * S2z +
                               Nu * (Complex(0, 3076) - 70 * S2z +
                                     385 * kappa2 * S2z))) +
                   PhiDOT * params.r2 * params.rDOT2 *
                       (Complex(0, 6) * (-62 + 219 * Nu) * S1z -
                        5 *
                            (189 - 756 * Nu + 388 * params.eta2 +
                             kappa1 * (98 - 266 * Nu + 276 * params.eta2)) *
                            params.S1z2 +
                        S2z *
                            (Complex(0, 372) + 945 * S2z + 490 * kappa2 * S2z +
                             20 * (97 + 69 * kappa2) * params.eta2 * S2z -
                             2 * Nu *
                                 (Complex(0, 657) +
                                  35 * (54 + 19 * kappa2) * S2z))) +
                   mass * PhiDOT * r *
                       (Complex(0, 8) * (-61 + 480 * Nu) * S1z +
                        5 *
                            (-392 + 448 * Nu + 474 * params.eta2 +
                             kappa1 * (-11 + 150 * Nu + 58 * params.eta2)) *
                            params.S1z2 -
                        S2z *
                            (Complex(0, -488) - 5 * (392 + 11 * kappa2) * S2z +
                             10 * (237 + 29 * kappa2) * params.eta2 * S2z +
                             10 * Nu *
                                 (Complex(0, 384) +
                                  (224 + 75 * kappa2) * S2z)))) +
              delta * mass *
                  (-240 * mass * PhiDOT * params.r3 *
                       ((197936 - 139360 * Nu - 367105 * params.eta2 +
                         253245 * params.eta3) *
                            params.PhiDOT4 * params.r4 +
                        Complex(0, 1) *
                            (279236 - 483940 * Nu - 2817805 * params.eta2 +
                             459180 * params.eta3) *
                            params.PhiDOT3 * params.r3 * rDOT -
                        6 *
                            (38627 + 89295 * Nu - 492740 * params.eta2 +
                             75975 * params.eta3) *
                            params.PhiDOT2 * params.r2 * params.rDOT2 -
                        Complex(0, 1) *
                            (-731008 + 2287930 * Nu + 981060 * params.eta2 +
                             10275 * params.eta3) *
                            PhiDOT * r * params.rDOT3 +
                        (-327667 + 436705 * Nu + 659790 * params.eta2 -
                         438255 * params.eta3) *
                            params.rDOT4) +
                   900 * PhiDOT * params.r4 *
                       (2 *
                            (-2594 + 27609 * Nu - 74032 * params.eta2 +
                             25974 * params.eta3) *
                            params.PhiDOT6 * params.r6 +
                        Complex(0, 4) *
                            (-5730 + 58833 * Nu - 137842 * params.eta2 +
                             17123 * params.eta3) *
                            params.PhiDOT5 * params.r5 * rDOT +
                        2 *
                            (-114 - 41622 * Nu + 147569 * params.eta2 +
                             4196 * params.eta3) *
                            params.PhiDOT4 * params.r4 * params.rDOT2 +
                        Complex(0, 4) *
                            (-9554 + 70788 * Nu - 156227 * params.eta2 +
                             5810 * params.eta3) *
                            params.PhiDOT3 * params.r3 * params.rDOT3 +
                        (17619 - 138450 * Nu + 322600 * params.eta2 -
                         80816 * params.eta3) *
                            params.PhiDOT2 * params.r2 * params.rDOT4 -
                        Complex(0, 2) *
                            (8793 - 52230 * Nu + 69340 * params.eta2 +
                             2536 * params.eta3) *
                            PhiDOT * r * params.rDOT5 +
                        2 *
                            (3957 - 24534 * Nu + 42584 * params.eta2 -
                             20800 * params.eta3) *
                            params.rDOT6) -
                   2 * params.Mtot3 *
                       (Complex(0, -23760) * rDOT *
                            (5 *
                                 (Nu * (-14 + 31 * Nu) +
                                  7 * kappa1 * (10 + 31 * Nu)) *
                                 params.S1z2 +
                             2 * S1z *
                                 (Complex(0, -156) + 155 * params.eta2 * S2z +
                                  2 * Nu * (Complex(0, 613) + 390 * S2z)) +
                             S2z * (Complex(0, -312) + 350 * kappa2 * S2z +
                                    155 * params.eta2 * S2z +
                                    Nu * (Complex(0, 2452) +
                                          35 * (-2 + 31 * kappa2) * S2z))) +
                        PhiDOT * r *
                            (8946400 * params.eta3 -
                             8 * (6991786 + Complex(0, 724680) * S1z +
                                  7425 * (392 + 11 * kappa1) * params.S1z2 +
                                  Complex(0, 724680) * S2z +
                                  7425 * (392 + 11 * kappa2) * params.S2z2) -
                             3600 * params.eta2 *
                                 (-628 +
                                  33 * (-19 + 92 * kappa1) * params.S1z2 -
                                  7326 * S1z * S2z +
                                  33 * (-19 + 92 * kappa2) * params.S2z2) +
                             15 * Nu *
                                 (994455 * M_PI2 +
                                  8 * (-2249485 +
                                       7920 * (-21 + 8 * kappa1) * params.S1z2 +
                                       Complex(0, 283536) * S2z +
                                       7920 * (-21 + 8 * kappa2) * params.S2z2 -
                                       1584 * S1z *
                                           (Complex(0, -179) + 170 * S2z))))) +
                   3 * params.Mtot2 * r *
                       (Complex(0, 31680) * params.rDOT3 *
                            (5 * (4 * params.eta2 + 7 * kappa1 * (-1 + 5 * Nu)) *
                                 params.S1z2 +
                             S1z * (Complex(0, 78) + 40 * params.eta2 * S2z +
                                    Nu * (Complex(0, 327) + 420 * S2z)) +
                             S2z * (Complex(0, 78) - 35 * kappa2 * S2z +
                                    20 * params.eta2 * S2z +
                                    Nu * (Complex(0, 327) +
                                          175 * kappa2 * S2z))) -
                        22 * PhiDOT * r * params.rDOT2 *
                            (2553200 * params.eta3 -
                             24 * (268267 + Complex(0, 5580) * S1z +
                                   525 * (27 + 14 * kappa1) * params.S1z2 +
                                   Complex(0, 5580) * S2z +
                                   525 * (27 + 14 * kappa2) * params.S2z2) -
                             200 * params.eta2 *
                                 (39445 +
                                  72 * (-4 + 21 * kappa1) * params.S1z2 -
                                  3600 * S1z * S2z +
                                  72 * (-4 + 21 * kappa2) * params.S2z2) +
                             25 * Nu *
                                 (23247 * M_PI2 +
                                  8 * (-69259 + Complex(0, 1026) * S1z +
                                       126 * (27 + 5 * kappa1) * params.S1z2 +
                                       Complex(0, 1026) * S2z +
                                       126 * (27 + 5 * kappa2) *
                                           params.S2z2))) +
                        params.PhiDOT3 * params.r3 *
                            (10071200 * params.eta3 +
                             96 * (-421183 - Complex(0, 34650) * S1z +
                                   825 * (-63 + kappa1) * params.S1z2 -
                                   Complex(0, 34650) * S2z +
                                   825 * (-63 + kappa2) * params.S2z2) -
                             400 * params.eta2 *
                                 (64177 +
                                  792 * (-5 + 6 * kappa1) * params.S1z2 -
                                  17424 * S1z * S2z +
                                  792 * (-5 + 6 * kappa2) * params.S2z2) +
                             15 * Nu *
                                 (426195 * M_PI2 +
                                  8 * (-509635 +
                                       330 * (210 + 83 * kappa1) * params.S1z2 +
                                       Complex(0, 29304) * S2z +
                                       330 * (210 + 83 * kappa2) * params.S2z2 -
                                       792 * S1z *
                                           (Complex(0, -37) + 70 * S2z)))) -
                        Complex(0, 2) * params.PhiDOT2 * params.r2 * rDOT *
                            (-8330400 * params.eta3 +
                             8 * (-2810116 - Complex(0, 415800) * S1z +
                                  1012275 * kappa1 * params.S1z2 -
                                  Complex(0, 415800) * S2z +
                                  1012275 * kappa2 * params.S2z2) +
                             4800 * params.eta2 *
                                 (13411 +
                                  33 * (19 + 12 * kappa1) * params.S1z2 +
                                  462 * S1z * S2z +
                                  33 * (19 + 12 * kappa2) * params.S2z2) +
                             5 * Nu *
                                 (1278585 * M_PI2 -
                                  8 * (5139685 +
                                       990 * (-21 + 139 * kappa1) *
                                           params.S1z2 -
                                       Complex(0, 313632) * S2z +
                                       990 * (-21 + 139 * kappa2) *
                                           params.S2z2 -
                                       3564 * S1z *
                                           (Complex(0, 88) + 185 * S2z)))))) -
              13559040 * delta * params.Mtot3 * PhiDOT * r *
                  (2 * mass - 3 * params.PhiDOT2 * params.r3 +
                   Complex(0, 6) * PhiDOT * params.r2 * rDOT +
                   6 * r * params.rDOT2) *
                  log(r / r0))) /
            params.r4) + /* Henry et al. QC spinning hereditary terms */
        ((-14 * (1 + delta) * S1z + (431 + 87 * delta) * Nu * S1z +
          (14 - 431 * Nu + delta * (-14 + 87 * Nu)) * S2z) *
         params.x4p5 * (Complex(0, 1) * M_PI + log(4))) /
            84.);
  }

  else {
    return 0;
  }
}

static COMPLEX16 hQC_2_m_1(REAL8 Nu, UINT4 vpnorder, REAL8 x, struct kepler_vars params) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  /* if(vpnorder == 1){
       return(Complex(0,0.6666666666666666)*delta*params.x1p5);
   } */

  /* else if(vpnorder == 3){
       return((Complex(0,-0.40476190476190477)*delta +
    Complex(0,0.47619047619047616)*delta*Nu)*params.x2p5);
   } */

  /* else */ if (vpnorder == 4) {
    return ((2 * delta * params.x3 * (Complex(0, 1) * M_PI + log(4))) / 3.);
  }

  /* else if(vpnorder == 5){
       return((Complex(0,-0.2275132275132275)*delta -
    Complex(0,2.693121693121693)*delta*Nu +
    Complex(0,0.3134920634920635)*delta*params.eta2)*params.x3p5);
   } */

  else if (vpnorder == 6) {
    return (
        M_PI * (Complex(0, -0.40476190476190477) * delta * params.x4 +
                Complex(0, 0.14285714285714285) * delta * Nu * params.x4) +
        ((-17 * delta * params.x4) / 21. + (2 * delta * Nu * params.x4) / 7.) *
            log(2));
  }

  else {

    return 0;
  }
}

static COMPLEX16 hl_2_m_1(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z, REAL8 x, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_2_m_1: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_2_m_1(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x, params) +
            hQC_2_m_1(Nu, vpnorder, x, params)) *
           cpolar(1, -1 * Phi);
  }
}

static COMPLEX16 hl_2_m_min1(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z, REAL8 x, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_2_m_min1: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           conj((hGO_2_m_1(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x, params) +
                 hQC_2_m_1(Nu, vpnorder, x, params))) *
           cpolar(1, 1 * Phi);
  }
}

// H33
static COMPLEX16 hGO_3_m_3(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z,
                           REAL8 x, struct kepler_vars params) {
  REAL8 delta = sqrt(1 - 4 * Nu);
  REAL8 kappa1 = 1.0;
  REAL8 kappa2 = 1.0;
  REAL8 r0 = 1.0;

  if (vpnorder == 1) {
    return ((sqrt(0.11904761904761904) * delta *
             (2 * r * pow(Complex(0, 1) * PhiDOT * r - rDOT, 3) +
              mass * (Complex(0, -7) * PhiDOT * r + 4 * rDOT))) /
            (2. * r));
  }

  else if (vpnorder == 3) {
    return ((sqrt(0.11904761904761904) * delta *
             (6 * (-5 + 19 * Nu) * params.r2 *
                  pow(PhiDOT * r + Complex(0, 1) * rDOT, 4) *
                  (Complex(0, 1) * PhiDOT * r + rDOT) +
              2 * params.Mtot2 *
                  (Complex(0, -3) * (-101 + 43 * Nu) * PhiDOT * r +
                   (-109 + 86 * Nu) * rDOT) +
              3 * mass * r *
                  (Complex(0, -12) * (1 + 4 * Nu) * params.PhiDOT3 * params.r3 +
                   6 * (14 + 31 * Nu) * params.PhiDOT2 * params.r2 * rDOT +
                   Complex(0, 3) * (33 + 62 * Nu) * PhiDOT * r * params.rDOT2 -
                   4 * (8 + 17 * Nu) * params.rDOT3))) /
            (36. * params.r2));
  }

  else if (vpnorder == 4) {
    return (
        (Complex(0, -0.125) * sqrt(0.11904761904761904) * params.Mtot2 *
         (4 * mass * (-1 + 5 * Nu) * ((1 + delta) * S1z + (-1 + delta) * S2z) +
          r * (2 * params.rDOT2 *
                   (6 * (1 + delta) * S1z - 5 * (5 + 3 * delta) * Nu * S1z +
                    (-6 + delta * (6 - 15 * Nu) + 25 * Nu) * S2z) +
               params.PhiDOT2 * params.r2 *
                   (-24 * (1 + delta) * S1z + (119 + 33 * delta) * Nu * S1z +
                    (24 - 119 * Nu + 3 * delta * (-8 + 11 * Nu)) * S2z) +
               Complex(0, 2) * PhiDOT * r * rDOT *
                   (-18 * (1 + delta) * S1z + (77 + 39 * delta) * Nu * S1z +
                    (18 - 77 * Nu + 3 * delta * (-6 + 13 * Nu)) * S2z)))) /
            params.r3 -
        (Complex(0, -0.375) * sqrt(1.0714285714285714) *
         (-4 * S1z + 19 * Nu * (S1z - S2z) + 4 * S2z - 4 * delta * (S1z + S2z) +
          5 * delta * Nu * (S1z + S2z)) *
         params.x3) +
        (Complex(0, -0.375) * sqrt(0.04285714285714286) * params.x3 *
         (5 * (-4 + 19 * Nu + delta * (-4 + 5 * Nu)) * S1z +
          5 * (4 - 19 * Nu + delta * (-4 + 5 * Nu)) * S2z)));
  }

  else if (vpnorder == 5) {
    return (
        (delta *
         (30 * (183 - 1579 * Nu + 3387 * params.eta2) * params.r3 *
              pow(PhiDOT * r - Complex(0, 1) * rDOT, 2) *
              pow(Complex(0, -1) * PhiDOT * r + rDOT, 5) +
          10 * params.Mtot3 *
              (Complex(0, -1) * (26473 - 27451 * Nu + 9921 * params.eta2) *
                   PhiDOT * r +
               4 * (623 - 732 * Nu + 1913 * params.eta2) * rDOT) +
          2 * params.Mtot2 * r *
              (Complex(0, -11) * (-5353 - 13493 * Nu + 4671 * params.eta2) *
                   params.PhiDOT3 * params.r3 +
               (-75243 - 142713 * Nu + 192821 * params.eta2) * params.PhiDOT2 *
                   params.r2 * rDOT +
               Complex(0, 220) * (-256 + 781 * Nu + 840 * params.eta2) * PhiDOT *
                   r * params.rDOT2 -
               10 * (-756 + 8238 * Nu + 7357 * params.eta2) * params.rDOT3) +
          3 * mass * params.r2 *
              (Complex(0, 2) * (-7633 + 9137 * Nu + 28911 * params.eta2) *
                   params.PhiDOT5 * params.r5 -
               4 * (-8149 + 1576 * Nu + 43533 * params.eta2) * params.PhiDOT4 *
                   params.r4 * rDOT -
               Complex(0, 2) * (-9297 - 19517 * Nu + 64839 * params.eta2) *
                   params.PhiDOT3 * params.r3 * params.rDOT2 -
               32 * (-1288 + 3667 * Nu + 4056 * params.eta2) * params.PhiDOT2 *
                   params.r2 * params.rDOT3 -
               Complex(0, 5) * (-9851 + 17954 * Nu + 40968 * params.eta2) *
                   PhiDOT * r * params.rDOT4 +
               20 * (-771 + 1126 * Nu + 3616 * params.eta2) * params.rDOT5))) /
            (1584. * sqrt(210) * params.r3)
        /* Henry et al. QC spin terms */
        /* +(Complex(0,1.125)*sqrt(1.0714285714285714)*(kappa1*(-1 - delta +
  2*(2 + delta)*Nu)*params.S1z2
  + S2z*(-4*delta*Nu*S1z + kappa2*(1 - delta + 2*(-2 +
  delta)*Nu)*S2z))*params.x3p5) */
        + /* Henry et al. ecc spin terms */ (
              (Complex(0, 0.125) * sqrt(1.0714285714285714) * params.Mtot3 *
               (7 * PhiDOT * r + Complex(0, 2) * rDOT) *
               (kappa1 * (-1 - delta + 2 * (2 + delta) * Nu) * params.S1z2 +
                S2z * (-4 * delta * Nu * S1z +
                       kappa2 * (1 - delta + 2 * (-2 + delta) * Nu) * S2z))) /
              params.r3));
  }

  else if (vpnorder == 6) {
    return (
        -(delta * params.Mtot2 * Nu *
          (668 * params.Mtot2 +
           2 * mass * r *
               (4081 * params.PhiDOT2 * params.r2 +
                Complex(0, 297) * PhiDOT * r * rDOT - 452 * params.rDOT2) +
           5 * params.r2 *
               (1329 * params.PhiDOT4 * params.r4 -
                Complex(0, 2926) * params.PhiDOT3 * params.r3 * rDOT -
                384 * params.PhiDOT2 * params.r2 * params.rDOT2 -
                Complex(0, 408) * PhiDOT * r * params.rDOT3 +
                200 * params.rDOT4))) /
            (36. * sqrt(210) * params.r4)
        /* Henry et al. QC spin terms */
        /* +(Complex(0,-0.125)*sqrt(0.04285714285714286)*(-(Nu*((279 +
  delta)*S1z + (-279 + delta)*S2z))
  + 10*(S1z - S2z + delta*(S1z + S2z)) + params.eta2*(407*(S1z - S2z) +
  241*delta*(S1z + S2z)))*params.x4) */
        + /* Henry et al. ecc spin terms*/ (
              (Complex(0, -0.006944444444444444) * params.Mtot2 *
               (10 * params.Mtot2 *
                    ((252 * (1 + delta) - (1277 + 1279 * delta) * Nu +
                      8 * (12 + 47 * delta) * params.eta2) *
                         S1z +
                     (252 * (-1 + delta) + (1277 - 1279 * delta) * Nu +
                      8 * (-12 + 47 * delta) * params.eta2) *
                         S2z) +
                2 * mass * r *
                    (2 * params.PhiDOT2 * params.r2 *
                         ((1320 * (1 + delta) - 2 * (4469 + 211 * delta) * Nu +
                           (8709 + 2777 * delta) * params.eta2) *
                              S1z +
                          (1320 * (-1 + delta) + 8938 * Nu - 422 * delta * Nu +
                           (-8709 + 2777 * delta) * params.eta2) *
                              S2z) +
                     Complex(0, 3) * PhiDOT * r * rDOT *
                         ((2000 * (1 + delta) - (9147 + 3173 * delta) * Nu +
                           (8911 + 5273 * delta) * params.eta2) *
                              S1z +
                          (2000 * (-1 + delta) + (9147 - 3173 * delta) * Nu +
                           (-8911 + 5273 * delta) * params.eta2) *
                              S2z) +
                     10 * params.rDOT2 *
                         ((-105 * (1 + delta) + (541 + 77 * delta) * Nu -
                           2 * (462 + 247 * delta) * params.eta2) *
                              S1z +
                          (105 + Nu * (-541 + 924 * Nu) +
                           delta * (-105 + (77 - 494 * Nu) * Nu)) *
                              S2z)) -
                3 * params.r2 *
                    (-3 * params.PhiDOT2 * params.r2 * params.rDOT2 *
                         ((480 * (1 + delta) - (1711 + 1889 * delta) * Nu +
                           2 * (-1161 + 757 * delta) * params.eta2) *
                              S1z +
                          (480 * (-1 + delta) + (1711 - 1889 * delta) * Nu +
                           2 * (1161 + 757 * delta) * params.eta2) *
                              S2z) +
                     2 * params.PhiDOT4 * params.r4 *
                         ((350 * (1 + delta) - 4 * (404 + 461 * delta) * Nu +
                           (883 + 769 * delta) * params.eta2) *
                              S1z +
                          (350 * (-1 + delta) + 4 * (404 - 461 * delta) * Nu +
                           (-883 + 769 * delta) * params.eta2) *
                              S2z) +
                     Complex(0, 2) * params.PhiDOT3 * params.r3 * rDOT *
                         ((660 * (1 + delta) - (4061 + 2899 * delta) * Nu +
                           (2643 + 4789 * delta) * params.eta2) *
                              S1z +
                          (660 * (-1 + delta) + (4061 - 2899 * delta) * Nu +
                           (-2643 + 4789 * delta) * params.eta2) *
                              S2z) +
                     10 * params.rDOT4 *
                         ((-30 * (1 + delta) + (187 + 101 * delta) * Nu -
                           2 * (159 + 61 * delta) * params.eta2) *
                              S1z +
                          (30 + Nu * (-187 + 318 * Nu) +
                           delta * (-30 + (101 - 122 * Nu) * Nu)) *
                              S2z) +
                     Complex(0, 2) * PhiDOT * r * params.rDOT3 *
                         ((90 + Nu * (-1321 + 5118 * Nu) +
                           delta * (90 + Nu * (-319 + 714 * Nu))) *
                              S1z +
                          (-90 + (1321 - 5118 * Nu) * Nu +
                           delta * (90 + Nu * (-319 + 714 * Nu))) *
                              S2z)))) /
              (sqrt(210) * params.r4)));
  }

  else if (vpnorder == 7) {

    return (
        /* Henry et al. QC spin terms */ /* (Complex(0,-0.020833333333333332)*params.x4p5*(-270*(-6*Nu*(-2
 + delta*(-2 + Nu) + 6*Nu)
  + kappa1*(4 - Nu*(13 + 8*Nu) + delta*(4 + Nu*(-5 + 12*Nu))))*params.S1z2 +
 S1z*(810*(-4 + 19*Nu + delta*(-4 + 5*Nu))*M_PI + 5*delta*Nu*(Complex(0,-541)
 + 216*(-10 + 9*Nu)*S2z) + Complex(0,61560)*Nu*atanh(1/5) +
 Complex(0,1)*(2349*delta
 - 14207*Nu + 81*(29 - 80*log(1.5))) + Complex(0,1620)*delta*(-4 +
 5*Nu)*log(1.5)) + S2z*(810*(4 - 19*Nu + delta*(-4 + 5*Nu))*M_PI +
 Complex(0,1)*(2349*delta + 14207*Nu
 + 81*(-29 + 80*log(1.5))) - 5*(54*kappa2*(-4 + Nu*(13 + 8*Nu) + delta*(4 +
 Nu*(-5 + 12*Nu)))*S2z + Nu*(Complex(0,541)*delta - 324*(2 + delta*(-2 + Nu) -
 6*Nu)*S2z - Complex(0,648)*(-19 + 5*delta)*atanh(1/5))
  - Complex(0,324)*delta*log(1.5802469135802468)) -
 Complex(0,486)*delta*log(1024))))/sqrt(210)
  + */
                                         /* Henry et al. ecc spin terms */
        ((Complex(0, 504504) * params.Mtot3 *
              (2 * mass *
                   (rDOT *
                        (S1z * (108 - 498 * Nu +
                                Complex(0, 5) *
                                    (-24 - 5 * kappa1 + 3 * (76 + kappa1) * Nu +
                                     4 * (-111 + 13 * kappa1) * params.eta2) *
                                    S1z) +
                         6 * (-18 + 83 * Nu) * S2z -
                         Complex(0, 5) *
                             (-24 - 5 * kappa2 + 3 * (76 + kappa2) * Nu +
                              4 * (-111 + 13 * kappa2) * params.eta2) *
                             params.S2z2) +
                    PhiDOT * r *
                        (S1z * (Complex(0, -3) * (-99 + 581 * Nu) +
                                5 *
                                    (-24 + 399 * kappa1 +
                                     48 * (7 - 19 * Nu) * Nu +
                                     kappa1 * Nu * (-1629 + 188 * Nu)) *
                                    S1z) +
                         Complex(0, 3) * (-99 + 581 * Nu) * S2z -
                         5 *
                             (-24 + 399 * kappa2 + 48 * (7 - 19 * Nu) * Nu +
                              kappa2 * Nu * (-1629 + 188 * Nu)) *
                             params.S2z2)) +
               r * (3 * PhiDOT * r * params.rDOT2 *
                        (S1z * (Complex(0, 216) + 545 * kappa1 * S1z +
                                40 * (45 + 8 * kappa1) * params.eta2 * S1z -
                                30 * Nu *
                                    (Complex(0, 50) + 20 * S1z +
                                     73 * kappa1 * S1z)) +
                         Complex(0, 12) * (-18 + 125 * Nu) * S2z -
                         5 *
                             (109 * kappa2 - 6 * (20 + 73 * kappa2) * Nu +
                              8 * (45 + 8 * kappa2) * params.eta2) *
                             params.S2z2) +
                    2 * params.rDOT3 *
                        (S1z * (-54 + Complex(0, 145) * kappa1 * S1z -
                                Complex(0, 30) * Nu *
                                    (Complex(0, 11) + 2 * Nu * S1z +
                                     kappa1 * (17 + 8 * Nu) * S1z)) +
                         6 * (9 - 55 * Nu) * S2z +
                         Complex(0, 5) *
                             (12 * params.eta2 +
                              kappa2 * (-29 + 6 * Nu * (17 + 8 * Nu))) *
                             params.S2z2) +
                    6 * params.PhiDOT2 * params.r2 * rDOT *
                        (S1z * (297 - 465 * Nu +
                                Complex(0, 5) *
                                    (6 * (20 - 87 * Nu) * Nu +
                                     kappa1 * (-50 + 3 * Nu * (47 + 76 * Nu))) *
                                    S1z) +
                         3 * (-99 + 155 * Nu) * S2z -
                         Complex(0, 5) *
                             (6 * (20 - 87 * Nu) * Nu +
                              kappa2 * (-50 + 3 * Nu * (47 + 76 * Nu))) *
                             params.S2z2) +
                    params.PhiDOT3 * params.r3 *
                        (Complex(0, 3) * (531 - 1295 * Nu) * S1z +
                         10 *
                             (-33 * kappa1 + 6 * (30 + 13 * kappa1) * Nu +
                              4 * (-96 + 67 * kappa1) * params.eta2) *
                             params.S1z2 +
                         S2z * (Complex(0, -1593) + 330 * kappa2 * S2z +
                                5 * Nu *
                                    (Complex(0, 777) -
                                     4 *
                                         (90 + 39 * kappa2 - 192 * Nu +
                                          134 * kappa2 * Nu) *
                                         S2z))))) +
          delta *
              (-17640 * (-4083 + Nu * (58311 + Nu * (-269240 + 405617 * Nu))) *
                   params.r4 * pow(PhiDOT * r + Complex(0, 1) * rDOT, 6) *
                   pow(Complex(0, 1) * PhiDOT * r + rDOT, 3) +
               168 * params.Mtot2 * params.r2 *
                   (Complex(0, 1) *
                        (-7508635 +
                         7 * Nu *
                             (-1318438 + Nu * (-10231834 + 9667755 * Nu))) *
                        params.PhiDOT5 * params.r5 +
                    7 *
                        (1235591 +
                         Nu * (884445 + (23935218 - 26913443 * Nu) * Nu)) *
                        params.PhiDOT4 * params.r4 * rDOT -
                    Complex(0, 1) *
                        (8961149 +
                         7 * Nu *
                             (-31755709 + Nu * (-11134798 + 22187331 * Nu))) *
                        params.PhiDOT3 * params.r3 * params.rDOT2 -
                    (-36806435 +
                     7 * Nu * (33178545 + Nu * (24565078 + 22873537 * Nu))) *
                        params.PhiDOT2 * params.r2 * params.rDOT3 -
                    Complex(0, 5) *
                        (-7761899 +
                         7 * Nu *
                             (2892563 + 5998602 * Nu + 7493619 * params.eta2)) *
                        PhiDOT * r * params.rDOT4 +
                    5 * (-2422057 + 7 * Nu * (501045 + Nu * (2033141 + 2771816 * Nu))) *
                        params.rDOT5) +
               1764 * mass * params.r3 *
                   (Complex(0, -2) *
                        (239087 +
                         Nu * (-1206515 + Nu * (422631 + 3979375 * Nu))) *
                        params.PhiDOT7 * params.r7 +
                    2 * (621284 + Nu * (-2279907 + 2 * Nu * (-1180187 + 5876531 * Nu))) *
                        params.PhiDOT6 * params.r6 * rDOT +
                    Complex(0, 2) *
                        (39270 +
                         Nu * (1235486 - 5319747 * Nu + 4406349 * params.eta2)) *
                        params.PhiDOT5 * params.r5 * params.rDOT2 +
                    8 * (349111 + 4 * Nu * (-519370 + 33 * Nu * (10939 + 42635 * Nu))) *
                        params.PhiDOT4 * params.r4 * params.rDOT3 +
                    Complex(0, 2) *
                        (1212607 +
                         3 * Nu *
                             (-2012698 - 67827 * Nu + 7955628 * params.eta2)) *
                        params.PhiDOT3 * params.r3 * params.rDOT4 +
                    4 * (201135 + 2 * Nu * (-773107 + Nu * (1214819 + 1157652 * Nu))) *
                        params.PhiDOT2 * params.r2 * params.rDOT5 +
                    Complex(0, 5) *
                        (333969 +
                         2 * Nu * (-981471 + 4 * Nu * (154039 + 750016 * Nu))) *
                        PhiDOT * r * params.rDOT6 -
                    40 *
                        (13245 +
                         2 * Nu * (-37005 + Nu * (14251 + 130160 * Nu))) *
                        params.rDOT7) +
               2 * params.Mtot4 *
                   (4 * rDOT *
                        (269279500 * params.eta3 +
                         2 * (-174108226 +
                              63063 * S1z *
                                  (Complex(0, 108) +
                                   5 * (24 + 5 * kappa1) * S1z) +
                              63063 * S2z *
                                  (Complex(0, 108) +
                                   5 * (24 + 5 * kappa2) * S2z)) -
                         21 * Nu *
                             (103100846 + 1846845 * M_PI2 +
                              Complex(0, 1693692) * S2z -
                              6006 * (S1z * (Complex(0, -282) +
                                             5 * (-180 + 7 * kappa1) * S1z) -
                                      980 * S1z * S2z +
                                      5 * (-180 + 7 * kappa2) * params.S2z2)) -
                         2940 * params.eta2 *
                             (-122855 +
                              4719 * ((-6 + 7 * kappa1) * params.S1z2 -
                                      26 * S1z * S2z +
                                      (-6 + 7 * kappa2) * params.S2z2))) +
                    Complex(0, 1) * PhiDOT * r *
                        (-1176172480 * params.eta3 +
                         8 * (-74084729 +
                              189189 * S1z *
                                  (Complex(0, 99) +
                                   5 * (-8 + 133 * kappa1) * S1z) +
                              189189 * S2z *
                                  (Complex(0, 99) +
                                   5 * (-8 + 133 * kappa2) * S2z)) +
                         35280 * params.eta2 *
                             (56255 +
                              429 * ((-22 + 65 * kappa1) * params.S1z2 -
                                     174 * S1z * S2z +
                                     (-22 + 65 * kappa2) * params.S2z2)) -
                         147 * Nu *
                             (-65012788 + 4485195 * M_PI2 +
                              Complex(0, 3943368) * S2z +
                              10296 *
                                  (S1z * (Complex(0, 383) +
                                          5 * (-96 + 277 * kappa1) * S1z) -
                                   3220 * S1z * S2z +
                                   5 * (-96 + 277 * kappa2) * params.S2z2)))) +
               params.Mtot3 * r *
                   (Complex(0, -12) * PhiDOT * r * params.rDOT2 *
                        (-1035895280 * params.eta3 -
                         2 * (-547993687 +
                              63063 * S1z *
                                  (Complex(0, 216) + 545 * kappa1 * S1z) +
                              63063 * S2z *
                                  (Complex(0, 216) + 545 * kappa2 * S2z)) +
                         77 * Nu *
                             (42451610 + 1511055 * M_PI2 +
                              Complex(0, 1749384) * S2z +
                              6552 * (S1z * (Complex(0, 267) +
                                             25 * (6 + 11 * kappa1) * S1z) -
                                      5 * S1z * S2z +
                                      25 * (6 + 11 * kappa2) * params.S2z2)) +
                         490 * params.eta2 *
                             (-5802767 +
                              5148 * ((-6 + 23 * kappa1) * params.S1z2 -
                                      58 * S1z * S2z +
                                      (-6 + 23 * kappa2) * params.S2z2))) +
                    4 * params.rDOT3 *
                        (-1359334480 * params.eta3 -
                         4 * (-150254558 +
                              63063 * S1z *
                                  (Complex(0, 54) + 145 * kappa1 * S1z) +
                              63063 * S2z *
                                  (Complex(0, 54) + 145 * kappa2 * S2z)) +
                         231 * Nu *
                             (8490448 + 503685 * M_PI2 +
                              Complex(0, 242424) * S2z +
                              2184 * (S1z * (Complex(0, 111) +
                                             110 * kappa1 * S1z) +
                                      70 * S1z * S2z +
                                      110 * kappa2 * params.S2z2)) +
                         11760 * params.eta2 *
                             (-312980 +
                              429 * ((3 + 25 * kappa1) * params.S1z2 -
                                     44 * S1z * S2z +
                                     (3 + 25 * kappa2) * params.S2z2))) +
                    6 * params.PhiDOT2 * params.r2 * rDOT *
                        (2368900688 * params.eta3 +
                         8 * (-812986529 +
                              63063 * S1z *
                                  (Complex(0, 297) + 250 * kappa1 * S1z) +
                              63063 * S2z *
                                  (Complex(0, 297) + 250 * kappa2 * S2z)) -
                         1176 * params.eta2 *
                             (2423171 +
                              4290 * ((-3 + 41 * kappa1) * params.S1z2 -
                                      88 * S1z * S2z +
                                      (-3 + 41 * kappa2) * params.S2z2)) +
                         539 * Nu *
                             (-24139772 + 647595 * M_PI2 +
                              Complex(0, 120744) * S2z -
                              936 * (S1z * (Complex(0, -129) + 600 * S1z +
                                            205 * kappa1 * S1z) -
                                     460 * S1z * S2z +
                                     5 * (120 + 41 * kappa2) * params.S2z2))) +
                    Complex(0, 1) * params.PhiDOT3 * params.r3 *
                        (-4538040136 * params.eta3 -
                         88 * (259018351 +
                               17199 * S1z *
                                   (Complex(0, -531) + 110 * kappa1 * S1z) +
                               17199 * S2z *
                                   (Complex(0, -531) + 110 * kappa2 * S2z)) +
                         2352 * params.eta2 *
                             (7332973 +
                              12870 * ((5 + 23 * kappa1) * params.S1z2 -
                                       36 * S1z * S2z +
                                       (5 + 23 * kappa2) * params.S2z2)) +
                         21 * Nu *
                             (49864815 * M_PI2 +
                              8 * (-88128538 - Complex(0, 2099097) * S2z +
                                   9009 *
                                       (S1z * (Complex(0, -233) +
                                               40 * (15 + kappa1) * S1z) +
                                        360 * S1z * S2z +
                                        40 * (15 + kappa2) * params.S2z2)))))) +
          74954880 * delta * params.Mtot3 *
              (Complex(0, 22) * mass * PhiDOT * r +
               Complex(0, 59) * params.PhiDOT3 * params.r4 + 8 * mass * rDOT +
               66 * params.PhiDOT2 * params.r3 * rDOT +
               Complex(0, 24) * PhiDOT * params.r2 * params.rDOT2 -
               4 * r * params.rDOT3) *
              log(r / r0)) /
         (2.4216192e7 * sqrt(210) *
          params.r4)) + /* Henry et al. QC spinning hereditary terms */
        ((9 * sqrt(0.04285714285714286) * params.x4p5 *
          (Complex(0, -5) * M_PI *
               ((-4 - 4 * delta + 19 * Nu + 5 * delta * Nu) * S1z +
                (4 - 4 * delta - 19 * Nu + 5 * delta * Nu) * S2z) +
           20 * Nu * (19 * S1z - 19 * S2z + 5 * delta * S2z) * atanh(1 / 5) +
           delta * S2z * (10 * log(1.5802469135802468) - 3 * log(1024)))) /
         8.));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hQC_3_m_3(REAL8 Nu, UINT4 vpnorder, REAL8 x, struct kepler_vars params) {
  REAL8 delta = sqrt(1 - 4 * Nu);
  double EulerGamma = 0.5772156649015329;

  /* if(vpnorder == 1){
       return(Complex(0,-1.5)*sqrt(1.0714285714285714)*delta*params.x1p5);
   } */

  /* else if(vpnorder == 3){
       return((Complex(0,3)*sqrt(4.285714285714286)*delta -
    Complex(0,3)*sqrt(1.0714285714285714)*delta*Nu)*params.x2p5);
   } */

  /* else */ if (vpnorder == 4) {
    return ((9 * sqrt(1.0714285714285714) * delta * params.x3 *
             (Complex(0, -1) * M_PI + log(2.25))) /
            2.);
  }

  /* else if(vpnorder == 5){
       return((Complex(0,-8.386363636363637)*sqrt(0.04285714285714286)*delta +
    Complex(0,83.54545454545455)*sqrt(0.04285714285714286)*delta*Nu -
    Complex(0,20.15909090909091)*sqrt(0.04285714285714286)*delta*
     params.eta2)*params.x3p5);
   } */

  else if (vpnorder == 6) {
    return (M_PI *
                (Complex(0, 9) * sqrt(4.285714285714286) * delta * params.x4 -
                 Complex(0, 6.75) * sqrt(1.0714285714285714) * delta * Nu *
                     params.x4) +
            (-18 * sqrt(4.285714285714286) * delta * params.x4 +
             (27 * sqrt(1.0714285714285714) * delta * Nu * params.x4) / 2.) *
                log(1.5));
  } else if (vpnorder == 7) {
    return (Complex(0, 1.1149564720993292e-6) * sqrt(0.04285714285714286) *
            delta * params.x4p5 *
            (-465315528 + 74954880 * EulerGamma + 13827800 * Nu +
             124985672 * params.eta2 - 19373424 * params.eta3 +
             Complex(0, 47279232) * M_PI - 10090080 * M_PI2 -
             4309305 * Nu * M_PI2 - 94558464 * log(1.5) -
             Complex(0, 121080960) * M_PI * log(1.5) + 37477440 * log(16 * x) +
             121080960 * log(1.5) * log(1.5)));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hl_3_m_3(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z, REAL8 x, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_3_m_3: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_3_m_3(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x, params) +
            hQC_3_m_3(Nu, vpnorder, x, params)) *
           cpolar(1, -3 * Phi);
  }
}

static COMPLEX16 hl_3_m_min3(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z, REAL8 x, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_3_m_min3: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((- 4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           conj(hGO_3_m_3(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x, params) +
                hQC_3_m_3(Nu, vpnorder, x, params)) *
           cpolar(1, 3 * Phi);
  }
}

// H32
static COMPLEX16 hGO_3_m_2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z,
                           REAL8 x, struct kepler_vars params) {
  REAL8 delta = sqrt(1 - 4 * Nu);
  REAL8 kappa1 = 1.0;
  REAL8 kappa2 = 1.0;

  if (vpnorder == 2) {
    return (-(sqrt(0.7142857142857143) * mass * (-1 + 3 * Nu) * PhiDOT *
              (4 * PhiDOT * r + Complex(0, 1) * rDOT)) /
            6.);
  }

  else if (vpnorder == 3) {
    return (
        (sqrt(0.7142857142857143) * params.Mtot2 * Nu *
         (4 * PhiDOT * r + Complex(0, 1) * rDOT) * (S1z + S2z)) /
            (3. * params.r2) -
        ((4 * sqrt(0.7142857142857143) * Nu * (S1z + S2z) * params.x2p5) / 3.) +
        ((4 * sqrt(0.7142857142857143) * Nu * (S1z + S2z) * params.x2p5) / 3.));
  }

  else if (vpnorder == 4) {
    return (-(mass * PhiDOT *
              (2 * mass *
                   ((167 - 925 * Nu + 1615 * params.eta2) * PhiDOT * r +
                    Complex(0, 5) * (-82 + 239 * Nu + 55 * params.eta2) * rDOT) -
               3 * r *
                   (2 * (-13 - 25 * Nu + 355 * params.eta2) * params.PhiDOT3 *
                        params.r3 -
                    Complex(0, 60) * (-8 + 25 * Nu + params.eta2) *
                        params.PhiDOT2 * params.r2 * rDOT +
                    12 * (-23 + 70 * Nu + 10 * params.eta2) * PhiDOT * r *
                        params.rDOT2 +
                    Complex(0, 5) * (-13 + 38 * Nu + 10 * params.eta2) *
                        params.rDOT3))) /
            (108. * sqrt(35) * r));
  }

  else if (vpnorder == 5) {
    return (
        (params.Mtot2 * Nu * PhiDOT *
         (Complex(0, 7) * mass +
          r * (Complex(0, 49) * params.PhiDOT2 * params.r2 +
               90 * PhiDOT * r * rDOT - Complex(0, 6) * params.rDOT2))) /
            (4. * sqrt(35) * params.r2)
        /* Henry et al. QC spin terms */
        /* +((sqrt(0.7142857142857143)*(43*delta*Nu*(S1z - S2z) - 3*Nu*(S1z +
  S2z)
  - 26*params.eta2*(S1z + S2z) - 8*(S1z + delta*S1z + S2z -
  delta*S2z))*params.x3p5)/9.) */
        +
        /* Henry et al. ecc spin terms */ (
            (sqrt(0.7142857142857143) * params.Mtot2 *
             (Complex(0, 2) * mass * rDOT *
                  ((-12 + Nu * (97 + 4 * Nu) + delta * (-12 + 5 * Nu)) * S1z +
                   (-12 + delta * (12 - 5 * Nu) + Nu * (97 + 4 * Nu)) * S2z) +
              4 * mass * PhiDOT * r *
                  (-((12 + delta * (12 - 23 * Nu) + Nu * (53 + 8 * Nu)) * S1z) -
                   (12 + Nu * (53 + 8 * Nu) + delta * (-12 + 23 * Nu)) * S2z) -
              3 * r *
                  (16 * params.eta2 * PhiDOT * r *
                       (4 * params.PhiDOT2 * params.r2 -
                        Complex(0, 2) * PhiDOT * r * rDOT + params.rDOT2) *
                       (S1z + S2z) +
                   Complex(0, 30) * params.PhiDOT2 * params.r2 * rDOT *
                       (S1z + delta * S1z + S2z - delta * S2z) +
                   Nu * (-4 * params.PhiDOT3 * params.r3 *
                             ((5 + 17 * delta) * S1z + (5 - 17 * delta) * S2z) -
                         Complex(0, 1) * params.PhiDOT2 * params.r2 * rDOT *
                             ((189 + 17 * delta) * S1z +
                              (189 - 17 * delta) * S2z) +
                         20 * PhiDOT * r * params.rDOT2 *
                             (-((-3 + delta) * S1z) + (3 + delta) * S2z) -
                         Complex(0, 4) * params.rDOT3 *
                             ((-4 + delta) * S1z - (4 + delta) * S2z))))) /
            (72. * params.r3)));
  }

  else if (vpnorder == 6) {
    return (
        -(mass * PhiDOT *
          (4 * params.Mtot2 *
               (2 *
                    (5377 + 6438 * Nu - 79866 * params.eta2 +
                     37348 * params.eta3) *
                    PhiDOT * r -
                Complex(0, 5) *
                    (-4115 + 18399 * Nu - 20276 * params.eta2 + 7 * params.eta3) *
                    rDOT) -
           4 * mass * r *
               ((4599 - 15737 * Nu + 36259 * params.eta2 + 108563 * params.eta3) *
                    params.PhiDOT3 * params.r3 -
                Complex(0, 1) *
                    (-34053 + 59698 * Nu + 192949 * params.eta2 +
                     16193 * params.eta3) *
                    params.PhiDOT2 * params.r2 * rDOT +
                (-59058 + 77983 * Nu + 322468 * params.eta2 -
                 4264 * params.eta3) *
                    PhiDOT * r * params.rDOT2 +
                Complex(0, 5) *
                    (-3387 + 8518 * Nu + 8968 * params.eta2 + 884 * params.eta3) *
                    params.rDOT3) +
           3 * params.r2 *
               (4 *
                    (-710 + 3892 * Nu - 10655 * params.eta2 +
                     24000 * params.eta3) *
                    params.PhiDOT5 * params.r5 +
                Complex(0, 11) *
                    (-1484 + 11693 * Nu - 25006 * params.eta2 +
                     428 * params.eta3) *
                    params.PhiDOT4 * params.r4 * rDOT +
                4 *
                    (4161 - 25618 * Nu + 29489 * params.eta2 +
                     22078 * params.eta3) *
                    params.PhiDOT3 * params.r3 * params.rDOT2 +
                Complex(0, 44) *
                    (-151 + 1067 * Nu - 2419 * params.eta2 + 57 * params.eta3) *
                    params.PhiDOT2 * params.r2 * params.rDOT3 +
                4 *
                    (2041 - 11680 * Nu + 19334 * params.eta2 +
                     3368 * params.eta3) *
                    PhiDOT * r * params.rDOT4 +
                Complex(0, 5) *
                    (477 - 2624 * Nu + 3862 * params.eta2 + 1160 * params.eta3) *
                    params.rDOT5))) /
            (4752. * sqrt(35) * params.r2)
        /* Henry et al. QC spin terms */
        /* +((2*sqrt(0.7142857142857143)*((2*Nu*(-5 - 5*delta + 4*Nu) +
  3*kappa1*(1 + delta
  - 2*(2 + delta)*Nu + 6*params.eta2))*params.S1z2 + S2z*(2*(4 +
  9*kappa2)*params.eta2*S2z
  - 3*(-1 + delta)*(Complex(0,2) + kappa2*S2z) + 2*Nu*(Complex(0,-15) + 6*M_PI
  + 5*(-1 + delta)*S2z + 3*(-2 + delta)*kappa2*S2z)) + 2*S1z*(Complex(0,3)
  + Complex(0,3)*delta + Nu*(Complex(0,-15) + 6*M_PI + 2*S2z -
  10*Nu*S2z)))*params.x4)/9.) */
        +
        /* Henry et al. ecc spin terms */ (
            (sqrt(0.7142857142857143) * params.Mtot3 *
             (2 * mass *
                  (Complex(0, 2) * (1 + delta - 2 * Nu) * S1z +
                   Nu * ((1 + delta) * (-6 + kappa1) + 12 * Nu) * params.S1z2 +
                   8 * Nu * (-1 + 3 * Nu) * S1z * S2z +
                   S2z * (Complex(0, -2) * (-1 + delta + 2 * Nu) +
                          Nu * (-6 - delta * (-6 + kappa2) + kappa2 + 12 * Nu) *
                              S2z)) +
              r * (4 * params.rDOT2 *
                       (Complex(0, 2) * (1 + delta - 2 * Nu) * S1z +
                        Nu * ((1 + delta) * (-2 + kappa1) + 4 * Nu) *
                            params.S1z2 +
                        8 * params.eta2 * S1z * S2z +
                        S2z * (Complex(0, -2) * (-1 + delta + 2 * Nu) +
                               Nu *
                                   (-2 - delta * (-2 + kappa2) + kappa2 +
                                    4 * Nu) *
                                   S2z)) +
                   2 * params.PhiDOT2 * params.r2 *
                       (Complex(0, 14) * (1 + delta - 2 * Nu) * S1z +
                        (6 * (1 + delta) * kappa1 -
                         (26 + 23 * kappa1 + delta * (26 + 11 * kappa1)) * Nu +
                         4 * (1 + 9 * kappa1) * params.eta2) *
                            params.S1z2 -
                        64 * params.eta2 * S1z * S2z +
                        S2z * (Complex(0, -14) * (-1 + delta + 2 * Nu) +
                               2 * Nu * (-13 + 13 * delta + 2 * Nu) * S2z +
                               kappa2 *
                                   (6 + delta * (-6 + 11 * Nu) +
                                    Nu * (-23 + 36 * Nu)) *
                                   S2z)) +
                   PhiDOT * r * rDOT *
                       (40 * (1 + delta - 2 * Nu) * S1z -
                        Complex(0, 1) *
                            (8 * Nu * (-2 - 2 * delta + Nu) +
                             kappa1 * (3 + delta * (3 + 11 * Nu) +
                                       Nu * (5 + 18 * Nu))) *
                            params.S1z2 +
                        Complex(0, 20) * (-3 + Nu) * Nu * S1z * S2z +
                        S2z * (-40 * (-1 + delta + 2 * Nu) -
                               Complex(0, 1) *
                                   (8 * Nu * (-2 + 2 * delta + Nu) +
                                    kappa2 * (3 - delta * (3 + 11 * Nu) +
                                              Nu * (5 + 18 * Nu))) *
                                   S2z))))) /
            (24. * params.r4)) + /* Henry et al. QC spinning hereditary terms */
        ((8 * sqrt(0.7142857142857143) * Nu * M_PI * (S1z + S2z) * params.x4) /
         3.));
  }

  else if (vpnorder == 7) {

    return (
        /* Henry et al. QC spin terms */ /* ((-26902*params.eta3*(S1z + S2z) -
 4664*(S1z + delta*S1z + S2z - delta*S2z) + Nu*(3960*(1 +
 delta)*kappa1*params.S1z3 + 3960*(1 + delta)*kappa1*params.S1z2*S2z + S1z*(28921
 - 18889*delta - 3960*(-1 + delta)*kappa2*params.S2z2) + S2z*(28921 + 18889*delta
 - 3960*(-1 + delta)*kappa2*params.S2z2)) - 2*params.eta2*(1351*delta*(S1z - S2z) +
 6*(S1z + S2z)*(6773 + 660*(kappa1*params.S1z2
 + S2z*(-2*S1z + kappa2*S2z)))))*params.x4p5)/(1188.*sqrt(35))
 + */
                                         /* Henry et al. ecc + spin terms */
        (-0.000014029180695847363 *
         (params.Mtot2 *
          (3 * params.r2 *
               (-120 * params.eta3 *
                    (3565 * params.PhiDOT5 * params.r5 +
                     Complex(0, 2321) * params.PhiDOT4 * params.r4 * rDOT +
                     8244 * params.PhiDOT3 * params.r3 * params.rDOT2 -
                     Complex(0, 869) * params.PhiDOT2 * params.r2 *
                         params.rDOT3 -
                     56 * PhiDOT * r * params.rDOT4 -
                     Complex(0, 120) * params.rDOT5) *
                    (S1z + S2z) +
                2475 * params.PhiDOT2 * params.r2 *
                    (6 * params.PhiDOT3 * params.r3 +
                     Complex(0, 77) * params.PhiDOT2 * params.r2 * rDOT -
                     72 * PhiDOT * r * params.rDOT2 +
                     Complex(0, 6) * params.rDOT3) *
                    (S1z + delta * S1z + S2z - delta * S2z) -
                3 * Nu *
                    (Complex(0, 22) * params.PhiDOT4 * params.r4 * rDOT *
                         (Complex(0, 36322) + 5 * (2993 + 3893 * delta) * S1z +
                          5 * (2993 - 3893 * delta) * S2z) -
                     Complex(0, 25) * params.rDOT5 *
                         ((1053 + 443 * delta) * S1z +
                          (1053 - 443 * delta) * S2z) +
                     Complex(0, 44) * params.PhiDOT2 * params.r2 *
                         params.rDOT3 *
                         (Complex(0, -5444) + 5 * (1424 + 849 * delta) * S1z +
                          7120 * S2z - 4245 * delta * S2z) -
                     20 * PhiDOT * r * params.rDOT4 *
                         (Complex(0, -1782) + (5963 + 2969 * delta) * S1z +
                          5963 * S2z - 2969 * delta * S2z) +
                     4 * params.PhiDOT5 * params.r5 *
                         (Complex(0, -86889) + 10 * (2063 + 225 * delta) * S1z +
                          20630 * S2z - 2250 * delta * S2z) +
                     4 * params.PhiDOT3 * params.r3 * params.rDOT2 *
                         (Complex(0, 234861) + 40 * (-1824 + 97 * delta) * S1z -
                          40 * (1824 + 97 * delta) * S2z)) +
                2 * params.eta2 *
                    (params.PhiDOT5 * params.r5 *
                         (Complex(0, -1549757) +
                          300 * (1448 + 1311 * delta) * S1z +
                          300 * (1448 - 1311 * delta) * S2z) +
                     Complex(0, 11) * params.PhiDOT4 * params.r4 * rDOT *
                         (Complex(0, 329548) +
                          15 * (4113 + 1411 * delta) * S1z + 61695 * S2z -
                          21165 * delta * S2z) +
                     Complex(0, 22) * params.PhiDOT2 * params.r2 *
                         params.rDOT3 *
                         (Complex(0, -23971) + 15 * (3829 + 243 * delta) * S1z +
                          57435 * S2z - 3645 * delta * S2z) +
                     Complex(0, 150) * params.rDOT5 *
                         ((-503 + 92 * delta) * S1z -
                          (503 + 92 * delta) * S2z) +
                     10 * PhiDOT * r * params.rDOT4 *
                         (Complex(0, 4565) + 6 * (-6327 + 991 * delta) * S1z -
                          6 * (6327 + 991 * delta) * S2z) +
                     21 * params.PhiDOT3 * params.r3 * params.rDOT2 *
                         (Complex(0, 161403) +
                          10 * (-1897 + 1471 * delta) * S1z -
                          10 * (1897 + 1471 * delta) * S2z))) -
           6 * mass * r *
               (60 * params.eta3 *
                    (2417 * params.PhiDOT3 * params.r3 +
                     Complex(0, 7258) * params.PhiDOT2 * params.r2 * rDOT -
                     4381 * PhiDOT * r * params.rDOT2 +
                     Complex(0, 480) * params.rDOT3) *
                    (S1z + S2z) -
                165 *
                    (1161 * params.PhiDOT3 * params.r3 -
                     Complex(0, 536) * params.PhiDOT2 * params.r2 * rDOT -
                     2412 * PhiDOT * r * params.rDOT2 -
                     Complex(0, 270) * params.rDOT3) *
                    (S1z + delta * S1z + S2z - delta * S2z) +
                2 * Nu *
                    (Complex(0, 1) * params.PhiDOT2 * params.r2 * rDOT *
                         (Complex(0, 1015784) +
                          5 * (30849 + 88721 * delta) * S1z +
                          5 * (30849 - 88721 * delta) * S2z) +
                     params.PhiDOT3 * params.r3 *
                         (Complex(0, -173371) +
                          5 * (61569 + 10789 * delta) * S1z + 307845 * S2z -
                          53945 * delta * S2z) -
                     5 * PhiDOT * r * params.rDOT2 *
                         (Complex(0, -115368) + (177417 + 52307 * delta) * S1z +
                          177417 * S2z - 52307 * delta * S2z) +
                     Complex(0, 100) * params.rDOT3 *
                         ((-1545 + 181 * delta) * S1z -
                          (1545 + 181 * delta) * S2z)) +
                params.eta2 *
                    (20 * PhiDOT * r * params.rDOT2 *
                         (Complex(0, -11187) - 48074 * S1z +
                          11057 * delta * S1z - 48074 * S2z -
                          11057 * delta * S2z) +
                     Complex(0, 725) * params.rDOT3 *
                         (-73 * S1z + 31 * delta * S1z - 73 * S2z -
                          31 * delta * S2z) +
                     params.PhiDOT3 * params.r3 *
                         (Complex(0, 603141) - 543040 * S1z +
                          404620 * delta * S1z -
                          20 * (27152 + 20231 * delta) * S2z) +
                     Complex(0, 1) * params.PhiDOT2 * params.r2 * rDOT *
                         (Complex(0, -1798104) - 648485 * S1z +
                          105755 * delta * S1z -
                          5 * (129697 + 21151 * delta) * S2z))) +
           10 * params.Mtot2 *
               (24 * params.eta3 *
                    (6981 * PhiDOT * r + Complex(0, 1600) * rDOT) *
                    (S1z + S2z) -
                66 * (2027 * PhiDOT * r + Complex(0, 380) * rDOT) *
                    (S1z + delta * S1z + S2z - delta * S2z) +
                Complex(0, 30) * Nu * rDOT *
                    (297 * (1 + delta) * kappa1 * params.S1z3 +
                     297 * (1 + delta) * kappa1 * params.S1z2 * S2z +
                     S1z * (17261 - 1641 * delta -
                            297 * (-1 + delta) * kappa2 * params.S2z2) +
                     S2z * (17261 + 1641 * delta -
                            297 * (-1 + delta) * kappa2 * params.S2z2)) +
                8 * Nu * PhiDOT * r *
                    (Complex(0, -7315) -
                     4455 * (1 + delta) * kappa1 * params.S1z3 +
                     3 * (3881 - 7757 * delta) * S2z -
                     4455 * (1 + delta) * kappa1 * params.S1z2 * S2z +
                     4455 * (-1 + delta) * kappa2 * params.S2z3 +
                     3 * S1z *
                         (3881 + 7757 * delta +
                          1485 * (-1 + delta) * kappa2 * params.S2z2)) +
                3 * params.eta2 *
                    (Complex(0, -5) * rDOT *
                         (S1z * (18793 + 223 * delta +
                                 1188 * kappa1 * params.S1z2) +
                          (18793 - 223 * delta +
                           1188 * (-2 + kappa1) * params.S1z2) *
                              S2z +
                          1188 * (-2 + kappa2) * S1z * params.S2z2 +
                          1188 * kappa2 * params.S2z3) +
                     4 * PhiDOT * r *
                         (Complex(0, -4939) - 23359 * S1z - 5563 * delta * S1z +
                          5940 * kappa1 * params.S1z3 +
                          (-23359 + 5563 * delta +
                           5940 * (-2 + kappa1) * params.S1z2) *
                              S2z +
                          5940 * (-2 + kappa2) * S1z * params.S2z2 +
                          5940 * kappa2 * params.S2z3))))) /
         (sqrt(35) * params.r4)));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hQC_3_m_2(REAL8 Nu, UINT4 vpnorder, REAL8 x, struct kepler_vars params) {

  /* if(vpnorder == 2){
       return(((2*sqrt(0.7142857142857143))/3. - 2*sqrt(0.7142857142857143)*Nu)*
       params.x2);
   }

   else if(vpnorder == 4){
       return((-193/(27.*sqrt(35)) + (145*sqrt(0.7142857142857143)*Nu)/27. -
       (73*sqrt(0.7142857142857143)*params.eta2)/27.)*params.x3);
   } */

  /* else */ if (vpnorder == 5) {
    return ((4 * sqrt(0.7142857142857143) * (1 - 3 * Nu) * M_PI * params.x3p5) /
            3.);
  }

  /* else if(vpnorder == 6){
       return((-1451/(1188.*sqrt(35)) - (17387*Nu)/(1188.*sqrt(35)) +
       (5557*params.eta2)/(66.*sqrt(35)) -
       (763*sqrt(1.4)*params.eta3)/396.)*params.x4);
   } */

  else {
    return 0;
  }
}

static COMPLEX16 hl_3_m_2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z, REAL8 x, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_3_m_2: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_3_m_2(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x, params) +
            hQC_3_m_2(Nu, vpnorder, x, params)) *
           cpolar(1, -2 * Phi);
  }
}

static COMPLEX16 hl_3_m_min2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z, REAL8 x, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_3_m_min2: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((- 4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           conj(hGO_3_m_2(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x, params) +
                hQC_3_m_2(Nu, vpnorder, x, params)) *
           cpolar(1, 2 * Phi);
  }
}

// H31
static COMPLEX16 hGO_3_m_1(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z,
                           REAL8 x, struct kepler_vars params) {
  REAL8 delta = sqrt(1 - 4 * Nu);
  REAL8 kappa1 = 1.0;
  REAL8 kappa2 = 1.0;
  REAL8 r0 = 1.0;
  REAL8 combination_a = (PhiDOT * r + Complex(0, 1) * rDOT);
  REAL8 combination_a2 = combination_a * combination_a;
  REAL8 combination_a3 = combination_a2 * combination_a;

  if (vpnorder == 1) {
    return (delta * (mass * (Complex(0, 7) * PhiDOT * r - 12 * rDOT) -
                     Complex(0, 6) * r * (PhiDOT * r - Complex(0, 1) * rDOT) *
                         pow(PhiDOT * r + Complex(0, 1) * rDOT, 2))) /
           (6. * sqrt(14) * r);
  }

  else if (vpnorder == 3) {
    return (delta *
            (Complex(0, 6) * (-5 + 19 * Nu) * params.r2 *
                 pow(PhiDOT * r - Complex(0, 1) * rDOT, 2) *
                 combination_a3 +
             2 * params.Mtot2 *
                 (Complex(0, 1) * (-101 + 43 * Nu) * PhiDOT * r +
                  (109 - 86 * Nu) * rDOT) +
             3 * mass * r *
                 (Complex(0, -4) * (-9 + 14 * Nu) * params.PhiDOT3 * params.r3 +
                  6 * (2 + 9 * Nu) * params.PhiDOT2 * params.r2 * rDOT -
                  Complex(0, 1) * (33 + 62 * Nu) * PhiDOT * r * params.rDOT2 +
                  4 * (8 + 17 * Nu) * params.rDOT3))) /
           (36. * sqrt(14) * params.r2);
  }

  else if (vpnorder == 4) {
    return (
        (Complex(0, 0.041666666666666664) * params.Mtot2 *
         (4 * mass * (-1 + 5 * Nu) * ((1 + delta) * S1z + (-1 + delta) * S2z) -
          r * (2 * params.rDOT2 *
                   (-6 * (1 + delta) * S1z + 5 * (5 + 3 * delta) * Nu * S1z +
                    6 * S2z - 6 * delta * S2z +
                    5 * (-5 + 3 * delta) * Nu * S2z) +
               params.PhiDOT2 * params.r2 *
                   ((24 + 24 * delta - 87 * Nu + 31 * delta * Nu) * S1z +
                    (-24 + 24 * delta + 87 * Nu + 31 * delta * Nu) * S2z) +
               Complex(0, 2) * PhiDOT * r * rDOT *
                   ((6 + 6 * delta - 31 * Nu + 35 * delta * Nu) * S1z +
                    (-6 + 6 * delta + 31 * Nu + 35 * delta * Nu) * S2z)))) /
            (sqrt(14) * params.r3) -
        ((Complex(0, 0.041666666666666664) *
          (-4 * S1z + 11 * Nu * (S1z - S2z) + 4 * S2z -
           4 * delta * (S1z + S2z) + 13 * delta * Nu * (S1z + S2z)) *
          params.x3) /
         sqrt(14)) +
        ((Complex(0, 0.041666666666666664) *
          (-4 * S1z + 11 * Nu * (S1z - S2z) + 4 * S2z -
           4 * delta * (S1z + S2z) + 13 * delta * Nu * (S1z + S2z)) *
          params.x3) /
         sqrt(14)));
  }

  else if (vpnorder == 5) {
    return (
        (delta *
         (Complex(0, -18) * (183 - 1579 * Nu + 3387 * params.eta2) * params.r3 *
              pow(PhiDOT * r - Complex(0, 1) * rDOT, 3) *
              pow(PhiDOT * r + Complex(0, 1) * rDOT, 4) +
          2 * params.Mtot3 *
              (Complex(0, 1) * (26473 - 27451 * Nu + 9921 * params.eta2) *
                   PhiDOT * r -
               12 * (623 - 732 * Nu + 1913 * params.eta2) * rDOT) +
          2 * params.Mtot2 * r *
              (Complex(0, -1) * (-8641 - 59189 * Nu + 31959 * params.eta2) *
                   params.PhiDOT3 * params.r3 +
               (-32635 - 29345 * Nu + 29541 * params.eta2) * params.PhiDOT2 *
                   params.r2 * rDOT -
               Complex(0, 44) * (-256 + 781 * Nu + 840 * params.eta2) * PhiDOT *
                   r * params.rDOT2 +
               6 * (-756 + 8238 * Nu + 7357 * params.eta2) * params.rDOT3) +
          3 * mass * params.r2 *
              (Complex(0, 2) * (-2479 - 4505 * Nu + 16785 * params.eta2) *
                   params.PhiDOT5 * params.r5 +
               4 * (817 + 1220 * Nu - 7449 * params.eta2) * params.PhiDOT4 *
                   params.r4 * rDOT +
               Complex(0, 6) * (-1679 + 1469 * Nu + 12233 * params.eta2) *
                   params.PhiDOT3 * params.r3 * params.rDOT2 -
               32 * (-460 + 421 * Nu + 2514 * params.eta2) * params.PhiDOT2 *
                   params.r2 * params.rDOT3 +
               Complex(0, 1) * (-9851 + 17954 * Nu + 40968 * params.eta2) *
                   PhiDOT * r * params.rDOT4 -
               12 * (-771 + 1126 * Nu + 3616 * params.eta2) * params.rDOT5))) /
            (4752. * sqrt(14) *
             params.r3) /* Henry et al. QC spin terms */ /* +((Complex(0,-0.041666666666666664)*(kappa1*(5
- 4*Nu + delta*(5 + 6*Nu))*params.S1z2 + S2z*(-12*delta*Nu*S1z + kappa2*(-5 +
5*delta + 4*Nu + 6*delta*Nu)*S2z))* params.x3p5)/sqrt(14)) */
        +
        /* Henry et al. ecc spin terms */ (
            (params.Mtot3 *
             (kappa1 *
                  (Complex(0, -1) *
                       (-13 - 13 * delta + 68 * Nu + 42 * delta * Nu) * PhiDOT *
                       r -
                   34 * (1 + delta) * rDOT + 4 * (26 + 9 * delta) * Nu * rDOT) *
                  params.S1z2 +
              S2z * (12 * delta * Nu * (Complex(0, 7) * PhiDOT * r - 6 * rDOT) *
                         S1z +
                     kappa2 *
                         (Complex(0, -1) *
                              (13 - 13 * delta - 68 * Nu + 42 * delta * Nu) *
                              PhiDOT * r +
                          2 * (17 - 17 * delta - 52 * Nu + 18 * delta * Nu) *
                              rDOT) *
                         S2z))) /
            (24. * sqrt(14) * params.r3)));
  }

  else if (vpnorder == 6) {
    return (
        (delta * params.Mtot2 * Nu *
         (668 * params.Mtot2 -
          2 * mass * r *
              (727 * params.PhiDOT2 * params.r2 -
               Complex(0, 99) * PhiDOT * r * rDOT + 452 * params.rDOT2) +
          params.r2 * (-499 * params.PhiDOT4 * params.r4 +
                       Complex(0, 1534) * params.PhiDOT3 * params.r3 * rDOT +
                       3072 * params.PhiDOT2 * params.r2 * params.rDOT2 -
                       Complex(0, 680) * PhiDOT * r * params.rDOT3 +
                       1000 * params.rDOT4))) /
            (180. * sqrt(14) *
             params.r4) /*Henry et al QC spin terms */ /* +((Complex(0,0.004629629629629629)*(70*(S1z
- S2z
+ delta*(S1z + S2z)) - params.eta2*(931*(S1z - S2z) + 45*delta*(S1z + S2z)) -
Nu*(59*(-S1z + S2z) + 99*delta*(S1z + S2z)))*params.x4)/sqrt(14)) */
        + /* Henry et al. ecc spin terms */ (
              (Complex(0, 0.0023148148148148147) * params.Mtot2 *
               (2 * params.Mtot2 *
                    ((252 * (1 + delta) - (1277 + 1279 * delta) * Nu +
                      8 * (12 + 47 * delta) * params.eta2) *
                         S1z +
                     (252 * (-1 + delta) + (1277 - 1279 * delta) * Nu +
                      8 * (-12 + 47 * delta) * params.eta2) *
                         S2z) +
                3 * params.r2 *
                    (2 * params.rDOT4 *
                         ((30 + Nu * (-187 + 318 * Nu) +
                           delta * (30 + Nu * (-101 + 122 * Nu))) *
                              S1z +
                          (-30 + (187 - 318 * Nu) * Nu +
                           delta * (30 + Nu * (-101 + 122 * Nu))) *
                              S2z) +
                     2 * params.PhiDOT4 * params.r4 *
                         ((90 - Nu * (28 + 579 * Nu) +
                           delta * (90 + Nu * (-800 + 551 * Nu))) *
                              S1z +
                          (-90 + Nu * (28 + 579 * Nu) +
                           delta * (90 + Nu * (-800 + 551 * Nu))) *
                              S2z) +
                     Complex(0, 2) * PhiDOT * r * params.rDOT3 *
                         ((186 - Nu * (745 + 354 * Nu) +
                           delta * (186 + Nu * (-191 + 554 * Nu))) *
                              S1z +
                          (-186 + Nu * (745 + 354 * Nu) +
                           delta * (186 + Nu * (-191 + 554 * Nu))) *
                              S2z) +
                     3 * params.PhiDOT2 * params.r2 * params.rDOT2 *
                         ((32 + Nu * (-451 + 230 * Nu) +
                           delta * (32 + Nu * (691 + 626 * Nu))) *
                              S1z +
                          (-32 + (451 - 230 * Nu) * Nu +
                           delta * (32 + Nu * (691 + 626 * Nu))) *
                              S2z) +
                     Complex(0, 2) * params.PhiDOT3 * params.r3 * rDOT *
                         ((-12 + Nu * (-341 + 315 * Nu) +
                           delta * (-12 + Nu * (-91 + 1213 * Nu))) *
                              S1z +
                          (12 + (341 - 315 * Nu) * Nu +
                           delta * (-12 + Nu * (-91 + 1213 * Nu))) *
                              S2z)) -
                2 * mass * r *
                    (2 * params.PhiDOT2 * params.r2 *
                         ((-312 * (1 + delta) + 2 * (827 - 923 * delta) * Nu +
                           5 * (-201 + 131 * delta) * params.eta2) *
                              S1z +
                          (-312 * (-1 + delta) - 2 * (827 + 923 * delta) * Nu +
                           5 * (201 + 131 * delta) * params.eta2) *
                              S2z) +
                     2 * params.rDOT2 *
                         ((105 + Nu * (-541 + 924 * Nu) +
                           delta * (105 + Nu * (-77 + 494 * Nu))) *
                              S1z +
                          (-105 + (541 - 924 * Nu) * Nu +
                           delta * (105 + Nu * (-77 + 494 * Nu))) *
                              S2z) +
                     Complex(0, 1) * PhiDOT * r * rDOT *
                         ((1104 - 7 * Nu * (439 + 597 * Nu) +
                           delta * (1104 + Nu * (-3071 + 3083 * Nu))) *
                              S1z +
                          (-1104 + 7 * Nu * (439 + 597 * Nu) +
                           delta * (1104 + Nu * (-3071 + 3083 * Nu))) *
                              S2z)))) /
              (sqrt(14) * params.r4)));
  }

  else if (vpnorder == 7) {

    return (/* (Complex(0,0.001388888888888889)*params.x4p5*(10*(32*(1 +
      delta)*kappa1
      - (36*(1 + delta) + (65 + delta)*kappa1)*Nu + 2*(150 + 33*delta + (76 +
      22*delta)*kappa1)*params.eta2)*params.S1z2 + S2z*(Complex(0,-315)*delta*Nu +
      30*(4 - 11*Nu + delta*(-4 + 13*Nu))*M_PI + 60*Nu*(6 - 50*Nu + delta*(-6 +
      11*Nu))*S2z + 10*kappa2*(-32 + (65 - 152*Nu)*Nu + delta*(32 + Nu*(-1 +
      44*Nu)))*S2z + Complex(0,3)*(-117 + 117*delta + 199*Nu - 80*log(2)) -
        Complex(0,60)*(-11*Nu + delta*(-4 + 13*Nu))*log(2)) + S1z*(30*(-4 +
      11*Nu + delta*(-4 + 13*Nu))*M_PI + 5*delta*Nu*(Complex(0,-63)
        + 8*(-34 + 11*Nu)*S2z) - Complex(0,60)*(11*Nu + delta*(-4 +
      13*Nu))*log(2) + Complex(0,3)*(117 + 117*delta - 199*Nu +
      80*log(2)))))/sqrt(14)
        + */
            /* Henry et al. ecc + spin terms */
            ((Complex(0, 1) *
                  (1513512 * params.Mtot3 *
                       (2 * mass *
                            (PhiDOT * r *
                                 (Complex(0, 1) * (-97 + 631 * Nu) * S1z +
                                  5 *
                                      (8 + 16 * Nu * (-7 + 15 * Nu) +
                                       3 * kappa1 *
                                           (-39 + Nu * (149 + 4 * Nu))) *
                                      params.S1z2 +
                                  S2z * (Complex(0, 97) - Complex(0, 631) * Nu -
                                         5 *
                                             (8 + 16 * Nu * (-7 + 15 * Nu) +
                                              3 * kappa2 *
                                                  (-39 + Nu * (149 + 4 * Nu))) *
                                             S2z)) +
                             rDOT *
                                 (2 * (-18 + 83 * Nu) * S1z -
                                  Complex(0, 5) *
                                      (-4 * (6 + Nu * (-25 + 7 * Nu)) +
                                       kappa1 *
                                           (155 + Nu * (-373 + 164 * Nu))) *
                                      params.S1z2 +
                                  S2z *
                                      (36 - 166 * Nu +
                                       Complex(0, 5) *
                                           (-4 * (6 + Nu * (-25 + 7 * Nu)) +
                                            kappa2 * (155 +
                                                      Nu * (-373 + 164 * Nu))) *
                                           S2z))) +
                        r * (2 * params.rDOT3 *
                                 (S1z *
                                      (18 - 110 * Nu -
                                       Complex(0, 5) *
                                           (69 * kappa1 - 214 * kappa1 * Nu +
                                            4 * (5 + 4 * kappa1) * params.eta2) *
                                           S1z) +
                                  2 * (-9 + 55 * Nu) * S2z +
                                  Complex(0, 5) *
                                      (69 * kappa2 - 214 * kappa2 * Nu +
                                       4 * (5 + 4 * kappa2) * params.eta2) *
                                      params.S2z2) +
                             params.PhiDOT3 * params.r3 *
                                 (S1z *
                                      (Complex(0, 255) - Complex(0, 1403) * Nu +
                                       10 *
                                           (28 * (3 - 8 * Nu) * Nu +
                                            kappa1 *
                                                (51 +
                                                 2 * Nu * (-91 + 118 * Nu))) *
                                           S1z) +
                                  Complex(0, 1) * (-255 + 1403 * Nu) * S2z -
                                  10 *
                                      (28 * (3 - 8 * Nu) * Nu +
                                       kappa2 *
                                           (51 + 2 * Nu * (-91 + 118 * Nu))) *
                                      params.S2z2) +
                             2 * params.PhiDOT2 * params.r2 * rDOT *
                                 (S1z * (255 - 1079 * Nu +
                                         Complex(0, 5) *
                                             (2 * (60 - 361 * Nu) * Nu +
                                              kappa1 *
                                                  (6 + Nu * (-47 + 164 * Nu))) *
                                             S1z) +
                                  (-255 + 1079 * Nu) * S2z -
                                  Complex(0, 5) *
                                      (2 * (60 - 361 * Nu) * Nu +
                                       kappa2 * (6 + Nu * (-47 + 164 * Nu))) *
                                      params.S2z2) +
                             PhiDOT * r * params.rDOT2 *
                                 (Complex(0, 4) * (-114 + 781 * Nu) * S1z +
                                  5 *
                                      (-213 * kappa1 - 72 * Nu +
                                       278 * kappa1 * Nu +
                                       8 * (7 + 44 * kappa1) * params.eta2) *
                                      params.S1z2 +
                                  S2z *
                                      (Complex(0, 456) + 1065 * kappa2 * S2z +
                                       2 * Nu *
                                           (Complex(0, -1562) - 5 *
                                                                    (-36 +
                                                                     139 *
                                                                         kappa2 +
                                                                     4 * (7 + 44 * kappa2) *
                                                                         Nu) *
                                                                    S2z))))) +
                   delta *
                       (52920 *
                            (-4083 +
                             Nu * (58311 + Nu * (-269240 + 405617 * Nu))) *
                            params.r4 *
                            pow(PhiDOT * r - Complex(0, 1) * rDOT, 4) *
                            pow(PhiDOT * r + Complex(0, 1) * rDOT, 5) +
                        840 * params.Mtot2 * params.r2 *
                            ((-2555489 +
                              7 * Nu *
                                  (820078 + Nu * (-6623390 + 4948497 * Nu))) *
                                 params.PhiDOT5 * params.r5 +
                             Complex(0, 1) *
                                 (3537631 +
                                  7 * Nu *
                                      (-2817653 +
                                       Nu * (-7052042 + 4017147 * Nu))) *
                                 params.PhiDOT4 * params.r4 * rDOT +
                             3 * (-1428997 + 7 * Nu * (-1230747 + Nu * (-237418 + 4061717 * Nu))) *
                                 params.PhiDOT3 * params.r3 * params.rDOT2 +
                             Complex(0, 1) *
                                 (-5153011 +
                                  7 * Nu *
                                      (-2375327 +
                                       9 * Nu * (218846 + 1640185 * Nu))) *
                                 params.PhiDOT2 * params.r2 * params.rDOT3 +
                             (-7761899 + 7 * Nu *
                                             (2892563 + 5998602 * Nu +
                                              7493619 * params.eta2)) *
                                 PhiDOT * r * params.rDOT4 +
                             Complex(0, 3) *
                                 (-2422057 +
                                  7 * Nu *
                                      (501045 +
                                       Nu * (2033141 + 2771816 * Nu))) *
                                 params.rDOT5) -
                        8820 * mass * params.r3 *
                            (2 *
                                 (111737 +
                                  Nu * (-366573 +
                                        Nu * (-618923 + 2278593 * Nu))) *
                                 params.PhiDOT7 * params.r7 +
                             Complex(0, 2) *
                                 (101844 + Nu * (-273675 - 871630 * Nu +
                                                 2069774 * params.eta2)) *
                                 params.PhiDOT6 * params.r6 * rDOT +
                             2 * (341322 + Nu * (-1429938 + Nu * (-1206083 + 7690681 * Nu))) *
                                 params.PhiDOT5 * params.r5 * params.rDOT2 +
                             Complex(0, 8) *
                                 (90241 + 2 * Nu *
                                              (-206022 +
                                               Nu * (-62113 + 1003558 * Nu))) *
                                 params.PhiDOT4 * params.r4 * params.rDOT3 +
                             2 * (410547 + Nu * (-2269686 + Nu * (762091 + 8400052 * Nu))) *
                                 params.PhiDOT3 * params.r3 * params.rDOT4 +
                             Complex(0, 4) *
                                 (217935 +
                                  2 * Nu *
                                      (-573699 +
                                       5 * Nu * (18671 + 445748 * Nu))) *
                                 params.PhiDOT2 * params.r2 * params.rDOT5 +
                             (333969 +
                              2 * Nu *
                                  (-981471 + 4 * Nu * (154039 + 750016 * Nu))) *
                                 PhiDOT * r * params.rDOT6 +
                             Complex(0, 24) *
                                 (13245 +
                                  2 * Nu *
                                      (-37005 + Nu * (14251 + 130160 * Nu))) *
                                 params.rDOT7) +
                        2 * params.Mtot4 *
                            (Complex(0, -4178597424) * rDOT +
                             Complex(0, 84) * rDOT *
                                 (38468500 * params.eta3 +
                                  Complex(0, 648648) * (S1z + S2z) -
                                  90090 * ((-24 + 155 * kappa1) * params.S1z2 +
                                           (-24 + 155 * kappa2) * params.S2z2) -
                                  420 * params.eta2 *
                                      (-122855 +
                                       3003 *
                                           ((2 + 11 * kappa1) * params.S1z2 -
                                            18 * S1z * S2z +
                                            (2 + 11 * kappa2) * params.S2z2)) +
                                  3 * Nu *
                                      (-103100846 - 1846845 * M_PI2 -
                                       Complex(0, 564564) * S2z +
                                       6006 * (S1z * (Complex(0, -94) +
                                                      5 * (-52 + 63 * kappa1) *
                                                          S1z) -
                                               20 * S1z * S2z +
                                               5 * (-52 + 63 * kappa2) *
                                                   params.S2z2))) +
                             PhiDOT * r *
                                 (1176172480 * params.eta3 +
                                  8 * (74084729 -
                                       189189 * S1z *
                                           (Complex(0, 97) +
                                            5 * (-8 + 117 * kappa1) * S1z) -
                                       189189 * S2z *
                                           (Complex(0, 97) +
                                            5 * (-8 + 117 * kappa2) * S2z)) -
                                  176400 * params.eta2 *
                                      (11251 +
                                       429 *
                                           ((2 + 13 * kappa1) * params.S1z2 -
                                            22 * S1z * S2z +
                                            (2 + 13 * kappa2) * params.S2z2)) +
                                  147 * Nu *
                                      (-65012788 + 4485195 * M_PI2 +
                                       Complex(0, 4499352) * S2z +
                                       10296 *
                                           (S1z * (Complex(0, 437) +
                                                   15 * (-32 + 71 * kappa1) *
                                                       S1z) -
                                            3860 * S1z * S2z +
                                            15 * (-32 + 71 * kappa2) *
                                                params.S2z2)))) -
                        3 * params.Mtot3 * r *
                            (Complex(0, -4) * params.rDOT3 *
                                 (601018232 - 1359334480 * params.eta3 -
                                  756756 * S1z *
                                      (Complex(0, 6) + 115 * kappa1 * S1z) -
                                  756756 * S2z *
                                      (Complex(0, 6) + 115 * kappa2 * S2z) +
                                  231 * Nu *
                                      (8490448 + 503685 * M_PI2 +
                                       Complex(0, 80808) * S2z +
                                       2184 * (S1z * (Complex(0, 37) +
                                                      190 * kappa1 * S1z) +
                                               70 * S1z * S2z +
                                               190 * kappa2 * params.S2z2)) +
                                  58800 * params.eta2 *
                                      (-62596 +
                                       429 *
                                           ((-1 + 5 * kappa1) * params.S1z2 -
                                            12 * S1z * S2z +
                                            (-1 + 5 * kappa2) * params.S2z2))) -
                             Complex(0, 14) * params.PhiDOT2 * params.r2 *
                                 rDOT *
                                 (-229522160 * params.eta3 +
                                  8 * (48303859 +
                                       135135 * S1z *
                                           (Complex(0, -17) +
                                            2 * kappa1 * S1z) +
                                       135135 * S2z *
                                           (Complex(0, -17) +
                                            2 * kappa2 * S2z)) +
                                  2520 * params.eta2 *
                                      (100913 +
                                       286 *
                                           ((-31 + 5 * kappa1) * params.S1z2 -
                                            72 * S1z * S2z +
                                            (-31 + 5 * kappa2) * params.S2z2)) +
                                  7 * Nu *
                                      (125038052 + 2374515 * M_PI2 +
                                       Complex(0, 5858424) * S2z -
                                       10296 * (S1z * (Complex(0, -569) +
                                                       25 * (-24 + 7 * kappa1) *
                                                           S1z) +
                                                700 * S1z * S2z +
                                                25 * (-24 + 7 * kappa2) *
                                                    params.S2z2))) +
                             4 * PhiDOT * r * params.rDOT2 *
                                 (-1095987374 + 1035895280 * params.eta3 +
                                  378378 * S1z *
                                      (Complex(0, 152) + 355 * kappa1 * S1z) +
                                  378378 * S2z *
                                      (Complex(0, 152) + 355 * kappa2 * S2z) -
                                  490 * params.eta2 *
                                      (-5802767 +
                                       5148 *
                                           ((2 + 23 * kappa1) * params.S1z2 -
                                            42 * S1z * S2z +
                                            (2 + 23 * kappa2) * params.S2z2)) -
                                  77 * Nu *
                                      (42451610 + 1511055 * M_PI2 +
                                       Complex(0, 3623256) * S2z -
                                       6552 * (S1z * (Complex(0, -553) +
                                                      5 * (18 + 37 * kappa1) *
                                                          S1z) +
                                               965 * S1z * S2z +
                                               5 * (18 + 37 * kappa2) *
                                                   params.S2z2))) +
                             7 * params.PhiDOT3 * params.r3 *
                                 (512893080 * params.eta3 -
                                  136 *
                                      (-2089567 +
                                       135135 * S1z *
                                           (Complex(0, 1) + 2 * kappa1 * S1z) +
                                       135135 * S2z *
                                           (Complex(0, 1) + 2 * kappa2 * S2z)) -
                                  560 * params.eta2 *
                                      (2457671 +
                                       2574 *
                                           ((11 + 53 * kappa1) * params.S1z2 -
                                            84 * S1z * S2z +
                                            (11 + 53 * kappa2) * params.S2z2)) +
                                  3 * Nu *
                                      (16621605 * M_PI2 +
                                       8 * (27468722 +
                                            Complex(0, 2681679) * S2z +
                                            3003 * (S1z * (Complex(0, 893) -
                                                           840 * S1z +
                                                           800 * kappa1 * S1z) -
                                                    3160 * S1z * S2z +
                                                    40 * (-21 + 20 * kappa2) *
                                                        params.S2z2))))))) +
              74954880 * delta * params.Mtot3 *
                  (mass * (Complex(0, -22) * PhiDOT * r - 24 * rDOT) +
                   3 * r *
                       (Complex(0, 7) * params.PhiDOT3 * params.r3 +
                        14 * params.PhiDOT2 * params.r2 * rDOT -
                        Complex(0, 8) * PhiDOT * r * params.rDOT2 +
                        4 * params.rDOT3)) *
                  log(r / r0)) /
             (3.6324288e8 * sqrt(14) *
              params.r4)) + /* Henry et al. QC spinning hereditary terms */
            ((((-4 + 11 * Nu + delta * (-4 + 13 * Nu)) * S1z +
               (4 - 11 * Nu + delta * (-4 + 13 * Nu)) * S2z) *
              params.x4p5 * (Complex(0, 1) * M_PI + log(4))) /
             (24. * sqrt(14))));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hQC_3_m_1(REAL8 Nu, UINT4 vpnorder, REAL8 x, struct kepler_vars params) {
  REAL8 delta = sqrt(1 - 4 * Nu);
  double EulerGamma = 0.5772156649015329;

  /* if(vpnorder == 1){
        return((Complex(0,0.16666666666666666)*delta*params.x1p5)/sqrt(14));
    } */

  /* else if(vpnorder == 3){
        return((Complex(0,-0.2222222222222222)*sqrt(0.2857142857142857)*delta -
        (Complex(0,0.1111111111111111)*delta*Nu)/sqrt(14))*params.x2p5);
    } */

  /* else */ if (vpnorder == 4) {
    return ((Complex(0, 0.16666666666666666) * delta * M_PI * params.x3) /
                sqrt(14) +
            (delta * params.x3 * log(2)) / (3. * sqrt(14)));
  }

  /* else if(vpnorder == 5){
        return(((Complex(0,0.5109427609427609)*delta)/sqrt(14) -
        Complex(0,0.11447811447811448)*sqrt(0.2857142857142857)*delta*Nu -
        (Complex(0,0.2079124579124579)*delta*params.eta2)/sqrt(14))*
        params.x3p5);
    } */

  else if (vpnorder == 6) {
    return ((Complex(0, -0.027777777777777776) * delta * (16 + 7 * Nu) *
             params.x4 * (M_PI - Complex(0, 2) * log(2))) /
            sqrt(14));
  }

  else if (vpnorder == 7) {
    return ((Complex(0, -2.752978943455134e-9) * delta * params.x4p5 *
             (-430135880 + 74954880 * EulerGamma + 681626456 * Nu -
              641035640 * params.eta2 + 68698000 * params.eta3 +
              Complex(0, 47279232) * M_PI - 10090080 * M_PI2 -
              38783745 * Nu * M_PI2 + 244468224 * log(2) +
              Complex(0, 121080960) * M_PI * log(2) +
              121080960 * log(2) * log(2) + 37477440 * log(x))) /
            sqrt(14));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hl_3_m_1(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z, REAL8 x, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_4_m_1: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_3_m_1(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x, params) +
            hQC_3_m_1(Nu, vpnorder, x, params)) *
           cpolar(1, -1 * Phi);
  }
}

static COMPLEX16 hl_3_m_min1(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z, REAL8 x, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_4_m_min1: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((- 4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           conj(hGO_3_m_1(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x, params) +
                hQC_3_m_1(Nu, vpnorder, x, params)) *
           cpolar(1, 1 * Phi);
  }
}

// H44
static COMPLEX16 hGO_4_m_4(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z, struct kepler_vars params) {
  REAL8 delta = sqrt(1 - 4 * Nu);
  REAL8 kappa1 = 1.0;
  REAL8 kappa2 = 1.0;

  if (vpnorder == 2) {
    return ((sqrt(0.7142857142857143) * (-1 + 3 * Nu) *
             (7 * params.Mtot2 +
              6 * params.r2 * pow(PhiDOT * r + Complex(0, 1) * rDOT, 4) +
              3 * mass * r *
                  (17 * params.PhiDOT2 * params.r2 +
                   Complex(0, 18) * PhiDOT * r * rDOT - 6 * params.rDOT2))) /
            (36. * params.r2));
  }

  else if (vpnorder == 4) {
    return ((40 * params.Mtot3 * (314 - 987 * Nu + 195 * params.eta2) -
             60 * (23 - 159 * Nu + 291 * params.eta2) * params.r3 *
                 (PhiDOT * r - Complex(0, 1) * rDOT) *
                 pow(PhiDOT * r + Complex(0, 1) * rDOT, 5) +
             params.Mtot2 * r *
                 ((53143 - 199660 * Nu + 127500 * params.eta2) * params.PhiDOT2 *
                      params.r2 +
                  Complex(0, 24) * (967 - 4615 * Nu + 5935 * params.eta2) *
                      PhiDOT * r * rDOT -
                  10 * (290 - 2033 * Nu + 4365 * params.eta2) * params.rDOT2) -
             3 * mass * params.r2 *
                 ((613 - 920 * Nu + 6420 * params.eta2) * params.PhiDOT4 *
                      params.r4 -
                  Complex(0, 8) * (-976 + 1745 * Nu + 3150 * params.eta2) *
                      params.PhiDOT3 * params.r3 * rDOT +
                  2 * (-6141 + 8980 * Nu + 31500 * params.eta2) *
                      params.PhiDOT2 * params.r2 * params.rDOT2 +
                  Complex(0, 4) * (-1853 + 1730 * Nu + 13230 * params.eta2) *
                      PhiDOT * r * params.rDOT3 -
                  20 * (-83 + 30 * Nu + 762 * params.eta2) * params.rDOT4)) /
            (1584. * sqrt(35) * params.r3));
  }

  else if (vpnorder == 5) {
    return ((params.Mtot2 * Nu *
             (6 * mass * (Complex(0, -43) * PhiDOT * r + 9 * rDOT) +
              r * (Complex(0, -734) * params.PhiDOT3 * params.r3 +
                   129 * params.PhiDOT2 * params.r2 * rDOT +
                   Complex(0, 156) * PhiDOT * r * params.rDOT2 -
                   26 * params.rDOT3))) /
                (24. * sqrt(35) * params.r3)
            /* Henry et al. QC spin terms */
            /* +((-32*(39*delta*Nu*(S1z - S2z) + 41*Nu*(S1z + S2z) -
      42*params.eta2*(S1z + S2z)
      - 10*(S1z + delta*S1z + S2z - delta*S2z))*params.x3p5)/(27.*sqrt(35))) */
            + /* Henry et al. ecc spin terms */ (
                  (params.Mtot2 *
                   (Complex(0, -3) * params.PhiDOT2 * params.r3 * rDOT *
                        ((-250 + 1221 * Nu - 1512 * params.eta2 +
                          delta * (-250 + 849 * Nu)) *
                             S1z +
                         (-250 + delta * (250 - 849 * Nu) + 1221 * Nu -
                          1512 * params.eta2) *
                             S2z) -
                    2 * mass * PhiDOT * r *
                        ((-130 + 757 * Nu - 1224 * params.eta2 +
                          delta * (-130 + 513 * Nu)) *
                             S1z +
                         (-130 + delta * (130 - 513 * Nu) + 757 * Nu -
                          1224 * params.eta2) *
                             S2z) -
                    Complex(0, 2) * mass * rDOT *
                        ((-100 + 577 * Nu - 864 * params.eta2 +
                          delta * (-100 + 333 * Nu)) *
                             S1z +
                         (-100 + delta * (100 - 333 * Nu) + 577 * Nu -
                          864 * params.eta2) *
                             S2z) -
                    6 * params.PhiDOT3 * params.r4 *
                        ((-65 + 263 * Nu - 291 * params.eta2 +
                          delta * (-65 + 282 * Nu)) *
                             S1z +
                         (-65 + delta * (65 - 282 * Nu) + 263 * Nu -
                          291 * params.eta2) *
                             S2z) +
                    12 * PhiDOT * params.r2 * params.rDOT2 *
                        ((-40 + 201 * Nu - 252 * params.eta2 +
                          delta * (-40 + 129 * Nu)) *
                             S1z +
                         (-40 + delta * (40 - 129 * Nu) + 201 * Nu -
                          252 * params.eta2) *
                             S2z) +
                    Complex(0, 6) * r * params.rDOT3 *
                        ((-20 + 107 * Nu - 144 * params.eta2 +
                          delta * (-20 + 63 * Nu)) *
                             S1z +
                         (-20 + delta * (20 - 63 * Nu) + 107 * Nu -
                          144 * params.eta2) *
                             S2z))) /
                  (72. * sqrt(35) * params.r3)));
  }

  else if (vpnorder == 6) {
    return (
        (10 * params.Mtot4 *
             (-4477296 + 12734393 * Nu - 6895 * params.eta2 +
              1043805 * params.eta3) +
         3150 * (-367 + 4337 * Nu - 17462 * params.eta2 + 23577 * params.eta3) *
             params.r4 * pow(PhiDOT * r - Complex(0, 1) * rDOT, 2) *
             pow(PhiDOT * r + Complex(0, 1) * rDOT, 6) +
         2 * params.Mtot3 * r *
             ((-36967579 + 245501977 * Nu - 459916170 * params.eta2 +
               150200680 * params.eta3) *
                  params.PhiDOT2 * params.r2 +
              Complex(0, 4) *
                  (7571073 - 10780154 * Nu - 56898800 * params.eta2 +
                   43665510 * params.eta3) *
                  PhiDOT * r * rDOT -
              10 *
                  (1283609 - 5800627 * Nu + 3725295 * params.eta2 +
                   4771935 * params.eta3) *
                  params.rDOT2) -
         params.Mtot2 * params.r2 *
             ((-28258134 + 3245207 * Nu + 144051250 * params.eta2 +
               136991820 * params.eta3) *
                  params.PhiDOT4 * params.r4 -
              Complex(0, 24) *
                  (2371982 - 7733376 * Nu - 7948185 * params.eta2 +
                   9074870 * params.eta3) *
                  params.PhiDOT3 * params.r3 * rDOT +
              7 *
                  (6557973 - 50558069 * Nu + 59901380 * params.eta2 +
                   104752320 * params.eta3) *
                  params.PhiDOT2 * params.r2 * params.rDOT2 +
              Complex(0, 168) *
                  (52044 - 1084807 * Nu + 1849450 * params.eta2 +
                   4171730 * params.eta3) *
                  PhiDOT * r * params.rDOT3 -
              35 *
                  (1083 - 1246819 * Nu + 2524240 * params.eta2 +
                   5995845 * params.eta3) *
                  params.rDOT4) -
         105 * mass * params.r3 *
             ((116396 - 551405 * Nu + 560658 * params.eta2 +
               293036 * params.eta3) *
                  params.PhiDOT6 * params.r6 +
              Complex(0, 2) *
                  (158192 - 670661 * Nu + 177718 * params.eta2 +
                   2163976 * params.eta3) *
                  params.PhiDOT5 * params.r5 * rDOT +
              (-393665 + 1322392 * Nu + 1589680 * params.eta2 -
               8622660 * params.eta3) *
                  params.PhiDOT4 * params.r4 * params.rDOT2 -
              Complex(0, 8) *
                  (-23048 + 209397 * Nu - 487057 * params.eta2 +
                   260396 * params.eta3) *
                  params.PhiDOT3 * params.r3 * params.rDOT3 -
              (630647 - 3391000 * Nu + 2501958 * params.eta2 +
               7664096 * params.eta3) *
                  params.PhiDOT2 * params.r2 * params.rDOT4 -
              Complex(0, 2) *
                  (218975 - 1037408 * Nu + 148970 * params.eta2 +
                   3699480 * params.eta3) *
                  PhiDOT * r * params.rDOT5 +
              10 *
                  (10233 - 44864 * Nu - 13050 * params.eta2 +
                   203280 * params.eta3) *
                  params.rDOT6)) /
            (1.44144e6 * sqrt(35) * params.r4)
        /* Henry et al. QC spin terms */ /* +((16*sqrt(0.7142857142857143)*(-1 +
  3*Nu)*(kappa1*(1 + delta - 2*Nu)*params.S1z2
  + S2z*(4*Nu*S1z - kappa2*(-1 + delta + 2*Nu)*S2z))*params.x4)/9.) */
        + /* Henry et al. ecc spin terms */ (
              (sqrt(0.7142857142857143) * params.Mtot3 * (-1 + 3 * Nu) *
               (12 * mass +
                r * (53 * params.PhiDOT2 * params.r2 +
                     Complex(0, 26) * PhiDOT * r * rDOT - 8 * params.rDOT2)) *
               (kappa1 * (1 + delta - 2 * Nu) * params.S1z2 +
                S2z * (4 * Nu * S1z + kappa2 * S2z - delta * kappa2 * S2z -
                       2 * kappa2 * Nu * S2z))) /
              (48. * params.r4)));
  }

  else if (vpnorder == 7) {

    return (
        /* Henry et al. QC spin terms */ /* (8*(6630*params.eta3*(S1z + S2z) -
  786*(S1z + delta*S1z + S2z - delta*S2z)
  + 21*Nu*(273*delta*(S1z - S2z) + 173*(S1z + S2z)) -
  2*params.eta2*(1557*delta*(S1z - S2z)
   + 2960*(S1z + S2z)))*params.x4p5)/(297.*sqrt(35))
   + */
        /* Henry et al. ecc+spin terms */ (
            (params.Mtot2 *
             (14 * params.Mtot2 *
                  (120 * params.eta3 *
                       (24635 * PhiDOT * r + Complex(0, 18657) * rDOT) *
                       (S1z + S2z) -
                   60 * (10039 * PhiDOT * r + Complex(0, 7706) * rDOT) *
                       (S1z + delta * S1z + S2z - delta * S2z) +
                   5 * Nu *
                       (PhiDOT * r *
                            (Complex(0, 448616) + 703833 * S1z +
                             505635 * delta * S1z + 703833 * S2z -
                             505635 * delta * S2z) +
                        Complex(0, 4) * rDOT *
                            (Complex(0, 4175) + 123114 * S1z +
                             76938 * delta * S1z + 123114 * S2z -
                             76938 * delta * S2z)) -
                   6 * params.eta2 *
                       (2 * PhiDOT * r *
                            (Complex(0, 217374) + 549175 * S1z +
                             61075 * delta * S1z + 549175 * S2z -
                             61075 * delta * S2z) +
                        Complex(0, 5) * rDOT *
                            (Complex(0, 9861) + 132241 * S1z +
                             18825 * delta * S1z + 132241 * S2z -
                             18825 * delta * S2z))) -
              3 * params.r2 *
                  (1680 * params.eta3 *
                       (2833 * params.PhiDOT5 * params.r5 +
                        Complex(0, 18796) * params.PhiDOT4 * params.r4 * rDOT -
                        13185 * params.PhiDOT3 * params.r3 * params.rDOT2 +
                        Complex(0, 1186) * params.PhiDOT2 * params.r2 *
                            params.rDOT3 -
                        5863 * PhiDOT * r * params.rDOT4 -
                        Complex(0, 2100) * params.rDOT5) *
                       (S1z + S2z) -
                   1050 *
                       (454 * params.PhiDOT5 * params.r5 +
                        Complex(0, 1195) * params.PhiDOT4 * params.r4 * rDOT -
                        1950 * params.PhiDOT3 * params.r3 * params.rDOT2 -
                        Complex(0, 442) * params.PhiDOT2 * params.r2 *
                            params.rDOT3 -
                        384 * PhiDOT * r * params.rDOT4 -
                        Complex(0, 184) * params.rDOT5) *
                       (S1z + delta * S1z + S2z - delta * S2z) -
                   6 * params.eta2 *
                       (2 * params.PhiDOT5 * params.r5 *
                            (Complex(0, -2459811) +
                             35 * (26517 + 10223 * delta) * S1z +
                             (928095 - 357805 * delta) * S2z) -
                        Complex(0, 10) * params.rDOT5 *
                            (Complex(0, 35291) +
                             28 * (2183 + 1155 * delta) * S1z +
                             (61124 - 32340 * delta) * S2z) +
                        Complex(0, 1) * params.PhiDOT2 * params.r2 *
                            params.rDOT3 *
                            (Complex(0, 917901) +
                             1120 * (-191 + 1616 * delta) * S1z -
                             1120 * (191 + 1616 * delta) * S2z) +
                        5 * params.PhiDOT3 * params.r3 * params.rDOT2 *
                            (Complex(0, 85426) +
                             7 * (-148363 + 4365 * delta) * S1z -
                             7 * (148363 + 4365 * delta) * S2z) -
                        4 * PhiDOT * r * params.rDOT4 *
                            (Complex(0, 280067) +
                             70 * (5844 + 4817 * delta) * S1z -
                             70 * (-5844 + 4817 * delta) * S2z) +
                        Complex(0, 1) * params.PhiDOT4 * params.r4 * rDOT *
                            (Complex(0, 10375501) +
                             70 * (65831 + 22871 * delta) * S1z -
                             70 * (-65831 + 22871 * delta) * S2z)) +
                   Nu * (Complex(0, -40) * params.rDOT5 *
                             (Complex(0, 12203) +
                              42 * (843 + 656 * delta) * S1z -
                              42 * (-843 + 656 * delta) * S2z) +
                         Complex(0, 4) * params.PhiDOT2 * params.r2 *
                             params.rDOT3 *
                             (Complex(0, -266071) +
                              210 * (-2279 + 1010 * delta) * S1z -
                              210 * (2279 + 1010 * delta) * S2z) -
                         16 * PhiDOT * r * params.rDOT4 *
                             (Complex(0, 58753) +
                              105 * (2030 + 1947 * delta) * S1z -
                              105 * (-2030 + 1947 * delta) * S2z) +
                         params.PhiDOT5 * params.r5 *
                             (Complex(0, -8997592) +
                              105 * (39959 + 15835 * delta) * S1z -
                              105 * (-39959 + 15835 * delta) * S2z) -
                         10 * params.PhiDOT3 * params.r3 * params.rDOT2 *
                             (Complex(0, 254228) +
                              21 * (65293 + 34551 * delta) * S1z -
                              21 * (-65293 + 34551 * delta) * S2z) +
                         Complex(0, 2) * params.PhiDOT4 * params.r4 * rDOT *
                             (Complex(0, 12351083) +
                              105 * (43193 + 37913 * delta) * S1z -
                              105 * (-43193 + 37913 * delta) * S2z))) +
              mass * r *
                  (2520 * params.eta3 *
                       (8036 * params.PhiDOT3 * params.r3 +
                        Complex(0, 41814) * params.PhiDOT2 * params.r2 * rDOT -
                        30537 * PhiDOT * r * params.rDOT2 -
                        Complex(0, 9064) * params.rDOT3) *
                       (S1z + S2z) -
                   210 *
                       (11849 * params.PhiDOT3 * params.r3 +
                        Complex(0, 31868) * params.PhiDOT2 * params.r2 * rDOT +
                        3572 * PhiDOT * r * params.rDOT2 +
                        Complex(0, 1508) * params.rDOT3) *
                       (S1z + delta * S1z + S2z - delta * S2z) -
                   12 * params.eta2 *
                       (Complex(0, 1) * params.PhiDOT2 * params.r2 * rDOT *
                            (Complex(0, -1150397) +
                             35 * (154765 + 157449 * delta) * S1z +
                             5416775 * S2z - 5510715 * delta * S2z) -
                        PhiDOT * r * params.rDOT2 *
                            (Complex(0, 2306552) +
                             35 * (4931 + 110733 * delta) * S1z + 172585 * S2z -
                             3875655 * delta * S2z) +
                        params.PhiDOT3 * params.r3 *
                            (Complex(0, 12461121) +
                             35 * (18331 + 87381 * delta) * S1z + 641585 * S2z -
                             3058335 * delta * S2z) -
                        Complex(0, 25) * params.rDOT3 *
                            (Complex(0, 36676) +
                             35 * (137 + 1053 * delta) * S1z + 4795 * S2z -
                             36855 * delta * S2z)) +
                   Nu * (Complex(0, -200) * params.rDOT3 *
                             (Complex(0, -2501) +
                              42 * (-308 + 51 * delta) * S1z -
                              42 * (308 + 51 * delta) * S2z) +
                         2 * PhiDOT * r * params.rDOT2 *
                             (Complex(0, -1951984) -
                              105 * (-37907 + 14661 * delta) * S1z +
                              105 * (37907 + 14661 * delta) * S2z) +
                         Complex(0, 8) * params.PhiDOT2 * params.r2 * rDOT *
                             (Complex(0, -10005028) +
                              105 * (40991 + 31689 * delta) * S1z -
                              105 * (-40991 + 31689 * delta) * S2z) +
                         params.PhiDOT3 * params.r3 *
                             (Complex(0, 88418488) +
                              105 * (57793 + 266391 * delta) * S1z -
                              105 * (-57793 + 266391 * delta) * S2z))))) /
            (332640. * sqrt(35) * params.r4)));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hQC_4_m_4(REAL8 Nu, UINT4 vpnorder, REAL8 x, struct kepler_vars params) {

  /* if(vpnorder == 2){
       return(((-16*sqrt(0.7142857142857143))/9. +
       (16*sqrt(0.7142857142857143)*Nu)/3.)*params.x2);
   } */

  /* else if(vpnorder == 4){
       return((4744/(99.*sqrt(35)) - (10184*sqrt(0.7142857142857143)*Nu)/297. +
       (200*sqrt(35)*params.eta2)/99.)*params.x3);
   } */

  /* else */ if (vpnorder == 5) {
    return ((64 * sqrt(0.7142857142857143) * (-1 + 3 * Nu) * params.x3p5 *
             (M_PI + Complex(0, 2) * log(2))) /
            9.);
  }

  /* else if(vpnorder == 6){
       return((-2137342/(45045.*sqrt(35)) + (2176238*Nu)/(6435.*sqrt(35)) -
       (587516*params.eta2)/(1053.*sqrt(35)) +
       (452194*params.eta3)/(3861.*sqrt(35)))*params.x4);
   } */

  else {
    return 0;
  }
}

static COMPLEX16 hl_4_m_4(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z, REAL8 x, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_4_m_4: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_4_m_4(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, params) +
            hQC_4_m_4(Nu, vpnorder, x, params)) *
           cpolar(1, -4 * Phi);
  }
}

static COMPLEX16 hl_4_m_min4(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z, REAL8 x, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_4_m_min4: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           conj(hGO_4_m_4(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, params) +
                hQC_4_m_4(Nu, vpnorder, x, params)) *
           cpolar(1, 4 * Phi);
  }
}

// H43
static COMPLEX16 hGO_4_m_3(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z,
                           REAL8 x, struct kepler_vars params) {
  REAL8 delta = sqrt(1 - 4 * Nu);
  REAL8 kappa1 = 1.0;
  REAL8 kappa2 = 1.0;

  if (vpnorder == 3) {
    return (Complex(0, 0.16666666666666666) * delta * mass * (-1 + 2 * Nu) *
            PhiDOT *
            (4 * mass +
             r * (23 * params.PhiDOT2 * params.r2 +
                  Complex(0, 10) * PhiDOT * r * rDOT - 2 * params.rDOT2))) /
           (sqrt(70) * r);
  }

  else if (vpnorder == 4) {
    return ((Complex(0, -0.041666666666666664) * sqrt(0.35714285714285715) *
             params.Mtot2 * Nu *
             (4 * mass +
              r * (23 * params.PhiDOT2 * params.r2 +
                   Complex(0, 10) * PhiDOT * r * rDOT - 2 * params.rDOT2)) *
             ((-1 + delta) * S1z + S2z + delta * S2z)) /
                params.r3 -
            (Complex(0, -1.125) * sqrt(0.35714285714285715) * Nu *
             ((-1 + delta) * S1z + S2z + delta * S2z) * params.x3) +
            (Complex(0, -1.125) * sqrt(0.35714285714285715) * Nu *
             (-S1z + S2z + delta * (S1z + S2z)) * params.x3));
  }

  else if (vpnorder == 5) {
    return (Complex(0, 0.0012626262626262627) * delta * mass * PhiDOT *
            (2 * params.Mtot2 * (972 - 2293 * Nu + 1398 * params.eta2) +
             2 * mass * r *
                 ((1788 - 9077 * Nu + 13416 * params.eta2) * params.PhiDOT2 *
                      params.r2 +
                  Complex(0, 3) * (-2796 + 5299 * Nu + 1622 * params.eta2) *
                      PhiDOT * r * rDOT -
                  2 * (-1200 + 2545 * Nu + 162 * params.eta2) * params.rDOT2) -
             3 * params.r2 *
                 ((-524 - 489 * Nu + 6392 * params.eta2) * params.PhiDOT4 *
                      params.r4 +
                  Complex(0, 4) * (796 - 1864 * Nu + 133 * params.eta2) *
                      params.PhiDOT3 * params.r3 * rDOT +
                  42 * (-51 + 94 * Nu + 56 * params.eta2) * params.PhiDOT2 *
                      params.r2 * params.rDOT2 +
                  Complex(0, 4) * (-229 + 366 * Nu + 358 * params.eta2) *
                      PhiDOT * r * params.rDOT3 -
                  4 * (-43 + 62 * Nu + 80 * params.eta2) * params.rDOT4))) /
           (sqrt(70) * params.r2);
  }

  else if (vpnorder == 6) {
    return (
        (delta * params.Mtot2 * Nu * PhiDOT *
         (6 * mass * (181 * PhiDOT * r - Complex(0, 89) * rDOT) +
          r * (4847 * params.PhiDOT3 * params.r3 -
               Complex(0, 7338) * params.PhiDOT2 * params.r2 * rDOT -
               408 * PhiDOT * r * params.rDOT2 +
               Complex(0, 112) * params.rDOT3))) /
            (180. * sqrt(70) *
             params.r2) /* Henry et al. QC spin terms */ /* +((Complex(0,0.03409090909090909)
*(220*(S1z - S2z + delta*(S1z + S2z)) + Nu*(-2403*(S1z - S2z) + 203*delta*(S1z +
S2z)) + params.eta2*(3359*(S1z - S2z) + 457*delta*(S1z +
S2z)))*params.x4)/sqrt(70)) */
        + /* Henry et al. ecc spin terms */ (
              (Complex(0, -0.0006313131313131314) * params.Mtot2 *
               (2 * params.Mtot2 *
                    ((-440 + 6801 * Nu - 1428 * params.eta2 +
                      delta * (-440 - 3193 * Nu + 300 * params.eta2)) *
                         S1z +
                     (440 - 6801 * Nu + 1428 * params.eta2 +
                      delta * (-440 - 3193 * Nu + 300 * params.eta2)) *
                         S2z) -
                2 * mass * r *
                    (Complex(0, -3) * PhiDOT * r * rDOT *
                         (delta * (-1320 + 9093 * Nu + 59 * params.eta2) * S1z -
                          5 * (264 - 311 * Nu + 823 * params.eta2) * S1z +
                          delta * (-1320 + 9093 * Nu + 59 * params.eta2) * S2z +
                          5 * (264 - 311 * Nu + 823 * params.eta2) * S2z) -
                     2 * params.rDOT2 *
                         ((220 + 1659 * Nu - 1512 * params.eta2 +
                           delta * (220 - 3067 * Nu + 240 * params.eta2)) *
                              S1z +
                          (-220 - 1659 * Nu + 1512 * params.eta2 +
                           delta * (220 - 3067 * Nu + 240 * params.eta2)) *
                              S2z) +
                     2 * params.PhiDOT2 * params.r2 *
                         ((1826 - 19530 * Nu + 20145 * params.eta2 +
                           delta * (1826 + 1534 * Nu + 567 * params.eta2)) *
                              S1z +
                          (-1826 + 19530 * Nu - 20145 * params.eta2 +
                           delta * (1826 + 1534 * Nu + 567 * params.eta2)) *
                              S2z)) -
                3 * params.r2 *
                    (Complex(0, 3080) * params.PhiDOT3 * params.r3 * rDOT *
                         ((1 + delta) * S1z + (-1 + delta) * S2z) +
                     2 * params.eta2 *
                         (129 * params.PhiDOT2 * params.r2 * params.rDOT2 *
                              (41 * S1z + 23 * delta * S1z - 41 * S2z +
                               23 * delta * S2z) -
                          2 * params.rDOT4 *
                              (149 * S1z + 75 * delta * S1z - 149 * S2z +
                               75 * delta * S2z) +
                          Complex(0, 2) * PhiDOT * r * params.rDOT3 *
                              (925 * S1z + 491 * delta * S1z - 925 * S2z +
                               491 * delta * S2z) -
                          Complex(0, 1) * params.PhiDOT3 * params.r3 * rDOT *
                              (-2105 * S1z + 753 * delta * S1z + 2105 * S2z +
                               753 * delta * S2z) +
                          params.PhiDOT4 * params.r4 *
                              (11617 * S1z + 4847 * delta * S1z - 11617 * S2z +
                               4847 * delta * S2z)) +
                     Nu * (16 * params.PhiDOT4 * params.r4 *
                               (-413 * S1z + 127 * delta * S1z + 413 * S2z +
                                127 * delta * S2z) -
                           2 * params.rDOT4 *
                               (-301 * S1z + 213 * delta * S1z + 301 * S2z +
                                213 * delta * S2z) +
                           Complex(0, 2) * PhiDOT * r * params.rDOT3 *
                               (-1625 * S1z + 1009 * delta * S1z + 1625 * S2z +
                                1009 * delta * S2z) +
                           3 * params.PhiDOT2 * params.r2 * params.rDOT2 *
                               (-2587 * S1z + 1267 * delta * S1z + 2587 * S2z +
                                1267 * delta * S2z) -
                           Complex(0, 2) * params.PhiDOT3 * params.r3 * rDOT *
                               (3635 * S1z + 7981 * delta * S1z - 3635 * S2z +
                                7981 * delta * S2z))))) /
              (sqrt(70) * params.r4)));
  }

  else if (vpnorder == 7) {

    return (
        /*Henry et al. QC spin terms */ /* (Complex(0,-0.020833333333333332)*params.x4p5*(54*(10*(4
  + delta)*params.eta2
  + kappa1*(5 + 5*delta*pow(1 - 2*Nu,2) - 30*Nu + 28*params.eta2))*params.S1z2 +
  S1z*(Complex(0,297) - Nu*(Complex(0,283) + 810*M_PI) + delta*(Complex(0,297) -
  5*Nu*(Complex(0,629) - 162*M_PI + 216*Nu*S2z)) + Complex(0,3240)*(-1 +
  delta)*Nu*atanh(1/5)) + S2z*(Complex(0,297)*(-1 + delta) + Nu*(Complex(0,283)
  + 810*M_PI) + 54*kappa2*(-5 + 5*delta*pow(1 - 2*Nu,2) + 30*Nu -
  28*params.eta2)*S2z + 5*Nu*(Complex(0,-629)*delta + 162*delta*M_PI - 432*Nu*S2z
  + 108*delta*Nu*S2z + Complex(0,648)*(1 + delta)*atanh(1/5)))))/sqrt(70)
  + */
                                        /* Henry et al. ecc+spin terms */
        ((Complex(0, -1.3875013875013875e-6) * mass *
          (5005 * params.Mtot2 *
               (-24 * PhiDOT * params.r2 * params.rDOT2 *
                    (Complex(0, 1) * (-11 + 48 * Nu) * S1z +
                     6 * (-5 + 4 * kappa1) * params.eta2 * params.S1z2 +
                     S2z * (Complex(0, 11) - Complex(0, 48) * Nu -
                            6 * (-5 + 4 * kappa2) * params.eta2 * S2z)) +
                4 * r * params.rDOT3 *
                    ((-11 + 93 * Nu) * S1z -
                     Complex(0, 6) * (-5 + 4 * kappa1) * params.eta2 *
                         params.S1z2 +
                     S2z * (11 - 93 * Nu +
                            Complex(0, 6) * (-5 + 4 * kappa2) * params.eta2 *
                                S2z)) +
                4 * mass * rDOT *
                    ((22 - 111 * Nu) * S1z +
                     Complex(0, 6) * (15 + 8 * kappa1) * params.eta2 *
                         params.S1z2 +
                     S2z * (-22 + 111 * Nu -
                            Complex(0, 6) * (15 + 8 * kappa2) * params.eta2 *
                                S2z)) +
                Complex(0, 6) * params.PhiDOT2 * params.r3 * rDOT *
                    (Complex(0, 1) * (-121 + 963 * Nu) * S1z +
                     6 *
                         (-55 * params.eta2 +
                          kappa1 * (-5 + 30 * Nu + 4 * params.eta2)) *
                         params.S1z2 +
                     S2z * (Complex(0, 121) + 30 * kappa2 * S2z -
                            6 * (-55 + 4 * kappa2) * params.eta2 * S2z -
                            9 * Nu * (Complex(0, 107) + 20 * kappa2 * S2z))) +
                2 * mass * PhiDOT * r *
                    (Complex(0, -1) * (-121 + 633 * Nu) * S1z +
                     6 *
                         (180 * params.eta2 +
                          kappa1 * (9 - 54 * Nu + 28 * params.eta2)) *
                         params.S1z2 -
                     S2z * (Complex(0, 121) + 54 * kappa2 * S2z +
                            24 * (45 + 7 * kappa2) * params.eta2 * S2z -
                            3 * Nu * (Complex(0, 211) + 108 * kappa2 * S2z))) +
                params.PhiDOT3 * params.r4 *
                    (Complex(0, -1) * (-649 + 4767 * Nu) * S1z +
                     6 *
                         (820 * params.eta2 +
                          kappa1 * (75 - 450 * Nu + 364 * params.eta2)) *
                         params.S1z2 -
                     S2z *
                         (Complex(0, 649) + 450 * kappa2 * S2z +
                          24 * (205 + 91 * kappa2) * params.eta2 * S2z -
                          3 * Nu * (Complex(0, 1589) + 900 * kappa2 * S2z)))) +
           delta *
               (2 * mass * PhiDOT * params.r3 *
                    (10 *
                         (234744 - 1010534 * Nu + 1024443 * params.eta2 +
                          5451096 * params.eta3) *
                         params.PhiDOT4 * params.r4 -
                     Complex(0, 3) *
                         (-2426804 + 1512854 * Nu + 4994115 * params.eta2 +
                          610960 * params.eta3) *
                         params.PhiDOT3 * params.r3 * rDOT +
                     (-30341028 + 23936528 * Nu + 89326545 * params.eta2 +
                      19329660 * params.eta3) *
                         params.PhiDOT2 * params.r2 * params.rDOT2 +
                     Complex(0, 21) *
                         (-668008 + 803028 * Nu + 1908955 * params.eta2 +
                          540370 * params.eta3) *
                         PhiDOT * r * params.rDOT3 -
                     14 *
                         (-172143 + 155683 * Nu + 680580 * params.eta2 +
                          111840 * params.eta3) *
                         params.rDOT4) -
                105 * PhiDOT * params.r4 *
                    ((-8280 + 24681 * Nu - 151973 * params.eta2 +
                      624074 * params.eta3) *
                         params.PhiDOT6 * params.r6 +
                     Complex(0, 2) *
                         (-32208 + 248485 * Nu - 524074 * params.eta2 +
                          24546 * params.eta3) *
                         params.PhiDOT5 * params.r5 * rDOT +
                     2 *
                         (48924 - 239802 * Nu + 137447 * params.eta2 +
                          358156 * params.eta3) *
                         params.PhiDOT4 * params.r4 * params.rDOT2 +
                     Complex(0, 4) *
                         (174 + 24488 * Nu - 102039 * params.eta2 +
                          44882 * params.eta3) *
                         params.PhiDOT3 * params.r3 * params.rDOT3 +
                     3 *
                         (10455 - 56490 * Nu + 84504 * params.eta2 +
                          54016 * params.eta3) *
                         params.PhiDOT2 * params.r2 * params.rDOT4 +
                     Complex(0, 2) *
                         (11175 - 52698 * Nu + 57436 * params.eta2 +
                          60808 * params.eta3) *
                         PhiDOT * r * params.rDOT5 -
                     6 *
                         (829 - 3726 * Nu + 3480 * params.eta2 +
                          4640 * params.eta3) *
                         params.rDOT6) +
                params.Mtot2 * r *
                    (20020 * params.rDOT3 * (S1z + S2z) *
                         (-11 + 71 * Nu +
                          Complex(0, 30) * params.eta2 * (S1z + S2z)) -
                     52 * PhiDOT * r * params.rDOT2 *
                         (129150 * params.eta3 +
                          7 * Nu *
                              (31961 + Complex(0, 8580) * S1z +
                               Complex(0, 8580) * S2z) -
                          Complex(0, 3) *
                              (Complex(0, 32671) + 8470 * S1z + 8470 * S2z) -
                          35 * params.eta2 *
                              (33313 + 1980 * params.S1z2 + 3960 * S1z * S2z +
                               1980 * params.S2z2)) +
                     params.PhiDOT3 * params.r3 *
                         (-10566168 - 70869960 * params.eta3 +
                          Complex(0, 3248245) * S1z +
                          2252250 * kappa1 * params.S1z2 +
                          Complex(0, 3248245) * S2z +
                          2252250 * kappa2 * params.S2z2 -
                          7 * Nu *
                              (4818166 + Complex(0, 2480335) * S1z +
                               1287000 * kappa1 * params.S1z2 +
                               Complex(0, 2480335) * S2z +
                               1287000 * kappa2 * params.S2z2) +
                          3850 * params.eta2 *
                              (38873 + 78 * (7 + 30 * kappa1) * params.S1z2 -
                               3588 * S1z * S2z +
                               78 * (7 + 30 * kappa2) * params.S2z2)) -
                     Complex(0, 2) * params.PhiDOT2 * params.r2 * rDOT *
                         (6531280 * params.eta3 +
                          3 * (5382288 + Complex(0, 605605) * S1z +
                               150150 * kappa1 * params.S1z2 +
                               Complex(0, 605605) * S2z +
                               150150 * kappa2 * params.S2z2) +
                          210 * params.eta2 *
                              (322144 + 2145 * (1 + 4 * kappa1) * params.S1z2 -
                               12870 * S1z * S2z +
                               2145 * (1 + 4 * kappa2) * params.S2z2) -
                          21 * Nu *
                              (2826484 + 85800 * kappa1 * params.S1z2 +
                               Complex(0, 515515) * S2z +
                               85800 * kappa2 * params.S2z2 -
                               5005 * S1z * (Complex(0, -103) + 60 * S2z)))) -
                26 * params.Mtot3 *
                    (770 * rDOT *
                         (Complex(0, -90) * params.eta2 * params.S1z2 +
                          S2z * (-22 + 67 * Nu -
                                 Complex(0, 90) * params.eta2 * S2z) +
                          S1z * (-22 + Nu * (67 + Complex(0, 300) * S2z) -
                                 Complex(0, 180) * params.eta2 * S2z)) +
                     PhiDOT * r *
                         (-38076 + 174720 * params.eta3 -
                          Complex(0, 46585) * S1z -
                          20790 * kappa1 * params.S1z2 -
                          Complex(0, 46585) * S2z -
                          20790 * kappa2 * params.S2z2 -
                          1260 * params.eta2 *
                              (158 + 33 * (5 + 2 * kappa1) * params.S1z2 +
                               198 * S1z * S2z +
                               33 * (5 + 2 * kappa2) * params.S2z2) +
                          7 * Nu *
                              (6188 + 11880 * kappa1 * params.S1z2 +
                               Complex(0, 21505) * S2z +
                               11880 * kappa2 * params.S2z2 +
                               55 * S1z * (Complex(0, 391) + 744 * S2z))))))) /
         (sqrt(70) *
          params.r4)) + /* Henry et al. QC spinning hereditary terms */
        ((27 * sqrt(0.35714285714285715) * Nu *
          ((-1 + delta) * S1z + S2z + delta * S2z) * params.x4p5 *
          (Complex(0, -1) * M_PI + 4 * atanh(1 / 5))) /
         8.));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hQC_4_m_3(REAL8 Nu, UINT4 vpnorder, REAL8 x, struct kepler_vars params) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  /* if(vpnorder == 3){
       return(((Complex(0,-4.5)*delta)/sqrt(70) +
   (Complex(0,9)*delta*Nu)/sqrt(70))* params.x2p5);
   } */

  /*  else if(vpnorder == 5){
       return(((Complex(0,15.954545454545455)*delta)/sqrt(70) -
        Complex(0,6.170454545454546)*sqrt(0.7)*delta*Nu +
        (Complex(0,17.863636363636363)*delta*params.eta2)/sqrt(70))*
        params.x3p5);
   } */

  /* else */ if (vpnorder == 6) {
    return ((Complex(0, 13.5) * delta * (-1 + 2 * Nu) * params.x4 *
             (M_PI + Complex(0, 2) * log(1.5))) /
            sqrt(70));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hl_4_m_3(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z, REAL8 x, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_4_m_3: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_4_m_3(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x, params) +
            hQC_4_m_3(Nu, vpnorder, x, params)) *
           cpolar(1, -3 * Phi);
  }
}

static COMPLEX16 hl_4_m_min3(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z, REAL8 x, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_4_m_min3: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           conj(hGO_4_m_3(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x, params) +
                hQC_4_m_3(Nu, vpnorder, x, params)) *
           cpolar(1, 3 * Phi);
  }
}

// H42
static COMPLEX16 hGO_4_m_2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z, struct kepler_vars params) {
  REAL8 delta = sqrt(1 - 4 * Nu);
  REAL8 kappa1 = 1.0;
  REAL8 kappa2 = 1.0;
  REAL8 combination_a = (PhiDOT * r + Complex(0, 1) * rDOT);
  REAL8 combination_a2 = combination_a * combination_a;
  REAL8 combination_a3 = combination_a2 * combination_a;

  if (vpnorder == 2) {
    return (-(sqrt(5) * (-1 + 3 * Nu) *
              (7 * params.Mtot2 -
               6 * params.r2 * (PhiDOT * r - Complex(0, 1) * rDOT) *
                   combination_a3 +
               3 * mass * r *
                   (params.PhiDOT2 * params.r2 +
                    Complex(0, 9) * PhiDOT * r * rDOT - 6 * params.rDOT2))) /
            (126. * params.r2));
  }

  else if (vpnorder == 4) {
    return (-(40 * params.Mtot3 * (314 - 987 * Nu + 195 * params.eta2) +
              60 * (23 - 159 * Nu + 291 * params.eta2) * params.r3 *
                  pow(PhiDOT * r - Complex(0, 1) * rDOT, 2) *
                  pow(PhiDOT * r + Complex(0, 1) * rDOT, 4) +
              params.Mtot2 * r *
                  ((1987 - 11200 * Nu + 12960 * params.eta2) * params.PhiDOT2 *
                       params.r2 +
                   Complex(0, 12) * (967 - 4615 * Nu + 5935 * params.eta2) *
                       PhiDOT * r * rDOT -
                   10 * (290 - 2033 * Nu + 4365 * params.eta2) * params.rDOT2) -
              3 * mass * params.r2 *
                  ((1577 - 7940 * Nu + 9920 * params.eta2) * params.PhiDOT4 *
                       params.r4 +
                   Complex(0, 4) * (-454 - 315 * Nu + 5980 * params.eta2) *
                       params.PhiDOT3 * params.r3 * rDOT -
                   2 * (549 - 2140 * Nu + 2140 * params.eta2) * params.PhiDOT2 *
                       params.r2 * params.rDOT2 +
                   Complex(0, 2) * (-1853 + 1730 * Nu + 13230 * params.eta2) *
                       PhiDOT * r * params.rDOT3 -
                   20 * (-83 + 30 * Nu + 762 * params.eta2) * params.rDOT4)) /
            (5544. * sqrt(5) * params.r3));
  }

  else if (vpnorder == 5) {
    return ((params.Mtot2 * Nu *
             (mass * (Complex(0, 129) * PhiDOT * r - 54 * rDOT) +
              r * (Complex(0, -73) * params.PhiDOT3 * params.r3 +
                   21 * params.PhiDOT2 * params.r2 * rDOT -
                   Complex(0, 78) * PhiDOT * r * params.rDOT2 +
                   26 * params.rDOT3))) /
                (84. * sqrt(5) * params.r3)
            /* Henry et al. QC spin terms */
            /* +((4*(21*delta*Nu*(S1z - S2z) + 59*Nu*(S1z + S2z) -
      78*params.eta2*(S1z + S2z)
      - 10*(S1z + delta*S1z + S2z - delta*S2z))*params.x3p5)/(189.*sqrt(5))) */
            + /* Henry et al. ecc spin terms */ (
                  (params.Mtot2 *
                   (Complex(0, 2) * mass * rDOT *
                        ((-100 - 100 * delta + 577 * Nu + 333 * delta * Nu -
                          864 * params.eta2) *
                             S1z +
                         (-100 + 100 * delta + 577 * Nu - 333 * delta * Nu -
                          864 * params.eta2) *
                             S2z) +
                    4 * mass * PhiDOT * r *
                        ((-70 - 70 * delta + 253 * Nu + 237 * delta * Nu -
                          156 * params.eta2) *
                             S1z +
                         (-70 + 70 * delta + 253 * Nu - 237 * delta * Nu -
                          156 * params.eta2) *
                             S2z) -
                    3 * r *
                        (Complex(0, 2) * params.rDOT3 *
                             ((-20 - 20 * delta + 107 * Nu + 63 * delta * Nu -
                               144 * params.eta2) *
                                  S1z +
                              (-20 + 20 * delta + 107 * Nu - 63 * delta * Nu -
                               144 * params.eta2) *
                                  S2z) +
                         4 * PhiDOT * r * params.rDOT2 *
                             ((-20 - 20 * delta + 63 * Nu + 67 * delta * Nu -
                               16 * params.eta2) *
                                  S1z +
                              (-20 + 20 * delta + 63 * Nu - 67 * delta * Nu -
                               16 * params.eta2) *
                                  S2z) +
                         4 * params.PhiDOT3 * params.r3 *
                             ((5 + 5 * delta - 7 * Nu + 2 * delta * Nu -
                               41 * params.eta2) *
                                  S1z -
                              (-5 + 5 * delta + 7 * Nu + 2 * delta * Nu +
                               41 * params.eta2) *
                                  S2z) -
                         Complex(0, 1) * params.PhiDOT2 * params.r2 * rDOT *
                             ((-10 - 10 * delta - 339 * Nu + 89 * delta * Nu +
                               1048 * params.eta2) *
                                  S1z +
                              (-10 + 10 * delta - 339 * Nu - 89 * delta * Nu +
                               1048 * params.eta2) *
                                  S2z)))) /
                  (504. * sqrt(5) * params.r3)));
  }

  else if (vpnorder == 6) {
    return (
        (-10 * params.Mtot4 *
             (-4477296 + 12734393 * Nu - 6895 * params.eta2 +
              1043805 * params.eta3) +
         3150 * (-367 + 4337 * Nu - 17462 * params.eta2 + 23577 * params.eta3) *
             params.r4 * pow(PhiDOT * r - Complex(0, 1) * rDOT, 3) *
             pow(PhiDOT * r + Complex(0, 1) * rDOT, 5) -
         2 * params.Mtot3 * r *
             (7 *
                  (-100473 + 3430399 * Nu - 9132990 * params.eta2 +
                   2885660 * params.eta3) *
                  params.PhiDOT2 * params.r2 +
              Complex(0, 2) *
                  (7571073 - 10780154 * Nu - 56898800 * params.eta2 +
                   43665510 * params.eta3) *
                  PhiDOT * r * rDOT -
              10 *
                  (1283609 - 5800627 * Nu + 3725295 * params.eta2 +
                   4771935 * params.eta3) *
                  params.rDOT2) +
         7 * params.Mtot2 * params.r2 *
             ((1071402 + 3846989 * Nu - 27339110 * params.eta2 +
               17538420 * params.eta3) *
                  params.PhiDOT4 * params.r4 +
              Complex(0, 12) *
                  (671714 - 1645932 * Nu - 1903365 * params.eta2 +
                   3346250 * params.eta3) *
                  params.PhiDOT3 * params.r3 * rDOT -
              (481563 + 4291961 * Nu - 17137220 * params.eta2 +
               9315720 * params.eta3) *
                  params.PhiDOT2 * params.r2 * params.rDOT2 +
              Complex(0, 12) *
                  (52044 - 1084807 * Nu + 1849450 * params.eta2 +
                   4171730 * params.eta3) *
                  PhiDOT * r * params.rDOT3 -
              5 *
                  (1083 - 1246819 * Nu + 2524240 * params.eta2 +
                   5995845 * params.eta3) *
                  params.rDOT4) -
         105 * mass * params.r3 *
             ((54272 - 58271 * Nu - 815454 * params.eta2 +
               1435572 * params.eta3) *
                  params.PhiDOT6 * params.r6 +
              Complex(0, 1) *
                  (73976 - 157355 * Nu - 811766 * params.eta2 +
                   2935488 * params.eta3) *
                  params.PhiDOT5 * params.r5 * rDOT +
              (72637 - 400832 * Nu + 282028 * params.eta2 +
               1063956 * params.eta3) *
                  params.PhiDOT4 * params.r4 * params.rDOT2 +
              Complex(0, 4) *
                  (90174 - 385167 * Nu - 126419 * params.eta2 +
                   1739072 * params.eta3) *
                  params.PhiDOT3 * params.r3 * params.rDOT3 +
              (-55817 + 58920 * Nu + 989942 * params.eta2 -
               2334016 * params.eta3) *
                  params.PhiDOT2 * params.r2 * params.rDOT4 +
              Complex(0, 1) *
                  (218975 - 1037408 * Nu + 148970 * params.eta2 +
                   3699480 * params.eta3) *
                  PhiDOT * r * params.rDOT5 +
              10 *
                  (-10233 + 44864 * Nu + 13050 * params.eta2 -
                   203280 * params.eta3) *
                  params.rDOT6)) /
            (5.04504e6 * sqrt(5) *
             params.r4) /* Henry et al. QC spin terms */ /* +((2*sqrt(5)*(kappa1*(1
+ delta - 2*Nu
+ 6*params.eta2)*params.S1z2 + S2z*(4*(1 - 3*Nu)*Nu*S1z - kappa2*(-1 + delta + 2*Nu
- 6*params.eta2)*S2z))*params.x4)/63.) */
        +
        /* Henry et al. ecc spin terms */ (
            (sqrt(5) * params.Mtot3 *
             (kappa1 *
                  (mass * (12 + delta * (12 - 34 * Nu) - 58 * Nu +
                           72 * params.eta2) +
                   r * ((5 - delta * (-5 + Nu) - 11 * Nu + 30 * params.eta2) *
                            params.PhiDOT2 * params.r2 +
                        Complex(0, 1) *
                            (13 + delta * (13 - 59 * Nu) - 85 * Nu +
                             78 * params.eta2) *
                            PhiDOT * r * rDOT +
                        4 *
                            (-2 + 11 * Nu - 12 * params.eta2 +
                             delta * (-2 + 7 * Nu)) *
                            params.rDOT2)) *
                  params.S1z2 +
              S2z *
                  (-2 * mass *
                       (-24 * Nu * S1z + 6 * (-1 + delta) * kappa2 * S2z +
                        (29 - 17 * delta) * kappa2 * Nu * S2z +
                        36 * params.eta2 * (2 * S1z - kappa2 * S2z)) +
                   r * (-((-1 + delta) * kappa2 *
                          (5 * params.PhiDOT2 * params.r2 +
                           Complex(0, 13) * PhiDOT * r * rDOT -
                           8 * params.rDOT2) *
                          S2z) -
                        6 * params.eta2 *
                            (5 * params.PhiDOT2 * params.r2 +
                             Complex(0, 13) * PhiDOT * r * rDOT -
                             8 * params.rDOT2) *
                            (2 * S1z - kappa2 * S2z) +
                        Nu * (params.PhiDOT2 * params.r2 *
                                  (20 * S1z + (-11 + delta) * kappa2 * S2z) -
                              4 * params.rDOT2 *
                                  (8 * S1z + (-11 + 7 * delta) * kappa2 * S2z) +
                              Complex(0, 1) * PhiDOT * r * rDOT *
                                  (52 * S1z +
                                   (-85 + 59 * delta) * kappa2 * S2z)))))) /
            (168. * params.r4)));
  }

  else if (vpnorder == 7) {

    return (
        /* Henry et al. QC spin terms */ /* ((-5958*delta*params.eta2*(S1z - S2z)
 + 45*delta*Nu*(-S1z + S2z)
  - 6009*Nu*(S1z + S2z) + 7528*params.eta2*(S1z + S2z) - 1134*params.eta3*(S1z +
 S2z) + 1104*(S1z + delta*S1z + S2z - delta*S2z))*params.x4p5)/(2079.*sqrt(5))
 + */
        /* Henry et al. ecc+spin terms */ (
            (params.Mtot2 *
             (-14 * params.Mtot2 *
                  (60 * params.eta3 *
                       (6830 * PhiDOT * r + Complex(0, 18657) * rDOT) *
                       (S1z + S2z) -
                   30 * (10617 * PhiDOT * r + Complex(0, 7706) * rDOT) *
                       (S1z + delta * S1z + S2z - delta * S2z) +
                   5 * Nu *
                       (PhiDOT * r *
                            (Complex(0, 224308) +
                             3 * (67379 + 79865 * delta) * S1z + 202137 * S2z -
                             239595 * delta * S2z) +
                        Complex(0, 4) * rDOT *
                            (Complex(0, 4175) + (61557 + 38469 * delta) * S1z +
                             61557 * S2z - 38469 * delta * S2z)) -
                   3 * params.eta2 *
                       (4 * PhiDOT * r *
                            (Complex(0, 108687) + 41225 * S1z +
                             29825 * delta * S1z + 41225 * S2z -
                             29825 * delta * S2z) +
                        Complex(0, 5) * rDOT *
                            (Complex(0, 19722) + 132241 * S1z +
                             18825 * delta * S1z + 132241 * S2z -
                             18825 * delta * S2z))) -
              3 * params.r2 *
                  (840 * params.eta3 *
                       (3159 * params.PhiDOT5 * params.r5 +
                        Complex(0, 18450) * params.PhiDOT4 * params.r4 * rDOT +
                        3678 * params.PhiDOT3 * params.r3 * params.rDOT2 +
                        Complex(0, 13970) * params.PhiDOT2 * params.r2 *
                            params.rDOT3 -
                        1546 * PhiDOT * r * params.rDOT4 +
                        Complex(0, 2100) * params.rDOT5) *
                       (S1z + S2z) -
                   525 *
                       (354 * params.PhiDOT5 * params.r5 +
                        Complex(0, 453) * params.PhiDOT4 * params.r4 * rDOT -
                        60 * params.PhiDOT3 * params.r3 * params.rDOT2 +
                        Complex(0, 1546) * params.PhiDOT2 * params.r2 *
                            params.rDOT3 -
                        336 * PhiDOT * r * params.rDOT4 +
                        Complex(0, 184) * params.rDOT5) *
                       (S1z + delta * S1z + S2z - delta * S2z) +
                   Nu * (Complex(0, 7) * params.PhiDOT4 * params.r4 * rDOT *
                             (Complex(0, 1062898) +
                              15 * (11275 + 23003 * delta) * S1z +
                              (169125 - 345045 * delta) * S2z) -
                         8 * PhiDOT * r * params.rDOT4 *
                             (Complex(0, -58753) +
                              420 * (145 + 6 * delta) * S1z -
                              420 * (-145 + 6 * delta) * S2z) +
                         Complex(0, 40) * params.rDOT5 *
                             (Complex(0, 12203) +
                              21 * (843 + 656 * delta) * S1z -
                              21 * (-843 + 656 * delta) * S2z) -
                         5 * params.PhiDOT5 * params.r5 *
                             (Complex(0, 240964) +
                              21 * (-20361 + 2747 * delta) * S1z -
                              21 * (20361 + 2747 * delta) * S2z) +
                         2 * params.PhiDOT3 * params.r3 * params.rDOT2 *
                             (Complex(0, 1832842) +
                              105 * (-13567 + 6603 * delta) * S1z -
                              105 * (13567 + 6603 * delta) * S2z) +
                         Complex(0, 4) * params.PhiDOT2 * params.r2 *
                             params.rDOT3 *
                             (Complex(0, 97999) +
                              105 * (9121 + 7144 * delta) * S1z -
                              105 * (-9121 + 7144 * delta) * S2z)) +
                   6 * params.eta2 *
                       (Complex(0, -10) * params.rDOT5 *
                            (Complex(0, 35291) +
                             14 * (2183 + 1155 * delta) * S1z +
                             (30562 - 16170 * delta) * S2z) -
                        Complex(0, 3) * params.PhiDOT2 * params.r2 *
                            params.rDOT3 *
                            (Complex(0, 120747) +
                             280 * (1571 + 165 * delta) * S1z -
                             280 * (-1571 + 165 * delta) * S2z) +
                        2 * PhiDOT * r * params.rDOT4 *
                            (Complex(0, -280067) -
                             140 * (-124 + 1293 * delta) * S1z +
                             140 * (124 + 1293 * delta) * S2z) +
                        5 * params.PhiDOT5 * params.r5 *
                            (Complex(0, 45393) +
                             14 * (-12899 + 1573 * delta) * S1z -
                             14 * (12899 + 1573 * delta) * S2z) -
                        Complex(0, 7) * params.PhiDOT4 * params.r4 * rDOT *
                            (Complex(0, 213613) +
                             5 * (31897 + 4025 * delta) * S1z -
                             5 * (-31897 + 4025 * delta) * S2z) +
                        params.PhiDOT3 * params.r3 * params.rDOT2 *
                            (Complex(0, -2032699) -
                             35 * (-33139 + 11613 * delta) * S1z +
                             35 * (33139 + 11613 * delta) * S2z))) +
              mass * r *
                  (2520 * params.eta3 *
                       (3772 * params.PhiDOT3 * params.r3 +
                        Complex(0, 22365) * params.PhiDOT2 * params.r2 * rDOT -
                        1427 * PhiDOT * r * params.rDOT2 +
                        Complex(0, 4532) * params.rDOT3) *
                       (S1z + S2z) +
                   420 *
                       (643 * params.PhiDOT3 * params.r3 -
                        Complex(0, 17143) * params.PhiDOT2 * params.r2 * rDOT +
                        4684 * PhiDOT * r * params.rDOT2 +
                        Complex(0, 377) * params.rDOT3) *
                       (S1z + delta * S1z + S2z - delta * S2z) -
                   6 * params.eta2 *
                       (Complex(0, 25) * params.rDOT3 *
                            (Complex(0, 73352) +
                             35 * (137 + 1053 * delta) * S1z +
                             (4795 - 36855 * delta) * S2z) +
                        params.PhiDOT3 * params.r3 *
                            (Complex(0, 3625463) +
                             70 * (45113 + 1515 * delta) * S1z -
                             70 * (-45113 + 1515 * delta) * S2z) +
                        4 * PhiDOT * r * params.rDOT2 *
                            (Complex(0, 576638) +
                             245 * (-884 + 2103 * delta) * S1z -
                             245 * (884 + 2103 * delta) * S2z) -
                        Complex(0, 1) * params.PhiDOT2 * params.r2 * rDOT *
                            (Complex(0, -3147986) +
                             35 * (-290963 + 34473 * delta) * S1z -
                             35 * (290963 + 34473 * delta) * S2z)) +
                   Nu * (Complex(0, 200) * params.rDOT3 *
                             (Complex(0, -2501) +
                              21 * (-308 + 51 * delta) * S1z -
                              21 * (308 + 51 * delta) * S2z) +
                         5 * params.PhiDOT3 * params.r3 *
                             (Complex(0, 1501652) -
                              21 * (-37535 + 8031 * delta) * S1z +
                              21 * (37535 + 8031 * delta) * S2z) -
                         2 * PhiDOT * r * params.rDOT2 *
                             (Complex(0, -975992) +
                              105 * (33593 + 24921 * delta) * S1z -
                              105 * (-33593 + 24921 * delta) * S2z) +
                         Complex(0, 4) * params.PhiDOT2 * params.r2 * rDOT *
                             (Complex(0, 4843364) +
                              105 * (88681 + 40215 * delta) * S1z -
                              105 * (-88681 + 40215 * delta) * S2z))))) /
            (1.16424e6 * sqrt(5) * params.r4)));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hQC_4_m_2(REAL8 Nu, UINT4 vpnorder, REAL8 x, struct kepler_vars params) {

  /* if(vpnorder == 2){
       return(((2*sqrt(5))/63. - (2*sqrt(5)*Nu)/21.)*params.x2);
   }

   else if(vpnorder == 4){
       return((-437/(693.*sqrt(5)) + (115*sqrt(5)*Nu)/297. -
       (19*sqrt(5)*params.eta2)/693.)*params.x3);
   } */

  /* else */ if (vpnorder == 5) {
    return (M_PI * ((4 * sqrt(5) * params.x3p5) / 63. -
                    (4 * sqrt(5) * Nu * params.x3p5) / 21.));
  }
  /* else if(vpnorder == 6){
       return((346013/(420420.*sqrt(5)) - (606751*Nu)/(180180.*sqrt(5)) +
       (400453*params.eta2)/(162162.*sqrt(5)) +
       (25783*params.eta3)/(108108.*sqrt(5)))*params.x4);
   } */
  else {
    return 0;
  }
}

static COMPLEX16 hl_4_m_2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z, REAL8 x, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_4_m_2: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_4_m_2(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, params) +
            hQC_4_m_2(Nu, vpnorder, x, params)) *
           cpolar(1, -2 * Phi);
  }
}

static COMPLEX16 hl_4_m_min2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z, REAL8 x, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_4_m_min2: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           conj(hGO_4_m_2(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, params) +
                hQC_4_m_2(Nu, vpnorder, x, params)) *
           cpolar(1, 2 * Phi);
  }
}

// H41
static COMPLEX16 hGO_4_m_1(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z,
                           REAL8 x, struct kepler_vars params) {
  REAL8 delta = sqrt(1 - 4 * Nu);
  REAL8 kappa1 = 1.0;
  REAL8 kappa2 = 1.0;

  if (vpnorder == 3) {
    return (delta * mass * (-1 + 2 * Nu) * PhiDOT *
            (Complex(0, -12) * mass +
             r * (Complex(0, 11) * params.PhiDOT2 * params.r2 +
                  10 * PhiDOT * r * rDOT + Complex(0, 6) * params.rDOT2))) /
           (42. * sqrt(10) * r);
  }

  else if (vpnorder == 4) {
    return ((Complex(0, 0.005952380952380952) * sqrt(2.5) * params.Mtot2 * Nu *
             (12 * mass +
              r * (-11 * params.PhiDOT2 * params.r2 +
                   Complex(0, 10) * PhiDOT * r * rDOT - 6 * params.rDOT2)) *
             ((-1 + delta) * S1z + S2z + delta * S2z)) /
                params.r3 -
            (Complex(0, 0.005952380952380952) * sqrt(2.5) * Nu *
             (-S1z + S2z + delta * (S1z + S2z)) * params.x3) +
            (Complex(0, 0.005952380952380952) * sqrt(2.5) * Nu *
             (-S1z + S2z + delta * (S1z + S2z)) * params.x3));
  }

  else if (vpnorder == 5) {
    return (Complex(0, -0.00018037518037518038) * delta * mass * PhiDOT *
            (6 * params.Mtot2 * (972 - 2293 * Nu + 1398 * params.eta2) -
             2 * mass * r *
                 ((-340 - 785 * Nu + 6504 * params.eta2) * params.PhiDOT2 *
                      params.r2 -
                  Complex(0, 3) * (-2796 + 5299 * Nu + 1622 * params.eta2) *
                      PhiDOT * r * rDOT +
                  6 * (-1200 + 2545 * Nu + 162 * params.eta2) * params.rDOT2) +
             3 * params.r2 *
                 ((-540 + 235 * Nu + 2648 * params.eta2) * params.PhiDOT4 *
                      params.r4 -
                  Complex(0, 4) * (-1764 + 3536 * Nu + 373 * params.eta2) *
                      params.PhiDOT3 * params.r3 * rDOT +
                  2 * (-723 + 1022 * Nu + 1384 * params.eta2) * params.PhiDOT2 *
                      params.r2 * params.rDOT2 -
                  Complex(0, 4) * (-229 + 366 * Nu + 358 * params.eta2) *
                      PhiDOT * r * params.rDOT3 +
                  12 * (-43 + 62 * Nu + 80 * params.eta2) * params.rDOT4))) /
           (sqrt(10) * params.r2);
  }

  else if (vpnorder == 6) {
    return (
        -(delta * params.Mtot2 * Nu * PhiDOT *
          (mass * (362 * PhiDOT * r - Complex(0, 534) * rDOT) +
           r * (149 * params.PhiDOT3 * params.r3 +
                Complex(0, 182) * params.PhiDOT2 * params.r2 * rDOT -
                136 * PhiDOT * r * params.rDOT2 +
                Complex(0, 112) * params.rDOT3))) /
            (420. * sqrt(10) * params.r2)
        /* Henry et al. QC spin terms */
        /* +((Complex(0,-0.00018037518037518038)*(220*(S1z - S2z + delta*(S1z +
  S2z))
  + Nu*(-2247*(S1z - S2z) + 47*delta*(S1z + S2z)) +
  params.eta2*(2891*(S1z - S2z) + 613*delta*(S1z + S2z)))*params.x4)/sqrt(10)) */
        + /*Henry et al. ecc spin terms*/ (
              (Complex(0, 0.00027056277056277056) * params.Mtot2 *
               (2 * params.Mtot2 *
                    ((-440 + 6801 * Nu - 1428 * params.eta2 +
                      delta * (-440 - 3193 * Nu + 300 * params.eta2)) *
                         S1z +
                     (440 - 6801 * Nu + 1428 * params.eta2 +
                      delta * (-440 - 3193 * Nu + 300 * params.eta2)) *
                         S2z) +
                2 * mass * r *
                    (Complex(0, 1) * PhiDOT * r * rDOT *
                         (delta * (-1320 + 9093 * Nu + 59 * params.eta2) * S1z -
                          5 * (264 - 311 * Nu + 823 * params.eta2) * S1z +
                          delta * (-1320 + 9093 * Nu + 59 * params.eta2) * S2z +
                          5 * (264 - 311 * Nu + 823 * params.eta2) * S2z) +
                     2 * params.rDOT2 *
                         ((220 + 1659 * Nu - 1512 * params.eta2 +
                           delta * (220 - 3067 * Nu + 240 * params.eta2)) *
                              S1z +
                          (-220 - 1659 * Nu + 1512 * params.eta2 +
                           delta * (220 - 3067 * Nu + 240 * params.eta2)) *
                              S2z) +
                     2 * params.PhiDOT2 * params.r2 *
                         ((-66 + 158 * Nu - 3389 * params.eta2 +
                           delta * (-66 + 238 * Nu + 301 * params.eta2)) *
                              S1z +
                          (66 - 158 * Nu + 3389 * params.eta2 +
                           delta * (-66 + 238 * Nu + 301 * params.eta2)) *
                              S2z)) +
                params.r2 *
                    (Complex(0, 3960) * params.PhiDOT3 * params.r3 * rDOT *
                         ((1 + delta) * S1z + (-1 + delta) * S2z) +
                     2 * params.eta2 *
                         (6 * params.rDOT4 *
                              (149 * S1z + 75 * delta * S1z - 149 * S2z +
                               75 * delta * S2z) +
                          7 * params.PhiDOT2 * params.r2 * params.rDOT2 *
                              (331 * S1z + 293 * delta * S1z - 331 * S2z +
                               293 * delta * S2z) -
                          Complex(0, 2) * PhiDOT * r * params.rDOT3 *
                              (925 * S1z + 491 * delta * S1z - 925 * S2z +
                               491 * delta * S2z) +
                          params.PhiDOT4 * params.r4 *
                              (2593 * S1z + 551 * delta * S1z - 2593 * S2z +
                               551 * delta * S2z) +
                          Complex(0, 1) * params.PhiDOT3 * params.r3 * rDOT *
                              (8815 * S1z + 2393 * delta * S1z - 8815 * S2z +
                               2393 * delta * S2z)) +
                     Nu * (6 * params.rDOT4 *
                               (-301 * S1z + 213 * delta * S1z + 301 * S2z +
                                213 * delta * S2z) +
                           8 * params.PhiDOT4 * params.r4 *
                               (-349 * S1z + 305 * delta * S1z + 349 * S2z +
                                305 * delta * S2z) -
                           Complex(0, 2) * PhiDOT * r * params.rDOT3 *
                               (-1625 * S1z + 1009 * delta * S1z + 1625 * S2z +
                                1009 * delta * S2z) +
                           params.PhiDOT2 * params.r2 * params.rDOT2 *
                               (-4973 * S1z + 2773 * delta * S1z + 4973 * S2z +
                                2773 * delta * S2z) -
                           Complex(0, 2) * params.PhiDOT3 * params.r3 * rDOT *
                               (1765 * S1z + 14779 * delta * S1z - 1765 * S2z +
                                14779 * delta * S2z))))) /
              (sqrt(10) * params.r4)));
  }

  else if (vpnorder == 7) {

    return (
        /* Henry et al. QC spiun terms */ /* (Complex(0,0.000992063492063492)*params.x4p5*(6*(10*(4
 + delta)*params.eta2
   + kappa1*(5 + 5*delta*pow(1 - 2*Nu,2) + 6*Nu*(-5 + 2*Nu)))*params.S1z2 +
 S2z*(Complex(0,11)*(-1 + delta) + Nu*(Complex(0,-31) + 30*M_PI) + 6*kappa2*(-5
 + 5*delta*pow(1 - 2*Nu,2)
 + 6*(5 - 2*Nu)*Nu)*S2z + 15*Nu*(-4*(4*Nu*S2z + Complex(0,1)*log(2)) +
 delta*(2*M_PI
 + 4*Nu*S2z - Complex(0,1)*(5 + log(16))))) + S1z*(Complex(0,11) +
 Nu*(Complex(0,31)
 - 30*M_PI + Complex(0,60)*log(2)) + delta*(Complex(0,11) - 15*Nu*(-2*M_PI +
 8*Nu*S2z + Complex(0,1)*(5 + log(16)))))))/sqrt(10)
 + */
                                          /* Henry et al. ecc+spin terms */
        ((Complex(0, 1.9821448392876965e-7) * mass *
          (-(delta *
             (14 * mass * PhiDOT * params.r3 *
                  (2 *
                       (243784 - 839354 * Nu + 840205 * params.eta2 +
                        1464600 * params.eta3) *
                       params.PhiDOT4 * params.r4 -
                   Complex(0, 3) *
                       (-14468 - 3185682 * Nu + 8381955 * params.eta2 +
                        30000 * params.eta3) *
                       params.PhiDOT3 * params.r3 * rDOT +
                   5 *
                       (-657932 + 496144 * Nu + 2150479 * params.eta2 +
                        542628 * params.eta3) *
                       params.PhiDOT2 * params.r2 * params.rDOT2 -
                   Complex(0, 3) *
                       (-668008 + 803028 * Nu + 1908955 * params.eta2 +
                        540370 * params.eta3) *
                       PhiDOT * r * params.rDOT3 +
                   6 *
                       (-172143 + 155683 * Nu + 680580 * params.eta2 +
                        111840 * params.eta3) *
                       params.rDOT4) -
              105 * PhiDOT * params.r4 *
                  ((1320 - 35155 * Nu + 36743 * params.eta2 +
                    241538 * params.eta3) *
                       params.PhiDOT6 * params.r6 +
                   Complex(0, 26) *
                       (-7536 + 50167 * Nu - 84782 * params.eta2 +
                        22 * params.eta3) *
                       params.PhiDOT5 * params.r5 * rDOT +
                   2 *
                       (46908 - 203266 * Nu + 59723 * params.eta2 +
                        272204 * params.eta3) *
                       params.PhiDOT4 * params.r4 * params.rDOT2 -
                   Complex(0, 4) *
                       (38798 - 217656 * Nu + 324761 * params.eta2 +
                        67186 * params.eta3) *
                       params.PhiDOT3 * params.r3 * params.rDOT3 +
                   (58177 - 260454 * Nu + 239688 * params.eta2 +
                    324160 * params.eta3) *
                       params.PhiDOT2 * params.r2 * params.rDOT4 -
                   Complex(0, 2) *
                       (11175 - 52698 * Nu + 57436 * params.eta2 +
                        60808 * params.eta3) *
                       PhiDOT * r * params.rDOT5 +
                   18 *
                       (829 - 3726 * Nu + 3480 * params.eta2 +
                        4640 * params.eta3) *
                       params.rDOT6) +
              26 * params.Mtot3 *
                  (2310 * rDOT *
                       (Complex(0, -30) * params.eta2 * params.S1z2 +
                        S2z * (-22 + 67 * Nu -
                               Complex(0, 30) * params.eta2 * S2z) +
                        S1z * (-22 + Nu * (67 + Complex(0, 100) * S2z) -
                               Complex(0, 60) * params.eta2 * S2z)) +
                   PhiDOT * r *
                       (-114228 + 524160 * params.eta3 -
                        Complex(0, 46585) * S1z - 62370 * kappa1 * params.S1z2 -
                        Complex(0, 46585) * S2z - 62370 * kappa2 * params.S2z2 -
                        420 * params.eta2 *
                            (1422 + 11 * (-25 + 54 * kappa1) * params.S1z2 -
                             1738 * S1z * S2z +
                             11 * (-25 + 54 * kappa2) * params.S2z2) +
                        7 * Nu *
                            (18564 + 35640 * kappa1 * params.S1z2 +
                             Complex(0, 21505) * S2z +
                             35640 * kappa2 * params.S2z2 +
                             55 * S1z * (Complex(0, 391) + 312 * S2z)))) +
              params.Mtot2 * r *
                  (60060 * params.rDOT3 * (S1z + S2z) *
                       (11 - 71 * Nu -
                        Complex(0, 10) * params.eta2 * (S1z + S2z)) +
                   156 * PhiDOT * r * params.rDOT2 *
                       (98013 + 129150 * params.eta3 - Complex(0, 8470) * S1z +
                        7 * Nu *
                            (31961 + Complex(0, 2860) * S1z +
                             Complex(0, 2860) * S2z) -
                        Complex(0, 8470) * S2z +
                        35 * params.eta2 *
                            (-33313 + 220 * params.S1z2 + 440 * S1z * S2z +
                             220 * params.S2z2)) -
                   5 * params.PhiDOT3 * params.r3 *
                       (3966680 + 5749576 * params.eta3 -
                        Complex(0, 231231) * S1z -
                        90090 * kappa1 * params.S1z2 -
                        Complex(0, 231231) * S2z -
                        90090 * kappa2 * params.S2z2 +
                        7 * Nu *
                            (-921242 + Complex(0, 110253) * S1z +
                             51480 * kappa1 * params.S1z2 +
                             Complex(0, 110253) * S2z +
                             51480 * kappa2 * params.S2z2) -
                        462 * params.eta2 *
                            (10483 + 130 * (-5 + 6 * kappa1) * params.S1z2 -
                             2860 * S1z * S2z +
                             130 * (-5 + 6 * kappa2) * params.S2z2)) +
                   Complex(0, 2) * params.PhiDOT2 * params.r2 * rDOT *
                       (6531280 * params.eta3 +
                        3 * (5382288 - Complex(0, 385385) * S1z +
                             150150 * kappa1 * params.S1z2 -
                             Complex(0, 385385) * S2z +
                             150150 * kappa2 * params.S2z2) +
                        210 * params.eta2 *
                            (322144 + 715 * (19 + 12 * kappa1) * params.S1z2 +
                             10010 * S1z * S2z +
                             715 * (19 + 12 * kappa2) * params.S2z2) -
                        21 * Nu *
                            (2826484 + 85800 * kappa1 * params.S1z2 -
                             Complex(0, 140855) * S2z +
                             85800 * kappa2 * params.S2z2 +
                             715 * S1z * (Complex(0, -197) + 540 * S2z)))))) +
           5005 * params.Mtot2 *
               (-3 * r *
                    (-8 * PhiDOT * r * params.rDOT2 *
                         (Complex(0, -1) * (-11 + 48 * Nu) * S1z +
                          2 * (-5 + 16 * kappa1) * params.eta2 * params.S1z2 +
                          S2z * (Complex(0, -11) + Complex(0, 48) * Nu +
                                 2 * (5 - 16 * kappa2) * params.eta2 * S2z)) +
                     4 * params.rDOT3 *
                         ((11 - 93 * Nu) * S1z +
                          Complex(0, 2) * (-5 + 4 * kappa1) * params.eta2 *
                              params.S1z2 +
                          S2z * (-11 + 93 * Nu -
                                 Complex(0, 2) * (-5 + 4 * kappa2) *
                                     params.eta2 * S2z)) +
                     Complex(0, 2) * params.PhiDOT2 * params.r2 * rDOT *
                         (Complex(0, 1) * (-77 + 351 * Nu) * S1z +
                          2 * (245 * params.eta2 + kappa1 * (15 - 90 * Nu + 68 * params.eta2)) *
                              params.S1z2 -
                          S2z * (Complex(0, -77) + 30 * kappa2 * S2z +
                                 2 * (245 + 68 * kappa2) * params.eta2 * S2z -
                                 9 * Nu *
                                     (Complex(0, -39) + 20 * kappa2 * S2z))) +
                     params.PhiDOT3 * params.r3 *
                         (Complex(0, -1) * (-77 + 411 * Nu) * S1z +
                          2 * (100 * params.eta2 + kappa1 * (15 - 90 * Nu + 172 * params.eta2)) *
                              params.S1z2 -
                          S2z * (Complex(0, 77) + 30 * kappa2 * S2z +
                                 8 * (25 + 43 * kappa2) * params.eta2 * S2z -
                                 3 * Nu *
                                     (Complex(0, 137) + 60 * kappa2 * S2z)))) +
                2 * mass *
                    (6 * rDOT *
                         ((22 - 111 * Nu) * S1z +
                          Complex(0, 2) * (15 + 8 * kappa1) * params.eta2 *
                              params.S1z2 +
                          S2z * (-22 + 111 * Nu -
                                 Complex(0, 2) * (15 + 8 * kappa2) *
                                     params.eta2 * S2z)) +
                     PhiDOT * r *
                         (Complex(0, -1) * (-121 + 633 * Nu) * S1z +
                          6 * (220 * params.eta2 + 3 * kappa1 * (9 - 54 * Nu + 76 * params.eta2)) *
                              params.S1z2 -
                          S2z * (Complex(0, 121) + 162 * kappa2 * S2z +
                                 24 * (55 + 57 * kappa2) * params.eta2 * S2z -
                                 3 * Nu *
                                     (Complex(0, 211) +
                                      324 * kappa2 * S2z))))))) /
         (sqrt(10) * params.r4)) +
        /* Henry et al. QC spinning hereditary terms */ (
            (sqrt(2.5) * Nu * ((-1 + delta) * S1z + S2z + delta * S2z) *
             params.x4p5 * (Complex(0, 1) * M_PI + log(4))) /
            168.));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hQC_4_m_1(REAL8 Nu, UINT4 vpnorder, REAL8 x, struct kepler_vars params) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  /* if(vpnorder == 3){
       return(((Complex(0,0.023809523809523808)*delta)/sqrt(10) -
       (Complex(0,0.047619047619047616)*delta*Nu)/sqrt(10))*params.x2p5);
   }
   else if(vpnorder == 5){
       return(((Complex(0,-0.07287157287157287)*delta)/sqrt(10) +
       (Complex(0,0.18235930735930736)*delta*Nu)/sqrt(10) -
       (Complex(0,0.05988455988455989)*delta*params.eta2)/sqrt(10))*
       params.x3p5);
   } */
  /* else */ if (vpnorder == 6) {
    return (
        M_PI *
            ((Complex(0, 0.023809523809523808) * delta * params.x4) / sqrt(10) -
             (Complex(0, 0.047619047619047616) * delta * Nu * params.x4) /
                 sqrt(10)) +
        ((delta * params.x4) / (21. * sqrt(10)) -
         (sqrt(0.4) * delta * Nu * params.x4) / 21.) *
            log(2));
  } else {
    return 0;
  }
}

static COMPLEX16 hl_4_m_1(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z, REAL8 x, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_4_m_1: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_4_m_1(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x, params) +
            hQC_4_m_1(Nu, vpnorder, x, params)) *
           cpolar(1, -1 * Phi);
  }
}

static COMPLEX16 hl_4_m_min1(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z, REAL8 x, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_4_m_min1: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           conj(hGO_4_m_1(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x, params) +
                hQC_4_m_1(Nu, vpnorder, x, params)) *
           cpolar(1, 1 * Phi);
  }
}

// H55
static COMPLEX16 hGO_5_m_5(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z, struct kepler_vars params) {
  REAL8 delta = sqrt(1 - 4 * Nu);
  REAL8 kappa1 = 1.0;
  REAL8 kappa2 = 1.0;

  if (vpnorder == 3) {
    return ((delta * (-1 + 2 * Nu) *
             (24 * params.r2 * pow(Complex(0, -1) * PhiDOT * r + rDOT, 5) +
              2 * params.Mtot2 * (Complex(0, -86) * PhiDOT * r + 41 * rDOT) +
              3 * mass * r *
                  (Complex(0, -143) * params.PhiDOT3 * params.r3 +
                   208 * params.PhiDOT2 * params.r2 * rDOT +
                   Complex(0, 132) * PhiDOT * r * params.rDOT2 -
                   32 * params.rDOT3))) /
            (48. * sqrt(66) * params.r2));
  }

  else if (vpnorder == 5) {
    return (
        (delta *
         (360 * (33 - 197 * Nu + 294 * params.eta2) * params.r3 *
              pow(PhiDOT * r + Complex(0, 1) * rDOT, 6) *
              (Complex(0, 1) * PhiDOT * r + rDOT) +
          2 * params.Mtot3 *
              (Complex(0, -5) * (53311 - 121906 * Nu + 42816 * params.eta2) *
                   PhiDOT * r +
               78 * (1141 - 2760 * Nu + 1420 * params.eta2) * rDOT) +
          2 * params.Mtot2 * r *
              (Complex(0, -10) * (40826 - 125981 * Nu + 87534 * params.eta2) *
                   params.PhiDOT3 * params.r3 +
               (112966 - 818425 * Nu + 1385970 * params.eta2) * params.PhiDOT2 *
                   params.r2 * rDOT +
               Complex(0, 20) * (-2636 - 11335 * Nu + 43962 * params.eta2) *
                   PhiDOT * r * params.rDOT2 -
               39 * (-639 - 735 * Nu + 5690 * params.eta2) * params.rDOT3) +
          15 * mass * params.r2 *
              (Complex(0, 4) * (139 - 1687 * Nu + 11376 * params.eta2) *
                   params.PhiDOT5 * params.r5 -
               2 * (11276 - 19559 * Nu + 5982 * params.eta2) * params.PhiDOT4 *
                   params.r4 * rDOT +
               Complex(0, 3) * (-14615 + 12440 * Nu + 37132 * params.eta2) *
                   params.PhiDOT3 * params.r3 * params.rDOT2 -
               8 * (-4666 + 139 * Nu + 22194 * params.eta2) * params.PhiDOT2 *
                   params.r2 * params.rDOT3 -
               Complex(0, 4) * (-3971 - 3226 * Nu + 27360 * params.eta2) *
                   PhiDOT * r * params.rDOT4 +
               48 * (-57 - 97 * Nu + 518 * params.eta2) * params.rDOT5))) /
        (18720. * sqrt(66) * params.r3));
  }

  else if (vpnorder == 6) {
    return (
        -(delta * params.Mtot2 * Nu *
          (3566 * params.Mtot2 +
           6 * mass * r *
               (11305 * params.PhiDOT2 * params.r2 +
                Complex(0, 3921) * PhiDOT * r * rDOT - 906 * params.rDOT2) +
           params.r2 * (104681 * params.PhiDOT4 * params.r4 +
                        Complex(0, 17192) * params.PhiDOT3 * params.r3 * rDOT -
                        27840 * params.PhiDOT2 * params.r2 * params.rDOT2 -
                        Complex(0, 9968) * PhiDOT * r * params.rDOT3 +
                        1424 * params.rDOT4))) /
            (1680. * sqrt(66) * params.r4)
        /* +((Complex(0,-21.70138888888889)*((1 + delta + 3*delta*(-1 + Nu)*Nu +
  Nu*(-7 + 11*Nu))*S1z
  + (-1 + delta + (7 - 11*Nu)*Nu + 3*delta*(-1 +
  Nu)*Nu)*S2z)*params.x4)/sqrt(66)) */
        +
        /* Henry et al. ecc spin terms */ (
            (Complex(0, -0.006944444444444444) * params.Mtot2 *
             (params.Mtot2 * ((90 - 622 * Nu + 966 * params.eta2 +
                               delta * (90 - 524 * Nu + 770 * params.eta2)) *
                                  S1z +
                              2 *
                                  (-45 + 311 * Nu - 483 * params.eta2 +
                                   delta * (45 - 262 * Nu + 385 * params.eta2)) *
                                  S2z) +
              3 * params.r2 *
                  (4 * params.rDOT4 *
                       ((15 - 92 * Nu + 126 * params.eta2 +
                         delta * (15 - 64 * Nu + 70 * params.eta2)) *
                            S1z +
                        (-15 + 92 * Nu - 126 * params.eta2 +
                         delta * (15 - 64 * Nu + 70 * params.eta2)) *
                            S2z) -
                   15 * params.PhiDOT2 * params.r2 * params.rDOT2 *
                       ((41 - 258 * Nu + 364 * params.eta2 +
                         delta * (41 - 154 * Nu + 156 * params.eta2)) *
                            S1z +
                        (-41 + 258 * Nu - 364 * params.eta2 +
                         delta * (41 - 154 * Nu + 156 * params.eta2)) *
                            S2z) -
                   Complex(0, 4) * PhiDOT * r * params.rDOT3 *
                       ((75 - 464 * Nu + 642 * params.eta2 +
                         delta * (75 - 298 * Nu + 310 * params.eta2)) *
                            S1z +
                        (-75 + 464 * Nu - 642 * params.eta2 +
                         delta * (75 - 298 * Nu + 310 * params.eta2)) *
                            S2z) +
                   params.PhiDOT4 * params.r4 *
                       ((300 - 2347 * Nu + 4041 * params.eta2 +
                         delta * (300 - 869 * Nu + 1085 * params.eta2)) *
                            S1z +
                        (-300 + 2347 * Nu - 4041 * params.eta2 +
                         delta * (300 - 869 * Nu + 1085 * params.eta2)) *
                            S2z) +
                   Complex(0, 1) * params.PhiDOT3 * params.r3 * rDOT *
                       ((675 - 4414 * Nu + 6492 * params.eta2 +
                         delta * (675 - 2558 * Nu + 2780 * params.eta2)) *
                            S1z +
                        (-675 + 4414 * Nu - 6492 * params.eta2 +
                         delta * (675 - 2558 * Nu + 2780 * params.eta2)) *
                            S2z)) +
              mass * r *
                  (-2 * params.rDOT2 *
                       ((255 - 1592 * Nu + 2226 * params.eta2 +
                         delta * (255 - 1144 * Nu + 1330 * params.eta2)) *
                            S1z +
                        (-255 + 1592 * Nu - 2226 * params.eta2 +
                         delta * (255 - 1144 * Nu + 1330 * params.eta2)) *
                            S2z) +
                   Complex(0, 2) * PhiDOT * r * rDOT *
                       ((870 - 5563 * Nu + 7989 * params.eta2 +
                         delta * (870 - 3743 * Nu + 4349 * params.eta2)) *
                            S1z +
                        (-870 + 5563 * Nu - 7989 * params.eta2 +
                         delta * (870 - 3743 * Nu + 4349 * params.eta2)) *
                            S2z) +
                   params.PhiDOT2 * params.r2 *
                       ((1329 - 9376 * Nu + 14838 * params.eta2 +
                         delta * (1329 - 5438 * Nu + 6962 * params.eta2)) *
                            S1z +
                        (-1329 + 9376 * Nu - 14838 * params.eta2 +
                         delta * (1329 - 5438 * Nu + 6962 * params.eta2)) *
                            S2z)))) /
            (sqrt(66) * params.r4)));
  }

  else if (vpnorder == 7) {

    return (/* (Complex(0,16.276041666666668)*(-1 + 2*Nu)*(kappa1*(-1 - delta +
      2*(2 + delta)*Nu)*params.S1z2
      + S2z*(-4*delta*Nu*S1z + kappa2*(1 - delta + 2*(-2 +
      delta)*Nu)*S2z))*params.x4p5)/sqrt(66) */
            /*  + */ /* Henry et al. ecc+spin terms */ (
                (Complex(0, 24570) * params.Mtot3 *
                     (1 - 6 * Nu + 8 * params.eta2) *
                     (mass * (298 * PhiDOT * r + Complex(0, 94) * rDOT) +
                      5 * r *
                          (95 * params.PhiDOT3 * params.r3 +
                           Complex(0, 58) * params.PhiDOT2 * params.r2 * rDOT -
                           34 * PhiDOT * r * params.rDOT2 -
                           Complex(0, 8) * params.rDOT3)) *
                     (kappa1 * params.S1z2 - kappa2 * params.S2z2) +
                 delta *
                     (7560 *
                          (-135 + 1467 * Nu - 5285 * params.eta2 +
                           6230 * params.eta3) *
                          params.r4 *
                          pow(Complex(0, 1) * PhiDOT * r - rDOT, 7) *
                          pow(PhiDOT * r - Complex(0, 1) * rDOT, 2) +
                      6 * params.Mtot2 * params.r2 *
                          (Complex(0, 2) *
                               (-2032864 - 4559381 * Nu - 3520590 * params.eta2 +
                                36950400 * params.eta3) *
                               params.PhiDOT5 * params.r5 -
                           8 *
                               (-1581827 + 6211490 * Nu + 607860 * params.eta2 +
                                6123006 * params.eta3) *
                               params.PhiDOT4 * params.r4 * rDOT +
                           Complex(0, 12) *
                               (738669 - 8183419 * Nu + 7682580 * params.eta2 +
                                9432160 * params.eta3) *
                               params.PhiDOT3 * params.r3 * params.rDOT2 -
                           2 *
                               (323812 - 40147318 * Nu + 43264365 * params.eta2 +
                                116543718 * params.eta3) *
                               params.PhiDOT2 * params.r2 * params.rDOT3 -
                           Complex(0, 14) *
                               (-19909 - 2880941 * Nu + 2929365 * params.eta2 +
                                11172930 * params.eta3) *
                               PhiDOT * r * params.rDOT4 +
                           21 *
                               (1285 - 417939 * Nu + 376375 * params.eta2 +
                                1794742 * params.eta3) *
                               params.rDOT5) +
                      45 * mass * params.r3 *
                          (Complex(0, -4) *
                               (-76218 + 328011 * Nu - 511853 * params.eta2 +
                                804238 * params.eta3) *
                               params.PhiDOT7 * params.r7 -
                           4 *
                               (241824 - 1017903 * Nu + 1022762 * params.eta2 +
                                918896 * params.eta3) *
                               params.PhiDOT6 * params.r6 * rDOT -
                           Complex(0, 4) *
                               (450537 - 1503933 * Nu - 642815 * params.eta2 +
                                6184642 * params.eta3) *
                               params.PhiDOT5 * params.r5 * params.rDOT2 +
                           8 *
                               (53226 + 197748 * Nu - 1673279 * params.eta2 +
                                3045358 * params.eta3) *
                               params.PhiDOT4 * params.r4 * params.rDOT3 -
                           Complex(0, 1) *
                               (1855503 - 9943746 * Nu + 10573448 * params.eta2 +
                                6270464 * params.eta3) *
                               params.PhiDOT3 * params.r3 * params.rDOT4 +
                           16 *
                               (147339 - 604776 * Nu + 17530 * params.eta2 +
                                1666996 * params.eta3) *
                               params.PhiDOT2 * params.r2 * params.rDOT5 +
                           Complex(0, 28) *
                               (41703 - 148362 * Nu - 104120 * params.eta2 +
                                627040 * params.eta3) *
                               PhiDOT * r * params.rDOT6 -
                           1344 *
                               (162 - 516 * Nu - 707 * params.eta2 +
                                2870 * params.eta3) *
                               params.rDOT7) +
                      28 * params.Mtot4 *
                          (-3 * rDOT *
                               (300238 - 535772 * params.eta3 +
                                27495 * kappa1 * params.S1z2 +
                                27495 * kappa2 * params.S2z2 +
                                180 * params.eta2 *
                                    (3335 + 611 * kappa1 * params.S1z2 -
                                     1222 * S1z * S2z +
                                     611 * kappa2 * params.S2z2) -
                                15 * Nu *
                                    (38797 + 7332 * kappa1 * params.S1z2 -
                                     7332 * S1z * S2z +
                                     7332 * kappa2 * params.S2z2)) -
                           Complex(0, 1) * PhiDOT * r *
                               (-5713898 + 2860180 * params.eta3 -
                                261495 * kappa1 * params.S1z2 -
                                261495 * kappa2 * params.S2z2 -
                                30 * params.eta2 *
                                    (278563 + 34866 * kappa1 * params.S1z2 -
                                     69732 * S1z * S2z +
                                     34866 * kappa2 * params.S2z2) +
                                12 * Nu *
                                    (1137419 + 87165 * kappa1 * params.S1z2 -
                                     87165 * S1z * S2z +
                                     87165 * kappa2 * params.S2z2))) +
                      6 * params.Mtot3 * r *
                          (28 * params.rDOT3 *
                               (-175718 - 773906 * params.eta3 +
                                5850 * kappa1 * params.S1z2 +
                                5850 * kappa2 * params.S2z2 +
                                45 * params.eta2 *
                                    (-11217 + 520 * kappa1 * params.S1z2 -
                                     1040 * S1z * S2z +
                                     520 * kappa2 * params.S2z2) -
                                3 * Nu *
                                    (-243517 + 7800 * kappa1 * params.S1z2 -
                                     7800 * S1z * S2z +
                                     7800 * kappa2 * params.S2z2)) +
                           Complex(0, 14) * PhiDOT * r * params.rDOT2 *
                               (1974194 + 5947390 * params.eta3 -
                                49725 * kappa1 * params.S1z2 -
                                49725 * kappa2 * params.S2z2 -
                                1105 * params.eta2 *
                                    (-307 + 180 * kappa1 * params.S1z2 -
                                     360 * S1z * S2z +
                                     180 * kappa2 * params.S2z2) +
                                Nu * (-5423729 + 198900 * kappa1 * params.S1z2 -
                                      198900 * S1z * S2z +
                                      198900 * kappa2 * params.S2z2)) +
                           2 * params.PhiDOT2 * params.r2 * rDOT *
                               (22105852 + 67721336 * params.eta3 -
                                593775 * kappa1 * params.S1z2 -
                                593775 * kappa2 * params.S2z2 -
                                60 * params.eta2 *
                                    (1466303 + 39585 * kappa1 * params.S1z2 -
                                     79170 * S1z * S2z +
                                     39585 * kappa2 * params.S2z2) +
                                4 * Nu *
                                    (-4950311 + 593775 * kappa1 * params.S1z2 -
                                     593775 * S1z * S2z +
                                     593775 * kappa2 * params.S2z2)) -
                           Complex(0, 1) * params.PhiDOT3 * params.r3 *
                               (82632040 * params.eta3 -
                                3 * (4583428 + 648375 * kappa1 * params.S1z2 +
                                     648375 * kappa2 * params.S2z2) -
                                10 * params.eta2 *
                                    (24739703 + 778050 * kappa1 * params.S1z2 -
                                     1556100 * S1z * S2z +
                                     778050 * kappa2 * params.S2z2) +
                                9 * Nu *
                                    (14885821 + 864500 * kappa1 * params.S1z2 -
                                     864500 * S1z * S2z +
                                     864500 * kappa2 * params.S2z2))))) /
                (1.57248e6 * sqrt(66) * params.r4)));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hQC_5_m_5(REAL8 Nu, UINT4 vpnorder, REAL8 x, struct kepler_vars params) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  /* if(vpnorder == 3){
       return(((Complex(0,13.020833333333334)*delta)/sqrt(66) -
       (Complex(0,26.041666666666668)*delta*Nu)/sqrt(66))*params.x2p5);
   }
   else if(vpnorder == 5){
       return(((Complex(0,-87.80715811965813)*delta)/sqrt(66) +
       (Complex(0,229.7008547008547)*delta*Nu)/sqrt(66) -
       Complex(0,42.73504273504273)*sqrt(0.06060606060606061)*delta*
       params.eta2)*params.x3p5);
   } */
  /* else */ if (vpnorder == 6) {
    return ((3125 * delta * (-1 + 2 * Nu) * params.x4 *
             (Complex(0, -1) * M_PI + log(6.25))) /
            (48. * sqrt(66)));
  } else {
    return 0;
  }
}

static COMPLEX16 hl_5_m_5(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z, REAL8 x, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_5_m_5: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_5_m_5(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, params) +
            hQC_5_m_5(Nu, vpnorder, x, params)) *
           cpolar(1, -5 * Phi);
  }
}

static COMPLEX16 hl_5_m_min5(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z, REAL8 x, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_5_m_min5: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((- 4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           conj(hGO_5_m_5(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, params) +
                hQC_5_m_5(Nu, vpnorder, x, params)) *
           cpolar(1, 5 * Phi);
  }
}

// H54
static COMPLEX16 hGO_5_m_4(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z, struct kepler_vars params) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  if (vpnorder == 4) {
    return (-(mass * (1 - 5 * Nu + 5 * params.eta2) * PhiDOT *
              (mass * (82 * PhiDOT * r + Complex(0, 22) * rDOT) +
               3 * r *
                   (58 * params.PhiDOT3 * params.r3 +
                    Complex(0, 33) * params.PhiDOT2 * params.r2 * rDOT -
                    12 * PhiDOT * r * params.rDOT2 -
                    Complex(0, 2) * params.rDOT3))) /
            (36. * sqrt(165) * r));
  }

  else if (vpnorder == 5) {

    return (/* (32*Nu*((-1 + delta + 2*Nu)*S1z - (1 + delta -
      2*Nu)*S2z)*params.x3p5)/(3.*sqrt(165))
      + */
            /* Henry et al. ecc spin terms */ (
                (params.Mtot2 * Nu *
                 (82 * mass * PhiDOT * r + 174 * params.PhiDOT3 * params.r4 +
                  Complex(0, 22) * mass * rDOT +
                  Complex(0, 99) * params.PhiDOT2 * params.r3 * rDOT -
                  36 * PhiDOT * params.r2 * params.rDOT2 -
                  Complex(0, 6) * r * params.rDOT3) *
                 ((-1 + delta + 2 * Nu) * S1z - (1 + delta - 2 * Nu) * S2z)) /
                (24. * sqrt(165) * params.r3)));
  }

  else if (vpnorder == 6) {
    return (mass * PhiDOT *
            (-4 * params.Mtot2 *
                 (26 *
                      (-5051 + 28623 * Nu - 46305 * params.eta2 +
                       29470 * params.eta3) *
                      PhiDOT * r +
                  Complex(0, 7) *
                      (5684 - 26697 * Nu + 14225 * params.eta2 +
                       25355 * params.eta3) *
                      rDOT) -
             4 * mass * r *
                 ((-149157 + 1133006 * Nu - 2731750 * params.eta2 +
                   2085685 * params.eta3) *
                      params.PhiDOT3 * params.r3 +
                  Complex(0, 7) *
                      (101118 - 491779 * Nu + 402185 * params.eta2 +
                       172105 * params.eta3) *
                      params.PhiDOT2 * params.r2 * rDOT -
                  (363357 - 1825691 * Nu + 1714720 * params.eta2 +
                   395255 * params.eta3) *
                      PhiDOT * r * params.rDOT2 -
                  Complex(0, 7) *
                      (9717 - 48896 * Nu + 46000 * params.eta2 +
                       10865 * params.eta3) *
                      params.rDOT3) +
             15 * params.r2 *
                 (4 *
                      (3449 - 6580 * Nu - 56728 * params.eta2 +
                       115269 * params.eta3) *
                      params.PhiDOT5 * params.r5 +
                  Complex(0, 7) *
                      (-8128 + 45859 * Nu - 62702 * params.eta2 +
                       13996 * params.eta3) *
                      params.PhiDOT4 * params.r4 * rDOT +
                  4 *
                      (10125 - 47852 * Nu + 26635 * params.eta2 +
                       44240 * params.eta3) *
                      params.PhiDOT3 * params.r3 * params.rDOT2 +
                  Complex(0, 182) *
                      (127 - 548 * Nu + 73 * params.eta2 + 816 * params.eta3) *
                      params.PhiDOT2 * params.r2 * params.rDOT3 -
                  8 *
                      (1009 - 4060 * Nu - 889 * params.eta2 +
                       7952 * params.eta3) *
                      PhiDOT * r * params.rDOT4 -
                  Complex(0, 28) *
                      (45 - 172 * Nu - 85 * params.eta2 + 400 * params.eta3) *
                      params.rDOT5))) /
           (65520. * sqrt(165) * params.r2);
  }

  else if (vpnorder == 7) {

    return (/* (16*(-530*params.eta3*(S1z + S2z) + 104*(S1z + delta*S1z + S2z -
     delta*S2z)
      + 2*params.eta2*(541*delta*(S1z - S2z) - 120*(S1z + S2z)) -
     Nu*(1139*delta*(S1z - S2z) + 109*(S1z + S2z)))*params.x4p5)/(117.*sqrt(165))
     + */
            /* Henry et al. ecc+spin terms */
            ((params.Mtot2 *
              (4 * params.Mtot2 *
                   (210 * params.eta3 *
                        (1950 * PhiDOT * r + Complex(0, 1373) * rDOT) *
                        (S1z + S2z) +
                    10920 * (20 * PhiDOT * r + Complex(0, 7) * rDOT) *
                        (S1z + delta * S1z + S2z - delta * S2z) +
                    21 * params.eta2 *
                        (Complex(0, 5) * rDOT *
                             ((8093 + 1355 * delta) * S1z +
                              (8093 - 1355 * delta) * S2z) +
                         26 * PhiDOT * r *
                             (Complex(0, 554) + 25 * (-157 + 89 * delta) * S1z -
                              25 * (157 + 89 * delta) * S2z)) -
                    5 * Nu *
                        (Complex(0, 42) * rDOT *
                             ((2479 + 1083 * delta) * S1z +
                              (2479 - 1083 * delta) * S2z) +
                         13 * PhiDOT * r *
                             (Complex(0, 1675) +
                              168 * (-74 + 295 * delta) * S1z -
                              168 * (74 + 295 * delta) * S2z))) +
               3 * params.r2 *
                   (-840 * params.eta3 *
                        (13072 * params.PhiDOT5 * params.r5 +
                         Complex(0, 679) * params.PhiDOT4 * params.r4 * rDOT +
                         8700 * params.PhiDOT3 * params.r3 * params.rDOT2 +
                         Complex(0, 4459) * params.PhiDOT2 * params.r2 *
                             params.rDOT3 -
                         1432 * PhiDOT * r * params.rDOT4 -
                         Complex(0, 210) * params.rDOT5) *
                        (S1z + S2z) +
                    Complex(0, 859950) * params.PhiDOT4 * params.r4 * rDOT *
                        (S1z + delta * S1z + S2z - delta * S2z) +
                    14 * params.eta2 *
                        (30 * params.PhiDOT3 * params.r3 * params.rDOT2 *
                             (Complex(0, -3796) + (77 + 4765 * delta) * S1z +
                              (77 - 4765 * delta) * S2z) -
                         Complex(0, 210) * params.rDOT5 *
                             ((-23 + 5 * delta) * S1z -
                              (23 + 5 * delta) * S2z) +
                         Complex(0, 13) * params.PhiDOT2 * params.r2 *
                             params.rDOT3 *
                             (Complex(0, 5296) +
                              105 * (-29 + 39 * delta) * S1z -
                              105 * (29 + 39 * delta) * S2z) +
                         Complex(0, 2) * params.PhiDOT4 * params.r4 * rDOT *
                             (Complex(0, -1058096) +
                              105 * (4821 + 131 * delta) * S1z -
                              105 * (-4821 + 131 * delta) * S2z) -
                         4 * PhiDOT * r * params.rDOT4 *
                             (Complex(0, 3796) +
                              15 * (-399 + 193 * delta) * S1z -
                              15 * (399 + 193 * delta) * S2z) +
                         4 * params.PhiDOT5 * params.r5 *
                             (Complex(0, 387296) +
                              15 * (-623 + 4413 * delta) * S1z -
                              15 * (623 + 4413 * delta) * S2z)) -
                    3 * Nu *
                        (Complex(0, -140) * params.rDOT5 *
                             ((-141 + 193 * delta) * S1z -
                              (141 + 193 * delta) * S2z) +
                         Complex(0, 182) * params.PhiDOT2 * params.r2 *
                             params.rDOT3 *
                             (Complex(0, 632) + 5 * (-319 + 555 * delta) * S1z -
                              5 * (319 + 555 * delta) * S2z) -
                         8 * PhiDOT * r * params.rDOT4 *
                             (Complex(0, 2574) +
                              35 * (-415 + 623 * delta) * S1z -
                              35 * (415 + 623 * delta) * S2z) +
                         2 * params.PhiDOT5 * params.r5 *
                             (Complex(0, 1106001) +
                              140 * (-1069 + 3045 * delta) * S1z -
                              140 * (1069 + 3045 * delta) * S2z) +
                         20 * params.PhiDOT3 * params.r3 * params.rDOT2 *
                             (Complex(0, -8502) +
                              7 * (-2883 + 6419 * delta) * S1z -
                              7 * (2883 + 6419 * delta) * S2z) +
                         Complex(0, 7) * params.PhiDOT4 * params.r4 * rDOT *
                             (Complex(0, -431288) +
                              5 * (72097 + 7957 * delta) * S1z -
                              5 * (-72097 + 7957 * delta) * S2z))) -
               6 * mass * r *
                   (420 * params.eta3 *
                        (2698 * params.PhiDOT3 * params.r3 +
                         Complex(0, 307) * params.PhiDOT2 * params.r2 * rDOT +
                         739 * PhiDOT * r * params.rDOT2 +
                         Complex(0, 183) * params.rDOT3) *
                        (S1z + S2z) -
                    5460 *
                        (86 * params.PhiDOT3 * params.r3 +
                         Complex(0, 137) * params.PhiDOT2 * params.r2 * rDOT -
                         30 * PhiDOT * r * params.rDOT2 -
                         Complex(0, 3) * params.rDOT3) *
                        (S1z + delta * S1z + S2z - delta * S2z) +
                    Nu * (Complex(0, 700) * params.rDOT3 *
                              ((-651 + 482 * delta) * S1z -
                               (651 + 482 * delta) * S2z) +
                          28 * params.PhiDOT2 * params.r2 * rDOT *
                              (34294 -
                               Complex(0, 5) * (-56523 + 1013 * delta) * S1z +
                               Complex(0, 5) * (56523 + 1013 * delta) * S2z) +
                          params.PhiDOT3 * params.r3 *
                              (Complex(0, 1423669) +
                               280 * (2589 + 18458 * delta) * S1z -
                               280 * (-2589 + 18458 * delta) * S2z) +
                          2 * PhiDOT * r * params.rDOT2 *
                              (Complex(0, 34203) +
                               35 * (-40527 + 19571 * delta) * S1z -
                               35 * (40527 + 19571 * delta) * S2z)) -
                    14 * params.eta2 *
                        (Complex(0, 50) * params.rDOT3 *
                             ((-1330 + 223 * delta) * S1z -
                              (1330 + 223 * delta) * S2z) +
                         2 * PhiDOT * r * params.rDOT2 *
                             (Complex(0, 7956) +
                              5 * (-39389 + 2528 * delta) * S1z -
                              5 * (39389 + 2528 * delta) * S2z) +
                         Complex(0, 1) * params.PhiDOT2 * params.r2 * rDOT *
                             (Complex(0, -213408) +
                              5 * (202025 + 36547 * delta) * S1z -
                              5 * (-202025 + 36547 * delta) * S2z) +
                         2 * params.PhiDOT3 * params.r3 *
                             (Complex(0, 160212) +
                              5 * (-2831 + 38849 * delta) * S1z -
                              5 * (2831 + 38849 * delta) * S2z))))) /
             (393120. * sqrt(165) * params.r4)));
  }

  else {
    return 0;
  }
}

/* static COMPLEX16 hQC_5_m_4(REAL8 Nu, UINT4 vpnorder, REAL8 x){

    if(vpnorder == 4){
        return((-64/(9.*sqrt(165)) + (64*sqrt(0.15151515151515152)*Nu)/9. -
        (64*sqrt(0.15151515151515152)*params.eta2)/9.)*params.x3);
    }
    else if(vpnorder == 6){
        return((142432/(4095.*sqrt(165)) -
(10528*sqrt(0.7333333333333333)*Nu)/585. + (33344*params.eta2)/(117.*sqrt(165)) -
        (3616*params.eta3)/(39.*sqrt(165)))*params.x4);
    }
    else{
        return 0;
    }
} */

static COMPLEX16 hl_5_m_4(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_5_m_4: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_5_m_4(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z,
                      S2z, params) /* +hQC_5_m_4(Nu,vpnorder,x) */) *
           cpolar(1, -4 * Phi);
  }
}

static COMPLEX16 hl_5_m_min4(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_5_m_min4: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((- 4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           conj(hGO_5_m_4(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z,
                          S2z, params) /* +hQC_5_m_4(Nu,vpnorder,x) */) *
           cpolar(1, 4 * Phi);
  }
}

// H53
static COMPLEX16 hGO_5_m_3(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z, struct kepler_vars params) {
  REAL8 delta = sqrt(1 - 4 * Nu);
  REAL8 kappa1 = 1.0;
  REAL8 kappa2 = 1.0;
  REAL8 Nu_factor = (1-2*Nu) * (1-2*Nu);

  if (vpnorder == 3) {
    return (delta * (-1 + 2 * Nu) *
            (2 * params.Mtot2 * (Complex(0, 258) * PhiDOT * r - 205 * rDOT) -
             Complex(0, 120) * params.r2 * (PhiDOT * r - Complex(0, 1) * rDOT) *
                 pow(PhiDOT * r + Complex(0, 1) * rDOT, 4) +
             3 * mass * r *
                 (Complex(0, -51) * params.PhiDOT3 * params.r3 -
                  240 * params.PhiDOT2 * params.r2 * rDOT -
                  Complex(0, 396) * PhiDOT * r * params.rDOT2 +
                  160 * params.rDOT3))) /
           (144. * sqrt(330) * params.r2);
  }

  else if (vpnorder == 5) {
    return (delta *
            (Complex(0, 120) * (33 - 197 * Nu + 294 * params.eta2) * params.r3 *
                 pow(PhiDOT * r - Complex(0, 1) * rDOT, 2) *
                 pow(PhiDOT * r + Complex(0, 1) * rDOT, 5) +
             2 * params.Mtot3 *
                 (Complex(0, 1) * (53311 - 121906 * Nu + 42816 * params.eta2) *
                      PhiDOT * r -
                  26 * (1141 - 2760 * Nu + 1420 * params.eta2) * rDOT) +
             2 * params.Mtot2 * r *
                 (Complex(0, -2) * (6350 - 12803 * Nu + 10314 * params.eta2) *
                      params.PhiDOT3 * params.r3 +
                  (-6546 + 64131 * Nu - 109702 * params.eta2) * params.PhiDOT2 *
                      params.r2 * rDOT -
                  Complex(0, 4) * (-2636 - 11335 * Nu + 43962 * params.eta2) *
                      PhiDOT * r * params.rDOT2 +
                  13 * (-639 - 735 * Nu + 5690 * params.eta2) * params.rDOT3) +
             3 * mass * params.r2 *
                 (Complex(0, -4) * (1223 - 4567 * Nu + 3396 * params.eta2) *
                      params.PhiDOT5 * params.r5 +
                  26 * (-412 - 437 * Nu + 3114 * params.eta2) * params.PhiDOT4 *
                      params.r4 * rDOT +
                  Complex(0, 1) * (-2331 - 31496 * Nu + 93276 * params.eta2) *
                      params.PhiDOT3 * params.r3 * params.rDOT2 +
                  8 * (-1994 + 1541 * Nu + 5718 * params.eta2) * params.PhiDOT2 *
                      params.r2 * params.rDOT3 +
                  Complex(0, 4) * (-3971 - 3226 * Nu + 27360 * params.eta2) *
                      PhiDOT * r * params.rDOT4 -
                  80 * (-57 - 97 * Nu + 518 * params.eta2) * params.rDOT5))) /
           (3744. * sqrt(330) * params.r3);
  }

  else if (vpnorder == 6) {
    return (
        (delta * params.Mtot2 * Nu *
         (17830 * params.Mtot2 +
          18 * mass * r *
              (5231 * params.PhiDOT2 * params.r2 +
               Complex(0, 3921) * PhiDOT * r * rDOT - 1510 * params.rDOT2) -
          params.r2 * (48579 * params.PhiDOT4 * params.r4 +
                       Complex(0, 31304) * params.PhiDOT3 * params.r3 * rDOT +
                       33024 * params.PhiDOT2 * params.r2 * params.rDOT2 +
                       Complex(0, 29904) * PhiDOT * r * params.rDOT3 -
                       7120 * params.rDOT4))) /
            (5040. * sqrt(330) * params.r4)
        /* +(Complex(0,0.1875)*sqrt(0.02727272727272727)*((5 + Nu*(-27 + 31*Nu)
   + delta*(5 + 23*(-1 + Nu)*Nu))*S1z
   + (-5 + (27 - 31*Nu)*Nu + delta*(5 + 23*(-1 + Nu)*Nu))*S2z)*params.x4) */
        +
        /* Henry et al. ecc spin terms */ (
            (Complex(0, 0.006944444444444444) * params.Mtot2 *
             (params.Mtot2 * ((90 - 622 * Nu + 966 * params.eta2 +
                               delta * (90 - 524 * Nu + 770 * params.eta2)) *
                                  S1z +
                              2 *
                                  (-45 + 311 * Nu - 483 * params.eta2 +
                                   delta * (45 - 262 * Nu + 385 * params.eta2)) *
                                  S2z) -
              3 * params.r2 *
                  (Complex(0, 4) * PhiDOT * r * params.rDOT3 *
                       (5 * (9 - 56 * Nu + 78 * params.eta2) * S1z +
                        delta * (45 - 142 * Nu + 114 * params.eta2) * S1z -
                        5 * (9 - 56 * Nu + 78 * params.eta2) * S2z +
                        delta * (45 - 142 * Nu + 114 * params.eta2) * S2z) +
                   Complex(0, 1) * params.PhiDOT3 * params.r3 * rDOT *
                       (25 * (3 - 14 * Nu + 12 * params.eta2) * S1z +
                        delta * (75 - 958 * Nu + 1516 * params.eta2) * S1z -
                        25 * (3 - 14 * Nu + 12 * params.eta2) * S2z +
                        delta * (75 - 958 * Nu + 1516 * params.eta2) * S2z) -
                   4 * params.rDOT4 *
                       ((15 - 92 * Nu + 126 * params.eta2 +
                         delta * (15 - 64 * Nu + 70 * params.eta2)) *
                            S1z +
                        (-15 + 92 * Nu - 126 * params.eta2 +
                         delta * (15 - 64 * Nu + 70 * params.eta2)) *
                            S2z) +
                   params.PhiDOT4 * params.r4 *
                       ((20 - 317 * Nu + 751 * params.eta2 +
                         delta * (20 + 13 * Nu + 91 * params.eta2)) *
                            S1z +
                        (-20 + 317 * Nu - 751 * params.eta2 +
                         delta * (20 + 13 * Nu + 91 * params.eta2)) *
                            S2z) -
                   3 * params.PhiDOT2 * params.r2 * params.rDOT2 *
                       ((-45 + 298 * Nu - 444 * params.eta2 +
                         delta * (-45 + 2 * Nu + 148 * params.eta2)) *
                            S1z +
                        (45 - 298 * Nu + 444 * params.eta2 +
                         delta * (-45 + 2 * Nu + 148 * params.eta2)) *
                            S2z)) +
              mass * r *
                  (Complex(0, 6) * PhiDOT * r * rDOT *
                       (5 * (38 - 243 * Nu + 349 * params.eta2) * S1z +
                        delta * (190 - 643 * Nu + 601 * params.eta2) * S1z -
                        5 * (38 - 243 * Nu + 349 * params.eta2) * S2z +
                        delta * (190 - 643 * Nu + 601 * params.eta2) * S2z) +
                   params.PhiDOT2 * params.r2 *
                       ((849 - 5360 * Nu + 7590 * params.eta2 +
                         delta * (849 - 1774 * Nu + 418 * params.eta2)) *
                            S1z +
                        (-849 + 5360 * Nu - 7590 * params.eta2 +
                         delta * (849 - 1774 * Nu + 418 * params.eta2)) *
                            S2z) -
                   2 * params.rDOT2 *
                       ((255 - 1592 * Nu + 2226 * params.eta2 +
                         delta * (255 - 1144 * Nu + 1330 * params.eta2)) *
                            S1z +
                        (-255 + 1592 * Nu - 2226 * params.eta2 +
                         delta * (255 - 1144 * Nu + 1330 * params.eta2)) *
                            S2z)))) /
            (sqrt(330) * params.r4)));
  }

  else if (vpnorder == 7) {

    return (/* Complex(0,-0.140625)*sqrt(0.02727272727272727)*(kappa1*(5 +
      5*delta*Nu_factor
      - 30*Nu + 8*params.eta2)*params.S1z2 + S2z*(20*delta*(1 - 2*Nu)*Nu*S1z
      + kappa2*(-5 + 5*delta*Nu_factor + 30*Nu -
      8*params.eta2)*S2z))*params.x4p5
      + */
            /* Henry et al. ecc+spin terms */
            ((8190 * params.Mtot3 *
                  (Complex(0, -2) * mass *
                       (447 - 2682 * Nu + 3224 * params.eta2) * PhiDOT * r +
                   2 * mass * (235 + 2 * Nu * (-705 + 1004 * Nu)) * rDOT +
                   r * (Complex(0, 1) * (15 - 90 * Nu + 2008 * params.eta2) *
                            params.PhiDOT3 * params.r3 +
                        66 * (5 + 6 * Nu * (-5 + 12 * Nu)) * params.PhiDOT2 *
                            params.r2 * rDOT +
                        Complex(0, 6) * (85 - 510 * Nu + 808 * params.eta2) *
                            PhiDOT * r * params.rDOT2 -
                        8 * (-5 + 12 * Nu) * (-5 + 18 * Nu) * params.rDOT3)) *
                  (kappa1 * params.S1z2 - kappa2 * params.S2z2) +
              delta *
                  (Complex(0, -12600) *
                       (-135 + Nu * (1467 + 35 * Nu * (-151 + 178 * Nu))) *
                       params.r4 * pow(PhiDOT * r - Complex(0, 1) * rDOT, 3) *
                       pow(PhiDOT * r + Complex(0, 1) * rDOT, 6) +
                   6 * params.Mtot2 * params.r2 *
                       (Complex(0, -10) *
                            (129184 +
                             Nu *
                                 (917839 + 2 * Nu * (-1601563 + 652752 * Nu))) *
                            params.PhiDOT5 * params.r5 +
                        8 *
                            (3087667 +
                             2 * Nu *
                                 (-3131631 +
                                  5 * Nu * (-385362 + 1212143 * Nu))) *
                            params.PhiDOT4 * params.r4 * rDOT +
                        Complex(0, 4) *
                            (1877161 +
                             Nu * (-1391951 +
                                   20 * Nu * (-594787 + 1654110 * Nu))) *
                            params.PhiDOT3 * params.r3 * params.rDOT2 +
                        2 *
                            (-919652 +
                             Nu * (-16365978 +
                                   5 * Nu * (5343599 + 5349122 * Nu))) *
                            params.PhiDOT2 * params.r2 * params.rDOT3 +
                        Complex(0, 14) *
                            (-19909 + Nu * (-2880941 +
                                            15 * Nu * (195291 + 744862 * Nu))) *
                            PhiDOT * r * params.rDOT4 -
                        35 *
                            (1285 +
                             Nu * (-417939 + Nu * (376375 + 1794742 * Nu))) *
                            params.rDOT5) +
                   45 * mass * params.r3 *
                       (Complex(0, 4) *
                            (54310 +
                             Nu * (-20821 + Nu * (-636477 + 820286 * Nu))) *
                            params.PhiDOT7 * params.r7 +
                        4 *
                            (-119296 + Nu * (251655 + 893222 * Nu -
                                             3122048 * params.eta2)) *
                            params.PhiDOT6 * params.r6 * rDOT -
                        Complex(0, 4) *
                            (42783 +
                             Nu * (123661 + 3 * Nu * (-447987 + 841994 * Nu))) *
                            params.PhiDOT5 * params.r5 * params.rDOT2 -
                        8 *
                            (237862 +
                             Nu * (-912196 + 3 * Nu * (-46779 + 933382 * Nu))) *
                            params.PhiDOT4 * params.r4 * params.rDOT3 -
                        Complex(0, 1) *
                            (1550225 +
                             6 * Nu *
                                 (-680405 +
                                  4 * Nu * (-402523 + 1277568 * Nu))) *
                            params.PhiDOT3 * params.r3 * params.rDOT4 -
                        16 *
                            (39855 +
                             2 * Nu * (-106258 + Nu * (113231 + 75438 * Nu))) *
                            params.PhiDOT2 * params.r2 * params.rDOT5 -
                        Complex(0, 28) *
                            (41703 +
                             2 * Nu *
                                 (-74181 + 20 * Nu * (-2603 + 15676 * Nu))) *
                            PhiDOT * r * params.rDOT6 +
                        2240 *
                            (162 + Nu * (-516 + 7 * Nu * (-101 + 410 * Nu))) *
                            params.rDOT7) +
                   28 * params.Mtot4 *
                       (Complex(0, 1) * PhiDOT * r *
                            (-5713898 - 261495 * kappa1 * params.S1z2 +
                             2 * Nu *
                                 (6824514 - 4178445 * Nu +
                                  1430090 * params.eta2 -
                                  522990 * kappa1 * (-1 + Nu) * params.S1z2) +
                             1045980 * Nu * (-1 + 2 * Nu) * S1z * S2z -
                             261495 * kappa2 * Nu_factor *
                                 params.S2z2) +
                        5 * rDOT *
                            (300238 + 27495 * kappa1 * params.S1z2 +
                             Nu * (-581955 + 600300 * Nu - 535772 * params.eta2 +
                                   109980 * kappa1 * (-1 + Nu) * params.S1z2) -
                             109980 * Nu * (-1 + 2 * Nu) * S1z * S2z +
                             27495 * kappa2 * Nu_factor *
                                 params.S2z2)) +
                   2 * params.Mtot3 * r *
                       (Complex(0, -1) * params.PhiDOT3 * params.r3 *
                            (-13 * (371212 + 4725 * kappa1 * params.S1z2) +
                             Nu * (-7716399 - 11038350 * Nu +
                                   21939560 * params.eta2 -
                                   245700 * kappa1 * (-1 + Nu) * params.S1z2) +
                             245700 * Nu * (-1 + 2 * Nu) * S1z * S2z -
                             61425 * kappa2 * Nu_factor *
                                 params.S2z2) -
                        Complex(0, 126) * PhiDOT * r * params.rDOT2 *
                            ((1974194 - 5423729 * Nu + 339235 * params.eta2 +
                              5947390 * params.eta3 -
                              49725 * kappa1 * Nu_factor *
                                  params.S1z2) /
                                 3. +
                             66300 * Nu * (-1 + 2 * Nu) * S1z * S2z -
                             16575 * kappa2 * Nu_factor *
                                 params.S2z2) +
                        140 * params.rDOT3 *
                            (175718 - 5850 * kappa1 * params.S1z2 +
                             Nu * (-730551 + 504765 * Nu + 773906 * params.eta2 -
                                   23400 * kappa1 * (-1 + Nu) * params.S1z2) +
                             23400 * Nu * (-1 + 2 * Nu) * S1z * S2z -
                             5850 * kappa2 * Nu_factor * params.S2z2) +
                        6 * params.PhiDOT2 * params.r2 * rDOT *
                            (-8099876 + 225225 * kappa1 * params.S1z2 +
                             4 * Nu *
                                 (1323649 + 9536535 * Nu -
                                  6899570 * params.eta2 +
                                  225225 * kappa1 * (-1 + Nu) * params.S1z2) -
                             900900 * Nu * (-1 + 2 * Nu) * S1z * S2z +
                             225225 * kappa2 * Nu_factor *
                                 params.S2z2)))) /
             (1.57248e6 * sqrt(330) * params.r4)));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hQC_5_m_3(REAL8 Nu, UINT4 vpnorder, REAL8 x, struct kepler_vars params) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  /* if(vpnorder == 3){
       return((Complex(0,-0.5625)*sqrt(0.02727272727272727)*delta +
       Complex(0,1.125)*sqrt(0.02727272727272727)*delta*Nu)*params.x2p5);
   }
   else if(vpnorder == 5){
       return((Complex(0,2.985576923076923)*sqrt(0.02727272727272727)*delta -
       Complex(0,6.6923076923076925)*sqrt(0.02727272727272727)*delta*Nu +
       Complex(0,0.11538461538461539)*sqrt(3.3)*delta*params.eta2)*
       params.x3p5);
   } */
  /* else */ if (vpnorder == 6) {
    return (Complex(0, 1.6875) * sqrt(0.02727272727272727) * delta *
            (-1 + 2 * Nu) * params.x4 * (M_PI + Complex(0, 2) * log(1.5)));
  } else {
    return 0;
  }
}

static COMPLEX16 hl_5_m_3(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z, REAL8 x, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_5_m_3: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_5_m_3(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, params) +
            hQC_5_m_3(Nu, vpnorder, x, params)) *
           cpolar(1, -3 * Phi);
  }
}

static COMPLEX16 hl_5_m_min3(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z, REAL8 x, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_5_m_min3: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((- 4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           conj(hGO_5_m_3(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, params) +
                hQC_5_m_3(Nu, vpnorder, x, params)) *
           cpolar(1, 3 * Phi);
  }
}

// H52
static COMPLEX16 hGO_5_m_2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z, struct kepler_vars params) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  if (vpnorder == 4) {
    return (mass * (1 - 5 * Nu + 5 * params.eta2) * PhiDOT *
            (mass * (41 * PhiDOT * r + Complex(0, 22) * rDOT) -
             3 * r *
                 (11 * params.PhiDOT3 * params.r3 -
                  Complex(0, 3) * params.PhiDOT2 * params.r2 * rDOT +
                  6 * PhiDOT * r * params.rDOT2 +
                  Complex(0, 2) * params.rDOT3))) /
           (54. * sqrt(55) * r);
  }

  else if (vpnorder == 5) {

    return (/* (2*Nu*(S1z - delta*S1z + S2z + delta*S2z - 2*Nu*(S1z +
      S2z))*params.x3p5)/(9.*sqrt(55))
      + */
            /* Henry et al. ecc spin terms */ (
                -0.027777777777777776 *
                (params.Mtot2 * Nu *
                 (mass * (41 * PhiDOT * r + Complex(0, 22) * rDOT) -
                  3 * r *
                      (11 * params.PhiDOT3 * params.r3 -
                       Complex(0, 3) * params.PhiDOT2 * params.r2 * rDOT +
                       6 * PhiDOT * r * params.rDOT2 +
                       Complex(0, 2) * params.rDOT3)) *
                 ((-1 + delta + 2 * Nu) * S1z - (1 + delta - 2 * Nu) * S2z)) /
                (sqrt(55) * params.r3)));
  }

  else if (vpnorder == 6) {
    return (mass * PhiDOT *
            (4 * params.Mtot2 *
                 (13 *
                      (-5051 + 28623 * Nu - 46305 * params.eta2 +
                       29470 * params.eta3) *
                      PhiDOT * r +
                  Complex(0, 7) *
                      (5684 - 26697 * Nu + 14225 * params.eta2 +
                       25355 * params.eta3) *
                      rDOT) -
             2 * mass * r *
                 ((23157 + 154 * Nu - 648410 * params.eta2 +
                   1133195 * params.eta3) *
                      params.PhiDOT3 * params.r3 -
                  Complex(0, 14) *
                      (6648 - 31729 * Nu + 23795 * params.eta2 +
                       13225 * params.eta3) *
                      params.PhiDOT2 * params.r2 * rDOT +
                  (363357 - 1825691 * Nu + 1714720 * params.eta2 +
                   395255 * params.eta3) *
                      PhiDOT * r * params.rDOT2 +
                  Complex(0, 14) *
                      (9717 - 48896 * Nu + 46000 * params.eta2 +
                       10865 * params.eta3) *
                      params.rDOT3) +
             15 * params.r2 *
                 (2 *
                      (2207 - 6076 * Nu - 18424 * params.eta2 +
                       38787 * params.eta3) *
                      params.PhiDOT5 * params.r5 -
                  Complex(0, 7) *
                      (4744 - 23965 * Nu + 23906 * params.eta2 +
                       1892 * params.eta3) *
                      params.PhiDOT4 * params.r4 * rDOT +
                  2 *
                      (5667 - 22652 * Nu - 5747 * params.eta2 +
                       45416 * params.eta3) *
                      params.PhiDOT3 * params.r3 * params.rDOT2 -
                  Complex(0, 14) *
                      (97 - 536 * Nu + 643 * params.eta2 + 36 * params.eta3) *
                      params.PhiDOT2 * params.r2 * params.rDOT3 +
                  4 *
                      (1009 - 4060 * Nu - 889 * params.eta2 +
                       7952 * params.eta3) *
                      PhiDOT * r * params.rDOT4 +
                  Complex(0, 28) *
                      (45 - 172 * Nu - 85 * params.eta2 + 400 * params.eta3) *
                      params.rDOT5))) /
           (98280. * sqrt(55) * params.r2);
  }

  else if (vpnorder == 7) {

    return (
        (/* (698*params.eta3*(S1z + S2z) - 104*(S1z + delta*S1z + S2z - delta*S2z)
   - 2*params.eta2*(457*delta*(S1z - S2z) + 48*(S1z + S2z)) +
  Nu*(1055*delta*(S1z - S2z) + 193*(S1z + S2z)))*params.x4p5)/(351.*sqrt(55))
  + */
         /* Henry et al. ecc+spin terms */
         (-1.6958350291683625e-6 *
          (params.Mtot2 *
           (4 * params.Mtot2 *
                (210 * params.eta3 *
                     (975 * PhiDOT * r + Complex(0, 1373) * rDOT) *
                     (S1z + S2z) +
                 10920 * (10 * PhiDOT * r + Complex(0, 7) * rDOT) *
                     (S1z + delta * S1z + S2z - delta * S2z) +
                 21 * params.eta2 *
                     (Complex(0, 5) * rDOT *
                          ((8093 + 1355 * delta) * S1z +
                           (8093 - 1355 * delta) * S2z) +
                      13 * PhiDOT * r *
                          (Complex(0, 1108) + 25 * (-157 + 89 * delta) * S1z -
                           25 * (157 + 89 * delta) * S2z)) -
                 5 * Nu *
                     (Complex(0, 42) * rDOT *
                          ((2479 + 1083 * delta) * S1z +
                           (2479 - 1083 * delta) * S2z) +
                      13 * PhiDOT * r *
                          (Complex(0, 1675) + 84 * (-74 + 295 * delta) * S1z -
                           84 * (74 + 295 * delta) * S2z))) +
            3 * params.r2 *
                (840 * params.eta3 *
                     (896 * params.PhiDOT5 * params.r5 +
                      Complex(0, 1390) * params.PhiDOT4 * params.r4 * rDOT +
                      2682 * params.PhiDOT3 * params.r3 * params.rDOT2 -
                      Complex(0, 235) * params.PhiDOT2 * params.r2 *
                          params.rDOT3 +
                      716 * PhiDOT * r * params.rDOT4 +
                      Complex(0, 210) * params.rDOT5) *
                     (S1z + S2z) -
                 Complex(0, 286650) * params.PhiDOT4 * params.r4 * rDOT *
                     (S1z + delta * S1z + S2z - delta * S2z) -
                 14 * params.eta2 *
                     (Complex(0, 2) * params.PhiDOT4 * params.r4 * rDOT *
                          (Complex(0, -41392) +
                           15 * (15219 + 1835 * delta) * S1z +
                           (228285 - 27525 * delta) * S2z) +
                      Complex(0, 210) * params.rDOT5 *
                          ((-23 + 5 * delta) * S1z - (23 + 5 * delta) * S2z) +
                      2 * PhiDOT * r * params.rDOT4 *
                          (Complex(0, 7592) + 15 * (-399 + 193 * delta) * S1z -
                           15 * (399 + 193 * delta) * S2z) +
                      10 * params.PhiDOT5 * params.r5 *
                          (Complex(0, 1508) + 3 * (-1889 + 483 * delta) * S1z -
                           3 * (1889 + 483 * delta) * S2z) +
                      params.PhiDOT2 * params.r2 * params.rDOT3 *
                          (34424 - Complex(0, 15) * (-161 + 555 * delta) * S1z +
                           Complex(0, 15) * (161 + 555 * delta) * S2z) +
                      3 * params.PhiDOT3 * params.r3 * params.rDOT2 *
                          (Complex(0, 15288) + 5 * (539 + 907 * delta) * S1z -
                           5 * (-539 + 907 * delta) * S2z)) +
                 3 * Nu *
                     (Complex(0, 140) * params.rDOT5 *
                          ((-141 + 193 * delta) * S1z -
                           (141 + 193 * delta) * S2z) +
                      14 * params.PhiDOT2 * params.r2 * params.rDOT3 *
                          (4108 - Complex(0, 5) * (-193 + 453 * delta) * S1z +
                           Complex(0, 5) * (193 + 453 * delta) * S2z) +
                      4 * PhiDOT * r * params.rDOT4 *
                          (Complex(0, 5148) + 35 * (-415 + 623 * delta) * S1z -
                           35 * (415 + 623 * delta) * S2z) +
                      10 * params.PhiDOT5 * params.r5 *
                          (Complex(0, 2223) + 14 * (-859 + 963 * delta) * S1z -
                           14 * (859 + 963 * delta) * S2z) +
                      7 * params.PhiDOT4 * params.r4 * rDOT *
                          (16276 -
                           Complex(0, 5) * (-31495 + 2453 * delta) * S1z +
                           Complex(0, 5) * (31495 + 2453 * delta) * S2z) +
                      2 * params.PhiDOT3 * params.r3 * params.rDOT2 *
                          (Complex(0, 39156) +
                           35 * (-1773 + 3437 * delta) * S1z -
                           35 * (1773 + 3437 * delta) * S2z))) -
            6 * mass * r *
                (210 * params.eta3 *
                     (946 * params.PhiDOT3 * params.r3 -
                      Complex(0, 1150) * params.PhiDOT2 * params.r2 * rDOT +
                      739 * PhiDOT * r * params.rDOT2 +
                      Complex(0, 366) * params.rDOT3) *
                     (S1z + S2z) -
                 5460 *
                     (3 * params.PhiDOT3 * params.r3 +
                      Complex(0, 17) * params.PhiDOT2 * params.r2 * rDOT -
                      15 * PhiDOT * r * params.rDOT2 -
                      Complex(0, 3) * params.rDOT3) *
                     (S1z + delta * S1z + S2z - delta * S2z) -
                 14 * params.eta2 *
                     (Complex(0, 1) * params.PhiDOT2 * params.r2 * rDOT *
                          (Complex(0, -106704) +
                           5 * (13991 + 6349 * delta) * S1z +
                           (69955 - 31745 * delta) * S2z) +
                      Complex(0, 50) * params.rDOT3 *
                          ((-1330 + 223 * delta) * S1z -
                           (1330 + 223 * delta) * S2z) +
                      PhiDOT * r * params.rDOT2 *
                          (Complex(0, 15912) +
                           5 * (-39389 + 2528 * delta) * S1z -
                           5 * (39389 + 2528 * delta) * S2z) +
                      params.PhiDOT3 * params.r3 *
                          (Complex(0, 62868) +
                           5 * (-5543 + 5105 * delta) * S1z -
                           5 * (5543 + 5105 * delta) * S2z)) +
                 Nu * (Complex(0, 28) * params.PhiDOT2 * params.r2 * rDOT *
                           (Complex(0, -17147) +
                            5 * (5007 + 1675 * delta) * S1z +
                            (25035 - 8375 * delta) * S2z) +
                       5 * params.PhiDOT3 * params.r3 *
                           (Complex(0, 57161) + 28 * (495 + 272 * delta) * S1z -
                            28 * (-495 + 272 * delta) * S2z) +
                       Complex(0, 700) * params.rDOT3 *
                           ((-651 + 482 * delta) * S1z -
                            (651 + 482 * delta) * S2z) +
                       PhiDOT * r * params.rDOT2 *
                           (Complex(0, 68406) +
                            35 * (-40527 + 19571 * delta) * S1z -
                            35 * (40527 + 19571 * delta) * S2z))))) /
          (sqrt(55) * params.r4))));
  } else {
    return 0;
  }
}

/* static COMPLEX16 hQC_5_m_2(REAL8 Nu, UINT4 vpnorder, REAL8 x){

    if(vpnorder == 4){
        return((4/(27.*sqrt(55)) - (4*sqrt(0.45454545454545453)*Nu)/27. +
        (4*sqrt(0.45454545454545453)*params.eta2)/27.)*params.x3);
    }
    else if(vpnorder == 6){
        return((-7822/(12285.*sqrt(55)) + (6158*Nu)/(1755.*sqrt(55)) -
        (1652*params.eta2)/(351.*sqrt(55)) +
        (14*sqrt(2.2)*params.eta3)/117.)*params.x4);
    }
    else{
        return 0;
    }
} */

static COMPLEX16 hl_5_m_2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_5_m_2: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_5_m_2(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z,
                      S2z, params) /* +hQC_5_m_2(Nu,vpnorder,x) */) *
           cpolar(1, -2 * Phi);
  }
}

static COMPLEX16 hl_5_m_min2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_5_m_min2: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((- 4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           conj(hGO_5_m_2(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z,
                          S2z, params) /* +hQC_5_m_2(Nu,vpnorder,x) */) *
           cpolar(1, 2 * Phi);
  }
}

// H51
static COMPLEX16 hGO_5_m_1(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z, struct kepler_vars params) {
  REAL8 delta = sqrt(1 - 4 * Nu);
  REAL8 kappa1 = 1.0;
  REAL8 kappa2 = 1.0;
  REAL8 Nu_factor = (1-2*Nu) * (1-2*Nu);

  if (vpnorder == 3) {
    return (delta * (-1 + 2 * Nu) *
            (120 * params.r2 * pow(Complex(0, 1) * PhiDOT * r - rDOT, 3) *
                 pow(PhiDOT * r - Complex(0, 1) * rDOT, 2) +
             2 * params.Mtot2 * (Complex(0, -86) * PhiDOT * r + 205 * rDOT) +
             Complex(0, 3) * mass * r *
                 (97 * params.PhiDOT3 * params.r3 +
                  Complex(0, 160) * params.PhiDOT2 * params.r2 * rDOT +
                  132 * PhiDOT * r * params.rDOT2 +
                  Complex(0, 160) * params.rDOT3))) /
           (144. * sqrt(385) * params.r2);
  }

  else if (vpnorder == 5) {
    return (delta *
            (-360 * (33 - 197 * Nu + 294 * params.eta2) * params.r3 *
                 pow(PhiDOT * r + Complex(0, 1) * rDOT, 4) *
                 pow(Complex(0, 1) * PhiDOT * r + rDOT, 3) +
             2 * params.Mtot3 *
                 (Complex(0, -1) * (53311 - 121906 * Nu + 42816 * params.eta2) *
                      PhiDOT * r +
                  78 * (1141 - 2760 * Nu + 1420 * params.eta2) * rDOT) +
             2 * params.Mtot2 * r *
                 (Complex(0, 2) * (29938 - 82195 * Nu + 59238 * params.eta2) *
                      params.PhiDOT3 * params.r3 +
                  (-27026 + 120623 * Nu - 199326 * params.eta2) *
                      params.PhiDOT2 * params.r2 * rDOT +
                  Complex(0, 4) * (-2636 - 11335 * Nu + 43962 * params.eta2) *
                      PhiDOT * r * params.rDOT2 -
                  39 * (-639 - 735 * Nu + 5690 * params.eta2) * params.rDOT3) +
             3 * mass * params.r2 *
                 (Complex(0, -4) * (3115 - 15385 * Nu + 22098 * params.eta2) *
                      params.PhiDOT5 * params.r5 +
                  2 * (-2108 - 17893 * Nu + 56466 * params.eta2) *
                      params.PhiDOT4 * params.r4 * rDOT -
                  Complex(0, 3) * (-8473 - 9528 * Nu + 65204 * params.eta2) *
                      params.PhiDOT3 * params.r3 * params.rDOT2 +
                  8 * (-2692 - 6587 * Nu + 29754 * params.eta2) *
                      params.PhiDOT2 * params.r2 * params.rDOT3 -
                  Complex(0, 4) * (-3971 - 3226 * Nu + 27360 * params.eta2) *
                      PhiDOT * r * params.rDOT4 +
                  240 * (-57 - 97 * Nu + 518 * params.eta2) * params.rDOT5))) /
           (11232. * sqrt(385) * params.r3);
  }

  else if (vpnorder == 6) {
    return (
        -(delta * params.Mtot2 * Nu *
          (17830 * params.Mtot2 -
           6 * mass * r *
               (4723 * params.PhiDOT2 * params.r2 -
                Complex(0, 3921) * PhiDOT * r * rDOT + 4530 * params.rDOT2) +
           params.r2 * (12629 * params.PhiDOT4 * params.r4 -
                        Complex(0, 24248) * params.PhiDOT3 * params.r3 * rDOT +
                        20064 * params.PhiDOT2 * params.r2 * params.rDOT2 -
                        Complex(0, 9968) * PhiDOT * r * params.rDOT3 +
                        7120 * params.rDOT4))) /
            (5040. * sqrt(385) * params.r4)
        /* +((Complex(0,-0.0023148148148148147)*((5 + Nu*(-23 + 19*Nu)
  + delta*(5 + 27*(-1 + Nu)*Nu))*S1z + (-5 + (23 - 19*Nu)*Nu
  + delta*(5 + 27*(-1 + Nu)*Nu))*S2z)*params.x4)/sqrt(385)) */
        +
        /* Henry et al. ecc spin terms */ (
            (Complex(0, -0.0023148148148148147) * params.Mtot2 *
             (params.Mtot2 * ((90 - 622 * Nu + 966 * params.eta2 +
                               delta * (90 - 524 * Nu + 770 * params.eta2)) *
                                  S1z +
                              2 *
                                  (-45 + 311 * Nu - 483 * params.eta2 +
                                   delta * (45 - 262 * Nu + 385 * params.eta2)) *
                                  S2z) +
              3 * params.r2 *
                  (Complex(0, 4) * PhiDOT * r * params.rDOT3 *
                       (-5 * (3 - 20 * Nu + 30 * params.eta2) * S1z +
                        delta * (-15 - 106 * Nu + 262 * params.eta2) * S1z +
                        5 * (3 - 20 * Nu + 30 * params.eta2) * S2z +
                        delta * (-15 - 106 * Nu + 262 * params.eta2) * S2z) +
                   Complex(0, 1) * params.PhiDOT3 * params.r3 * rDOT *
                       (-35 * (3 - 22 * Nu + 36 * params.eta2) * S1z +
                        delta * (-105 - 494 * Nu + 1268 * params.eta2) * S1z +
                        35 * (3 - 22 * Nu + 36 * params.eta2) * S2z +
                        delta * (-105 - 494 * Nu + 1268 * params.eta2) * S2z) +
                   4 * params.rDOT4 *
                       ((15 - 92 * Nu + 126 * params.eta2 +
                         delta * (15 - 64 * Nu + 70 * params.eta2)) *
                            S1z +
                        (-15 + 92 * Nu - 126 * params.eta2 +
                         delta * (15 - 64 * Nu + 70 * params.eta2)) *
                            S2z) +
                   params.PhiDOT4 * params.r4 *
                       ((-180 + 1009 * Nu - 1227 * params.eta2 +
                         delta * (-180 + 95 * Nu + 601 * params.eta2)) *
                            S1z +
                        (180 - 1009 * Nu + 1227 * params.eta2 +
                         delta * (-180 + 95 * Nu + 601 * params.eta2)) *
                            S2z) +
                   3 * params.PhiDOT2 * params.r2 * params.rDOT2 *
                       ((35 - 198 * Nu + 244 * params.eta2 +
                         delta * (35 - 382 * Nu + 612 * params.eta2)) *
                            S1z +
                        (-35 + 198 * Nu - 244 * params.eta2 +
                         delta * (35 - 382 * Nu + 612 * params.eta2)) *
                            S2z)) +
              mass * r *
                  (Complex(0, -2) * PhiDOT * r * rDOT *
                       (-5 * (78 - 499 * Nu + 717 * params.eta2) * S1z +
                        delta * (-390 - 677 * Nu + 2759 * params.eta2) * S1z +
                        5 * (78 - 499 * Nu + 717 * params.eta2) * S2z +
                        delta * (-390 - 677 * Nu + 2759 * params.eta2) * S2z) +
                   params.PhiDOT2 * params.r2 *
                       ((609 - 3352 * Nu + 3966 * params.eta2 +
                         delta * (609 + 58 * Nu - 2854 * params.eta2)) *
                            S1z +
                        (-609 + 3352 * Nu - 3966 * params.eta2 +
                         delta * (609 + 58 * Nu - 2854 * params.eta2)) *
                            S2z) -
                   2 * params.rDOT2 *
                       ((255 - 1592 * Nu + 2226 * params.eta2 +
                         delta * (255 - 1144 * Nu + 1330 * params.eta2)) *
                            S1z +
                        (-255 + 1592 * Nu - 2226 * params.eta2 +
                         delta * (255 - 1144 * Nu + 1330 * params.eta2)) *
                            S2z)))) /
            (sqrt(385) * params.r4)));
  }

  else if (vpnorder == 7) {

    return (/* (Complex(0,0.001736111111111111)*(kappa1*(5 + 5*delta*Nu_factor
      - 2*Nu*(15 + 4*Nu))*params.S1z2 + S2z*(20*delta*(1 - 2*Nu)*Nu*S1z
      + kappa2*(-5 + 5*delta*Nu_factor + 30*Nu +
      8*params.eta2)*S2z))*params.x4p5)/sqrt(385)
      + */
            /* Henry et al. ecc+spin terms */
            ((3510 * params.Mtot3 *
                  (Complex(0, 2) * mass * (149 + 2 * Nu * (-447 + 508 * Nu)) *
                       PhiDOT * r -
                   2 * mass * (235 + 2 * Nu * (-705 + 1036 * Nu)) * rDOT +
                   r * (Complex(0, -7) * (35 - 210 * Nu + 232 * params.eta2) *
                            params.PhiDOT3 * params.r3 +
                        2 * (115 + 2 * Nu * (-345 + 628 * Nu)) *
                            params.PhiDOT2 * params.r2 * rDOT -
                        Complex(0, 2) * (85 - 510 * Nu + 872 * params.eta2) *
                            PhiDOT * r * params.rDOT2 +
                        8 * (-5 + 14 * Nu) * (-5 + 16 * Nu) * params.rDOT3)) *
                  (kappa1 * params.S1z2 - kappa2 * params.S2z2) +
              delta *
                  (5400 * (-135 + Nu * (1467 + 35 * Nu * (-151 + 178 * Nu))) *
                       params.r4 * pow(PhiDOT * r - Complex(0, 1) * rDOT, 4) *
                       pow(Complex(0, -1) * PhiDOT * r + rDOT, 5) +
                   6 * params.Mtot2 * params.r2 *
                       (Complex(0, -2) *
                            (556064 +
                             7 * Nu *
                                 (259253 + 130 * Nu * (-10297 + 8136 * Nu))) *
                            params.PhiDOT5 * params.r5 +
                        8 *
                            (349973 +
                             4 * Nu *
                                 (-94417 + 5 * Nu * (-65429 + 126618 * Nu))) *
                            params.PhiDOT4 * params.r4 * rDOT -
                        Complex(0, 60) *
                            (37369 + Nu * (-136791 - 60188 * Nu +
                                           607348 * params.eta2)) *
                            params.PhiDOT3 * params.r3 * params.rDOT2 +
                        2 *
                            (706852 +
                             Nu * (-3817342 +
                                   15 * Nu * (-114953 + 1628610 * Nu))) *
                            params.PhiDOT2 * params.r2 * params.rDOT3 -
                        Complex(0, 2) *
                            (-19909 + Nu * (-2880941 +
                                            15 * Nu * (195291 + 744862 * Nu))) *
                            PhiDOT * r * params.rDOT4 +
                        15 *
                            (1285 +
                             Nu * (-417939 + Nu * (376375 + 1794742 * Nu))) *
                            params.rDOT5) +
                   45 * mass * params.r3 *
                       (Complex(0, 4) *
                            (11270 +
                             Nu * (40035 + Nu * (-338333 + 490662 * Nu))) *
                            params.PhiDOT7 * params.r7 -
                        4 *
                            (12816 +
                             Nu * (1923 - 262738 * Nu + 582488 * params.eta2)) *
                            params.PhiDOT6 * params.r6 * rDOT +
                        Complex(0, 4) *
                            (76389 +
                             Nu * (-196893 + Nu * (-529811 + 1606402 * Nu))) *
                            params.PhiDOT5 * params.r5 * params.rDOT2 +
                        8 *
                            (-45746 +
                             Nu * (115572 + (313303 - 946222 * Nu) * Nu)) *
                            params.PhiDOT4 * params.r4 * params.rDOT3 +
                        Complex(0, 1) *
                            (464727 +
                             2 * Nu *
                                 (-792537 + 4 * Nu * (-164359 + 877280 * Nu))) *
                            params.PhiDOT3 * params.r3 * params.rDOT4 -
                        32 *
                            (13500 +
                             Nu * (-39687 + Nu * (-69661 + 249182 * Nu))) *
                            params.PhiDOT2 * params.r2 * params.rDOT5 +
                        Complex(0, 4) *
                            (41703 +
                             2 * Nu *
                                 (-74181 + 20 * Nu * (-2603 + 15676 * Nu))) *
                            PhiDOT * r * params.rDOT6 -
                        960 * (162 + Nu * (-516 + 7 * Nu * (-101 + 410 * Nu))) *
                            params.rDOT7) +
                   4 * params.Mtot4 *
                       (Complex(0, -1) * PhiDOT * r *
                            (-5713898 - 261495 * kappa1 * params.S1z2 +
                             2 * Nu *
                                 (6824514 - 4178445 * Nu +
                                  1430090 * params.eta2 -
                                  522990 * kappa1 * (-1 + Nu) * params.S1z2) +
                             1045980 * Nu * (-1 + 2 * Nu) * S1z * S2z -
                             261495 * kappa2 * Nu_factor *
                                 params.S2z2) +
                        15 * rDOT *
                            (-300238 - 27495 * kappa1 * params.S1z2 +
                             Nu * (581955 - 600300 * Nu + 535772 * params.eta2 -
                                   109980 * kappa1 * (-1 + Nu) * params.S1z2) +
                             109980 * Nu * (-1 + 2 * Nu) * S1z * S2z -
                             27495 * kappa2 * Nu_factor *
                                 params.S2z2)) +
                   6 * params.Mtot3 * r *
                       (Complex(0, 2) * PhiDOT * r * params.rDOT2 *
                            (1974194 - 49725 * kappa1 * params.S1z2 +
                             Nu * (-5423729 + 339235 * Nu +
                                   5947390 * params.eta2 -
                                   198900 * kappa1 * (-1 + Nu) * params.S1z2) +
                             198900 * Nu * (-1 + 2 * Nu) * S1z * S2z -
                             49725 * kappa2 * Nu_factor *
                                 params.S2z2) +
                        Complex(0, 25) * params.PhiDOT3 * params.r3 *
                            ((-265372 +
                              Nu * (1803657 +
                                    2 * Nu * (-1845967 + 746940 * Nu))) /
                                 5. -
                             5733 * kappa1 * Nu_factor * params.S1z2 +
                             22932 * Nu * (-1 + 2 * Nu) * S1z * S2z -
                             5733 * kappa2 * Nu_factor * params.S2z2) +
                        20 * params.rDOT3 *
                            (-175718 + 5850 * kappa1 * params.S1z2 +
                             Nu * (730551 - 504765 * Nu - 773906 * params.eta2 +
                                   23400 * kappa1 * (-1 + Nu) * params.S1z2) -
                             23400 * Nu * (-1 + 2 * Nu) * S1z * S2z +
                             5850 * kappa2 * Nu_factor * params.S2z2) +
                        2 * params.PhiDOT2 * params.r2 * rDOT *
                            (-2687884 + 67275 * kappa1 * params.S1z2 +
                             4 * Nu *
                                 (917051 + 1724565 * Nu - 1611110 * params.eta2 +
                                  67275 * kappa1 * (-1 + Nu) * params.S1z2) -
                             269100 * Nu * (-1 + 2 * Nu) * S1z * S2z +
                             67275 * kappa2 * Nu_factor *
                                 params.S2z2)))) /
             (673920. * sqrt(385) * params.r4)));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hQC_5_m_1(REAL8 Nu, UINT4 vpnorder, REAL8 x, struct kepler_vars params) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  /* if(vpnorder == 3){
       return(((Complex(0,0.006944444444444444)*delta)/sqrt(385) -
       (Complex(0,0.013888888888888888)*delta*Nu)/sqrt(385))*params.x2p5);
   }
   else if(vpnorder == 5){
       return(((Complex(0,-0.03187321937321937)*delta)/sqrt(385) +
        Complex(0,0.005698005698005698)*sqrt(0.3142857142857143)*delta*
        Nu - (Complex(0,0.0007122507122507123)*delta*params.eta2)/
        sqrt(385))*params.x3p5);
   } */
  /* else */ if (vpnorder == 6) {
    return ((Complex(0, -0.006944444444444444) * delta * (-1 + 2 * Nu) *
             params.x4 * (M_PI - Complex(0, 2) * log(2))) /
            sqrt(385));
  } else {
    return 0;
  }
}

static COMPLEX16 hl_5_m_1(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z, REAL8 x, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_5_m_1: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_5_m_1(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, params) +
            hQC_5_m_1(Nu, vpnorder, x, params)) *
           cpolar(1, -1 * Phi);
  }
}

static COMPLEX16 hl_5_m_min1(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z, REAL8 x, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_5_m_min1: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((- 4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           conj(hGO_5_m_1(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, params) +
                hQC_5_m_1(Nu, vpnorder, x, params)) *
           cpolar(1, 1 * Phi);
  }
}

// 66
static COMPLEX16 hGO_6_m_6(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z, struct kepler_vars params) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  if (vpnorder == 4) {
    return (((1 - 5 * Nu + 5 * params.eta2) *
             (172 * params.Mtot3 +
              120 * params.r3 * pow(PhiDOT * r + Complex(0, 1) * rDOT, 6) +
              params.Mtot2 * r *
                  (3269 * params.PhiDOT2 * params.r2 +
                   Complex(0, 2920) * PhiDOT * r * rDOT - 806 * params.rDOT2) +
              15 * mass * params.r2 *
                  (281 * params.PhiDOT4 * params.r4 +
                   Complex(0, 494) * params.PhiDOT3 * params.r3 * rDOT -
                   444 * params.PhiDOT2 * params.r2 * params.rDOT2 -
                   Complex(0, 208) * PhiDOT * r * params.rDOT3 +
                   40 * params.rDOT4))) /
            (360. * sqrt(143) * params.r3));
  }

  else if (vpnorder == 6) {
    return (
        (14 * params.Mtot4 *
             (-5740 + 29361 * Nu - 33348 * params.eta2 + 7334 * params.eta3) -
         7560 * (-1 + 9 * Nu - 26 * params.eta2 + 23 * params.eta3) * params.r4 *
             (PhiDOT * r - Complex(0, 1) * rDOT) *
             pow(PhiDOT * r + Complex(0, 1) * rDOT, 7) +
         2 * params.Mtot3 * r *
             ((-539645 + 2950311 * Nu - 4086684 * params.eta2 +
               1644517 * params.eta3) *
                  params.PhiDOT2 * params.r2 +
              Complex(0, 6) *
                  (-47984 + 275121 * Nu - 442540 * params.eta2 +
                   255850 * params.eta3) *
                  PhiDOT * r * rDOT +
              14 *
                  (3614 - 21621 * Nu + 39684 * params.eta2 -
                   29332 * params.eta3) *
                  params.rDOT2) +
         3 * params.Mtot2 * params.r2 *
             ((-311847 + 1966993 * Nu - 3502751 * params.eta2 +
               1752968 * params.eta3) *
                  params.PhiDOT4 * params.r4 +
              Complex(0, 4) *
                  (629 + 160412 * Nu - 846370 * params.eta2 +
                   912975 * params.eta3) *
                  params.PhiDOT3 * params.r3 * rDOT -
              3 *
                  (65519 - 144403 * Nu - 684796 * params.eta2 +
                   1205253 * params.eta3) *
                  params.PhiDOT2 * params.r2 * params.rDOT2 -
              Complex(0, 32) *
                  (3921 - 11389 * Nu - 27265 * params.eta2 +
                   58450 * params.eta3) *
                  PhiDOT * r * params.rDOT3 +
              14 *
                  (1867 - 5501 * Nu - 12824 * params.eta2 + 28137 * params.eta3) *
                  params.rDOT4) -
         45 * mass * params.r3 *
             ((195 + 3619 * Nu - 36617 * params.eta2 + 66836 * params.eta3) *
                  params.PhiDOT6 * params.r6 +
              Complex(0, 4) *
                  (-1878 + 10969 * Nu - 20741 * params.eta2 +
                   18263 * params.eta3) *
                  params.PhiDOT5 * params.r5 * rDOT +
              (17169 - 75446 * Nu + 35497 * params.eta2 + 47054 * params.eta3) *
                  params.PhiDOT4 * params.r4 * params.rDOT2 +
              Complex(0, 2) *
                  (9183 - 30296 * Nu - 37835 * params.eta2 +
                   95060 * params.eta3) *
                  params.PhiDOT3 * params.r3 * params.rDOT3 -
              4 * (2781 - 6062 * Nu - 28595 * params.eta2 + 49070 * params.eta3) *
                  params.PhiDOT2 * params.r2 * params.rDOT4 -
              Complex(0, 16) *
                  (228 - 217 * Nu - 3871 * params.eta2 + 5803 * params.eta3) *
                  PhiDOT * r * params.rDOT5 +
              56 * (9 + 4 * Nu - 221 * params.eta2 + 308 * params.eta3) *
                  params.rDOT6)) /
        (15120. * sqrt(143) * params.r4));
  }

  else if (vpnorder == 7) {

    return (/* (108*(110*params.eta3*(S1z + S2z) - 14*(S1z + delta*S1z + S2z -
      delta*S2z)
      - 50*params.eta2*(2*delta*(S1z - S2z) + 3*(S1z + S2z)) +  Nu*(85*delta*(S1z
      - S2z)
      + 83*(S1z + S2z)))*params.x4p5)/(35.*sqrt(143))
      + */
            /* Henry et al. ecc+spin terms */
            ((params.Mtot2 *
              (8 * params.Mtot2 *
                   (30 * params.eta3 *
                        (9071 * PhiDOT * r + Complex(0, 4976) * rDOT) *
                        (S1z + S2z) -
                    21 * (761 * PhiDOT * r + Complex(0, 476) * rDOT) *
                        (S1z + delta * S1z + S2z - delta * S2z) +
                    Nu * (PhiDOT * r *
                              (Complex(0, 61450) + 125439 * S1z +
                               94725 * delta * S1z + 125439 * S2z -
                               94725 * delta * S2z) +
                          Complex(0, 1) * rDOT *
                              (Complex(0, 14017) + 77370 * S1z +
                               54090 * delta * S1z + 77370 * S2z -
                               54090 * delta * S2z)) -
                    3 * params.eta2 *
                        (Complex(0, 1) * rDOT *
                             (Complex(0, 15107) +
                              40 * (1592 + 485 * delta) * S1z + 63680 * S2z -
                              19400 * delta * S2z) +
                         5 * PhiDOT * r *
                             (Complex(0, 13198) + (21493 + 7303 * delta) * S1z +
                              21493 * S2z - 7303 * delta * S2z))) -
               3 * params.r2 *
                   (-900 * params.eta3 *
                        (975 * params.PhiDOT5 * params.r5 +
                         Complex(0, 2382) * params.PhiDOT4 * params.r4 * rDOT -
                         2208 * params.PhiDOT3 * params.r3 * params.rDOT2 -
                         Complex(0, 1592) * params.PhiDOT2 * params.r2 *
                             params.rDOT3 +
                         688 * PhiDOT * r * params.rDOT4 +
                         Complex(0, 128) * params.rDOT5) *
                        (S1z + S2z) +
                    1260 *
                        (61 * params.PhiDOT5 * params.r5 +
                         Complex(0, 151) * params.PhiDOT4 * params.r4 * rDOT -
                         172 * params.PhiDOT3 * params.r3 * params.rDOT2 -
                         Complex(0, 122) * params.PhiDOT2 * params.r2 *
                             params.rDOT3 +
                         48 * PhiDOT * r * params.rDOT4 +
                         Complex(0, 8) * params.rDOT5) *
                        (S1z + delta * S1z + S2z - delta * S2z) -
                    3 * Nu *
                        (Complex(0, 1) * params.PhiDOT4 * params.r4 * rDOT *
                             (Complex(0, 56511) +
                              10 * (43269 + 36265 * delta) * S1z +
                              432690 * S2z - 362650 * delta * S2z) -
                         Complex(0, 4) * params.PhiDOT2 * params.r2 *
                             params.rDOT3 *
                             (Complex(0, 14259) +
                              40 * (2154 + 1675 * delta) * S1z + 86160 * S2z -
                              67000 * delta * S2z) +
                         4 * params.PhiDOT5 * params.r5 *
                             (Complex(0, 94823) +
                              25 * (1494 + 1901 * delta) * S1z + 37350 * S2z -
                              47525 * delta * S2z) -
                         16 * params.PhiDOT3 * params.r3 * params.rDOT2 *
                             (Complex(0, 6984) +
                              5 * (5967 + 4855 * delta) * S1z + 29835 * S2z -
                              24275 * delta * S2z) +
                         Complex(0, 8) * params.rDOT5 *
                             (Complex(0, 251) + 10 * (303 + 215 * delta) * S1z +
                              3030 * S2z - 2150 * delta * S2z) +
                         64 * PhiDOT * r * params.rDOT4 *
                             (Complex(0, 251) + 5 * (438 + 325 * delta) * S1z +
                              2190 * S2z - 1625 * delta * S2z)) +
                    params.eta2 *
                        (Complex(0, 8) * params.rDOT5 *
                             (Complex(0, 2449) + 600 * (34 + 11 * delta) * S1z -
                              600 * (-34 + 11 * delta) * S2z) +
                         64 * PhiDOT * r * params.rDOT4 *
                             (Complex(0, 2449) + 75 * (188 + 67 * delta) * S1z -
                              75 * (-188 + 67 * delta) * S2z) -
                         16 * params.PhiDOT3 * params.r3 * params.rDOT2 *
                             (Complex(0, 66866) +
                              75 * (2423 + 1039 * delta) * S1z -
                              75 * (-2423 + 1039 * delta) * S2z) -
                         Complex(0, 4) * params.PhiDOT2 * params.r2 *
                             params.rDOT3 *
                             (Complex(0, 139541) +
                              150 * (3551 + 1399 * delta) * S1z -
                              150 * (-3551 + 1399 * delta) * S2z) +
                         Complex(0, 1) * params.PhiDOT4 * params.r4 * rDOT *
                             (Complex(0, 515089) +
                              600 * (4703 + 2041 * delta) * S1z -
                              600 * (-4703 + 2041 * delta) * S2z) +
                         4 * params.PhiDOT5 * params.r5 *
                             (Complex(0, 925502) +
                              75 * (2966 + 2521 * delta) * S1z -
                              75 * (-2966 + 2521 * delta) * S2z))) +
               12 * mass * r *
                   (30 * params.eta3 *
                        (20149 * params.PhiDOT3 * params.r3 +
                         Complex(0, 33176) * params.PhiDOT2 * params.r2 * rDOT -
                         18519 * PhiDOT * r * params.rDOT2 -
                         Complex(0, 4092) * params.rDOT3) *
                        (S1z + S2z) -
                    42 *
                        (1167 * params.PhiDOT3 * params.r3 +
                         Complex(0, 2042) * params.PhiDOT2 * params.r2 * rDOT -
                         1131 * PhiDOT * r * params.rDOT2 -
                         Complex(0, 236) * params.rDOT3) *
                        (S1z + delta * S1z + S2z - delta * S2z) +
                    params.eta2 *
                        (-2 * params.PhiDOT3 * params.r3 *
                             (Complex(0, 452503) +
                              5 * (72902 + 36559 * delta) * S1z + 364510 * S2z -
                              182795 * delta * S2z) +
                         3 * PhiDOT * r * params.rDOT2 *
                             (Complex(0, 53821) +
                              30 * (8510 + 2941 * delta) * S1z + 255300 * S2z -
                              88230 * delta * S2z) -
                         Complex(0, 6) * params.PhiDOT2 * params.r2 * rDOT *
                             (Complex(0, 65461) +
                              30 * (7463 + 2831 * delta) * S1z + 223890 * S2z -
                              84930 * delta * S2z) +
                         Complex(0, 1) * params.rDOT3 *
                             (Complex(0, 26693) +
                              20 * (8419 + 2675 * delta) * S1z + 168380 * S2z -
                              53500 * delta * S2z)) +
                    Nu * (params.PhiDOT3 * params.r3 *
                              (Complex(0, 280937) +
                               (329808 + 305330 * delta) * S1z +
                               (329808 - 305330 * delta) * S2z) -
                          Complex(0, 1) * params.rDOT3 *
                              (Complex(0, 8281) +
                               10 * (7293 + 5153 * delta) * S1z +
                               (72930 - 51530 * delta) * S2z) -
                          3 * PhiDOT * r * params.rDOT2 *
                              (Complex(0, 16877) +
                               90 * (1261 + 930 * delta) * S1z -
                               90 * (-1261 + 930 * delta) * S2z) +
                          Complex(0, 2) * params.PhiDOT2 * params.r2 * rDOT *
                              (Complex(0, 60821) +
                               24 * (12576 + 9775 * delta) * S1z -
                               24 * (-12576 + 9775 * delta) * S2z))))) /
             (30240. * sqrt(143) * params.r4)));
  }

  else {
    return 0;
  }
}

/* static COMPLEX16 hQC_6_m_6(REAL8 Nu, UINT4 vpnorder, REAL8 x, struct kepler_vars params){

    if(vpnorder == 4){
        return((108/(5.*sqrt(143)) - (108*Nu)/sqrt(143) +
(108*params.eta2)/sqrt(143))* params.x3);
    }
    else if(vpnorder == 6){
        return((-6102/(35.*sqrt(143)) + (378*sqrt(1.1818181818181819)*Nu)/5. -
        (6912*params.eta2)/(5.*sqrt(143)) +
        (162*sqrt(1.1818181818181819)*params.eta3)/5.)*params.x4);
    }
    else{
        return 0;
    }
} */

static COMPLEX16 hl_6_m_6(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_6_m_6: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_6_m_6(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z,
                      S2z, params) /* +hQC_6_m_6(Nu,vpnorder,x,params) */) *
           cpolar(1, -6 * Phi);
  }
}

static COMPLEX16 hl_6_m_min6(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_6_m_min6: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           conj(hGO_6_m_6(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z,
                          S2z, params) /* +hQC_6_m_6(Nu,vpnorder,x,params) */) *
           cpolar(1, 6 * Phi);
  }
}

// 65
static COMPLEX16 hGO_6_m_5(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z, struct kepler_vars params) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  if (vpnorder == 5) {
    return (Complex(0, 0.003968253968253968) * delta * mass *
            (1 - 4 * Nu + 3 * params.eta2) * PhiDOT *
            (82 * params.Mtot2 +
             2 * mass * r *
                 (701 * params.PhiDOT2 * params.r2 +
                  Complex(0, 343) * PhiDOT * r * rDOT - 62 * params.rDOT2) +
             3 * params.r2 *
                 (547 * params.PhiDOT4 * params.r4 +
                  Complex(0, 364) * params.PhiDOT3 * params.r3 * rDOT -
                  180 * params.PhiDOT2 * params.r2 * params.rDOT2 -
                  Complex(0, 56) * PhiDOT * r * params.rDOT3 +
                  8 * params.rDOT4))) /
           (sqrt(429) * params.r2);
  }

  else if (vpnorder == 6) {

    return (/* (Complex(0,-21.70138888888889)*Nu*(S1z + delta*(-1 + Nu)*S1z
      - 3*Nu*S1z - (1 + delta)*S2z + (3 + delta)*Nu*S2z)*params.x4)/sqrt(429)
      + */
            /* Henry et al. ecc spin terms */ (
                (Complex(0, -0.006944444444444444) * params.Mtot2 * Nu *
                 (82 * params.Mtot2 +
                  2 * mass * r *
                      (701 * params.PhiDOT2 * params.r2 +
                       Complex(0, 343) * PhiDOT * r * rDOT -
                       62 * params.rDOT2) +
                  3 * params.r2 *
                      (547 * params.PhiDOT4 * params.r4 +
                       Complex(0, 364) * params.PhiDOT3 * params.r3 * rDOT -
                       180 * params.PhiDOT2 * params.r2 * params.rDOT2 -
                       Complex(0, 56) * PhiDOT * r * params.rDOT3 +
                       8 * params.rDOT4)) *
                 ((1 + delta * (-1 + Nu) - 3 * Nu) * S1z +
                  (-1 + delta * (-1 + Nu) + 3 * Nu) * S2z)) /
                (sqrt(429) * params.r4)));
  }

  else if (vpnorder == 7) {

    return (
        (Complex(0, 0.00003306878306878307) * delta * mass * PhiDOT *
         (16 * params.Mtot3 *
              (-3893 + 16386 * Nu - 17175 * params.eta2 + 6922 * params.eta3) +
          24 * params.Mtot2 * r *
              ((-25472 + 127239 * Nu - 198605 * params.eta2 +
                117623 * params.eta3) *
                   params.PhiDOT2 * params.r2 +
               Complex(0, 10) *
                   (1954 - 7090 * Nu + 1431 * params.eta2 + 5232 * params.eta3) *
                   PhiDOT * r * rDOT -
               (6603 - 26446 * Nu + 16325 * params.eta2 + 7138 * params.eta3) *
                   params.rDOT2) +
          12 * mass * params.r2 *
              (2 *
                   (-23971 + 144172 * Nu - 256320 * params.eta2 +
                    127374 * params.eta3) *
                   params.PhiDOT4 * params.r4 +
               Complex(0, 5) *
                   (40015 - 154144 * Nu + 84408 * params.eta2 +
                    41694 * params.eta3) *
                   params.PhiDOT3 * params.r3 * rDOT -
               3 *
                   (46471 - 183392 * Nu + 109660 * params.eta2 +
                    47046 * params.eta3) *
                   params.PhiDOT2 * params.r2 * params.rDOT2 -
               Complex(0, 5) *
                   (9713 - 38048 * Nu + 21348 * params.eta2 +
                    11562 * params.eta3) *
                   PhiDOT * r * params.rDOT3 +
               4 * (1811 - 7007 * Nu + 3585 * params.eta2 + 2511 * params.eta3) *
                   params.rDOT4) -
          45 * params.r3 *
              ((3177 - 1626 * Nu - 54862 * params.eta2 + 73376 * params.eta3) *
                   params.PhiDOT6 * params.r6 +
               Complex(0, 2) *
                   (-5337 + 26646 * Nu - 34862 * params.eta2 +
                    11212 * params.eta3) *
                   params.PhiDOT5 * params.r5 * rDOT +
               4 * (1917 - 7098 * Nu + 1187 * params.eta2 + 6278 * params.eta3) *
                   params.PhiDOT4 * params.r4 * params.rDOT2 +
               Complex(0, 8) *
                   (639 - 1986 * Nu - 1156 * params.eta2 + 3296 * params.eta3) *
                   params.PhiDOT3 * params.r3 * params.rDOT3 -
               16 * (153 - 402 * Nu - 577 * params.eta2 + 1022 * params.eta3) *
                   params.PhiDOT2 * params.r2 * params.rDOT4 -
               Complex(0, 16) *
                   (45 - 102 * Nu - 236 * params.eta2 + 352 * params.eta3) *
                   PhiDOT * r * params.rDOT5 +
               32 * (3 - 6 * Nu - 19 * params.eta2 + 26 * params.eta3) *
                   params.rDOT6))) /
        (sqrt(429) * params.r4));
  }

  else {
    return 0;
  }
}

/* static COMPLEX16 hQC_6_m_5(REAL8 Nu, UINT4 vpnorder, REAL8 x, struct kepler_vars params){
    REAL8 delta = sqrt(1-4*Nu);

    if(vpnorder == 5){
        return(((Complex(0,12.40079365079365)*delta)/sqrt(429) -
        (Complex(0,49.6031746031746)*delta*Nu)/sqrt(429) +
        (Complex(0,37.20238095238095)*delta*params.eta2)/sqrt(429))*
        params.x3p5);
    }
    else{
        return 0;
    }
} */

static COMPLEX16 hl_6_m_5(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_6_m_m5: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_6_m_5(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z,
                      S2z, params) /* +hQC_6_m_5(Nu,vpnorder,x, params) */) *
           cpolar(1, -5 * Phi);
  }
}

static COMPLEX16 hl_6_m_min5(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_6_m_min5: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           conj(hGO_6_m_5(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z,
                          S2z, params) /* +hQC_6_m_5(Nu,vpnorder,x, params) */) *
           cpolar(1, 5 * Phi);
  }
}

// 64
static COMPLEX16 hGO_6_m_4(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z, struct kepler_vars params) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  if (vpnorder == 4) {
    return (((1 - 5 * Nu + 5 * params.eta2) *
             (-516 * params.Mtot3 +
              360 * params.r3 * (PhiDOT * r - Complex(0, 1) * rDOT) *
                  pow(PhiDOT * r + Complex(0, 1) * rDOT, 5) +
              params.Mtot2 * r *
                  (-3587 * params.PhiDOT2 * params.r2 -
                   Complex(0, 5840) * PhiDOT * r * rDOT + 2418 * params.rDOT2) +
              15 * mass * params.r2 *
                  (113 * params.PhiDOT4 * params.r4 -
                   Complex(0, 108) * params.PhiDOT3 * params.r3 * rDOT +
                   468 * params.PhiDOT2 * params.r2 * params.rDOT2 +
                   Complex(0, 416) * PhiDOT * r * params.rDOT3 -
                   120 * params.rDOT4))) /
            (1980. * sqrt(78) * params.r3));
  }

  else if (vpnorder == 6) {
    return (
        -(14 * params.Mtot4 *
              (-5740 + 29361 * Nu - 33348 * params.eta2 + 7334 * params.eta3) +
          7560 * (-1 + 9 * Nu - 26 * params.eta2 + 23 * params.eta3) * params.r4 *
              pow(PhiDOT * r - Complex(0, 1) * rDOT, 2) *
              pow(PhiDOT * r + Complex(0, 1) * rDOT, 6) +
          2 * params.Mtot3 * r *
              ((-196625 + 1082991 * Nu - 1522164 * params.eta2 +
                618457 * params.eta3) *
                   params.PhiDOT2 * params.r2 +
               Complex(0, 4) *
                   (-47984 + 275121 * Nu - 442540 * params.eta2 +
                    255850 * params.eta3) *
                   PhiDOT * r * rDOT +
               14 *
                   (3614 - 21621 * Nu + 39684 * params.eta2 -
                    29332 * params.eta3) *
                   params.rDOT2) +
          params.Mtot2 * params.r2 *
              ((133599 - 779681 * Nu + 1417087 * params.eta2 -
                1130416 * params.eta3) *
                   params.PhiDOT4 * params.r4 +
               Complex(0, 8) *
                   (3849 + 4172 * Nu - 80290 * params.eta2 +
                    64435 * params.eta3) *
                   params.PhiDOT3 * params.r3 * rDOT +
               (-226971 + 551047 * Nu + 2049124 * params.eta2 -
                3713857 * params.eta3) *
                   params.PhiDOT2 * params.r2 * params.rDOT2 -
               Complex(0, 64) *
                   (3921 - 11389 * Nu - 27265 * params.eta2 +
                    58450 * params.eta3) *
                   PhiDOT * r * params.rDOT3 +
               42 *
                   (1867 - 5501 * Nu - 12824 * params.eta2 +
                    28137 * params.eta3) *
                   params.rDOT4) +
          15 * mass * params.r3 *
              ((2267 - 12733 * Nu + 13895 * params.eta2 + 6300 * params.eta3) *
                   params.PhiDOT6 * params.r6 -
               Complex(0, 8) *
                   (908 - 2597 * Nu - 5873 * params.eta2 + 11809 * params.eta3) *
                   params.PhiDOT5 * params.r5 * rDOT +
               (5241 + 10066 * Nu - 173159 * params.eta2 + 235382 * params.eta3) *
                   params.PhiDOT4 * params.r4 * params.rDOT2 +
               Complex(0, 4) *
                   (-1651 + 11312 * Nu - 25417 * params.eta2 +
                    20916 * params.eta3) *
                   params.PhiDOT3 * params.r3 * params.rDOT3 +
               4 *
                   (3127 - 8386 * Nu - 23569 * params.eta2 +
                    45122 * params.eta3) *
                   params.PhiDOT2 * params.r2 * params.rDOT4 +
               Complex(0, 32) *
                   (228 - 217 * Nu - 3871 * params.eta2 + 5803 * params.eta3) *
                   PhiDOT * r * params.rDOT5 -
               168 * (9 + 4 * Nu - 221 * params.eta2 + 308 * params.eta3) *
                   params.rDOT6)) /
        (27720. * sqrt(78) * params.r4));
  }

  else if (vpnorder == 7) {

    return (/* (256*sqrt(0.05128205128205128)*(-150*params.eta3*(S1z + S2z)
      + 14*(S1z + delta*S1z + S2z - delta*S2z) + 10*params.eta2*(6*delta*(S1z -
      S2z) + 23*(S1z + S2z)) - Nu*(65*delta*(S1z - S2z) + 103*(S1z +
      S2z)))*params.x4p5)/3465.
      + */
            /* Henry et al. ecc+spin terms */
            ((params.Mtot2 *
              (-8 * params.Mtot2 *
                   (30 * params.eta3 *
                        (10641 * PhiDOT * r + Complex(0, 9952) * rDOT) *
                        (S1z + S2z) -
                    21 * (1283 * PhiDOT * r + Complex(0, 952) * rDOT) *
                        (S1z + delta * S1z + S2z - delta * S2z) +
                    Nu * (PhiDOT * r *
                              (Complex(0, 122900) +
                               3 * (62389 + 51235 * delta) * S1z +
                               187167 * S2z - 153705 * delta * S2z) +
                          Complex(0, 3) * rDOT *
                              (Complex(0, 14017) +
                               20 * (2579 + 1803 * delta) * S1z + 51580 * S2z -
                               36060 * delta * S2z)) -
                    3 * params.eta2 *
                        (5 * PhiDOT * r *
                             (Complex(0, 26396) +
                              (27731 + 11513 * delta) * S1z +
                              (27731 - 11513 * delta) * S2z) +
                         Complex(0, 1) * rDOT *
                             (Complex(0, 45321) +
                              80 * (1592 + 485 * delta) * S1z -
                              80 * (-1592 + 485 * delta) * S2z))) -
               3 * params.r2 *
                   (-60 * params.eta3 *
                        (401 * params.PhiDOT5 * params.r5 +
                         Complex(0, 30572) * params.PhiDOT4 * params.r4 * rDOT -
                         18000 * params.PhiDOT3 * params.r3 * params.rDOT2 +
                         Complex(0, 3632) * params.PhiDOT2 * params.r2 *
                             params.rDOT3 -
                         10736 * PhiDOT * r * params.rDOT4 -
                         Complex(0, 3840) * params.rDOT5) *
                        (S1z + S2z) +
                    1260 *
                        (17 * params.PhiDOT5 * params.r5 +
                         Complex(0, 58) * params.PhiDOT4 * params.r4 * rDOT +
                         16 * params.PhiDOT3 * params.r3 * params.rDOT2 +
                         Complex(0, 84) * params.PhiDOT2 * params.r2 *
                             params.rDOT3 -
                         64 * PhiDOT * r * params.rDOT4 -
                         Complex(0, 16) * params.rDOT5) *
                        (S1z + delta * S1z + S2z - delta * S2z) -
                    3 * Nu *
                        (Complex(0, 20) * params.PhiDOT2 * params.r2 *
                             params.rDOT3 *
                             (Complex(0, 3109) + (9208 + 9384 * delta) * S1z +
                              (9208 - 9384 * delta) * S2z) -
                         32 * PhiDOT * r * params.rDOT4 *
                             (Complex(0, 1004) +
                              5 * (1091 + 869 * delta) * S1z +
                              (5455 - 4345 * delta) * S2z) -
                         Complex(0, 8) * params.rDOT5 *
                             (Complex(0, 753) + 20 * (303 + 215 * delta) * S1z +
                              (6060 - 4300 * delta) * S2z) +
                         20 * params.PhiDOT5 * params.r5 *
                             (Complex(0, 4690) + (652 + 3807 * delta) * S1z +
                              (652 - 3807 * delta) * S2z) +
                         16 * params.PhiDOT3 * params.r3 * params.rDOT2 *
                             (Complex(0, 1872) +
                              5 * (-591 + 605 * delta) * S1z -
                              5 * (591 + 605 * delta) * S2z) +
                         Complex(0, 1) * params.PhiDOT4 * params.r4 * rDOT *
                             (Complex(0, 78819) +
                              20 * (12535 + 5539 * delta) * S1z -
                              20 * (-12535 + 5539 * delta) * S2z)) +
                    params.eta2 *
                        (Complex(0, 20) * params.PhiDOT2 * params.r2 *
                             params.rDOT3 *
                             (Complex(0, 30451) +
                              12 * (2831 + 2487 * delta) * S1z +
                              (33972 - 29844 * delta) * S2z) -
                         64 * PhiDOT * r * params.rDOT4 *
                             (Complex(0, 4898) +
                              15 * (1062 + 449 * delta) * S1z +
                              (15930 - 6735 * delta) * S2z) -
                         Complex(0, 24) * params.rDOT5 *
                             (Complex(0, 2449) + 400 * (34 + 11 * delta) * S1z -
                              400 * (-34 + 11 * delta) * S2z) +
                         96 * params.PhiDOT3 * params.r3 * params.rDOT2 *
                             (Complex(0, 2838) +
                              25 * (-445 + 79 * delta) * S1z -
                              25 * (445 + 79 * delta) * S2z) +
                         Complex(0, 1) * params.PhiDOT4 * params.r4 * rDOT *
                             (Complex(0, 761581) +
                              240 * (9441 + 1247 * delta) * S1z -
                              240 * (-9441 + 1247 * delta) * S2z) +
                         20 * params.PhiDOT5 * params.r5 *
                             (Complex(0, 44900) +
                              3 * (-1858 + 5829 * delta) * S1z -
                              3 * (1858 + 5829 * delta) * S2z))) +
               12 * mass * r *
                   (30 * params.eta3 *
                        (73 * params.PhiDOT3 * params.r3 -
                         Complex(0, 11024) * params.PhiDOT2 * params.r2 * rDOT +
                         20149 * PhiDOT * r * params.rDOT2 +
                         Complex(0, 8184) * params.rDOT3) *
                        (S1z + S2z) +
                    42 *
                        (871 * params.PhiDOT3 * params.r3 +
                         Complex(0, 1684) * params.PhiDOT2 * params.r2 * rDOT -
                         1543 * PhiDOT * r * params.rDOT2 -
                         Complex(0, 472) * params.rDOT3) *
                        (S1z + delta * S1z + S2z - delta * S2z) +
                    Nu * (-2 * params.PhiDOT3 * params.r3 *
                              (Complex(0, 69725) +
                               (89007 + 95170 * delta) * S1z +
                               (89007 - 95170 * delta) * S2z) +
                          Complex(0, 1) * params.rDOT3 *
                              (Complex(0, 24843) +
                               20 * (7293 + 5153 * delta) * S1z -
                               20 * (-7293 + 5153 * delta) * S2z) -
                          Complex(0, 2) * params.PhiDOT2 * params.r2 * rDOT *
                              (Complex(0, 66847) +
                               8 * (24951 + 24245 * delta) * S1z -
                               8 * (-24951 + 24245 * delta) * S2z) +
                          2 * PhiDOT * r * params.rDOT2 *
                              (Complex(0, 50631) +
                               5 * (43338 + 34327 * delta) * S1z -
                               5 * (-43338 + 34327 * delta) * S2z)) +
                    params.eta2 *
                        (10 * params.PhiDOT3 * params.r3 *
                             (Complex(0, 44862) +
                              (16586 + 19777 * delta) * S1z +
                              (16586 - 19777 * delta) * S2z) +
                         Complex(0, 2) * params.PhiDOT2 * params.r2 * rDOT *
                             (Complex(0, 216321) +
                              100 * (3079 + 2111 * delta) * S1z -
                              100 * (-3079 + 2111 * delta) * S2z) -
                         Complex(0, 1) * params.rDOT3 *
                             (Complex(0, 80079) +
                              40 * (8419 + 2675 * delta) * S1z -
                              40 * (-8419 + 2675 * delta) * S2z) -
                         2 * PhiDOT * r * params.rDOT2 *
                             (Complex(0, 161463) +
                              5 * (89002 + 36251 * delta) * S1z -
                              5 * (-89002 + 36251 * delta) * S2z))))) /
             (166320. * sqrt(78) * params.r4)));
  }

  else {
    return 0;
  }
}

/* static COMPLEX16 hQC_6_m_4(REAL8 Nu, UINT4 vpnorder, REAL8 x, struct kepler_vars params){

    if(vpnorder == 4){
        return(((-256*sqrt(0.05128205128205128))/495. +
        (256*sqrt(0.05128205128205128)*Nu)/99. -
        (256*sqrt(0.05128205128205128)*params.eta2)/99.)*params.x3);
    }
    else if(vpnorder == 6){
        return(((3968*sqrt(0.05128205128205128))/1155. -
        (9088*sqrt(0.05128205128205128)*Nu)/495. +
        (1024*sqrt(0.05128205128205128)*params.eta2)/45. -
        (2432*sqrt(0.05128205128205128)*params.eta3)/495.)*params.x4);
    }
    else{
        return 0;
    }

} */

static COMPLEX16 hl_6_m_4(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_6_m_4: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_6_m_4(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z,
                      S2z, params) /* +hQC_6_m_4(Nu,vpnorder,x, params) */) *
           cpolar(1, -4 * Phi);
  }
}

static COMPLEX16 hl_6_m_min4(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_6_m_min4: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           conj(hGO_6_m_4(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z,
                          S2z, params) /* +hQC_6_m_4(Nu,vpnorder,x, params) */) *
           cpolar(1, 4 * Phi);
  }
}

// 63
static COMPLEX16 hGO_6_m_3(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z, struct kepler_vars params) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  if (vpnorder == 5) {
    return (Complex(0, -0.00036075036075036075) * delta * mass *
            (1 - 4 * Nu + 3 * params.eta2) * PhiDOT *
            (410 * params.Mtot2 +
             2 * mass * r *
                 (929 * params.PhiDOT2 * params.r2 +
                  Complex(0, 1029) * PhiDOT * r * rDOT - 310 * params.rDOT2) -
             3 * params.r2 *
                 (513 * params.PhiDOT4 * params.r4 +
                  Complex(0, 28) * params.PhiDOT3 * params.r3 * rDOT +
                  228 * params.PhiDOT2 * params.r2 * params.rDOT2 +
                  Complex(0, 168) * PhiDOT * r * params.rDOT3 -
                  40 * params.rDOT4))) /
           (sqrt(65) * params.r2);
  }

  else if (vpnorder == 6) {

    return (/* (Complex(0,0.4602272727272727)*Nu*(S1z + delta*(-1 + Nu)*S1z -
      3*Nu*S1z
      - (1 + delta)*S2z + (3 + delta)*Nu*S2z)*params.x4)/sqrt(65)
      + */
            /* Henry et al. ecc spin terms */ (
                (Complex(0, 0.0006313131313131314) * params.Mtot2 * Nu *
                 (410 * params.Mtot2 +
                  2 * mass * r *
                      (929 * params.PhiDOT2 * params.r2 +
                       Complex(0, 1029) * PhiDOT * r * rDOT -
                       310 * params.rDOT2) -
                  3 * params.r2 *
                      (513 * params.PhiDOT4 * params.r4 +
                       Complex(0, 28) * params.PhiDOT3 * params.r3 * rDOT +
                       228 * params.PhiDOT2 * params.r2 * params.rDOT2 +
                       Complex(0, 168) * PhiDOT * r * params.rDOT3 -
                       40 * params.rDOT4)) *
                 ((1 + delta * (-1 + Nu) - 3 * Nu) * S1z +
                  (-1 + delta * (-1 + Nu) + 3 * Nu) * S2z)) /
                (sqrt(65) * params.r4)));
  }

  else if (vpnorder == 7) {

    return (
        (Complex(0, -0.000015031265031265032) * delta * mass * PhiDOT *
         (16 * params.Mtot3 *
              (-3893 + 16386 * Nu - 17175 * params.eta2 + 6922 * params.eta3) +
          8 * params.Mtot2 * r *
              ((-22256 + 111525 * Nu - 170247 * params.eta2 +
                94453 * params.eta3) *
                   params.PhiDOT2 * params.r2 +
               Complex(0, 18) *
                   (1954 - 7090 * Nu + 1431 * params.eta2 + 5232 * params.eta3) *
                   PhiDOT * r * rDOT -
               3 *
                   (6603 - 26446 * Nu + 16325 * params.eta2 +
                    7138 * params.eta3) *
                   params.rDOT2) -
          12 * mass * params.r2 *
              (2 * (771 + 6004 * Nu - 44896 * params.eta2 + 48978 * params.eta3) *
                   params.PhiDOT4 * params.r4 +
               Complex(0, 1) *
                   (10771 - 40224 * Nu + 14456 * params.eta2 +
                    21414 * params.eta3) *
                   params.PhiDOT3 * params.r3 * rDOT +
               (34437 - 136928 * Nu + 85844 * params.eta2 + 30834 * params.eta3) *
                   params.PhiDOT2 * params.r2 * params.rDOT2 +
               Complex(0, 3) *
                   (9713 - 38048 * Nu + 21348 * params.eta2 +
                    11562 * params.eta3) *
                   PhiDOT * r * params.rDOT3 -
               4 * (1811 - 7007 * Nu + 3585 * params.eta2 + 2511 * params.eta3) *
                   params.rDOT4) +
          9 * params.r3 *
              ((4339 - 4350 * Nu - 50298 * params.eta2 + 61600 * params.eta3) *
                   params.PhiDOT6 * params.r6 +
               Complex(0, 2) *
                   (-12661 + 53118 * Nu - 45494 * params.eta2 +
                    2652 * params.eta3) *
                   params.PhiDOT5 * params.r5 * rDOT +
               12 * (885 - 2346 * Nu - 3253 * params.eta2 + 5846 * params.eta3) *
                   params.PhiDOT4 * params.r4 * params.rDOT2 +
               Complex(0, 8) *
                   (211 - 90 * Nu - 2692 * params.eta2 + 2880 * params.eta3) *
                   params.PhiDOT3 * params.r3 * params.rDOT3 +
               16 * (181 - 522 * Nu - 493 * params.eta2 + 1062 * params.eta3) *
                   params.PhiDOT2 * params.r2 * params.rDOT4 +
               Complex(0, 48) *
                   (45 - 102 * Nu - 236 * params.eta2 + 352 * params.eta3) *
                   PhiDOT * r * params.rDOT5 -
               160 * (3 - 6 * Nu - 19 * params.eta2 + 26 * params.eta3) *
                   params.rDOT6))) /
        (sqrt(65) * params.r3));
  }

  else {
    return 0;
  }
}

/* static COMPLEX16 hQC_6_m_3(REAL8 Nu, UINT4 vpnorder, REAL8 x, struct kepler_vars params){
    REAL8 delta = sqrt(1-4*Nu);

    if(vpnorder == 5){
        return((Complex(0,-0.262987012987013)*delta*(1 - 4*Nu + 3*params.eta2)*
        params.x3p5)/sqrt(65));
    }
    else{
        return 0;
    }
} */

static COMPLEX16 hl_6_m_3(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_6_m_3: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_6_m_3(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z,
                      S2z, params) /* +hQC_6_m_3(Nu,vpnorder,x, params) */) *
           cpolar(1, -3 * Phi);
  }
}

static COMPLEX16 hl_6_m_min3(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_6_m_min3: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           conj(hGO_6_m_3(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z,
                          S2z, params) /* +hQC_6_m_3(Nu,vpnorder,x, params) */) *
           cpolar(1, 3 * Phi);
  }
}

// 62
static COMPLEX16 hGO_6_m_2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z, struct kepler_vars params) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  if (vpnorder == 4) {
    return ((1 - 5 * Nu + 5 * params.eta2) *
            (516 * params.Mtot3 +
             360 * params.r3 * pow(PhiDOT * r - Complex(0, 1) * rDOT, 2) *
                 pow(PhiDOT * r + Complex(0, 1) * rDOT, 4) +
             params.Mtot2 * r *
                 (-145 * params.PhiDOT2 * params.r2 +
                  Complex(0, 2920) * PhiDOT * r * rDOT - 2418 * params.rDOT2) -
             3 * mass * params.r2 *
                 (233 * params.PhiDOT4 * params.r4 +
                  Complex(0, 1050) * params.PhiDOT3 * params.r3 * rDOT -
                  252 * params.PhiDOT2 * params.r2 * params.rDOT2 +
                  Complex(0, 1040) * PhiDOT * r * params.rDOT3 -
                  600 * params.rDOT4))) /
           (2376. * sqrt(65) * params.r3);
  }

  else if (vpnorder == 6) {
    return (14 * params.Mtot4 *
                (-5740 + 29361 * Nu - 33348 * params.eta2 + 7334 * params.eta3) -
            7560 * (-1 + 9 * Nu - 26 * params.eta2 + 23 * params.eta3) *
                params.r4 * pow(PhiDOT * r - Complex(0, 1) * rDOT, 3) *
                pow(PhiDOT * r + Complex(0, 1) * rDOT, 5) +
            2 * params.Mtot3 * r *
                ((9187 - 37401 * Nu + 16548 * params.eta2 + 2821 * params.eta3) *
                     params.PhiDOT2 * params.r2 +
                 Complex(0, 2) *
                     (-47984 + 275121 * Nu - 442540 * params.eta2 +
                      255850 * params.eta3) *
                     PhiDOT * r * rDOT +
                 14 *
                     (3614 - 21621 * Nu + 39684 * params.eta2 -
                      29332 * params.eta3) *
                     params.rDOT2) +
            params.Mtot2 * params.r2 *
                ((54699 - 336749 * Nu + 596995 * params.eta2 -
                  337960 * params.eta3) *
                     params.PhiDOT4 * params.r4 -
                 Complex(0, 4) *
                     (-5781 + 89572 * Nu - 379358 * params.eta2 +
                      444689 * params.eta3) *
                     params.PhiDOT3 * params.r3 * rDOT +
                 (-9351 + 101899 * Nu - 419300 * params.eta2 +
                  566195 * params.eta3) *
                     params.PhiDOT2 * params.r2 * params.rDOT2 -
                 Complex(0, 32) *
                     (3921 - 11389 * Nu - 27265 * params.eta2 +
                      58450 * params.eta3) *
                     PhiDOT * r * params.rDOT3 +
                 42 *
                     (1867 - 5501 * Nu - 12824 * params.eta2 +
                      28137 * params.eta3) *
                     params.rDOT4) +
            3 * mass * params.r3 *
                ((-7885 + 64211 * Nu - 170905 * params.eta2 +
                  146580 * params.eta3) *
                     params.PhiDOT6 * params.r6 +
                 Complex(0, 4) *
                     (1438 + 9779 * Nu - 86023 * params.eta2 +
                      109949 * params.eta3) *
                    params.PhiDOT5 * params.r5 * rDOT +
                 (16353 - 68054 * Nu + 10297 * params.eta2 +
                  77294 * params.eta3) *
                     params.PhiDOT4 * params.r4 * params.rDOT2 +
                 Complex(0, 2) *
                     (14341 - 392 * Nu - 316841 * params.eta2 +
                      452508 * params.eta3) *
                     params.PhiDOT3 * params.r3 * params.rDOT3 -
                 4 *
                     (13 + 12530 * Nu - 68803 * params.eta2 +
                      80654 * params.eta3) *
                     params.PhiDOT2 * params.r2 * params.rDOT4 +
                 Complex(0, 80) *
                     (228 - 217 * Nu - 3871 * params.eta2 + 5803 * params.eta3) *
                     PhiDOT * r * params.rDOT5 -
                 840 * (9 + 4 * Nu - 221 * params.eta2 + 308 * params.eta3) *
                     params.rDOT6)) /
           (33264. * sqrt(65) * params.r4);
  }

  else if (vpnorder == 7) {

    return (/* (4*(174*params.eta3*(S1z + S2z) - 14*(S1z + delta*S1z + S2z -
      delta*S2z)
      + Nu*(53*delta*(S1z - S2z) + 115*(S1z + S2z)) - 2*params.eta2*(18*delta*(S1z
      - S2z)
      + 139*(S1z + S2z)))*params.x4p5)/(2079.*sqrt(65))
      + */
            /*Henry et al. ecc+spin terms */
            ((params.Mtot2 *
              (40 * params.Mtot2 *
                   (6 * params.eta3 *
                        (3489 * PhiDOT * r + Complex(0, 24880) * rDOT) *
                        (S1z + S2z) -
                    21 * (683 * PhiDOT * r + Complex(0, 476) * rDOT) *
                        (S1z + delta * S1z + S2z - delta * S2z) -
                    3 * params.eta2 *
                        (PhiDOT * r *
                             (Complex(0, 65990) +
                              (28411 + 26377 * delta) * S1z + 28411 * S2z -
                              26377 * delta * S2z) +
                         Complex(0, 1) * rDOT *
                             (Complex(0, 45321) +
                              40 * (1592 + 485 * delta) * S1z + 63680 * S2z -
                              19400 * delta * S2z)) +
                    Nu * (PhiDOT * r *
                              (Complex(0, 61450) + 73677 * S1z +
                               75423 * delta * S1z + 73677 * S2z -
                               75423 * delta * S2z) +
                          Complex(0, 3) * rDOT *
                              (Complex(0, 14017) + 25790 * S1z +
                               18030 * delta * S1z + 25790 * S2z -
                               18030 * delta * S2z))) -
               12 * mass * r *
                   (30 * params.eta3 *
                        (9893 * params.PhiDOT3 * params.r3 +
                         Complex(0, 55432) * params.PhiDOT2 * params.r2 * rDOT -
                         5479 * PhiDOT * r * params.rDOT2 +
                         Complex(0, 20460) * params.rDOT3) *
                        (S1z + S2z) -
                    210 *
                        (67 * params.PhiDOT3 * params.r3 -
                         Complex(0, 122) * params.PhiDOT2 * params.r2 * rDOT +
                         433 * PhiDOT * r * params.rDOT2 +
                         Complex(0, 236) * params.rDOT3) *
                        (S1z + delta * S1z + S2z - delta * S2z) +
                    Nu * (params.PhiDOT3 * params.r3 *
                              (Complex(0, 285011) +
                               10 * (10608 + 10973 * delta) * S1z +
                               (106080 - 109730 * delta) * S2z) +
                          5 * PhiDOT * r * params.rDOT2 *
                              (Complex(0, 50631) +
                               (80562 + 97252 * delta) * S1z +
                               (80562 - 97252 * delta) * S2z) +
                          Complex(0, 5) * params.rDOT3 *
                              (Complex(0, 24843) +
                               10 * (7293 + 5153 * delta) * S1z +
                               (72930 - 51530 * delta) * S2z) -
                          Complex(0, 2) * params.PhiDOT2 * params.r2 * rDOT *
                              (Complex(0, -12613) +
                               40 * (-2676 + 1801 * delta) * S1z -
                               40 * (2676 + 1801 * delta) * S2z)) -
                    params.eta2 *
                        (5 * PhiDOT * r * params.rDOT2 *
                             (Complex(0, 161463) +
                              2 * (22706 + 51787 * delta) * S1z +
                              (45412 - 103574 * delta) * S2z) +
                         2 * params.PhiDOT3 * params.r3 *
                             (Complex(0, 460269) +
                              5 * (28838 + 14911 * delta) * S1z +
                              (144190 - 74555 * delta) * S2z) +
                         Complex(0, 5) * params.rDOT3 *
                             (Complex(0, 80079) +
                              20 * (8419 + 2675 * delta) * S1z -
                              20 * (-8419 + 2675 * delta) * S2z) -
                         Complex(0, 2) * params.PhiDOT2 * params.r2 * rDOT *
                             (Complex(0, -36879) +
                              10 * (-78341 + 8003 * delta) * S1z -
                              10 * (78341 + 8003 * delta) * S2z))) +
               3 * params.r2 *
                   (60 * params.eta3 *
                        (11215 *params.PhiDOT5 * params.r5 +
                         Complex(0, 57242) * params.PhiDOT4 * params.r4 * rDOT +
                         28128 * params.PhiDOT3 * params.r3 * params.rDOT2 +
                         Complex(0, 57112) * params.PhiDOT2 * params.r2 *
                             params.rDOT3 -
                         6992 * PhiDOT * r * params.rDOT4 +
                         Complex(0, 9600) * params.rDOT5) *
                        (S1z + S2z) +
                    6300 *
                        (9 *params.PhiDOT5 * params.r5 +
                         Complex(0, 9) * params.PhiDOT4 * params.r4 * rDOT -
                         28 * params.PhiDOT3 * params.r3 * params.rDOT2 -
                         Complex(0, 6) * params.PhiDOT2 * params.r2 *
                             params.rDOT3 -
                         16 * PhiDOT * r * params.rDOT4 -
                         Complex(0, 8) * params.rDOT5) *
                        (S1z + delta * S1z + S2z - delta * S2z) +
                    params.eta2 *
                        (Complex(0, -4) * params.PhiDOT2 * params.r2 *
                             params.rDOT3 *
                             (Complex(0, 37829) +
                              30 * (30617 + 1089 * delta) * S1z +
                              (918510 - 32670 * delta) * S2z) -
                         48 * params.PhiDOT3 * params.r3 * params.rDOT2 *
                             (Complex(0, 44162) +
                              5 * (10459 + 3923 * delta) * S1z +
                              (52295 - 19615 * delta) * S2z) -
                         Complex(0, 120) * params.rDOT5 *
                             (Complex(0, 2449) + 200 * (34 + 11 * delta) * S1z -
                              200 * (-34 + 11 * delta) * S2z) +
                         4 *params.PhiDOT5 * params.r5 *
                             (Complex(0, -354934) +
                              75 * (-650 + 297 * delta) * S1z -
                              75 * (650 + 297 * delta) * S2z) -
                         320 * PhiDOT * r * params.rDOT4 *
                             (Complex(0, 2449) + 3 * (36 + 577 * delta) * S1z -
                              3 * (-36 + 577 * delta) * S2z) +
                         Complex(0, 1) * params.PhiDOT4 * params.r4 * rDOT *
                             (Complex(0, 880129) +
                              120 * (-28647 + 4751 * delta) * S1z -
                              120 * (28647 + 4751 * delta) * S2z)) -
                    3 * Nu *
                        (Complex(0, -40) * params.rDOT5 *
                             (Complex(0, 753) + 10 * (303 + 215 * delta) * S1z +
                              (3030 - 2150 * delta) * S2z) -
                         320 * PhiDOT * r * params.rDOT4 *
                             (Complex(0, 251) + (422 + 551 * delta) * S1z +
                              (422 - 551 * delta) * S2z) -
                         Complex(0, 4) * params.PhiDOT2 * params.r2 *
                             params.rDOT3 *
                             (Complex(0, 3971) +
                              40 * (1858 + 333 * delta) * S1z -
                              40 * (-1858 + 333 * delta) * S2z) +
                         4 *params.PhiDOT5 * params.r5 *
                             (Complex(0, -37091) +
                              5 * (3454 + 3105 * delta) * S1z -
                              5 * (-3454 + 3105 * delta) * S2z) -
                         16 * params.PhiDOT3 * params.r3 * params.rDOT2 *
                             (Complex(0, 13464) +
                              5 * (5007 + 3799 * delta) * S1z -
                              5 * (-5007 + 3799 * delta) * S2z) +
                         Complex(0, 1) * params.PhiDOT4 * params.r4 * rDOT *
                             (Complex(0, 94671) +
                              10 * (-16313 + 14227 * delta) * S1z -
                              10 * (16313 + 14227 * delta) * S2z))))) /
             (997920. * sqrt(65) * params.r4)));
  }

  else {
    return 0;
  }
}

/* static COMPLEX16 hQC_6_m_2(REAL8 Nu, UINT4 vpnorder, REAL8 x, struct kepler_vars params){

    if(vpnorder == 4){
        return((4/(297.*sqrt(65)) - (4*sqrt(0.38461538461538464)*Nu)/297. +
        (4*sqrt(0.38461538461538464)*params.eta2)/297.)*params.x3);
    }
    else if(vpnorder == 6){
        return((-6/(77.*sqrt(65)) + (118*Nu)/(297.*sqrt(65)) -
        (128*params.eta2)/(297.*sqrt(65)) +
        (14*params.eta3)/(297.*sqrt(65)))*params.x4);
    }
    else{
        return 0;
    }
} */

static COMPLEX16 hl_6_m_2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_6_m_1: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_6_m_2(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z,
                      S2z, params) /* +hQC_6_m_2(Nu,vpnorder,x,params) */) *
           cpolar(1, -2 * Phi);
  }
}

static COMPLEX16 hl_6_m_min2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_6_m_min2: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           conj(hGO_6_m_2(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z,
                          S2z, params) /* +hQC_6_m_2(Nu,vpnorder,x,params) */) *
           cpolar(1, 2 * Phi);
  }
}

// 61
static COMPLEX16 hGO_6_m_1(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z, struct kepler_vars params) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  if (vpnorder == 5) {
    return (Complex(0, 0.0002405002405002405) * delta * mass *
            (1 - 4 * Nu + 3 * params.eta2) * PhiDOT *
            (410 * params.Mtot2 -
             2 * mass * r *
                 (359 * params.PhiDOT2 * params.r2 -
                  Complex(0, 343) * PhiDOT * r * rDOT + 310 * params.rDOT2) +
             3 * params.r2 *
                 (103 * params.PhiDOT4 * params.r4 -
                  Complex(0, 196) * params.PhiDOT3 * params.r3 * rDOT +
                  108 * params.PhiDOT2 * params.r2 * params.rDOT2 -
                  Complex(0, 56) * PhiDOT * r * params.rDOT3 +
                  40 * params.rDOT4))) /
           (sqrt(26) * params.r2);
  }

  else if (vpnorder == 6) {

    return (/* (Complex(0,-0.00042087542087542086)*Nu*(S1z + delta*(-1 + Nu)*S1z
      - 3*Nu*S1z
      - (1 + delta)*S2z + (3 + delta)*Nu*S2z)*params.x4)/sqrt(26)
      + */
            /* Henry et al. ecc spin terms */ (
                (Complex(0, -0.00042087542087542086) * params.Mtot2 * Nu *
                 (410 * params.Mtot2 -
                  2 * mass * r *
                      (359 * params.PhiDOT2 * params.r2 -
                       Complex(0, 343) * PhiDOT * r * rDOT +
                       310 * params.rDOT2) +
                  3 * params.r2 *
                      (103 * params.PhiDOT4 * params.r4 -
                       Complex(0, 196) * params.PhiDOT3 * params.r3 * rDOT +
                       108 * params.PhiDOT2 * params.r2 * params.rDOT2 -
                       Complex(0, 56) * PhiDOT * r * params.rDOT3 +
                       40 * params.rDOT4)) *
                 ((1 + delta * (-1 + Nu) - 3 * Nu) * S1z +
                  (-1 + delta * (-1 + Nu) + 3 * Nu) * S2z)) /
                (sqrt(26) * params.r4)));
  }

  else if (vpnorder == 7) {

    return ((delta * mass * PhiDOT *
             (Complex(0, 16) * params.Mtot3 *
                  (-3893 + Nu * (16386 + Nu * (-17175 + 6922 * Nu))) +
              24 * params.Mtot2 * r *
                  (Complex(0, -1) *
                       (-1608 + Nu * (7857 + Nu * (-14179 + 11585 * Nu))) *
                       params.PhiDOT2 * params.r2 -
                   2 * (1954 + Nu * (-7090 + 3 * Nu * (477 + 1744 * Nu))) *
                       PhiDOT * r * rDOT -
                   Complex(0, 1) *
                       (6603 + Nu * (-26446 + Nu * (16325 + 7138 * Nu))) *
                       params.rDOT2) +
              12 * mass * params.r2 *
                  (Complex(0, 2) *
                       (973 + 2 * Nu * (-898 + Nu * (-4168 + 6015 * Nu))) *
                       params.PhiDOT4 * params.r4 +
                   (25393 + 2 * Nu * (-48592 + Nu * (24716 + 15777 * Nu))) *
                       params.PhiDOT3 * params.r3 * rDOT +
                   Complex(0, 3) *
                       (6017 + 2 * Nu * (-11616 + Nu * (5954 + 4053 * Nu))) *
                       params.PhiDOT2 * params.r2 * params.rDOT2 +
                   (9713 + 2 * Nu * (-19024 + 3 * Nu * (3558 + 1927 * Nu))) *
                       PhiDOT * r * params.rDOT3 +
                   Complex(0, 4) *
                       (1811 + Nu * (-7007 + 3 * Nu * (1195 + 837 * Nu))) *
                       params.rDOT4) +
              9 * params.r3 *
                  (Complex(0, -1) *
                       (781 + 2 * Nu * (-81 + Nu * (-5699 + 6336 * Nu))) *
                       params.PhiDOT6 * params.r6 -
                   2 * (8025 + 2 * Nu * (-15003 + Nu * (6791 + 5258 * Nu))) *
                       params.PhiDOT5 * params.r5 * rDOT -
                   Complex(0, 844) * (3 + Nu * (-6 + Nu * (-19 + 26 * Nu))) *
                       params.PhiDOT4 * params.r4 * params.rDOT2 +
                   8 * (-425 + 2 * Nu * (519 + 2 * (481 - 772 * Nu) * Nu)) *
                       params.PhiDOT3 * params.r3 * params.rDOT3 -
                   Complex(0, 592) * (3 + Nu * (-6 + Nu * (-19 + 26 * Nu))) *
                       params.PhiDOT2 * params.r2 * params.rDOT4 +
                   16 * (-45 + 2 * Nu * (51 + 2 * (59 - 88 * Nu) * Nu)) *
                       PhiDOT * r * params.rDOT5 -
                   Complex(0, 160) * (3 + Nu * (-6 + Nu * (-19 + 26 * Nu))) *
                       params.rDOT6))) /
            (99792. * sqrt(26) * params.r3));
  }

  else {
    return 0;
  }
}

/* static COMPLEX16 hQC_6_m_1(REAL8 Nu, UINT4 vpnorder, REAL8 x, struct kepler_vars params){
    REAL8 delta = sqrt(1-4*Nu);

    if(vpnorder == 5){
        return((Complex(0,0.0002405002405002405)*delta*(1 - 4*Nu + 3*params.eta2)*
        params.x3p5)/sqrt(26));
    }
    else{
        return 0;
    }
} */

static COMPLEX16 hl_6_m_1(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_6_m_1: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_6_m_1(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z,
                      S2z, params) /* +hQC_6_m_1(Nu,vpnorder,x,params) */) *
           cpolar(1, -1 * Phi);
  }
}

static COMPLEX16 hl_6_m_min1(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_6_m_min1: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           conj(hGO_6_m_1(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z,
                          S2z, params) /* +hQC_6_m_1(Nu,vpnorder,x,params) */) *
           cpolar(1, 1 * Phi);
  }
}

// 77
static COMPLEX16 hGO_7_m_7(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, struct kepler_vars params) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  if (vpnorder == 5) {
    return (delta * (1 - 4 * Nu + 3 * params.eta2) *
            (720 * params.r3 * pow(Complex(0, 1) * PhiDOT * r - rDOT, 7) +
             2 * params.Mtot3 * (Complex(0, -4559) * PhiDOT * r + 1976 * rDOT) +
             18 * params.Mtot2 * r *
                 (Complex(0, -3317) * params.PhiDOT3 * params.r3 +
                  4143 * params.PhiDOT2 * params.r2 * rDOT +
                  Complex(0, 2178) * PhiDOT * r * params.rDOT2 -
                  442 * params.rDOT3) +
             45 * mass * params.r2 *
                 (Complex(0, -1069) * params.PhiDOT5 * params.r5 +
                  2112 * params.PhiDOT4 * params.r4 * rDOT +
                  Complex(0, 2370) * params.PhiDOT3 * params.r3 * params.rDOT2 -
                  1600 * params.PhiDOT2 * params.r2 * params.rDOT3 -
                  Complex(0, 600) * PhiDOT * r * params.rDOT4 +
                  96 * params.rDOT5))) /
           (720. * sqrt(6006) * params.r3);
  }

  else if (vpnorder == 7) {

    return (
        (delta *
         (5040 * (-59 + 473 * Nu - 1185 * params.eta2 + 831 * params.eta3) *
              params.r4 * pow(PhiDOT * r + Complex(0, 1) * rDOT, 8) *
              (Complex(0, 1) * PhiDOT * r + rDOT) +
          4 * params.Mtot4 *
              (Complex(0, -7) *
                   (-845315 + 3567934 * Nu - 3505793 * params.eta2 +
                    1006326 * params.eta3) *
                   PhiDOT * r +
               4 *
                   (-482641 + 2069935 * Nu - 2198091 * params.eta2 +
                    803481 * params.eta3) *
                   rDOT) +
          2 * params.Mtot3 * r *
              (Complex(0, -7) *
                   (-8046605 + 37481621 * Nu - 46964439 * params.eta2 +
                    19173243 * params.eta3) *
                   params.PhiDOT3 * params.r3 +
               2 *
                   (-17292661 + 90383188 * Nu - 148920723 * params.eta2 +
                    88022760 * params.eta3) *
                   params.PhiDOT2 * params.r2 * rDOT +
               Complex(0, 7) *
                   (-1175695 + 7490560 * Nu - 17116047 * params.eta2 +
                    13239024 * params.eta3) *
                   PhiDOT * r * params.rDOT2 -
               8 *
                   (-74093 + 755216 * Nu - 2568447 * params.eta2 +
                    2398116 * params.eta3) *
                   params.rDOT3) +
          18 * params.Mtot2 * params.r2 *
              (Complex(0, -14) *
                   (-270811 + 1455597 * Nu - 2122103 * params.eta2 +
                    757575 * params.eta3) *
                   params.PhiDOT5 * params.r5 +
               (1113940 + 4956889 * Nu - 34515698 * params.eta2 +
                28651791 * params.eta3) *
                   params.PhiDOT4 * params.r4 * rDOT +
               Complex(0, 7) *
                   (715061 - 1359787 * Nu - 4470709 * params.eta2 +
                    5729499 * params.eta3) *
                   params.PhiDOT3 * params.r3 * params.rDOT2 +
               (-4243543 + 9349627 * Nu + 22353219 * params.eta2 -
                32044971 * params.eta3) *
                   params.PhiDOT2 * params.r2 * params.rDOT3 -
               Complex(0, 14) *
                   (119115 - 252749 * Nu - 689509 * params.eta2 +
                    975153 * params.eta3) *
                   PhiDOT * r * params.rDOT4 +
               2 *
                   (131615 - 252847 * Nu - 892507 * params.eta2 +
                    1206639 * params.eta3) *
                   params.rDOT5) +
          315 * mass * params.r3 *
              (Complex(0, 1) *
                   (5005 + 44124 * Nu - 401470 * params.eta2 +
                    512250 * params.eta3) *
                   params.PhiDOT7 * params.r7 +
               2 *
                   (32339 - 197757 * Nu + 459526 * params.eta2 -
                    383013 * params.eta3) *
                   params.PhiDOT6 * params.r6 * rDOT -
               Complex(0, 1) *
                   (-168255 + 685408 * Nu - 635707 * params.eta2 +
                    199944 * params.eta3) *
                   params.PhiDOT5 * params.r5 * params.rDOT2 +
               4 *
                   (-53106 + 146273 * Nu + 123916 * params.eta2 -
                    235713 * params.eta3) *
                   params.PhiDOT4 * params.r4 * params.rDOT3 -
               Complex(0, 10) *
                   (16287 - 22820 * Nu - 136891 * params.eta2 +
                    159864 * params.eta3) *
                   params.PhiDOT3 * params.r3 * params.rDOT4 +
               16 *
                   (4761 + 725 * Nu - 72818 * params.eta2 + 75357 * params.eta3) *
                   params.PhiDOT2 * params.r2 * params.rDOT5 +
               Complex(0, 40) *
                   (499 + 1018 * Nu - 11789 * params.eta2 + 11502 * params.eta3) *
                   PhiDOT * r * params.rDOT6 -
               32 * (70 + 311 * Nu - 2394 * params.eta2 + 2253 * params.eta3) *
                   params.rDOT7))) /
        (171360. * sqrt(6006) * params.r4));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hl_7_m_7(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_7_m_7: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           hGO_7_m_7(mass, Nu, r, rDOT, PhiDOT, vpnorder, params) * cpolar(1, -7 * Phi);
  }
}

static COMPLEX16 hl_7_m_min7(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_7_m_min7: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((- 4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           conj(hGO_7_m_7(mass, Nu, r, rDOT, PhiDOT, vpnorder, params)) *
           cpolar(1, 7 * Phi);
  }
}

// 75
static COMPLEX16 hGO_7_m_5(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, struct kepler_vars params) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  if (vpnorder == 5) {
    return (
        (delta * (1 - 4 * Nu + 3 * params.eta2) *
         (2 * params.Mtot3 * (Complex(0, 22795) * PhiDOT * r - 13832 * rDOT) -
          Complex(0, 5040) * params.r3 * (PhiDOT * r - Complex(0, 1) * rDOT) *
              pow(PhiDOT * r + Complex(0, 1) * rDOT, 6) +
          18 * params.Mtot2 * r *
              (Complex(0, 5105) * params.PhiDOT3 * params.r3 -
               13041 * params.PhiDOT2 * params.r2 * rDOT -
               Complex(0, 10890) * PhiDOT * r * params.rDOT2 +
               3094 * params.rDOT3) +
          45 * mass * params.r2 *
              (Complex(0, -1207) * params.PhiDOT5 * params.r5 +
               336 * params.PhiDOT4 * params.r4 * rDOT -
               Complex(0, 3114) * params.PhiDOT3 * params.r3 * params.rDOT2 +
               4928 * params.PhiDOT2 * params.r2 * params.rDOT3 +
               Complex(0, 3000) * PhiDOT * r * params.rDOT4 -
               672 * params.rDOT5))) /
        (65520. * sqrt(66) * params.r3));
  }

  else if (vpnorder == 7) {

    return (
        (delta *
         (Complex(0, 5040) *
              (-59 + 473 * Nu - 1185 * params.eta2 + 831 * params.eta3) *
              params.r4 * pow(PhiDOT * r - Complex(0, 1) * rDOT, 2) *
              pow(PhiDOT * r + Complex(0, 1) * rDOT, 7) +
          4 * params.Mtot4 *
              (Complex(0, 5) *
                   (-845315 + 3567934 * Nu - 3505793 * params.eta2 +
                    1006326 * params.eta3) *
                   PhiDOT * r +
               4 *
                   (482641 - 2069935 * Nu + 2198091 * params.eta2 -
                    803481 * params.eta3) *
                   rDOT) +
          2 * params.Mtot3 * r *
              (Complex(0, 5) *
                   (-2480333 + 11838029 * Nu - 15223575 * params.eta2 +
                    5981667 * params.eta3) *
                   params.PhiDOT3 * params.r3 -
               2 *
                   (-7732189 + 40782964 * Nu - 67735947 * params.eta2 +
                    39807720 * params.eta3) *
                   params.PhiDOT2 * params.r2 * rDOT -
               Complex(0, 5) *
                   (-1175695 + 7490560 * Nu - 17116047 * params.eta2 +
                    13239024 * params.eta3) *
                   PhiDOT * r * params.rDOT2 +
               8 *
                   (-74093 + 755216 * Nu - 2568447 * params.eta2 +
                    2398116 * params.eta3) *
                   params.rDOT3) +
          18 * params.Mtot2 * params.r2 *
              (Complex(0, -10) *
                   (-58221 + 312011 * Nu - 571313 * params.eta2 +
                    397665 * params.eta3) *
                   params.PhiDOT5 * params.r5 +
               (-117884 + 723063 * Nu - 2639174 * params.eta2 +
                3313409 * params.eta3) *
                   params.PhiDOT4 * params.r4 * rDOT -
               Complex(0, 5) *
                   (236029 - 577139 * Nu - 804413 * params.eta2 +
                    1190115 * params.eta3) *
                   params.PhiDOT3 * params.r3 * params.rDOT2 +
               (1940751 - 4452163 * Nu - 9331755 * params.eta2 +
                13753811 * params.eta3) *
                   params.PhiDOT2 * params.r2 * params.rDOT3 +
               Complex(0, 10) *
                   (119115 - 252749 * Nu - 689509 * params.eta2 +
                    975153 * params.eta3) *
                   PhiDOT * r * params.rDOT4 +
               2 *
                   (-131615 + 252847 * Nu + 892507 * params.eta2 -
                    1206639 * params.eta3) *
                   params.rDOT5) +
          45 * mass * params.r3 *
              (Complex(0, 1) *
                   (39847 - 121612 * Nu - 163210 * params.eta2 +
                    376622 * params.eta3) *
                   params.PhiDOT7 * params.r7 +
               2 *
                   (79987 - 238709 * Nu - 93106 * params.eta2 +
                    259939 * params.eta3) *
                   params.PhiDOT6 * params.r6 * rDOT +
               Complex(0, 1) *
                   (173997 + 166656 * Nu - 3315983 * params.eta2 +
                    3362728 * params.eta3) *
                   params.PhiDOT5 * params.r5 * params.rDOT2 -
               4 *
                   (-13074 + 266583 * Nu - 998488 * params.eta2 +
                    847097 * params.eta3) *
                   params.PhiDOT4 * params.r4 * params.rDOT3 +
               Complex(0, 2) *
                   (135447 - 440484 * Nu - 25811 * params.eta2 +
                    357784 * params.eta3) *
                   params.PhiDOT3 * params.r3 * params.rDOT4 -
               16 *
                   (15423 - 4597 * Nu - 205150 * params.eta2 +
                    217363 * params.eta3) *
                   params.PhiDOT2 * params.r2 * params.rDOT5 -
               Complex(0, 200) *
                   (499 + 1018 * Nu - 11789 * params.eta2 + 11502 * params.eta3) *
                   PhiDOT * r * params.rDOT6 +
               224 * (70 + 311 * Nu - 2394 * params.eta2 + 2253 * params.eta3) *
                   params.rDOT7))) /
        (2.22768e6 * sqrt(66) * params.r4));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hl_7_m_5(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_7_m_5: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           hGO_7_m_5(mass, Nu, r, rDOT, PhiDOT, vpnorder, params) * cpolar(1, -5 * Phi);
  }
}

static COMPLEX16 hl_7_m_min5(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_7_m_min5: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((- 4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           conj(hGO_7_m_5(mass, Nu, r, rDOT, PhiDOT, vpnorder, params)) *
           cpolar(1, 5 * Phi);
  }
}

// 73
static COMPLEX16 hGO_7_m_3(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, struct kepler_vars params) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  if (vpnorder == 5) {
    return (delta * (1 - 4 * Nu + 3 * params.eta2) *
            (5040 * params.r3 * pow(PhiDOT * r - Complex(0, 1) * rDOT, 2) *
                 pow(Complex(0, -1) * PhiDOT * r + rDOT, 5) +
             2 * params.Mtot3 *
                 (Complex(0, -13677) * PhiDOT * r + 13832 * rDOT) +
             18 * params.Mtot2 * r *
                 (Complex(0, 1529) * params.PhiDOT3 * params.r3 +
                  2401 * params.PhiDOT2 * params.r2 * rDOT +
                  Complex(0, 6534) * PhiDOT * r * params.rDOT2 -
                  3094 * params.rDOT3) +
             15 * mass * params.r2 *
                 (Complex(0, 179) * params.PhiDOT5 * params.r5 -
                  4368 * params.PhiDOT4 * params.r4 * rDOT -
                  Complex(0, 4878) * params.PhiDOT3 * params.r3 * params.rDOT2 -
                  2240 * params.PhiDOT2 * params.r2 * params.rDOT3 -
                  Complex(0, 5400) * PhiDOT * r * params.rDOT4 +
                  2016 * params.rDOT5))) /
           (240240. * sqrt(6) * params.r3);
  }

  else if (vpnorder == 7) {
    return (
        (delta *
         (-5040 * (-59 + 473 * Nu - 1185 * params.eta2 + 831 * params.eta3) *
              params.r4 * pow(PhiDOT * r + Complex(0, 1) * rDOT, 6) *
              pow(Complex(0, 1) * PhiDOT * r + rDOT, 3) +
          4 * params.Mtot4 *
              (Complex(0, -3) *
                   (-845315 + 3567934 * Nu - 3505793 * params.eta2 +
                    1006326 * params.eta3) *
                   PhiDOT * r +
               4 *
                   (-482641 + 2069935 * Nu - 2198091 * params.eta2 +
                    803481 * params.eta3) *
                   rDOT) +
          2 * params.Mtot3 * r *
              (Complex(0, 3) *
                   (-1230515 + 5257699 * Nu - 5937001 * params.eta2 +
                    2812717 * params.eta3) *
                   params.PhiDOT3 * params.r3 +
               2 *
                   (-1358541 + 7716148 * Nu - 13612763 * params.eta2 +
                    7664360 * params.eta3) *
                   params.PhiDOT2 * params.r2 * rDOT +
               Complex(0, 3) *
                   (-1175695 + 7490560 * Nu - 17116047 * params.eta2 +
                    13239024 * params.eta3) *
                   PhiDOT * r * params.rDOT2 -
               8 *
                   (-74093 + 755216 * Nu - 2568447 * params.eta2 +
                    2398116 * params.eta3) *
                   params.rDOT3) +
          6 * params.Mtot2 * params.r2 *
              (Complex(0, -2) *
                   (76741 - 370227 * Nu + 248553 * params.eta2 +
                    279655 * params.eta3) *
                   params.PhiDOT5 * params.r5 -
               3 *
                   (65300 + 667351 * Nu - 4038422 * params.eta2 +
                    3825889 * params.eta3) *
                   params.PhiDOT4 * params.r4 * rDOT -
               Complex(0, 3) *
                   (249977 + 166121 * Nu - 4919353 * params.eta2 +
                    5508423 * params.eta3) *
                   params.PhiDOT3 * params.r3 * params.rDOT2 +
               (-1216669 + 3561561 * Nu + 1952337 * params.eta2 -
                4679113 * params.eta3) *
                   params.PhiDOT2 * params.r2 * params.rDOT3 -
               Complex(0, 18) *
                   (119115 - 252749 * Nu - 689509 * params.eta2 +
                    975153 * params.eta3) *
                   PhiDOT * r * params.rDOT4 +
               6 *
                   (131615 - 252847 * Nu - 892507 * params.eta2 +
                    1206639 * params.eta3) *
                   params.rDOT5) +
          15 * mass * params.r3 *
              (Complex(0, -17) *
                   (-3693 + 24516 * Nu - 50930 * params.eta2 +
                    30982 * params.eta3) *
                   params.PhiDOT7 * params.r7 +
               6 *
                   (14141 + 89013 * Nu - 604958 * params.eta2 +
                    566877 * params.eta3) *
                   params.PhiDOT6 * params.r6 * rDOT +
               Complex(0, 1) *
                   (-125697 + 1209408 * Nu - 3554741 * params.eta2 +
                    2822200 * params.eta3) *
                   params.PhiDOT5 * params.r5 * params.rDOT2 +
               4 *
                   (90786 - 8859 * Nu - 1290784 * params.eta2 +
                    1354859 * params.eta3) *
                   params.PhiDOT4 * params.r4 * params.rDOT3 +
               Complex(0, 6) *
                   (27423 + 212284 * Nu - 1343099 * params.eta2 +
                    1240856 * params.eta3) *
                   params.PhiDOT3 * params.r3 * params.rDOT4 +
               16 *
                   (10461 - 33135 * Nu - 6298 * params.eta2 +
                    31817 * params.eta3) *
                   params.PhiDOT2 * params.r2 * params.rDOT5 +
               Complex(0, 360) *
                   (499 + 1018 * Nu - 11789 * params.eta2 + 11502 * params.eta3) *
                   PhiDOT * r * params.rDOT6 -
               672 * (70 + 311 * Nu - 2394 * params.eta2 + 2253 * params.eta3) *
                   params.rDOT7))) /
        (8.16816e6 * sqrt(6) * params.r4));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hl_7_m_3(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_7_m_3: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           hGO_7_m_3(mass, Nu, r, rDOT, PhiDOT, vpnorder, params) * cpolar(1, -3 * Phi);
  }
}

static COMPLEX16 hl_7_m_min3(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_7_m_min3: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((- 4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           conj(hGO_7_m_3(mass, Nu, r, rDOT, PhiDOT, vpnorder, params)) *
           cpolar(1, 3 * Phi);
  }
}

// 71
static COMPLEX16 hGO_7_m_1(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, struct kepler_vars params) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  if (vpnorder == 5) {
    return (delta * (1 - 4 * Nu + 3 * params.eta2) *
            (2 * params.Mtot3 * (Complex(0, 4559) * PhiDOT * r - 13832 * rDOT) -
             Complex(0, 5040) * params.r3 *
                 pow(PhiDOT * r - Complex(0, 1) * rDOT, 3) *
                 pow(PhiDOT * r + Complex(0, 1) * rDOT, 4) +
             18 * params.Mtot2 * r *
                 (Complex(0, -1275) * params.PhiDOT3 * params.r3 +
                  2919 * params.PhiDOT2 * params.r2 * rDOT -
                  Complex(0, 2178) * PhiDOT * r * params.rDOT2 +
                  3094 * params.rDOT3) +
             Complex(0, 27) * mass * params.r2 *
                 (699 * params.PhiDOT5 * params.r5 +
                  Complex(0, 1120) * params.PhiDOT4 * params.r4 * rDOT +
                  1874 * params.PhiDOT3 * params.r3 * params.rDOT2 +
                  Complex(0, 2240) * params.PhiDOT2 * params.r2 * params.rDOT3 +
                  1000 * PhiDOT * r * params.rDOT4 +
                  Complex(0, 1120) * params.rDOT5))) /
           (432432. * sqrt(2) * params.r3);
  }

  else if (vpnorder == 7) {

    return (
        (delta *
         (Complex(0, 5040) *
              (-59 + 473 * Nu - 1185 * params.eta2 + 831 * params.eta3) *
              params.r4 * pow(PhiDOT * r - Complex(0, 1) * rDOT, 4) *
              pow(PhiDOT * r + Complex(0, 1) * rDOT, 5) +
          4 * params.Mtot4 *
              (Complex(0, 1) *
                   (-845315 + 3567934 * Nu - 3505793 * params.eta2 +
                    1006326 * params.eta3) *
                   PhiDOT * r +
               4 *
                   (482641 - 2069935 * Nu + 2198091 * params.eta2 -
                    803481 * params.eta3) *
                   rDOT) +
          2 * params.Mtot3 * r *
              (Complex(0, -1) *
                   (-3085939 + 13805563 * Nu - 16517289 * params.eta2 +
                    7209909 * params.eta3) *
                   params.PhiDOT3 * params.r3 +
               2 *
                   (-1828283 + 8817260 * Nu - 13448829 * params.eta2 +
                    8407320 * params.eta3) *
                   params.PhiDOT2 * params.r2 * rDOT -
               Complex(0, 1) *
                   (-1175695 + 7490560 * Nu - 17116047 * params.eta2 +
                    13239024 * params.eta3) *
                   PhiDOT * r * params.rDOT2 +
               8 *
                   (-74093 + 755216 * Nu - 2568447 * params.eta2 +
                    2398116 * params.eta3) *
                   params.rDOT3) +
          18 * params.Mtot2 * params.r2 *
              (Complex(0, 2) *
                   (-97035 + 529085 * Nu - 946023 * params.eta2 +
                    605111 * params.eta3) *
                   params.PhiDOT5 * params.r5 +
               (12636 - 513209 * Nu + 2273154 * params.eta2 -
                2157167 * params.eta3) *
                   params.PhiDOT4 * params.r4 * rDOT +
               Complex(0, 3) *
                   (81001 - 68503 * Nu - 953961 * params.eta2 +
                    1116423 * params.eta3) *
                   params.PhiDOT3 * params.r3 * params.rDOT2 +
               (-362041 + 445301 * Nu + 3689709 * params.eta2 -
                4537349 * params.eta3) *
                   params.PhiDOT2 * params.r2 * params.rDOT3 +
               Complex(0, 2) *
                   (119115 - 252749 * Nu - 689509 * params.eta2 +
                    975153 * params.eta3) *
                   PhiDOT * r * params.rDOT4 +
               2 *
                   (-131615 + 252847 * Nu + 892507 * params.eta2 -
                    1206639 * params.eta3) *
                   params.rDOT5) -
          Complex(0, 9) * mass * params.r3 *
              ((-107407 + 883020 * Nu - 2319286 * params.eta2 +
                1727170 * params.eta3) *
                   params.PhiDOT7 * params.r7 +
               Complex(0, 2) *
                   (-937 + 307287 * Nu - 1356450 * params.eta2 +
                    1189583 * params.eta3) *
                   params.PhiDOT6 * params.r6 * rDOT +
               (216555 + 690656 * Nu - 6235017 * params.eta2 +
                5984984 * params.eta3) *
                   params.PhiDOT5 * params.r5 * params.rDOT2 +
               Complex(0, 4) *
                   (34014 + 335581 * Nu - 1987500 * params.eta2 +
                    1820899 * params.eta3) *
                   params.PhiDOT4 * params.r4 * params.rDOT3 +
               2 *
                   (136281 + 310468 * Nu - 3370653 * params.eta2 +
                    3281032 * params.eta3) *
                   params.PhiDOT3 * params.r3 * params.rDOT4 +
               Complex(0, 80) *
                   (2481 + 14269 * Nu - 99426 * params.eta2 +
                    92773 * params.eta3) *
                   params.PhiDOT2 * params.r2 * params.rDOT5 +
               200 *
                   (499 + 1018 * Nu - 11789 * params.eta2 + 11502 * params.eta3) *
                   PhiDOT * r * params.rDOT6 +
               Complex(0, 1120) *
                   (70 + 311 * Nu - 2394 * params.eta2 + 2253 * params.eta3) *
                   params.rDOT7))) /
        (1.4702688e7 * sqrt(2) * params.r4));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hl_7_m_1(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_7_m_1: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           hGO_7_m_1(mass, Nu, r, rDOT, PhiDOT, vpnorder, params) * cpolar(1, -1 * Phi);
  }
}

static COMPLEX16 hl_7_m_min1(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_7_m_min1: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((- 4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           conj(hGO_7_m_1(mass, Nu, r, rDOT, PhiDOT, vpnorder, params)) *
           cpolar(1, 1 * Phi);
  }
}

// 72
static COMPLEX16 hGO_7_m_2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z, struct kepler_vars params) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  if (vpnorder == 6) {
    return (-(mass * (-1 + 7 * Nu - 14 * params.eta2 + 7 * params.eta3) * PhiDOT *
              (8 * params.Mtot2 * (494 * PhiDOT * r + Complex(0, 411) * rDOT) -
               12 * mass * r *
                   (530 * params.PhiDOT3 * params.r3 -
                    Complex(0, 6) * params.PhiDOT2 * params.r2 * rDOT +
                    453 * PhiDOT * r * params.rDOT2 +
                    Complex(0, 197) * params.rDOT3) +
               3 * params.r2 *
                   (824 * params.PhiDOT5 * params.r5 -
                    Complex(0, 671) * params.PhiDOT4 * params.r4 * rDOT +
                    864 * params.PhiDOT3 * params.r3 * params.rDOT2 +
                    Complex(0, 44) * params.PhiDOT2 * params.r2 * params.rDOT3 +
                    320 * PhiDOT * r * params.rDOT4 +
                    Complex(0, 120) * params.rDOT5))) /
            (96096. * sqrt(3) * params.r2));
  }

  else if (vpnorder == 7) {

    return (/* (4*Nu*(S1z - delta*S1z + 2*(-2 + delta)*Nu*S1z + S2z + delta*S2z
      - 2*(2 + delta)*Nu*S2z + 2*params.eta2*(S1z +
      S2z))*params.x4p5)/(3003.*sqrt(3))
      + */
            /* Henry et al. ecc+spin terms */ (
                (params.Mtot2 * Nu *
                 (8 * params.Mtot2 *
                      (494 * PhiDOT * r + Complex(0, 411) * rDOT) -
                  12 * mass * r *
                      (530 * params.PhiDOT3 * params.r3 -
                       Complex(0, 6) * params.PhiDOT2 * params.r2 * rDOT +
                       453 * PhiDOT * r * params.rDOT2 +
                       Complex(0, 197) * params.rDOT3) +
                  3 * params.r2 *
                      (824 * params.PhiDOT5 * params.r5 -
                       Complex(0, 671) * params.PhiDOT4 * params.r4 * rDOT +
                       864 * params.PhiDOT3 * params.r3 * params.rDOT2 +
                       Complex(0, 44) * params.PhiDOT2 * params.r2 *
                           params.rDOT3 +
                       320 * PhiDOT * r * params.rDOT4 +
                       Complex(0, 120) * params.rDOT5)) *
                 ((1 - 4 * Nu + 2 * params.eta2 + delta * (-1 + 2 * Nu)) * S1z +
                  (1 + delta - 4 * Nu - 2 * delta * Nu + 2 * params.eta2) *
                      S2z)) /
                (48048. * sqrt(3) * params.r4)));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hl_7_m_2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_7_m_2: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           hGO_7_m_2(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, params) *
           cpolar(1, -2 * Phi);
  }
}

static COMPLEX16 hl_7_m_min2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_7_m_min2: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((- 4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           conj(hGO_7_m_2(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, params)) *
           cpolar(1, 2 * Phi);
  }
}

// 74
static COMPLEX16 hGO_7_m_4(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z, struct kepler_vars params) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  if (vpnorder == 6) {
    return ((mass * (-1 + 7 * Nu - 14 * params.eta2 + 7 * params.eta3) * PhiDOT *
             (8 * params.Mtot2 * (988 * PhiDOT * r + Complex(0, 411) * rDOT) +
              12 * mass * r *
                  (844 * params.PhiDOT3 * params.r3 +
                   Complex(0, 1518) * params.PhiDOT2 * params.r2 * rDOT -
                   906 * PhiDOT * r * params.rDOT2 -
                   Complex(0, 197) * params.rDOT3) -
              15 * params.r2 *
                  (656 * params.PhiDOT5 * params.r5 +
                   Complex(0, 179) * params.PhiDOT4 * params.r4 * rDOT +
                   192 * params.PhiDOT3 * params.r3 * params.rDOT2 +
                   Complex(0, 260) * params.PhiDOT2 * params.r2 * params.rDOT3 -
                   128 * PhiDOT * r * params.rDOT4 -
                   Complex(0, 24) * params.rDOT5))) /
            (21840. * sqrt(66) * params.r2));
  }

  else if (vpnorder == 7) {

    return (/* (-512*sqrt(0.06060606060606061)*Nu*(S1z - delta*S1z + 2*(-2 +
      delta)*Nu*S1z
      + S2z + delta*S2z - 2*(2 + delta)*Nu*S2z + 2*params.eta2*(S1z +
      S2z))*params.x4p5)/1365.
      + */
            /* Henry et al. ecc+spin terms */ (
                -0.00009157509157509158 *
                (params.Mtot2 * Nu *
                 (8 * params.Mtot2 *
                      (988 * PhiDOT * r + Complex(0, 411) * rDOT) +
                  12 * mass * r *
                      (844 * params.PhiDOT3 * params.r3 +
                       Complex(0, 1518) * params.PhiDOT2 * params.r2 * rDOT -
                       906 * PhiDOT * r * params.rDOT2 -
                       Complex(0, 197) * params.rDOT3) -
                  15 * params.r2 *
                      (656 * params.PhiDOT5 * params.r5 +
                       Complex(0, 179) * params.PhiDOT4 * params.r4 * rDOT +
                       192 * params.PhiDOT3 * params.r3 * params.rDOT2 +
                       Complex(0, 260) * params.PhiDOT2 * params.r2 *
                           params.rDOT3 -
                       128 * PhiDOT * r * params.rDOT4 -
                       Complex(0, 24) * params.rDOT5)) *
                 ((1 - 4 * Nu + 2 * params.eta2 + delta * (-1 + 2 * Nu)) * S1z +
                  (1 + delta - 4 * Nu - 2 * delta * Nu + 2 * params.eta2) *
                      S2z)) /
                (sqrt(66) * params.r4)));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hl_7_m_4(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_7_m_4: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           hGO_7_m_4(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, params) *
           cpolar(1, -4 * Phi);
  }
}

static COMPLEX16 hl_7_m_min4(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_7_m_min4: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((- 4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           conj(hGO_7_m_4(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, params)) *
           cpolar(1, 4 * Phi);
  }
}

// 76
static COMPLEX16 hGO_7_m_6(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z, struct kepler_vars params) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  if (vpnorder == 6) {
    return (
        -(mass * (-1 + 7 * Nu - 14 * params.eta2 + 7 * params.eta3) * PhiDOT *
          (8 * params.Mtot2 * (494 * PhiDOT * r + Complex(0, 137) * rDOT) +
           4 * mass * r *
               (6026 * params.PhiDOT3 * params.r3 +
                Complex(0, 4038) * params.PhiDOT2 * params.r2 * rDOT -
                1359 * PhiDOT * r * params.rDOT2 -
                Complex(0, 197) * params.rDOT3) +
           15 * params.r2 *
               (1240 * params.PhiDOT5 * params.r5 +
                Complex(0, 911) * params.PhiDOT4 * params.r4 * rDOT -
                544 * params.PhiDOT3 * params.r3 * params.rDOT2 -
                Complex(0, 236) * params.PhiDOT2 * params.r2 * params.rDOT3 +
                64 * PhiDOT * r * params.rDOT4 +
                Complex(0, 8) * params.rDOT5))) /
        (3360. * sqrt(429) * params.r2));
  } else if (vpnorder == 7) {

    return (/* (324*sqrt(0.02097902097902098)*Nu*(S1z - delta*S1z
      + 2*(-2 + delta)*Nu*S1z + S2z + delta*S2z - 2*(2 + delta)*Nu*S2z
      + 2*params.eta2*(S1z + S2z))*params.x4p5)/35.
      + */
            /* Henry et al. ecc+spin terms */ (
                (params.Mtot2 * Nu *
                 (8 * params.Mtot2 *
                      (494 * PhiDOT * r + Complex(0, 137) * rDOT) +
                  4 * mass * r *
                      (6026 * params.PhiDOT3 * params.r3 +
                       Complex(0, 4038) * params.PhiDOT2 * params.r2 * rDOT -
                       1359 * PhiDOT * r * params.rDOT2 -
                       Complex(0, 197) * params.rDOT3) +
                  15 * params.r2 *
                      (1240 * params.PhiDOT5 * params.r5 +
                       Complex(0, 911) * params.PhiDOT4 * params.r4 * rDOT -
                       544 * params.PhiDOT3 * params.r3 * params.rDOT2 -
                       Complex(0, 236) * params.PhiDOT2 * params.r2 *
                           params.rDOT3 +
                       64 * PhiDOT * r * params.rDOT4 +
                       Complex(0, 8) * params.rDOT5)) *
                 ((1 - 4 * Nu + 2 * params.eta2 + delta * (-1 + 2 * Nu)) * S1z +
                  (1 + delta - 4 * Nu - 2 * delta * Nu + 2 * params.eta2) *
                      S2z)) /
                (1680. * sqrt(429) * params.r4)));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hl_7_m_6(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_7_m_6: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           hGO_7_m_6(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, params) *
           cpolar(1, -6 * Phi);
  }
}

static COMPLEX16 hl_7_m_min6(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_7_m_min6: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((- 4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           conj(hGO_7_m_6(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, params)) *
           cpolar(1, 6 * Phi);
  }
}

// 88
static COMPLEX16 hGO_8_m_8(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder,struct kepler_vars params) {

  if (vpnorder == 6) {
    return (
        ((-1 + 7 * Nu - 14 * params.eta2 + 7 * params.eta3) *
         (9118 * params.Mtot4 +
          5040 * params.r4 * pow(PhiDOT * r + Complex(0, 1) * rDOT, 8) +
          4 * params.Mtot3 * r *
              (82543 * params.PhiDOT2 * params.r2 +
               Complex(0, 67760) * PhiDOT * r * rDOT - 16717 * params.rDOT2) +
          9 * params.Mtot2 * params.r2 *
              (124583 * params.PhiDOT4 * params.r4 +
               Complex(0, 192640) * params.PhiDOT3 * params.r3 * rDOT -
               144024 * params.PhiDOT2 * params.r2 * params.rDOT2 -
               Complex(0, 56280) * PhiDOT * r * params.rDOT3 +
               9188 * params.rDOT4) +
          315 * mass * params.r3 *
              (2005 * params.PhiDOT6 * params.r6 +
               Complex(0, 4250) * params.PhiDOT5 * params.r5 * rDOT -
               5538 * params.PhiDOT4 * params.r4 * params.rDOT2 -
               Complex(0, 4760) * params.PhiDOT3 * params.r3 * params.rDOT3 +
               2600 * params.PhiDOT2 * params.r2 * params.rDOT4 +
               Complex(0, 816) * PhiDOT * r * params.rDOT5 -
               112 * params.rDOT6))) /
        (2016. * sqrt(170170) * params.r4));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hl_8_m_8(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_8_m_8: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           hGO_8_m_8(mass, Nu, r, rDOT, PhiDOT, vpnorder, params) * cpolar(1, -8 * Phi);
  }
}

static COMPLEX16 hl_8_m_min8(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_8_m_min8: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           conj(hGO_8_m_8(mass, Nu, r, rDOT, PhiDOT, vpnorder, params)) *
           cpolar(1, 8 * Phi);
  }
}

// 86
static COMPLEX16 hGO_8_m_6(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder,struct kepler_vars params) {

  if (vpnorder == 6) {
    return (
        -((-1 + 7 * Nu - 14 * params.eta2 + 7 * params.eta3) *
          (18236 * params.Mtot4 -
           10080 * params.r4 * (PhiDOT * r - Complex(0, 1) * rDOT) *
               pow(PhiDOT * r + Complex(0, 1) * rDOT, 7) +
           8 * params.Mtot3 * r *
               (42923 * params.PhiDOT2 * params.r2 +
                Complex(0, 50820) * PhiDOT * r * rDOT - 16717 * params.rDOT2) +
           18 * params.Mtot2 * params.r2 *
               (15663 * params.PhiDOT4 * params.r4 +
                Complex(0, 58240) * params.PhiDOT3 * params.r3 * rDOT -
                74094 * params.PhiDOT2 * params.r2 * params.rDOT2 -
                Complex(0, 42210) * PhiDOT * r * params.rDOT3 +
                9188 * params.rDOT4) -
           315 * mass * params.r3 *
               (678 * params.PhiDOT6 * params.r6 +
                Complex(0, 681) * params.PhiDOT5 * params.r5 * rDOT +
                852 * params.PhiDOT4 * params.r4 * params.rDOT2 +
                Complex(0, 2660) * params.PhiDOT3 * params.r3 * params.rDOT3 -
                2640 * params.PhiDOT2 * params.r2 * params.rDOT4 -
                Complex(0, 1224) * PhiDOT * r * params.rDOT5 +
                224 * params.rDOT6))) /
        (10080. * sqrt(51051) * params.r4));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hl_8_m_6(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_8_m_6: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           hGO_8_m_6(mass, Nu, r, rDOT, PhiDOT, vpnorder, params) * cpolar(1, -6 * Phi);
  }
}

static COMPLEX16 hl_8_m_min6(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_8_m_min6: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           conj(hGO_8_m_6(mass, Nu, r, rDOT, PhiDOT, vpnorder, params)) *
           cpolar(1, 6 * Phi);
  }
}

// 84
static COMPLEX16 hGO_8_m_4(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, struct kepler_vars params) {

  if (vpnorder == 6) {
    return (
        ((-1 + 7 * Nu - 14 * params.eta2 + 7 *params.eta3) *
         (9118 * params.Mtot4 +
          5040 * params.r4 * pow(PhiDOT * r - Complex(0, 1) * rDOT, 2) *
              pow(PhiDOT * r + Complex(0, 1) * rDOT, 6) +
          4 * params.Mtot3 * r *
              (14623 * params.PhiDOT2 * params.r2 +
               Complex(0, 33880) * PhiDOT * r * rDOT - 16717 * params.rDOT2) -
          9 * params.Mtot2 * params.r2 *
              (8377 * params.PhiDOT4 * params.r4 +
               Complex(0, 2240) * params.PhiDOT3 * params.r3 * rDOT +
               24144 * params.PhiDOT2 * params.r2 * params.rDOT2 +
               Complex(0, 28140) * PhiDOT * r * params.rDOT3 -
               9188 * params.rDOT4) +
          45 * mass * params.r3 *
              (243 * params.PhiDOT6 * params.r6 -
               Complex(0, 1701) * params.PhiDOT5 * params.r5 * rDOT +
               3762 * params.PhiDOT4 * params.r4 * params.rDOT2 +
               Complex(0, 1260) * params.PhiDOT3 * params.r3 * params.rDOT3 +
               2840 * params.PhiDOT2 * params.r2 * params.rDOT4 +
               Complex(0, 2856) * PhiDOT * r * params.rDOT5 -
               784 * params.rDOT6))) /
        (65520. * sqrt(374) * params.r4));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hl_8_m_4(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_8_m_4: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           hGO_8_m_4(mass, Nu, r, rDOT, PhiDOT, vpnorder, params) * cpolar(1, -4 * Phi);
  }
}

static COMPLEX16 hl_8_m_min4(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_8_m_min4: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           conj(hGO_8_m_4(mass, Nu, r, rDOT, PhiDOT, vpnorder, params)) *
           cpolar(1, 4 * Phi);
  }
}

// 82
static COMPLEX16 hGO_8_m_2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, struct kepler_vars params) {

  if (vpnorder == 6) {
    return (
        ((-1 + 7 * Nu - 14 * params.eta2 + 7 * params.eta3) *
         (-18236 * params.Mtot4 +
          10080 * params.r4 * pow(PhiDOT * r - Complex(0, 1) * rDOT, 3) *
              pow(PhiDOT * r + Complex(0, 1) * rDOT, 5) +
          8 * params.Mtot3 * r *
              (2357 * params.PhiDOT2 * params.r2 -
               Complex(0, 16940) * PhiDOT * r * rDOT + 16717 * params.rDOT2) +
          18 * params.Mtot2 * params.r2 *
              (1297 * params.PhiDOT4 * params.r4 +
               Complex(0, 13440) * params.PhiDOT3 * params.r3 * rDOT -
               5826 * params.PhiDOT2 * params.r2 * params.rDOT2 +
               Complex(0, 14070) * PhiDOT * r * params.rDOT3 -
               9188 * params.rDOT4) -
          45 * mass * params.r3 *
              (758 * params.PhiDOT6 * params.r6 +
               Complex(0, 2891) * params.PhiDOT5 * params.r5 * rDOT +
               564 * params.PhiDOT4 * params.r4 * params.rDOT2 +
               Complex(0, 5740) * params.PhiDOT3 * params.r3 * params.rDOT3 -
               2000 * params.PhiDOT2 * params.r2 * params.rDOT4 +
               Complex(0, 2856) * PhiDOT * r * params.rDOT5 -
               1568 * params.rDOT6))) /
        (288288. * sqrt(85) * params.r4));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hl_8_m_2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_8_m_2: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           hGO_8_m_2(mass, Nu, r, rDOT, PhiDOT, vpnorder, params) * cpolar(1, -2 * Phi);
  }
}

static COMPLEX16 hl_8_m_min2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, struct kepler_vars params) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_8_m_min2: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           conj(hGO_8_m_2(mass, Nu, r, rDOT, PhiDOT, vpnorder, params)) *
           cpolar(1, 2 * Phi);
  }
}

static COMPLEX16 hlmGOresult(UINT4 l, INT4 m, REAL8 mass, REAL8 Nu, REAL8 r,
                             REAL8 rDOT, REAL8 Phi, REAL8 PhiDOT, REAL8 R,
                             UINT4 vpnorder, REAL8 S1z, REAL8 S2z, REAL8 x,
                             struct kepler_vars orbital_vars) {
  if (vpnorder > 8) {
    XLAL_ERROR(
        XLAL_EINVAL,
        "Error in hlmGOresult: Input PN order should be between [0, 8].");
  }

  COMPLEX16 hlm = 0;

  switch (l) {
  case 2:
    switch (m) {
    case 2:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_2_m_2(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x, orbital_vars);
      }
      return (hlm);

    case 1:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_2_m_1(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x, orbital_vars);
      }
      return (hlm);
    case 0:
      return (0);

    case -1:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_2_m_min1(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x, orbital_vars);
      }
      return (hlm);
    case -2:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_2_m_min2(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x, orbital_vars);
      }
      return (hlm);
    default:
      XLAL_ERROR(XLAL_EINVAL, "Error in m for l=2 of hlmGOresult: Input m "
                              "value should be between -2 to +2.");
      break;
    }
    break;

  case 3:
    switch (m) {
    case 3:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_3_m_3(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x, orbital_vars);
      }
      return (hlm);
    case 2:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_3_m_2(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x, orbital_vars);
      }
      return (hlm);
    case 1:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_3_m_1(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x, orbital_vars);
      }
      return (hlm);
    case 0:
      return (0);

    case -1:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_3_m_min1(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x, orbital_vars);
      }
      return (hlm);
    case -2:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_3_m_min2(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x, orbital_vars);
      }
      return (hlm);
    case -3:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_3_m_min3(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x, orbital_vars);
      }
      return (hlm);
    default:
      XLAL_ERROR(XLAL_EINVAL, "Error in m for l=3 of hlmGOresult: Input m "
                              "value should be between -3 to +3.");
      break;
    }
    break;

  case 4:
    switch (m) {
    case 4:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_4_m_4(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x, orbital_vars);
      }
      return (hlm);
    case 3:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_4_m_3(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x, orbital_vars);
      }
      return (hlm);
    case 2:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_4_m_2(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x, orbital_vars);
      }
      return (hlm);
    case 1:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_4_m_1(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x, orbital_vars);
      }
      return (hlm);
    case 0:
      return (0);

    case -1:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_4_m_min1(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x, orbital_vars);
      }
      return (hlm);
    case -2:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_4_m_min2(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x, orbital_vars);
      }
      return (hlm);
    case -3:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_4_m_min3(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x, orbital_vars);
      }
      return (hlm);
    case -4:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_4_m_min4(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x, orbital_vars);
      }
      return (hlm);
    default:
      XLAL_ERROR(XLAL_EINVAL, "Error in m for l=4 of hlmGOresult: Input m "
                              "value should be between -4 to +4.");
      break;
    }
    break;

  case 5:
    switch (m) {
    case 5:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_5_m_5(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x, orbital_vars);
      }
      return (hlm);
    case 4:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_5_m_4(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, orbital_vars);
      }
      return (hlm);
    case 3:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_5_m_3(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x, orbital_vars);
      }
      return (hlm);
    case 2:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_5_m_2(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, orbital_vars);
      }
      return (hlm);
    case 1:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_5_m_1(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x, orbital_vars);
      }
      return (hlm);
    case 0:
      return (0);
    case -1:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_5_m_min1(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x, orbital_vars);
      }
      return (hlm);
    case -2:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_5_m_min2(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, orbital_vars);
      }
      return (hlm);
    case -3:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_5_m_min3(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x, orbital_vars);
      }
      return (hlm);
    case -4:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_5_m_min4(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, orbital_vars);
      }
      return (hlm);
    case -5:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_5_m_min5(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x, orbital_vars);
      }
      return (hlm);
    default:
      XLAL_ERROR(XLAL_EINVAL, "Error in m for l=5 of hlmGOresult: Input m "
                              "value should be between -5 to +5.");
      break;
    }
    break;

  case 6:
    switch (m) {
    case 6:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_6_m_6(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, orbital_vars);
      }
      return (hlm);
    case 5:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_6_m_5(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, orbital_vars);
      }
      return (hlm);
    case 4:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_6_m_4(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, orbital_vars);
      }
      return (hlm);
    case 3:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_6_m_3(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, orbital_vars);
      }
      return (hlm);
    case 2:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_6_m_2(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, orbital_vars);
      }
      return (hlm);
    case 1:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_6_m_1(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, orbital_vars);
      }
      return (hlm);
    case 0:
      return (0);
    case -1:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_6_m_min1(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, orbital_vars);
      }
      return (hlm);
    case -2:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_6_m_min2(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, orbital_vars);
      }
      return (hlm);
    case -3:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_6_m_min3(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, orbital_vars);
      }
      return (hlm);
    case -4:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_6_m_min4(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, orbital_vars);
      }
      return (hlm);
    case -5:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_6_m_min5(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, orbital_vars);
      }
      return (hlm);
    case -6:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_6_m_min6(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, orbital_vars);
      }
      return (hlm);
    default:
      XLAL_ERROR(XLAL_EINVAL, "Error in m for l=6 of hlmGOresult: Input m "
                              "value should be between -6 to +6.");
      break;
    }
    break;

  case 7:
    switch (m) {
    case 7:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_7_m_7(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, orbital_vars);
      }
      return (hlm);
    case 6:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_7_m_6(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, orbital_vars);
      }
      return (hlm);
    case 5:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_7_m_5(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, orbital_vars);
      }
      return (hlm);
    case 4:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_7_m_4(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, orbital_vars);
      }
      return (hlm);
    case 3:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_7_m_3(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, orbital_vars);
      }
      return (hlm);
    case 2:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_7_m_2(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, orbital_vars);
      }
      return (hlm);
    case 1:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_7_m_1(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, orbital_vars);
      }
      return (hlm);
    case 0:
      return (0);
    case -1:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_7_m_min1(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, orbital_vars);
      }
      return (hlm);
    case -2:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_7_m_min2(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, orbital_vars);
      }
      return (hlm);
    case -3:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_7_m_min3(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, orbital_vars);
      }
      return (hlm);
    case -4:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_7_m_min4(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, orbital_vars);
      }
      return (hlm);
    case -5:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_7_m_min5(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, orbital_vars);
      }
      return (hlm);
    case -6:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_7_m_min6(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, orbital_vars);
      }
      return (hlm);
    case -7:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_7_m_min7(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, orbital_vars);
      }
      return (hlm);
    default:
      XLAL_ERROR(XLAL_EINVAL, "Error in m for l=7 of hlmGOresult: Input m "
                              "value should be between -7 to +7.");
      break;
    }
    break;

  case 8:
    switch (m) {
    case 8:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_8_m_8(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, orbital_vars);
      }
      return (hlm);
    case 7:
      return (0);
    case 6:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_8_m_6(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, orbital_vars);
      }
      return (hlm);
    case 5:
      return (0);
    case 4:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_8_m_4(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, orbital_vars);
      }
      return (hlm);
    case 3:
      return (0);
    case 2:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_8_m_2(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, orbital_vars);
      }
      return (hlm);
    case 1:
      return (0);
    case 0:
      return (0);
    case -1:
      return (0);
    case -2:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_8_m_min2(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, orbital_vars);
      }
      return (hlm);
    case -3:
      return (0);
    case -4:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_8_m_min4(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, orbital_vars);
      }
      return (hlm);
    case -5:
      return (0);
    case -6:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_8_m_min6(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, orbital_vars);
      }
      return (hlm);
    case -7:
      return (0);
    case -8:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_8_m_min8(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, orbital_vars);
      }
      return (hlm);
    default:
      XLAL_ERROR(XLAL_EINVAL, "Error in m for l=8 of hlmGOresult: Input m "
                              "value should be between -8 to +8.");
      break;
    }
    break;
  default:
    XLAL_ERROR(XLAL_EINVAL, "Error in input value of l in hlmGOresult: Input l "
                            "value should be between [2, 8].");
    break;
  }
}

// static COMPLEX16 h05PNGOresult(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
// REAL8 Phi, REAL8 PhiDOT, REAL8 inc, REAL8 euler_beta, REAL8 R, UINT4
// vpnorder, REAL8 S1z, REAL8 S2z, REAL8 x)
// {

//    //0 PN

//    if (vpnorder == 0)
//    {
//       COMPLEX16 hcombined = XLALSpinWeightedSphericalHarmonic(inc,
//       euler_beta, -2, 2, 2) * hl_2_m_2(mass, Nu, r, rDOT, Phi, PhiDOT, R, 0,
//       S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 2, -2) * hl_2_m_min2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 0, S1z, S2z, x);
//       return hcombined;
//    }

//    //0.5PN

//    else if (vpnorder == 1)
//    {
//       COMPLEX16 hcombined = XLALSpinWeightedSphericalHarmonic(inc,
//       euler_beta, -2, 2, 1) * hl_2_m_1(mass, Nu, r, rDOT, Phi, PhiDOT, R, 1,
//       S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 2, -1) * hl_2_m_min1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 1, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 3, 1) * hl_3_m_1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 1, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 3, -1) * hl_3_m_min1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 1, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 3, 3) * hl_3_m_3(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 1, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 3, -3) * hl_3_m_min3(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 1, S1z, S2z, x);
//       return hcombined;
//    }

//    //1PN
//    else if (vpnorder == 2)
//    {
//       COMPLEX16 hcombined = XLALSpinWeightedSphericalHarmonic(inc,
//       euler_beta, -2, 2, 2) * hl_2_m_2(mass, Nu, r, rDOT, Phi, PhiDOT, R, 2,
//       S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 2, -2) * hl_2_m_min2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 2, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 2, 1) * hl_2_m_1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 2, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 2, -1) * hl_2_m_min1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 2, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 3, 2) * hl_3_m_2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 2, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 3, -2) * hl_3_m_min2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 2, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, 2) * hl_4_m_2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 2, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, -2) * hl_4_m_min2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 2, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, 4) * hl_4_m_4(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 2, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, -4) * hl_4_m_min4(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 2, S1z, S2z, x);
//       return hcombined;
//    }

//    //1.5PN
//    else if (vpnorder == 3)
//    {
//       COMPLEX16 hcombined = XLALSpinWeightedSphericalHarmonic(inc,
//       euler_beta, -2, 2, 2) * hl_2_m_2(mass, Nu, r, rDOT, Phi, PhiDOT, R, 3,
//       S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 2, -2) * hl_2_m_min2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 3, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 2, 1) * hl_2_m_1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 3, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 2, -1) * hl_2_m_min1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 3, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 3, 3) * hl_3_m_3(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 3, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 3, -3) * hl_3_m_min3(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 3, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 5, 5) * hl_5_m_5(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 3, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 5, -5) * hl_5_m_min5(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 3, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 3, 2) * hl_3_m_2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 3, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 3, -2) * hl_3_m_min2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 3, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 3, 1) * hl_3_m_1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 3, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 3, -1) * hl_3_m_min1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 3, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, 3) * hl_4_m_3(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 3, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, -3) * hl_4_m_min3(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 3, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, 1) * hl_4_m_1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 3, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, -1) * hl_4_m_min1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 3, S1z, S2z, x);
//       return hcombined;
//    }

//    //2PN
//    else if (vpnorder == 4)
//    {
//       COMPLEX16 hcombined = XLALSpinWeightedSphericalHarmonic(inc,
//       euler_beta, -2, 2, 2) * hl_2_m_2(mass, Nu, r, rDOT, Phi, PhiDOT, R, 4,
//       S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 2, -2) * hl_2_m_min2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 4, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 2, 1) * hl_2_m_1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 4, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 2, -1) * hl_2_m_min1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 4, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 3, 3) * hl_3_m_3(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 4, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 3, -3) * hl_3_m_min3(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 4, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 3, 2) * hl_3_m_2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 4, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 3, -2) * hl_3_m_min2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 4, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 3, 1) * hl_3_m_1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 4, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 3, -1) * hl_3_m_min1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 4, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, 3) * hl_4_m_3(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 4, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, -3) * hl_4_m_min3(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 4, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, 2) * hl_4_m_2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 4, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, -2) * hl_4_m_min2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 4, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, 1) * hl_4_m_1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 4, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, -1) * hl_4_m_min1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 4, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 5, 2) * hl_5_m_2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 4, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 5, -2) * hl_5_m_min2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 4, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 6, 2) * hl_6_m_2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 4, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 6, -2) * hl_6_m_min2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 4, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, 4) * hl_4_m_4(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 4, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, -4) * hl_4_m_min4(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 4, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 5, 4) * hl_5_m_4(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 4, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 5, -4) * hl_5_m_min4(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 4, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 6, 4) * hl_6_m_4(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 4, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 6, -4) * hl_6_m_min4(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 4, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 6, 6) * hl_6_m_6(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 4, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 6, -6) * hl_6_m_min6(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 4, S1z, S2z);

//       return hcombined;
//    }

//    //2.5PN
//    else if (vpnorder == 5)
//    {
//       COMPLEX16 hcombined = XLALSpinWeightedSphericalHarmonic(inc,
//       euler_beta, -2, 6, 5) * hl_6_m_5(mass, Nu, r, rDOT, Phi, PhiDOT, R, 5,
//       S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 6, -5) * hl_6_m_min5(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 6, 3) * hl_6_m_3(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 6, -3) * hl_6_m_min3(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 6, 1) * hl_6_m_1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 6, -1) * hl_6_m_min1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 5, 5) * hl_5_m_5(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 5, -5) * hl_5_m_min5(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 5, 4) * hl_5_m_4(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 5, -4) * hl_5_m_min4(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 5, 1) * hl_5_m_1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 5, -1) * hl_5_m_min1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 5, 3) * hl_5_m_3(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 5, -3) * hl_5_m_min3(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 5, 2) * hl_5_m_2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 5, -2) * hl_5_m_min2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 7, 7) * hl_7_m_7(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 7, -7) * hl_7_m_min7(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 7, 5) * hl_7_m_5(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 7, -5) * hl_7_m_min5(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 7, 3) * hl_7_m_3(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 7, -3) * hl_7_m_min3(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 7, 1) * hl_7_m_1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 7, -1) * hl_7_m_min1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, 1) * hl_4_m_1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, -1) * hl_4_m_min1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, 2) * hl_4_m_2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, -2) * hl_4_m_min2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, 3) * hl_4_m_3(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, -3) * hl_4_m_min3(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, 4) * hl_4_m_4(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, -4) * hl_4_m_min4(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 3, 3) * hl_3_m_3(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 3, -3) * hl_3_m_min3(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 3, 2) * hl_3_m_2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 3, -2) * hl_3_m_min2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 3, 1) * hl_3_m_1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 3, -1) * hl_3_m_min1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 2, 1) * hl_2_m_1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 2, -1) * hl_2_m_min1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 2, 2) * hl_2_m_2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 2, -2) * hl_2_m_min2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 5, S1z, S2z, x);
//       return hcombined;
//    }
//    //3PN
//    else if (vpnorder == 6)
//    {
//       COMPLEX16 hcombined = XLALSpinWeightedSphericalHarmonic(inc,
//       euler_beta, -2, 6, 2) * hl_6_m_2(mass, Nu, r, rDOT, Phi, PhiDOT, R, 6,
//       S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 6, -2) * hl_6_m_min2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 6, 3) * hl_6_m_3(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 6, -3) * hl_6_m_min3(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 6, 1) * hl_6_m_1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 6, -1) * hl_6_m_min1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 6, 5) * hl_6_m_5(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 6, -5) * hl_6_m_min5(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 6, 4) * hl_6_m_4(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 6, -4) * hl_6_m_min4(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 6, 6) * hl_6_m_6(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 6, -6) * hl_6_m_min6(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 5, 5) * hl_5_m_5(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 5, -5) * hl_5_m_min5(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 5, 1) * hl_5_m_1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 5, -1) * hl_5_m_min1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 5, 3) * hl_5_m_3(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 5, -3) * hl_5_m_min3(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 5, 2) * hl_5_m_2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 5, -2) * hl_5_m_min2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 5, 4) * hl_5_m_4(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 5, -4) * hl_5_m_min4(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 7, 6) * hl_7_m_6(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 7, -6) * hl_7_m_min6(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 7, 4) * hl_7_m_4(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 7, -4) * hl_7_m_min4(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 7, 2) * hl_7_m_2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 7, -2) * hl_7_m_min2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, 1) * hl_4_m_1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, -1) * hl_4_m_min1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, 2) * hl_4_m_2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, -2) * hl_4_m_min2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, 3) * hl_4_m_3(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, -3) * hl_4_m_min3(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, 4) * hl_4_m_4(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, -4) * hl_4_m_min4(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 3, 3) * hl_3_m_3(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 3, -3) * hl_3_m_min3(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 3, 2) * hl_3_m_2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 3, -2) * hl_3_m_min2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 3, 1) * hl_3_m_1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 3, -1) * hl_3_m_min1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 2, 1) * hl_2_m_1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 2, -1) * hl_2_m_min1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 2, 2) * hl_2_m_2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 2, -2) * hl_2_m_min2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 8, 8) * hl_8_m_8(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 8, -8) * hl_8_m_min8(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 8, 6) * hl_8_m_6(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 8, -6) * hl_8_m_min6(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 8, 4) * hl_8_m_4(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 8, -4) * hl_8_m_min4(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 8, 2) * hl_8_m_2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 8, -2) * hl_8_m_min2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 6);

//       return hcombined;
//    }

//    else if (vpnorder == 7)
//    {
//       COMPLEX16 hcombined = XLALSpinWeightedSphericalHarmonic(inc,
//       euler_beta, -2, 6, 1) * hl_6_m_1(mass, Nu, r, rDOT, Phi, PhiDOT, R, 7,
//       S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 6, -1) * hl_6_m_min1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 6, 2) * hl_6_m_2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 6, -2) * hl_6_m_min2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 6, 3) * hl_6_m_3(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 6, -3) * hl_6_m_min3(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 6, 4) * hl_6_m_4(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 6, -4) * hl_6_m_min4(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 6, 5) * hl_6_m_5(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 6, -5) * hl_6_m_min5(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 6, 6) * hl_6_m_6(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 6, -6) * hl_6_m_min6(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 5, 5) * hl_5_m_5(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 5, -5) * hl_5_m_min5(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 5, 1) * hl_5_m_1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 5, -1) * hl_5_m_min1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 5, 3) * hl_5_m_3(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 5, -3) * hl_5_m_min3(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 5, 2) * hl_5_m_2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 5, -2) * hl_5_m_min2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 5, 4) * hl_5_m_4(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 5, -4) * hl_5_m_min4(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 7, 7) * hl_7_m_7(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 7, -7) * hl_7_m_min7(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 7, 6) * hl_7_m_6(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 7, -6) * hl_7_m_min6(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 7, 5) * hl_7_m_5(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 7, -5) * hl_7_m_min5(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 7, 4) * hl_7_m_4(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 7, -4) * hl_7_m_min4(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 7, 3) * hl_7_m_3(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 7, -3) * hl_7_m_min3(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 7, 2) * hl_7_m_2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 7, -2) * hl_7_m_min2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 7, 1) * hl_7_m_1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 7, -1) * hl_7_m_min1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, 1) * hl_4_m_1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, -1) * hl_4_m_min1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, 2) * hl_4_m_2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, -2) * hl_4_m_min2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, 3) * hl_4_m_3(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, -3) * hl_4_m_min3(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, 4) * hl_4_m_4(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 4, -4) * hl_4_m_min4(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 3, 3) * hl_3_m_3(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 3, -3) * hl_3_m_min3(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 3, 2) * hl_3_m_2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 3, -2) * hl_3_m_min2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 3, 1) * hl_3_m_1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 3, -1) * hl_3_m_min1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 2, 1) * hl_2_m_1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 2, -1) * hl_2_m_min1(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 2, 2) * hl_2_m_2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 2, -2) * hl_2_m_min2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 7, S1z, S2z, x);

//       return hcombined;
//    }
//    //4PN quasi-circular (2,2) mode terms
//    else
//    {
//       COMPLEX16 hcombined = XLALSpinWeightedSphericalHarmonic(inc,
//       euler_beta, -2, 2, 2) * hl_2_m_2(mass, Nu, r, rDOT, Phi, PhiDOT, R, 8,
//       S1z, S2z, x) +
//                             XLALSpinWeightedSphericalHarmonic(inc,
//                             euler_beta, -2, 2, -2) * hl_2_m_min2(mass, Nu, r,
//                             rDOT, Phi, PhiDOT, R, 8, S1z, S2z, x);
//       return hcombined;
//    }
// }

// static REAL8 hplusGO(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
// REAL8 PhiDOT, REAL8 inc, REAL8 euler_beta, REAL8 R, UINT4 vpnorder, REAL8
// S1z, REAL8 S2z, REAL8 x)
// {
//    // nu = symmetric mass ratio m1m2/M^2
//    REAL8 hplus = creal(h05PNGOresult(mass, Nu, r, rDOT, Phi, PhiDOT, inc,
//    euler_beta, R, vpnorder, S1z, S2z, x)); return hplus;
// }

// static REAL8 hcrossGO(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
// REAL8 PhiDOT, REAL8 inc, REAL8 euler_beta, REAL8 R, UINT4 vpnorder, REAL8
// S1z, REAL8 S2z, REAL8 x)
// {
//    // nu = symmetric mass ratio m1m2/M^2
//    REAL8 hcross = -1.0 * cimag(h05PNGOresult(mass, Nu, r, rDOT, Phi, PhiDOT,
//    inc, euler_beta, R, vpnorder, S1z, S2z, x)); return hcross;
// }
// End of this file