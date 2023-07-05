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
                           REAL8 x) {
  REAL8 kappa1 = 1.0; /*for black holes kappa and lambda is 1*/
  REAL8 kappa2 = 1.0;
  REAL8 lambda1 = 1.0;
  REAL8 lambda2 = 1.0;
  REAL8 delta = sqrt(1 - 4 * Nu);

  if (vpnorder == 0) {
    return (mass / r + pow(PhiDOT, 2) * pow(r, 2) +
            Complex(0, 2) * PhiDOT * r * rDOT - pow(rDOT, 2));
  }

  else if (vpnorder == 2) {
    return (
        (21 * pow(mass, 2) * (-10 + Nu) -
         27 * (-1 + 3 * Nu) * pow(r, 2) * (PhiDOT * r - Complex(0, 1) * rDOT) *
             pow(PhiDOT * r + Complex(0, 1) * rDOT, 3) +
         mass * r *
             ((11 + 156 * Nu) * pow(PhiDOT, 2) * pow(r, 2) +
              Complex(0, 10) * (5 + 27 * Nu) * PhiDOT * r * rDOT -
              3 * (15 + 32 * Nu) * pow(rDOT, 2))) /
        (42. * pow(r, 2)));
  }

  else if (vpnorder == 3) {
    return (
        (pow(mass, 2) * (Complex(0, -1) * rDOT *
                             ((3 + 3 * delta - 8 * Nu) * S1z +
                              (3 - 3 * delta - 8 * Nu) * S2z) +
                         PhiDOT * r *
                             ((-3 - 3 * delta + 5 * Nu) * S1z +
                              (-3 + 3 * delta + 5 * Nu) * S2z))) /
            (3. * pow(r, 2)) /* (<--This is the general orbit term) */
        /* (This is the quasi-circular limit of the general orbit term-->) */
        - ((-4 * ((1 + delta - Nu) * S1z + S2z - (delta + Nu) * S2z) *
            pow(x, 2.5)) /
           3.) +
        ((-4 * (S1z + delta * S1z + S2z - delta * S2z - Nu * (S1z + S2z)) *
          pow(x, 2.5)) /
         3.)); /* (<--This is Quentins quasi-circular term) */
  }

  else if (vpnorder == 4) {
    return ((6 * pow(mass, 3) * (3028 + 1267 * Nu + 158 * pow(Nu, 2)) +
             9 * (83 - 589 * Nu + 1111 * pow(Nu, 2)) * pow(r, 3) *
                 pow(PhiDOT * r - Complex(0, 1) * rDOT, 2) *
                 pow(PhiDOT * r + Complex(0, 1) * rDOT, 4) +
             pow(mass, 2) * r *
                 ((-11891 - 36575 * Nu + 13133 * pow(Nu, 2)) * pow(PhiDOT, 2) *
                      pow(r, 2) +
                  Complex(0, 8) * (-773 - 3767 * Nu + 2852 * pow(Nu, 2)) *
                      PhiDOT * r * rDOT -
                  6 * (-619 + 2789 * Nu + 934 * pow(Nu, 2)) * pow(rDOT, 2)) -
             3 * mass * pow(r, 2) *
                 (2 * (-835 - 19 * Nu + 2995 * pow(Nu, 2)) * pow(PhiDOT, 4) *
                      pow(r, 4) +
                  Complex(0, 6) * (-433 - 721 * Nu + 1703 * pow(Nu, 2)) *
                      pow(PhiDOT, 3) * pow(r, 3) * rDOT +
                  6 * (-33 + 1014 * Nu + 232 * pow(Nu, 2)) * pow(PhiDOT, 2) *
                      pow(r, 2) * pow(rDOT, 2) +
                  Complex(0, 4) * (-863 + 1462 * Nu + 2954 * pow(Nu, 2)) *
                      PhiDOT * r * pow(rDOT, 3) -
                  3 * (-557 + 664 * Nu + 1712 * pow(Nu, 2)) * pow(rDOT, 4))) /
                (1512. * pow(r, 3)) +
            (3 * pow(mass, 3) *
             (S1z * (4 * Nu * S2z + (1 + delta - 2 * Nu) * S1z * kappa1) -
              (-1 + delta + 2 * Nu) * pow(S2z, 2) * kappa2)) /
                (4. * pow(r, 3)) -
            ((kappa1 * (1 + delta - 2 * Nu) * pow(S1z, 2) +
              S2z * (4 * Nu * S1z - kappa2 * (-1 + delta + 2 * Nu) * S2z)) *
             pow(x, 3)) +
            ((kappa1 * (1 + delta - 2 * Nu) * pow(S1z, 2) +
              S2z * (4 * Nu * S1z - kappa2 * (-1 + delta + 2 * Nu) * S2z)) *
             pow(x, 3)));
  }

  else if (vpnorder == 5) {
    return (
        (pow(mass, 2) * Nu *
         (2 * mass * (Complex(0, -702) * PhiDOT * r + rDOT) +
          3 * r *
              (Complex(0, -316) * pow(PhiDOT, 3) * pow(r, 3) -
               847 * pow(PhiDOT, 2) * pow(r, 2) * rDOT +
               Complex(0, 184) * PhiDOT * r * pow(rDOT, 2) -
               122 * pow(rDOT, 3)))) /
            (105. * pow(r, 3)) /* Henry et al. QC spin terms */ /* +
((2*(56*delta*Nu*(-S1z + S2z)
+ 101*Nu*(S1z + S2z) + 132*pow(Nu,2)*(S1z + S2z) - 80*(S1z + delta*S1z + S2z -
delta*S2z))*pow(x,3.5))/63.) */
        +
        /* Henry et al. ecc spin terms */ (
            (pow(mass, 2) *
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
              r * (PhiDOT * r * pow(rDOT, 2) *
                       (-((18 * (1 + delta) + 5 * (-63 + 55 * delta) * Nu +
                           188 * pow(Nu, 2)) *
                          S1z) +
                        (18 * (-1 + delta) + 5 * (63 + 55 * delta) * Nu -
                         188 * pow(Nu, 2)) *
                            S2z) -
                   Complex(0, 2) * pow(rDOT, 3) *
                       ((-27 * (1 + delta) + 6 * (5 + 7 * delta) * Nu -
                         4 * pow(Nu, 2)) *
                            S1z +
                        (-27 + 27 * delta + 30 * Nu - 42 * delta * Nu -
                         4 * pow(Nu, 2)) *
                            S2z) +
                   Complex(0, 2) * pow(PhiDOT, 2) * pow(r, 2) * rDOT *
                       ((51 + 88 * Nu * (-3 + 5 * Nu) +
                         delta * (51 + 62 * Nu)) *
                            S1z +
                        (51 + 88 * Nu * (-3 + 5 * Nu) -
                         delta * (51 + 62 * Nu)) *
                            S2z) +
                   pow(PhiDOT, 3) * pow(r, 3) *
                       ((120 * (1 + delta) + (-483 + 83 * delta) * Nu +
                         234 * pow(Nu, 2)) *
                            S1z +
                        (120 + 3 * Nu * (-161 + 78 * Nu) -
                         delta * (120 + 83 * Nu)) *
                            S2z)))) /
            (84. * pow(r, 3))));
  }

  else if (vpnorder == 6) {
    return (
        (4 * pow(mass, 4) *
             (-8203424 + 2180250 * pow(Nu, 2) + 592600 * pow(Nu, 3) +
              15 * Nu * (-5503804 + 142065 * pow(M_PI, 2))) -
         2700 * (-507 + 6101 * Nu - 25050 * pow(Nu, 2) + 34525 * pow(Nu, 3)) *
             pow(r, 4) * pow(PhiDOT * r - Complex(0, 1) * rDOT, 3) *
             pow(PhiDOT * r + Complex(0, 1) * rDOT, 5) +
         pow(mass, 3) * r *
             (pow(PhiDOT, 2) *
                  (337510808 - 198882000 * pow(Nu, 2) + 56294600 * pow(Nu, 3) +
                   Nu * (183074880 - 6392925 * pow(M_PI, 2))) *
                  pow(r, 2) +
              Complex(0, 110) * PhiDOT *
                  (-5498800 - 785120 * pow(Nu, 2) + 909200 * pow(Nu, 3) +
                   3 * Nu * (-1849216 + 38745 * pow(M_PI, 2))) *
                  r * rDOT +
              2 *
                  (51172744 - 94929000 * pow(Nu, 2) - 5092400 * pow(Nu, 3) +
                   45 * Nu * (2794864 + 142065 * pow(M_PI, 2))) *
                  pow(rDOT, 2)) -
         20 * pow(mass, 2) * pow(r, 2) *
             ((-986439 + 1873255 * Nu - 9961400 * pow(Nu, 2) +
               6704345 * pow(Nu, 3)) *
                  pow(PhiDOT, 4) * pow(r, 4) +
              Complex(0, 4) *
                  (-273687 - 978610 * Nu - 4599055 * pow(Nu, 2) +
                   2783005 * pow(Nu, 3)) *
                  pow(PhiDOT, 3) * pow(r, 3) * rDOT +
              (-181719 + 19395325 * Nu + 8237980 * pow(Nu, 2) +
               2612735 * pow(Nu, 3)) *
                  pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) +
              Complex(0, 8) *
                  (-234312 + 1541140 * Nu + 1230325 * pow(Nu, 2) +
                   1828625 * pow(Nu, 3)) *
                  PhiDOT * r * pow(rDOT, 3) -
              3 *
                  (-370268 + 1085140 * Nu + 2004715 * pow(Nu, 2) +
                   1810425 * pow(Nu, 3)) *
                  pow(rDOT, 4)) +
         300 * mass * pow(r, 3) *
             (4 *
                  (12203 - 36427 * Nu - 27334 * pow(Nu, 2) +
                   149187 * pow(Nu, 3)) *
                  pow(PhiDOT, 6) * pow(r, 6) +
              Complex(0, 2) *
                  (44093 - 68279 * Nu - 295346 * pow(Nu, 2) +
                   541693 * pow(Nu, 3)) *
                  pow(PhiDOT, 5) * pow(r, 5) * rDOT +
              2 *
                  (27432 - 202474 * Nu + 247505 * pow(Nu, 2) +
                   394771 * pow(Nu, 3)) *
                  pow(PhiDOT, 4) * pow(r, 4) * pow(rDOT, 2) +
              Complex(0, 2) *
                  (97069 - 383990 * Nu - 8741 * pow(Nu, 2) +
                   1264800 * pow(Nu, 3)) *
                  pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 3) +
              (-42811 + 53992 * Nu + 309136 * pow(Nu, 2) -
               470840 * pow(Nu, 3)) *
                  pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 4) +
              Complex(0, 2) *
                  (51699 - 252256 * Nu + 131150 * pow(Nu, 2) +
                   681160 * pow(Nu, 3)) *
                  PhiDOT * r * pow(rDOT, 5) -
              3 *
                  (16743 - 75104 * Nu + 26920 * pow(Nu, 2) +
                   207200 * pow(Nu, 3)) *
                  pow(rDOT, 6))) /
            (3.3264e6 * pow(r, 4)) /* Henry et al. QC spin terms */ /* + (((4*(1
+ delta)*(-7 + 9*kappa1)
- 7*(9 + 17*delta)*Nu - 9*(15 + 7*delta)*kappa1*Nu + 12*(7 -
17*kappa1)*pow(Nu,2))*pow(S1z,2) + 2*S1z*(Complex(0,-42)*(1 + delta - 2*Nu) -
84*(1 + delta - Nu)*M_PI + Nu*(-271 + 288*Nu)*S2z) + S2z*(12*(7 -
17*kappa2)*pow(Nu,2)*S2z + 4*(-1 + delta)*(Complex(0,21) + 42*M_PI + 7*S2z -
9*kappa2*S2z) + Nu*(168*(Complex(0,1) + M_PI) + 7*delta*(17 + 9*kappa2)*S2z -
9*(7 + 15*kappa2)*S2z)))*pow(x,4))/63. */
        + /* Henry et al. ecc spin terms */ (
              -0.005952380952380952 *
              (pow(mass, 3) *
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
                         pow(S2z, 2)) +
                r * (pow(rDOT, 2) *
                         (Complex(0, 56) * (1 + delta - 2 * Nu) * S1z +
                          (-56 * (1 + delta - 2 * Nu) * Nu +
                           kappa1 * (291 * (1 + delta) -
                                     2 * (445 + 154 * delta) * Nu +
                                     24 * pow(Nu, 2))) *
                              pow(S1z, 2) +
                          4 * Nu * (-3 + 44 * Nu) * S1z * S2z +
                          S2z * (Complex(0, -56) * (-1 + delta + 2 * Nu) +
                                 56 * Nu * (-1 + delta + 2 * Nu) * S2z +
                                 kappa2 *
                                     (291 - 890 * Nu + 24 * pow(Nu, 2) +
                                      delta * (-291 + 308 * Nu)) *
                                     S2z)) +
                     pow(PhiDOT, 2) * pow(r, 2) *
                         (Complex(0, 196) * (1 + delta - 2 * Nu) * S1z +
                          (56 * Nu * (7 + 7 * delta + Nu) +
                           kappa1 * (-153 * (1 + delta) -
                                     2 * (62 + 215 * delta) * Nu +
                                     804 * pow(Nu, 2))) *
                              pow(S1z, 2) +
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
                               2 * (14 - 195 * kappa1) * pow(Nu, 2)) *
                              pow(S1z, 2) +
                          4 * S1z *
                              (35 * (1 + delta - 2 * Nu) +
                               Complex(0, 1) * (9 - 209 * Nu) * Nu * S2z) +
                          S2z * (-140 * (-1 + delta + 2 * Nu) +
                                 Complex(0, 1) *
                                     (-14 * Nu * (-31 + 31 * delta + 2 * Nu) +
                                      kappa2 * (117 * (-1 + delta) +
                                                (23 + 211 * delta) * Nu +
                                                390 * pow(Nu, 2))) *
                                     S2z))))) /
              pow(r, 4)) + /* Henry et al. QC spinning hereditary terms  */
        (((-8 * M_PI * ((1 + delta - Nu) * S1z + S2z - (delta + Nu) * S2z) *
           pow(x, 4)) /
          3.)));
  } else if (vpnorder == 7) {

    return (
        /* Henry et al QC spin terms */ /* ((3318*pow(Nu,3)*(S1z + S2z) +
    Nu*(-504*((7 + delta)*kappa1 - 3*(3 + delta)*lambda1)*pow(S1z,3) -
    1008*pow(S1z,2)*(3*kappa1*M_PI - 3*(1 + delta)*S2z + 2*(1 +
    delta)*kappa1*S2z) + S1z*(17387 + 20761*delta + 1008*S2z*(6*M_PI + (-1 +
    delta)*(-3 + 2*kappa2)*S2z)) + S2z*(17387 - 20761*delta +
    504*S2z*(-6*kappa2*M_PI + (-7 + delta)*kappa2*S2z - 3*(-3 +
    delta)*lambda2*S2z))) + 2*(2809*(1 + delta)*S1z + 756*(1 +
    delta)*kappa1*M_PI*pow(S1z,2) + 756*(1 + delta)*(kappa1 -
    lambda1)*pow(S1z,3) -
    (-1 + delta)*S2z*(2809 + 756*S2z*(-(lambda2*S2z) + kappa2*(M_PI + S2z)))) -
    2*pow(Nu,2)*(708*delta*(-S1z + S2z) + (S1z + S2z)*(4427 +
    1008*(kappa1*pow(S1z,2) + S2z*(-2*S1z + kappa2*S2z)))))*pow(x,4.5))/756. */
        /* + */ /* Henry et al. ecc+spin terms */ (
            (pow(mass, 2) *
             (-3 * mass * r *
                  (Complex(0, -16) * pow(rDOT, 3) *
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
                   4 * PhiDOT * r * pow(rDOT, 2) *
                       (-338520 * pow(Nu, 3) * (S1z + S2z) +
                        48930 * (S1z + delta * S1z + S2z - delta * S2z) +
                        Nu * (Complex(0, 3420696) -
                              35 * (14154 + 21167 * delta) * S1z -
                              495390 * S2z + 740845 * delta * S2z) +
                        pow(Nu, 2) * (Complex(0, -612336) +
                                      245 * (3566 - 1565 * delta) * S1z +
                                      245 * (3566 + 1565 * delta) * S2z)) +
                   pow(PhiDOT, 3) * pow(r, 3) *
                       (2515380 * pow(Nu, 3) * (S1z + S2z) -
                        5 * Nu *
                            (Complex(0, 1859936) +
                             7 * (82329 + 37061 * delta) * S1z +
                             7 * (82329 - 37061 * delta) * S2z) -
                        128100 * (S1z + delta * S1z + S2z - delta * S2z) +
                        4 * pow(Nu, 2) *
                            (Complex(0, 381348) +
                             35 * (-18505 + 1777 * delta) * S1z -
                             35 * (18505 + 1777 * delta) * S2z)) +
                   Complex(0, 8) * pow(PhiDOT, 2) * pow(r, 2) * rDOT *
                       (779100 * pow(Nu, 3) * (S1z + S2z) +
                        5 * Nu *
                            (Complex(0, 828806) +
                             7 * (4839 + 5971 * delta) * S1z +
                             7 * (4839 - 5971 * delta) * S2z) -
                        62475 * (S1z + delta * S1z + S2z - delta * S2z) +
                        pow(Nu, 2) * (Complex(0, -976002) +
                                      35 * (-29599 + 3109 * delta) * S1z -
                                      35 * (29599 + 3109 * delta) * S2z))) +
              3 * pow(r, 2) *
                  (Complex(0, 4) * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 3) *
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
                   PhiDOT * r * pow(rDOT, 4) *
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
                   Complex(0, 6) * pow(rDOT, 5) *
                       (-11760 * pow(Nu, 3) * (S1z + S2z) +
                        4 * pow(Nu, 2) *
                            (Complex(0, 57338) +
                             35 * (569 + 301 * delta) * S1z +
                             35 * (569 - 301 * delta) * S2z) -
                        16 * Nu *
                            (Complex(0, -2517) + 35 * (72 + 71 * delta) * S1z +
                             35 * (72 - 71 * delta) * S2z) +
                        8715 * (S1z + delta * S1z + S2z - delta * S2z)) +
                   pow(PhiDOT, 5) * pow(r, 5) *
                       (2263380 * pow(Nu, 3) * (S1z + S2z) +
                        3 * Nu *
                            (Complex(0, 653432) +
                             35 * (13283 + 7839 * delta) * S1z +
                             35 * (13283 - 7839 * delta) * S2z) -
                        219240 * (S1z + delta * S1z + S2z - delta * S2z) +
                        4 * pow(Nu, 2) *
                            (Complex(0, -291268) +
                             105 * (-7669 + 165 * delta) * S1z -
                             105 * (7669 + 165 * delta) * S2z)) +
                   Complex(0, 14) * pow(PhiDOT, 4) * pow(r, 4) * rDOT *
                       (385170 * pow(Nu, 3) * (S1z + S2z) +
                        15 * Nu *
                            (Complex(0, 16914) + (8633 + 2267 * delta) * S1z +
                             8633 * S2z - 2267 * delta * S2z) -
                        18630 * (S1z + delta * S1z + S2z - delta * S2z) +
                        2 * pow(Nu, 2) *
                            (Complex(0, -67904) +
                             15 * (-13932 + 679 * delta) * S1z -
                             15 * (13932 + 679 * delta) * S2z)) +
                   6 * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) *
                       (-6720 * pow(Nu, 3) * (S1z + S2z) +
                        5 * Nu *
                            (Complex(0, -241704) +
                             7 * (4377 + 3083 * delta) * S1z +
                             7 * (4377 - 3083 * delta) * S2z) -
                        36960 * (S1z + delta * S1z + S2z - delta * S2z) +
                        pow(Nu, 2) * (Complex(0, -364384) +
                                      35 * (5059 - 4753 * delta) * S1z +
                                      35 * (5059 + 4753 * delta) * S2z))) +
              14 * pow(mass, 2) *
                  (Complex(0, 30) * rDOT *
                       (15280 * pow(Nu, 3) * (S1z + S2z) -
                        4 * pow(Nu, 2) *
                            (Complex(0, -111) + 3562 * S1z + 682 * delta * S1z +
                             945 * kappa1 * pow(S1z, 3) +
                             (3562 - 682 * delta +
                              945 * (-2 + kappa1) * pow(S1z, 2)) *
                                 S2z +
                             945 * (-2 + kappa2) * S1z * pow(S2z, 2) +
                             945 * kappa2 * pow(S2z, 3)) +
                        Nu * (Complex(0, -27520) +
                              378 *
                                  (5 * (1 + delta) * kappa1 +
                                   3 * (3 + delta) * lambda1) *
                                  pow(S1z, 3) +
                              (29749 - 13605 * delta) * S2z -
                              1512 * (1 + delta) * kappa1 * pow(S1z, 2) * S2z -
                              378 *
                                  (5 * (-1 + delta) * kappa2 +
                                   3 * (-3 + delta) * lambda2) *
                                  pow(S2z, 3) +
                              S1z * (29749 + 13605 * delta +
                                     1512 * (-1 + delta) * kappa2 *
                                         pow(S2z, 2))) +
                        2 * (-8009 * (1 + delta) * S1z -
                             567 * (1 + delta) * lambda1 * pow(S1z, 3) +
                             (-1 + delta) * S2z *
                                 (8009 + 567 * lambda2 * pow(S2z, 2)))) +
                   PhiDOT * r *
                       (285840 * pow(Nu, 3) * (S1z + S2z) -
                        30 * pow(Nu, 2) *
                            (Complex(0, 11504) + 1890 * kappa1 * pow(S1z, 3) +
                             23090 * S2z - 1823 * delta * S2z +
                             1890 * (-2 + kappa1) * pow(S1z, 2) * S2z +
                             1890 * kappa2 * pow(S2z, 3) +
                             S1z * (23090 + 1823 * delta +
                                    1890 * (-2 + kappa2) * pow(S2z, 2))) +
                        30 * (689 * (1 + delta) * S1z -
                              1134 * (1 + delta) * lambda1 * pow(S1z, 3) +
                              (-1 + delta) * S2z *
                                  (-689 + 1134 * lambda2 * pow(S2z, 2))) +
                        2 * Nu *
                            (Complex(0, 415432) - 66840 * S2z +
                             15 * (8 * (-557 + 1342 * delta) * S1z +
                                   189 *
                                       (5 * (1 + delta) * kappa1 +
                                        6 * (3 + delta) * lambda1) *
                                       pow(S1z, 3) -
                                   10736 * delta * S2z -
                                   2457 * (1 + delta) * kappa1 * pow(S1z, 2) *
                                       S2z +
                                   2457 * (-1 + delta) * kappa2 * S1z *
                                       pow(S2z, 2) -
                                   189 *
                                       (5 * (-1 + delta) * kappa2 +
                                        6 * (-3 + delta) * lambda2) *
                                       pow(S2z, 3))))))) /
            (317520. *
             pow(r, 4))) + /* Henry et al. QC spinning hereditary terms */
        (2 * M_PI *
         (kappa1 * (1 + delta - 2 * Nu) * pow(S1z, 2) +
          S2z * (4 * Nu * S1z - kappa2 * (-1 + delta + 2 * Nu) * S2z)) *
         pow(x, 4.5)));
  }

  else {
    return 0;
  }
}
// hQC_l_m() functions contains only the non-spinning hereditary terms at
// particular PN order.

static COMPLEX16 hQC_2_m_2(REAL8 Nu, UINT4 vpnorder, REAL8 x) {
  double EulerGamma = 0.5772156649015329;
  // keeping only the hereditary terms
  // if(vpnorder == 0){
  //   return(2*x);
  // }

  /* else if(vpnorder == 2){
     return((-5.095238095238095 + (55*Nu)/21.)*pow(x,2));
   } */

  /* else */ if (vpnorder == 3) {
    return (4 * M_PI * pow(x, 2.5));
  }

  /* else if(vpnorder == 4){
      return((-2.874338624338624 - (1069*Nu)/108. + (2047*pow(Nu,2))/756.)*
      pow(x,3));
   } */

  else if (vpnorder == 5) {
    // return((Complex(0,-48)*Nu - (214*M_PI)/21. +
    // (68*Nu*M_PI)/21.)*pow(x,3.5));
    return ((2 * (-107 + 34 * Nu) * M_PI * pow(x, 3.5)) / 21.);
  }

  else if (vpnorder == 6) {
    return (
        (pow(x, 4) * (-27392 * EulerGamma +
                      M_PI * (Complex(0, 13696) + 35 * (64 + 41 * Nu) * M_PI) -
                      13696 * log(16 * x))) /
        1680.);
  }

  else if (vpnorder == 7) {
    return (((-2173 - 4990 * Nu + 1120 * pow(Nu, 2)) * M_PI * pow(x, 4.5)) /
            378.);
  }

  // 4PN non-spinning quasi-circular (2,2) mode has been obtained from Blanchet
  // et al. arXiv:2304.11185

  else if (vpnorder == 8) {
    return (
        (pow(x, 5) *
         (276756480 * EulerGamma * (11449 + 19105 * Nu) -
          12 * (846557506853 + 1008017482431 * Nu) +
          35 * (28 * pow(Nu, 2) *
                    (5385456111 + 5 * Nu * (-163158374 + 26251249 * Nu)) -
                Complex(0, 3953664) * (11449 + 109657 * Nu) * M_PI -
                135135 * (54784 + 5 * Nu * (1951 + 6560 * Nu)) * pow(M_PI, 2)) +
          138378240 * (11449 + 19105 * Nu) * log(16 * x))) /
        7.62810048e10);
  }

  else {
    return 0;
  }
}

static COMPLEX16 hl_2_m_2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z, REAL8 x) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_2_m_2: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_2_m_2(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x) +
            hQC_2_m_2(Nu, vpnorder, x)) *
           cpolar(1, -2 * Phi);
  }
}

static COMPLEX16 hl_2_m_min2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z, REAL8 x) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_2_m_min2: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) * pow(-1, 2) *
           conj((hGO_2_m_2(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x) +
                 hQC_2_m_2(Nu, vpnorder, x))) *
           cpolar(1, 2 * Phi);
  }
}

// H21

static COMPLEX16 hGO_2_m_1(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z,
                           REAL8 x) {
  REAL8 delta = sqrt(1 - 4 * Nu);
  REAL8 kappa1 = 1.0;
  REAL8 kappa2 = 1.0;
  REAL8 r0 = 1.0;

  if (vpnorder == 1) {
    return (Complex(0, 0.6666666666666666) * delta * mass * PhiDOT);
  }

  else if (vpnorder == 2) {
    return ((Complex(0, -0.5) * pow(mass, 2) *
             ((1 + delta) * S1z + (-1 + delta) * S2z)) /
                pow(r, 2) -
            (Complex(0, -0.5) * (S1z - S2z + delta * (S1z + S2z)) * pow(x, 2)) +
            (Complex(0, -0.5) * (S1z - S2z + delta * (S1z + S2z)) * pow(x, 2)));
  }

  else if (vpnorder == 3) {
    return ((Complex(0, 0.023809523809523808) * delta * mass * PhiDOT *
             (4 * mass * (-9 + 11 * Nu) +
              r * ((19 - 24 * Nu) * pow(PhiDOT, 2) * pow(r, 2) +
                   Complex(0, 2) * (83 + 2 * Nu) * PhiDOT * r * rDOT +
                   2 * (-33 + 10 * Nu) * pow(rDOT, 2)))) /
            r);
  }

  else if (vpnorder == 4) {
    return (
        (Complex(0, 0.011904761904761904) * pow(mass, 2) *
         (2 * mass *
              ((77 + 59 * Nu + 11 * delta * (7 + Nu)) * S1z +
               (-77 - 59 * Nu + 11 * delta * (7 + Nu)) * S2z) +
          r * (Complex(0, -2) * PhiDOT * r * rDOT *
                   (147 * (1 + delta) * S1z + (-83 + 13 * delta) * Nu * S1z +
                    147 * (-1 + delta) * S2z + (83 + 13 * delta) * Nu * S2z) +
               pow(rDOT, 2) *
                   ((105 * (1 + delta) - 4 * (13 + 15 * delta) * Nu) * S1z +
                    (-105 + 15 * delta * (7 - 4 * Nu) + 52 * Nu) * S2z) +
               4 * pow(PhiDOT, 2) * pow(r, 2) *
                   ((-21 - 21 * delta + 66 * Nu + 4 * delta * Nu) * S1z +
                    (21 - 21 * delta - 66 * Nu + 4 * delta * Nu) * S2z)))) /
            pow(r, 3) -
        (Complex(0, 0.023809523809523808) *
         ((-7 + 205 * Nu + delta * (-7 + 33 * Nu)) * S1z +
          (7 - 205 * Nu + delta * (-7 + 33 * Nu)) * S2z) *
         pow(x, 3)) +
        (Complex(0, 0.023809523809523808) * pow(x, 3) *
         ((-7 + 205 * Nu + delta * (-7 + 33 * Nu)) * S1z +
          (7 - 205 * Nu + delta * (-7 + 33 * Nu)) * S2z)));
  }

  else if (vpnorder == 5) {
    return (
        (Complex(0, 0.0013227513227513227) * delta * mass * PhiDOT *
         (10 * pow(mass, 2) * (31 - 205 * Nu + 111 * pow(Nu, 2)) -
          2 * mass * r *
              ((-197 + 5 * Nu + 660 * pow(Nu, 2)) * pow(PhiDOT, 2) * pow(r, 2) +
               Complex(0, 1) * (-3167 - 5278 * Nu + 201 * pow(Nu, 2)) * PhiDOT *
                   r * rDOT +
               8 * (202 + 587 * Nu - 177 * pow(Nu, 2)) * pow(rDOT, 2)) +
          3 * pow(r, 2) *
              ((152 - 692 * Nu + 333 * pow(Nu, 2)) * pow(PhiDOT, 4) *
                   pow(r, 4) +
               Complex(0, 2) * (308 - 1607 * Nu + 111 * pow(Nu, 2)) *
                   pow(PhiDOT, 3) * pow(r, 3) * rDOT -
               3 * (75 - 560 * Nu + 68 * pow(Nu, 2)) * pow(PhiDOT, 2) *
                   pow(r, 2) * pow(rDOT, 2) -
               Complex(0, 2) * (-265 + 526 * Nu + 18 * pow(Nu, 2)) * PhiDOT *
                   r * pow(rDOT, 3) +
               (-241 + 550 * Nu - 264 * pow(Nu, 2)) * pow(rDOT, 4)))) /
            pow(r, 2)
        /* Henry et al QC spin terms */ /* +
  (Complex(0,-0.08333333333333333)*pow(x,3.5)*(2*(-((1 + delta)*(5 + kappa1))
  + 2*(6 + delta + (4 + 3*delta)*kappa1)*Nu)*pow(S1z,2) + S1z*(Complex(0,-3) -
  Complex(0,3)*delta + 6*(1 + delta)*M_PI
  - 16*delta*Nu*S2z - Complex(0,3)*(1 + delta)*log(16)) + S2z*(6*(-1 +
  delta)*M_PI - 2*(-1 + delta)*(5 + kappa2)*S2z
  + 4*(-6 + delta - 4*kappa2 + 3*delta*kappa2)*Nu*S2z - Complex(0,3)*(-1 +
  delta)*(1 + log(16))))) */
        +
        /* Henry et al. ecc spin terms */ (
            (pow(mass, 3) *
             (-4 * rDOT *
                  (kappa1 * (1 + delta - 2 * Nu) * pow(S1z, 2) +
                   kappa2 * (-1 + delta + 2 * Nu) * pow(S2z, 2)) -
              Complex(0, 1) * PhiDOT * r *
                  ((-((1 + delta) * (9 + kappa1)) +
                    2 * (9 + (4 + 3 * delta) * kappa1) * Nu) *
                       pow(S1z, 2) -
                   12 * delta * Nu * S1z * S2z +
                   (9 + kappa2 - 2 * (9 + 4 * kappa2) * Nu +
                    delta * (-9 - kappa2 + 6 * kappa2 * Nu)) *
                       pow(S2z, 2)))) /
            (6. * pow(r, 3))) + /* Henry et al. QC spinning hereditary terms */
        (Complex(0, -0.5) * ((1 + delta) * S1z + (-1 + delta) * S2z) *
         pow(x, 3.5) * (M_PI - Complex(0, 2) * log(2))));
  }

  else if (vpnorder == 6) {
    return (
        (delta * pow(mass, 2) * Nu * PhiDOT *
         (mass * (195 * PhiDOT * r - Complex(0, 946) * rDOT) +
          9 * r *
              (270 * pow(PhiDOT, 3) * pow(r, 3) -
               Complex(0, 483) * pow(PhiDOT, 2) * pow(r, 2) * rDOT -
               580 * PhiDOT * r * pow(rDOT, 2) +
               Complex(0, 42) * pow(rDOT, 3)))) /
            (315. * pow(r, 2))
        /* Henry et al. QC spin terms */
        /* +(Complex(0,-0.0006613756613756613)*(-4*pow(Nu,2)*(-1365*(S1z - S2z)
   + 179*delta*(S1z + S2z)) + 6*(208*(1 + delta)*S1z + 63*(1 +
delta)*kappa1*pow(S1z,3) + (-1 + delta)*S2z*(208 + 63*kappa2*pow(S2z,2))) -
Nu*(378*(3 + delta)*kappa1*pow(S1z,3) + 378*(1 + delta)*(-2 +
kappa1)*pow(S1z,2)*S2z + S1z*(8351 + 7027*delta + 378*(-1 + delta)*(-2 +
kappa2)*pow(S2z,2)) + S2z*(-8351 + 7027*delta + 378*(-3 +
delta)*kappa2*pow(S2z,2))))*pow(x,4)) */
        +
        /* Henry et al. ecc spin terms */ (
            Complex(0, 0.00033068783068783067) * pow(mass, 2) *
            (3 * pow(r, 2) *
                 (Complex(0, 4) * PhiDOT * r * pow(rDOT, 3) *
                      ((-315 * (1 + delta) + 2 * (251 + 463 * delta) * Nu +
                        4 * (15 + delta) * pow(Nu, 2)) *
                           S1z +
                       (315 - 315 * delta - 502 * Nu + 926 * delta * Nu +
                        4 * (-15 + delta) * pow(Nu, 2)) *
                           S2z) +
                  12 * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) *
                      ((189 * (1 + delta) - 2 * (521 + 293 * delta) * Nu +
                        7 * (55 + 23 * delta) * pow(Nu, 2)) *
                           S1z +
                       (189 * (-1 + delta) + 2 * (521 - 293 * delta) * Nu +
                        7 * (-55 + 23 * delta) * pow(Nu, 2)) *
                           S2z) +
                  pow(rDOT, 4) *
                      ((567 * (1 + delta) - 16 * (77 + 64 * delta) * Nu +
                        8 * (177 + 173 * delta) * pow(Nu, 2)) *
                           S1z +
                       (567 * (-1 + delta) + 16 * (77 - 64 * delta) * Nu +
                        8 * (-177 + 173 * delta) * pow(Nu, 2)) *
                           S2z) -
                  Complex(0, 4) * pow(PhiDOT, 3) * pow(r, 3) * rDOT *
                      ((936 * (1 + delta) - 5 * (979 + 215 * delta) * Nu +
                        2 * (1353 + 293 * delta) * pow(Nu, 2)) *
                           S1z +
                       (936 * (-1 + delta) + 5 * (979 - 215 * delta) * Nu +
                        2 * (-1353 + 293 * delta) * pow(Nu, 2)) *
                           S2z) +
                  4 * pow(PhiDOT, 4) * pow(r, 4) *
                      ((-252 * (1 + delta) + (1315 + 857 * delta) * Nu +
                        4 * (-285 + 43 * delta) * pow(Nu, 2)) *
                           S1z +
                       (252 + 5 * Nu * (-263 + 228 * Nu) +
                        delta * (-252 + Nu * (857 + 172 * Nu))) *
                           S2z)) -
             2 * mass * r *
                 (Complex(0, -1) * PhiDOT * r * rDOT *
                      ((2043 * (1 + delta) + (37 + 2597 * delta) * Nu +
                        (10635 + 139 * delta) * pow(Nu, 2)) *
                           S1z +
                       (2043 * (-1 + delta) + (-37 + 2597 * delta) * Nu +
                        (-10635 + 139 * delta) * pow(Nu, 2)) *
                           S2z) +
                  pow(PhiDOT, 2) * pow(r, 2) *
                      ((-765 - Nu * (667 + 7773 * Nu) +
                        delta * (-765 + 7 * Nu * (-533 + 245 * Nu))) *
                           S1z +
                       (765 + Nu * (667 + 7773 * Nu) +
                        delta * (-765 + 7 * Nu * (-533 + 245 * Nu))) *
                           S2z) +
                  4 * pow(rDOT, 2) *
                      ((-234 * (1 + delta) - 4 * (560 + 901 * delta) * Nu +
                        (483 + 1111 * delta) * pow(Nu, 2)) *
                           S1z +
                       (234 + 7 * (320 - 69 * Nu) * Nu +
                        delta * (-234 + Nu * (-3604 + 1111 * Nu))) *
                           S2z)) +
             2 * pow(mass, 2) *
                 (1134 * kappa1 * (-1 - delta + (3 + delta) * Nu) *
                      pow(S1z, 3) +
                  1134 * (1 + delta) * (-2 + kappa1) * Nu * pow(S1z, 2) * S2z +
                  S1z *
                      (-5661 - 5661 * delta - 17156 * Nu - 9172 * delta * Nu +
                       231 * pow(Nu, 2) + 775 * delta * pow(Nu, 2) +
                       1134 * (-1 + delta) * (-2 + kappa2) * Nu * pow(S2z, 2)) +
                  S2z * (5661 - 5661 * delta + 17156 * Nu - 9172 * delta * Nu -
                         231 * pow(Nu, 2) + 775 * delta * pow(Nu, 2) +
                         1134 * kappa2 * (1 - delta + (-3 + delta) * Nu) *
                             pow(S2z, 2))))) /
            pow(r, 4));
  }

  else if (vpnorder == 7) {

    return (
        /* Henry et al. QC spin terms */ /* Complex(0,0.0003968253968253968)*pow(x,4.5)*(10*(6*(1
  + delta)*(63 + 26*kappa1) - (84*(26 + 17*delta) + (1069 +
  757*delta)*kappa1)*Nu + (2844 + 312*delta + (1409 +
  36*delta)*kappa1)*pow(Nu,2))*pow(S1z,2) + S2z*(30*(14 - 431*Nu + delta*(-14 +
  87*Nu))*M_PI + 10*(6*(-1 + delta)*(63 + 26*kappa2) + (2184 - 1428*delta +
  1069*kappa2 - 757*delta*kappa2)*Nu +
     (-2844 + 312*delta - 1409*kappa2 + 36*delta*kappa2)*pow(Nu,2))*S2z +
  Complex(0,3)*(66 - 66*delta - 2957*Nu + 4405*delta*Nu - 20*(14 - 14*delta -
  431*Nu + 87*delta*Nu)*log(2))) + 3*S1z*(Complex(0,-1)*(66 + 66*delta - 2957*Nu
  - 280*log(2)) + 5*
   (-28*(1 + delta)*M_PI + 2*(431 + 87*delta)*Nu*M_PI + delta*Nu*(Complex(0,881)
  + 16*(7 + 23*Nu)*S2z) - Complex(0,4)*(431*Nu + delta*(-14 + 87*Nu))*log(2))))
   + */
        /* Henry et al. ecc+spin terms */ (
            (Complex(0, -1.0020843354176688e-7) *
             (-23760 * pow(mass, 3) *
                  (pow(PhiDOT, 3) * pow(r, 4) *
                       (Complex(0, -12) * (-35 + 107 * Nu) * S1z +
                        5 *
                            (126 - 462 * Nu + 80 * pow(Nu, 2) +
                             kappa1 * (-2 - 79 * Nu + 153 * pow(Nu, 2))) *
                            pow(S1z, 2) +
                        S2z * (Complex(0, -420) + Complex(0, 1284) * Nu +
                               10 * (-63 + kappa2) * S2z +
                               5 * (462 + 79 * kappa2) * Nu * S2z -
                               5 * (80 + 153 * kappa2) * pow(Nu, 2) * S2z)) -
                   Complex(0, 1) * pow(PhiDOT, 2) * pow(r, 3) * rDOT *
                       (Complex(0, -24) * (-35 + 202 * Nu) * S1z +
                        5 *
                            (-14 * Nu * (3 + 58 * Nu) +
                             kappa1 * (-409 + 1096 * Nu + 6 * pow(Nu, 2))) *
                            pow(S1z, 2) +
                        S2z * (Complex(0, -840) + 2045 * kappa2 * S2z +
                               10 * (406 - 3 * kappa2) * pow(Nu, 2) * S2z +
                               Nu * (Complex(0, 4848) + 210 * S2z -
                                     5480 * kappa2 * S2z))) -
                   Complex(0, 4) * r * pow(rDOT, 3) *
                       (Complex(0, 3) * (26 + 57 * Nu) * S1z +
                        5 *
                            (4 * pow(Nu, 2) +
                             kappa1 * (-7 + 49 * Nu + 15 * pow(Nu, 2))) *
                            pow(S1z, 2) -
                        S2z * (Complex(0, 78) - 35 * kappa2 * S2z +
                               5 * (4 + 15 * kappa2) * pow(Nu, 2) * S2z +
                               Nu * (Complex(0, 171) + 245 * kappa2 * S2z))) +
                   Complex(0, 2) * mass * rDOT *
                       (Complex(0, -4) * (-78 + 769 * Nu) * S1z +
                        5 *
                            ((14 - 59 * Nu) * Nu +
                             kappa1 * (-70 - 77 * Nu + 18 * pow(Nu, 2))) *
                            pow(S1z, 2) +
                        S2z * (Complex(0, -312) + 350 * kappa2 * S2z +
                               5 * (59 - 18 * kappa2) * pow(Nu, 2) * S2z +
                               Nu * (Complex(0, 3076) - 70 * S2z +
                                     385 * kappa2 * S2z))) +
                   PhiDOT * pow(r, 2) * pow(rDOT, 2) *
                       (Complex(0, 6) * (-62 + 219 * Nu) * S1z -
                        5 *
                            (189 - 756 * Nu + 388 * pow(Nu, 2) +
                             kappa1 * (98 - 266 * Nu + 276 * pow(Nu, 2))) *
                            pow(S1z, 2) +
                        S2z *
                            (Complex(0, 372) + 945 * S2z + 490 * kappa2 * S2z +
                             20 * (97 + 69 * kappa2) * pow(Nu, 2) * S2z -
                             2 * Nu *
                                 (Complex(0, 657) +
                                  35 * (54 + 19 * kappa2) * S2z))) +
                   mass * PhiDOT * r *
                       (Complex(0, 8) * (-61 + 480 * Nu) * S1z +
                        5 *
                            (-392 + 448 * Nu + 474 * pow(Nu, 2) +
                             kappa1 * (-11 + 150 * Nu + 58 * pow(Nu, 2))) *
                            pow(S1z, 2) -
                        S2z *
                            (Complex(0, -488) - 5 * (392 + 11 * kappa2) * S2z +
                             10 * (237 + 29 * kappa2) * pow(Nu, 2) * S2z +
                             10 * Nu *
                                 (Complex(0, 384) +
                                  (224 + 75 * kappa2) * S2z)))) +
              delta * mass *
                  (-240 * mass * PhiDOT * pow(r, 3) *
                       ((197936 - 139360 * Nu - 367105 * pow(Nu, 2) +
                         253245 * pow(Nu, 3)) *
                            pow(PhiDOT, 4) * pow(r, 4) +
                        Complex(0, 1) *
                            (279236 - 483940 * Nu - 2817805 * pow(Nu, 2) +
                             459180 * pow(Nu, 3)) *
                            pow(PhiDOT, 3) * pow(r, 3) * rDOT -
                        6 *
                            (38627 + 89295 * Nu - 492740 * pow(Nu, 2) +
                             75975 * pow(Nu, 3)) *
                            pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) -
                        Complex(0, 1) *
                            (-731008 + 2287930 * Nu + 981060 * pow(Nu, 2) +
                             10275 * pow(Nu, 3)) *
                            PhiDOT * r * pow(rDOT, 3) +
                        (-327667 + 436705 * Nu + 659790 * pow(Nu, 2) -
                         438255 * pow(Nu, 3)) *
                            pow(rDOT, 4)) +
                   900 * PhiDOT * pow(r, 4) *
                       (2 *
                            (-2594 + 27609 * Nu - 74032 * pow(Nu, 2) +
                             25974 * pow(Nu, 3)) *
                            pow(PhiDOT, 6) * pow(r, 6) +
                        Complex(0, 4) *
                            (-5730 + 58833 * Nu - 137842 * pow(Nu, 2) +
                             17123 * pow(Nu, 3)) *
                            pow(PhiDOT, 5) * pow(r, 5) * rDOT +
                        2 *
                            (-114 - 41622 * Nu + 147569 * pow(Nu, 2) +
                             4196 * pow(Nu, 3)) *
                            pow(PhiDOT, 4) * pow(r, 4) * pow(rDOT, 2) +
                        Complex(0, 4) *
                            (-9554 + 70788 * Nu - 156227 * pow(Nu, 2) +
                             5810 * pow(Nu, 3)) *
                            pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 3) +
                        (17619 - 138450 * Nu + 322600 * pow(Nu, 2) -
                         80816 * pow(Nu, 3)) *
                            pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 4) -
                        Complex(0, 2) *
                            (8793 - 52230 * Nu + 69340 * pow(Nu, 2) +
                             2536 * pow(Nu, 3)) *
                            PhiDOT * r * pow(rDOT, 5) +
                        2 *
                            (3957 - 24534 * Nu + 42584 * pow(Nu, 2) -
                             20800 * pow(Nu, 3)) *
                            pow(rDOT, 6)) -
                   2 * pow(mass, 3) *
                       (Complex(0, -23760) * rDOT *
                            (5 *
                                 (Nu * (-14 + 31 * Nu) +
                                  7 * kappa1 * (10 + 31 * Nu)) *
                                 pow(S1z, 2) +
                             2 * S1z *
                                 (Complex(0, -156) + 155 * pow(Nu, 2) * S2z +
                                  2 * Nu * (Complex(0, 613) + 390 * S2z)) +
                             S2z * (Complex(0, -312) + 350 * kappa2 * S2z +
                                    155 * pow(Nu, 2) * S2z +
                                    Nu * (Complex(0, 2452) +
                                          35 * (-2 + 31 * kappa2) * S2z))) +
                        PhiDOT * r *
                            (8946400 * pow(Nu, 3) -
                             8 * (6991786 + Complex(0, 724680) * S1z +
                                  7425 * (392 + 11 * kappa1) * pow(S1z, 2) +
                                  Complex(0, 724680) * S2z +
                                  7425 * (392 + 11 * kappa2) * pow(S2z, 2)) -
                             3600 * pow(Nu, 2) *
                                 (-628 +
                                  33 * (-19 + 92 * kappa1) * pow(S1z, 2) -
                                  7326 * S1z * S2z +
                                  33 * (-19 + 92 * kappa2) * pow(S2z, 2)) +
                             15 * Nu *
                                 (994455 * pow(M_PI, 2) +
                                  8 * (-2249485 +
                                       7920 * (-21 + 8 * kappa1) * pow(S1z, 2) +
                                       Complex(0, 283536) * S2z +
                                       7920 * (-21 + 8 * kappa2) * pow(S2z, 2) -
                                       1584 * S1z *
                                           (Complex(0, -179) + 170 * S2z))))) +
                   3 * pow(mass, 2) * r *
                       (Complex(0, 31680) * pow(rDOT, 3) *
                            (5 * (4 * pow(Nu, 2) + 7 * kappa1 * (-1 + 5 * Nu)) *
                                 pow(S1z, 2) +
                             S1z * (Complex(0, 78) + 40 * pow(Nu, 2) * S2z +
                                    Nu * (Complex(0, 327) + 420 * S2z)) +
                             S2z * (Complex(0, 78) - 35 * kappa2 * S2z +
                                    20 * pow(Nu, 2) * S2z +
                                    Nu * (Complex(0, 327) +
                                          175 * kappa2 * S2z))) -
                        22 * PhiDOT * r * pow(rDOT, 2) *
                            (2553200 * pow(Nu, 3) -
                             24 * (268267 + Complex(0, 5580) * S1z +
                                   525 * (27 + 14 * kappa1) * pow(S1z, 2) +
                                   Complex(0, 5580) * S2z +
                                   525 * (27 + 14 * kappa2) * pow(S2z, 2)) -
                             200 * pow(Nu, 2) *
                                 (39445 +
                                  72 * (-4 + 21 * kappa1) * pow(S1z, 2) -
                                  3600 * S1z * S2z +
                                  72 * (-4 + 21 * kappa2) * pow(S2z, 2)) +
                             25 * Nu *
                                 (23247 * pow(M_PI, 2) +
                                  8 * (-69259 + Complex(0, 1026) * S1z +
                                       126 * (27 + 5 * kappa1) * pow(S1z, 2) +
                                       Complex(0, 1026) * S2z +
                                       126 * (27 + 5 * kappa2) *
                                           pow(S2z, 2)))) +
                        pow(PhiDOT, 3) * pow(r, 3) *
                            (10071200 * pow(Nu, 3) +
                             96 * (-421183 - Complex(0, 34650) * S1z +
                                   825 * (-63 + kappa1) * pow(S1z, 2) -
                                   Complex(0, 34650) * S2z +
                                   825 * (-63 + kappa2) * pow(S2z, 2)) -
                             400 * pow(Nu, 2) *
                                 (64177 +
                                  792 * (-5 + 6 * kappa1) * pow(S1z, 2) -
                                  17424 * S1z * S2z +
                                  792 * (-5 + 6 * kappa2) * pow(S2z, 2)) +
                             15 * Nu *
                                 (426195 * pow(M_PI, 2) +
                                  8 * (-509635 +
                                       330 * (210 + 83 * kappa1) * pow(S1z, 2) +
                                       Complex(0, 29304) * S2z +
                                       330 * (210 + 83 * kappa2) * pow(S2z, 2) -
                                       792 * S1z *
                                           (Complex(0, -37) + 70 * S2z)))) -
                        Complex(0, 2) * pow(PhiDOT, 2) * pow(r, 2) * rDOT *
                            (-8330400 * pow(Nu, 3) +
                             8 * (-2810116 - Complex(0, 415800) * S1z +
                                  1012275 * kappa1 * pow(S1z, 2) -
                                  Complex(0, 415800) * S2z +
                                  1012275 * kappa2 * pow(S2z, 2)) +
                             4800 * pow(Nu, 2) *
                                 (13411 +
                                  33 * (19 + 12 * kappa1) * pow(S1z, 2) +
                                  462 * S1z * S2z +
                                  33 * (19 + 12 * kappa2) * pow(S2z, 2)) +
                             5 * Nu *
                                 (1278585 * pow(M_PI, 2) -
                                  8 * (5139685 +
                                       990 * (-21 + 139 * kappa1) *
                                           pow(S1z, 2) -
                                       Complex(0, 313632) * S2z +
                                       990 * (-21 + 139 * kappa2) *
                                           pow(S2z, 2) -
                                       3564 * S1z *
                                           (Complex(0, 88) + 185 * S2z)))))) -
              13559040 * delta * pow(mass, 3) * PhiDOT * r *
                  (2 * mass - 3 * pow(PhiDOT, 2) * pow(r, 3) +
                   Complex(0, 6) * PhiDOT * pow(r, 2) * rDOT +
                   6 * r * pow(rDOT, 2)) *
                  log(r / r0))) /
            pow(r, 4)) + /* Henry et al. QC spinning hereditary terms */
        ((-14 * (1 + delta) * S1z + (431 + 87 * delta) * Nu * S1z +
          (14 - 431 * Nu + delta * (-14 + 87 * Nu)) * S2z) *
         pow(x, 4.5) * (Complex(0, 1) * M_PI + log(4))) /
            84.);
  }

  else {
    return 0;
  }
}

static COMPLEX16 hQC_2_m_1(REAL8 Nu, UINT4 vpnorder, REAL8 x) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  /* if(vpnorder == 1){
       return(Complex(0,0.6666666666666666)*delta*pow(x,1.5));
   } */

  /* else if(vpnorder == 3){
       return((Complex(0,-0.40476190476190477)*delta +
    Complex(0,0.47619047619047616)*delta*Nu)*pow(x,2.5));
   } */

  /* else */ if (vpnorder == 4) {
    return ((2 * delta * pow(x, 3) * (Complex(0, 1) * M_PI + log(4))) / 3.);
  }

  /* else if(vpnorder == 5){
       return((Complex(0,-0.2275132275132275)*delta -
    Complex(0,2.693121693121693)*delta*Nu +
    Complex(0,0.3134920634920635)*delta*pow(Nu,2))*pow(x,3.5));
   } */

  else if (vpnorder == 6) {
    return (
        M_PI * (Complex(0, -0.40476190476190477) * delta * pow(x, 4) +
                Complex(0, 0.14285714285714285) * delta * Nu * pow(x, 4)) +
        ((-17 * delta * pow(x, 4)) / 21. + (2 * delta * Nu * pow(x, 4)) / 7.) *
            log(2));
  }

  else {

    return 0;
  }
}

static COMPLEX16 hl_2_m_1(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z, REAL8 x) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_2_m_1: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_2_m_1(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x) +
            hQC_2_m_1(Nu, vpnorder, x)) *
           cpolar(1, -1 * Phi);
  }
}

static COMPLEX16 hl_2_m_min1(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z, REAL8 x) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_2_m_min1: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) * pow(-1, 2) *
           conj((hGO_2_m_1(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x) +
                 hQC_2_m_1(Nu, vpnorder, x))) *
           cpolar(1, 1 * Phi);
  }
}

// H33
static COMPLEX16 hGO_3_m_3(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z,
                           REAL8 x) {
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
             (6 * (-5 + 19 * Nu) * pow(r, 2) *
                  pow(PhiDOT * r + Complex(0, 1) * rDOT, 4) *
                  (Complex(0, 1) * PhiDOT * r + rDOT) +
              2 * pow(mass, 2) *
                  (Complex(0, -3) * (-101 + 43 * Nu) * PhiDOT * r +
                   (-109 + 86 * Nu) * rDOT) +
              3 * mass * r *
                  (Complex(0, -12) * (1 + 4 * Nu) * pow(PhiDOT, 3) * pow(r, 3) +
                   6 * (14 + 31 * Nu) * pow(PhiDOT, 2) * pow(r, 2) * rDOT +
                   Complex(0, 3) * (33 + 62 * Nu) * PhiDOT * r * pow(rDOT, 2) -
                   4 * (8 + 17 * Nu) * pow(rDOT, 3)))) /
            (36. * pow(r, 2)));
  }

  else if (vpnorder == 4) {
    return (
        (Complex(0, -0.125) * sqrt(0.11904761904761904) * pow(mass, 2) *
         (4 * mass * (-1 + 5 * Nu) * ((1 + delta) * S1z + (-1 + delta) * S2z) +
          r * (2 * pow(rDOT, 2) *
                   (6 * (1 + delta) * S1z - 5 * (5 + 3 * delta) * Nu * S1z +
                    (-6 + delta * (6 - 15 * Nu) + 25 * Nu) * S2z) +
               pow(PhiDOT, 2) * pow(r, 2) *
                   (-24 * (1 + delta) * S1z + (119 + 33 * delta) * Nu * S1z +
                    (24 - 119 * Nu + 3 * delta * (-8 + 11 * Nu)) * S2z) +
               Complex(0, 2) * PhiDOT * r * rDOT *
                   (-18 * (1 + delta) * S1z + (77 + 39 * delta) * Nu * S1z +
                    (18 - 77 * Nu + 3 * delta * (-6 + 13 * Nu)) * S2z)))) /
            pow(r, 3) -
        (Complex(0, -0.375) * sqrt(1.0714285714285714) *
         (-4 * S1z + 19 * Nu * (S1z - S2z) + 4 * S2z - 4 * delta * (S1z + S2z) +
          5 * delta * Nu * (S1z + S2z)) *
         pow(x, 3)) +
        (Complex(0, -0.375) * sqrt(0.04285714285714286) * pow(x, 3) *
         (5 * (-4 + 19 * Nu + delta * (-4 + 5 * Nu)) * S1z +
          5 * (4 - 19 * Nu + delta * (-4 + 5 * Nu)) * S2z)));
  }

  else if (vpnorder == 5) {
    return (
        (delta *
         (30 * (183 - 1579 * Nu + 3387 * pow(Nu, 2)) * pow(r, 3) *
              pow(PhiDOT * r - Complex(0, 1) * rDOT, 2) *
              pow(Complex(0, -1) * PhiDOT * r + rDOT, 5) +
          10 * pow(mass, 3) *
              (Complex(0, -1) * (26473 - 27451 * Nu + 9921 * pow(Nu, 2)) *
                   PhiDOT * r +
               4 * (623 - 732 * Nu + 1913 * pow(Nu, 2)) * rDOT) +
          2 * pow(mass, 2) * r *
              (Complex(0, -11) * (-5353 - 13493 * Nu + 4671 * pow(Nu, 2)) *
                   pow(PhiDOT, 3) * pow(r, 3) +
               (-75243 - 142713 * Nu + 192821 * pow(Nu, 2)) * pow(PhiDOT, 2) *
                   pow(r, 2) * rDOT +
               Complex(0, 220) * (-256 + 781 * Nu + 840 * pow(Nu, 2)) * PhiDOT *
                   r * pow(rDOT, 2) -
               10 * (-756 + 8238 * Nu + 7357 * pow(Nu, 2)) * pow(rDOT, 3)) +
          3 * mass * pow(r, 2) *
              (Complex(0, 2) * (-7633 + 9137 * Nu + 28911 * pow(Nu, 2)) *
                   pow(PhiDOT, 5) * pow(r, 5) -
               4 * (-8149 + 1576 * Nu + 43533 * pow(Nu, 2)) * pow(PhiDOT, 4) *
                   pow(r, 4) * rDOT -
               Complex(0, 2) * (-9297 - 19517 * Nu + 64839 * pow(Nu, 2)) *
                   pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) -
               32 * (-1288 + 3667 * Nu + 4056 * pow(Nu, 2)) * pow(PhiDOT, 2) *
                   pow(r, 2) * pow(rDOT, 3) -
               Complex(0, 5) * (-9851 + 17954 * Nu + 40968 * pow(Nu, 2)) *
                   PhiDOT * r * pow(rDOT, 4) +
               20 * (-771 + 1126 * Nu + 3616 * pow(Nu, 2)) * pow(rDOT, 5)))) /
            (1584. * sqrt(210) * pow(r, 3))
        /* Henry et al. QC spin terms */
        /* +(Complex(0,1.125)*sqrt(1.0714285714285714)*(kappa1*(-1 - delta +
  2*(2 + delta)*Nu)*pow(S1z,2)
  + S2z*(-4*delta*Nu*S1z + kappa2*(1 - delta + 2*(-2 +
  delta)*Nu)*S2z))*pow(x,3.5)) */
        + /* Henry et al. ecc spin terms */ (
              (Complex(0, 0.125) * sqrt(1.0714285714285714) * pow(mass, 3) *
               (7 * PhiDOT * r + Complex(0, 2) * rDOT) *
               (kappa1 * (-1 - delta + 2 * (2 + delta) * Nu) * pow(S1z, 2) +
                S2z * (-4 * delta * Nu * S1z +
                       kappa2 * (1 - delta + 2 * (-2 + delta) * Nu) * S2z))) /
              pow(r, 3)));
  }

  else if (vpnorder == 6) {
    return (
        -(delta * pow(mass, 2) * Nu *
          (668 * pow(mass, 2) +
           2 * mass * r *
               (4081 * pow(PhiDOT, 2) * pow(r, 2) +
                Complex(0, 297) * PhiDOT * r * rDOT - 452 * pow(rDOT, 2)) +
           5 * pow(r, 2) *
               (1329 * pow(PhiDOT, 4) * pow(r, 4) -
                Complex(0, 2926) * pow(PhiDOT, 3) * pow(r, 3) * rDOT -
                384 * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) -
                Complex(0, 408) * PhiDOT * r * pow(rDOT, 3) +
                200 * pow(rDOT, 4)))) /
            (36. * sqrt(210) * pow(r, 4))
        /* Henry et al. QC spin terms */
        /* +(Complex(0,-0.125)*sqrt(0.04285714285714286)*(-(Nu*((279 +
  delta)*S1z + (-279 + delta)*S2z))
  + 10*(S1z - S2z + delta*(S1z + S2z)) + pow(Nu,2)*(407*(S1z - S2z) +
  241*delta*(S1z + S2z)))*pow(x,4)) */
        + /* Henry et al. ecc spin terms*/ (
              (Complex(0, -0.006944444444444444) * pow(mass, 2) *
               (10 * pow(mass, 2) *
                    ((252 * (1 + delta) - (1277 + 1279 * delta) * Nu +
                      8 * (12 + 47 * delta) * pow(Nu, 2)) *
                         S1z +
                     (252 * (-1 + delta) + (1277 - 1279 * delta) * Nu +
                      8 * (-12 + 47 * delta) * pow(Nu, 2)) *
                         S2z) +
                2 * mass * r *
                    (2 * pow(PhiDOT, 2) * pow(r, 2) *
                         ((1320 * (1 + delta) - 2 * (4469 + 211 * delta) * Nu +
                           (8709 + 2777 * delta) * pow(Nu, 2)) *
                              S1z +
                          (1320 * (-1 + delta) + 8938 * Nu - 422 * delta * Nu +
                           (-8709 + 2777 * delta) * pow(Nu, 2)) *
                              S2z) +
                     Complex(0, 3) * PhiDOT * r * rDOT *
                         ((2000 * (1 + delta) - (9147 + 3173 * delta) * Nu +
                           (8911 + 5273 * delta) * pow(Nu, 2)) *
                              S1z +
                          (2000 * (-1 + delta) + (9147 - 3173 * delta) * Nu +
                           (-8911 + 5273 * delta) * pow(Nu, 2)) *
                              S2z) +
                     10 * pow(rDOT, 2) *
                         ((-105 * (1 + delta) + (541 + 77 * delta) * Nu -
                           2 * (462 + 247 * delta) * pow(Nu, 2)) *
                              S1z +
                          (105 + Nu * (-541 + 924 * Nu) +
                           delta * (-105 + (77 - 494 * Nu) * Nu)) *
                              S2z)) -
                3 * pow(r, 2) *
                    (-3 * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) *
                         ((480 * (1 + delta) - (1711 + 1889 * delta) * Nu +
                           2 * (-1161 + 757 * delta) * pow(Nu, 2)) *
                              S1z +
                          (480 * (-1 + delta) + (1711 - 1889 * delta) * Nu +
                           2 * (1161 + 757 * delta) * pow(Nu, 2)) *
                              S2z) +
                     2 * pow(PhiDOT, 4) * pow(r, 4) *
                         ((350 * (1 + delta) - 4 * (404 + 461 * delta) * Nu +
                           (883 + 769 * delta) * pow(Nu, 2)) *
                              S1z +
                          (350 * (-1 + delta) + 4 * (404 - 461 * delta) * Nu +
                           (-883 + 769 * delta) * pow(Nu, 2)) *
                              S2z) +
                     Complex(0, 2) * pow(PhiDOT, 3) * pow(r, 3) * rDOT *
                         ((660 * (1 + delta) - (4061 + 2899 * delta) * Nu +
                           (2643 + 4789 * delta) * pow(Nu, 2)) *
                              S1z +
                          (660 * (-1 + delta) + (4061 - 2899 * delta) * Nu +
                           (-2643 + 4789 * delta) * pow(Nu, 2)) *
                              S2z) +
                     10 * pow(rDOT, 4) *
                         ((-30 * (1 + delta) + (187 + 101 * delta) * Nu -
                           2 * (159 + 61 * delta) * pow(Nu, 2)) *
                              S1z +
                          (30 + Nu * (-187 + 318 * Nu) +
                           delta * (-30 + (101 - 122 * Nu) * Nu)) *
                              S2z) +
                     Complex(0, 2) * PhiDOT * r * pow(rDOT, 3) *
                         ((90 + Nu * (-1321 + 5118 * Nu) +
                           delta * (90 + Nu * (-319 + 714 * Nu))) *
                              S1z +
                          (-90 + (1321 - 5118 * Nu) * Nu +
                           delta * (90 + Nu * (-319 + 714 * Nu))) *
                              S2z)))) /
              (sqrt(210) * pow(r, 4))));
  }

  else if (vpnorder == 7) {

    return (
        /* Henry et al. QC spin terms */ /* (Complex(0,-0.020833333333333332)*pow(x,4.5)*(-270*(-6*Nu*(-2
 + delta*(-2 + Nu) + 6*Nu)
  + kappa1*(4 - Nu*(13 + 8*Nu) + delta*(4 + Nu*(-5 + 12*Nu))))*pow(S1z,2) +
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
        ((Complex(0, 504504) * pow(mass, 3) *
              (2 * mass *
                   (rDOT *
                        (S1z * (108 - 498 * Nu +
                                Complex(0, 5) *
                                    (-24 - 5 * kappa1 + 3 * (76 + kappa1) * Nu +
                                     4 * (-111 + 13 * kappa1) * pow(Nu, 2)) *
                                    S1z) +
                         6 * (-18 + 83 * Nu) * S2z -
                         Complex(0, 5) *
                             (-24 - 5 * kappa2 + 3 * (76 + kappa2) * Nu +
                              4 * (-111 + 13 * kappa2) * pow(Nu, 2)) *
                             pow(S2z, 2)) +
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
                             pow(S2z, 2))) +
               r * (3 * PhiDOT * r * pow(rDOT, 2) *
                        (S1z * (Complex(0, 216) + 545 * kappa1 * S1z +
                                40 * (45 + 8 * kappa1) * pow(Nu, 2) * S1z -
                                30 * Nu *
                                    (Complex(0, 50) + 20 * S1z +
                                     73 * kappa1 * S1z)) +
                         Complex(0, 12) * (-18 + 125 * Nu) * S2z -
                         5 *
                             (109 * kappa2 - 6 * (20 + 73 * kappa2) * Nu +
                              8 * (45 + 8 * kappa2) * pow(Nu, 2)) *
                             pow(S2z, 2)) +
                    2 * pow(rDOT, 3) *
                        (S1z * (-54 + Complex(0, 145) * kappa1 * S1z -
                                Complex(0, 30) * Nu *
                                    (Complex(0, 11) + 2 * Nu * S1z +
                                     kappa1 * (17 + 8 * Nu) * S1z)) +
                         6 * (9 - 55 * Nu) * S2z +
                         Complex(0, 5) *
                             (12 * pow(Nu, 2) +
                              kappa2 * (-29 + 6 * Nu * (17 + 8 * Nu))) *
                             pow(S2z, 2)) +
                    6 * pow(PhiDOT, 2) * pow(r, 2) * rDOT *
                        (S1z * (297 - 465 * Nu +
                                Complex(0, 5) *
                                    (6 * (20 - 87 * Nu) * Nu +
                                     kappa1 * (-50 + 3 * Nu * (47 + 76 * Nu))) *
                                    S1z) +
                         3 * (-99 + 155 * Nu) * S2z -
                         Complex(0, 5) *
                             (6 * (20 - 87 * Nu) * Nu +
                              kappa2 * (-50 + 3 * Nu * (47 + 76 * Nu))) *
                             pow(S2z, 2)) +
                    pow(PhiDOT, 3) * pow(r, 3) *
                        (Complex(0, 3) * (531 - 1295 * Nu) * S1z +
                         10 *
                             (-33 * kappa1 + 6 * (30 + 13 * kappa1) * Nu +
                              4 * (-96 + 67 * kappa1) * pow(Nu, 2)) *
                             pow(S1z, 2) +
                         S2z * (Complex(0, -1593) + 330 * kappa2 * S2z +
                                5 * Nu *
                                    (Complex(0, 777) -
                                     4 *
                                         (90 + 39 * kappa2 - 192 * Nu +
                                          134 * kappa2 * Nu) *
                                         S2z))))) +
          delta *
              (-17640 * (-4083 + Nu * (58311 + Nu * (-269240 + 405617 * Nu))) *
                   pow(r, 4) * pow(PhiDOT * r + Complex(0, 1) * rDOT, 6) *
                   pow(Complex(0, 1) * PhiDOT * r + rDOT, 3) +
               168 * pow(mass, 2) * pow(r, 2) *
                   (Complex(0, 1) *
                        (-7508635 +
                         7 * Nu *
                             (-1318438 + Nu * (-10231834 + 9667755 * Nu))) *
                        pow(PhiDOT, 5) * pow(r, 5) +
                    7 *
                        (1235591 +
                         Nu * (884445 + (23935218 - 26913443 * Nu) * Nu)) *
                        pow(PhiDOT, 4) * pow(r, 4) * rDOT -
                    Complex(0, 1) *
                        (8961149 +
                         7 * Nu *
                             (-31755709 + Nu * (-11134798 + 22187331 * Nu))) *
                        pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) -
                    (-36806435 +
                     7 * Nu * (33178545 + Nu * (24565078 + 22873537 * Nu))) *
                        pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 3) -
                    Complex(0, 5) *
                        (-7761899 +
                         7 * Nu *
                             (2892563 + 5998602 * Nu + 7493619 * pow(Nu, 2))) *
                        PhiDOT * r * pow(rDOT, 4) +
                    5 * (-2422057 + 7 * Nu * (501045 + Nu * (2033141 + 2771816 * Nu))) *
                        pow(rDOT, 5)) +
               1764 * mass * pow(r, 3) *
                   (Complex(0, -2) *
                        (239087 +
                         Nu * (-1206515 + Nu * (422631 + 3979375 * Nu))) *
                        pow(PhiDOT, 7) * pow(r, 7) +
                    2 * (621284 + Nu * (-2279907 + 2 * Nu * (-1180187 + 5876531 * Nu))) *
                        pow(PhiDOT, 6) * pow(r, 6) * rDOT +
                    Complex(0, 2) *
                        (39270 +
                         Nu * (1235486 - 5319747 * Nu + 4406349 * pow(Nu, 2))) *
                        pow(PhiDOT, 5) * pow(r, 5) * pow(rDOT, 2) +
                    8 * (349111 + 4 * Nu * (-519370 + 33 * Nu * (10939 + 42635 * Nu))) *
                        pow(PhiDOT, 4) * pow(r, 4) * pow(rDOT, 3) +
                    Complex(0, 2) *
                        (1212607 +
                         3 * Nu *
                             (-2012698 - 67827 * Nu + 7955628 * pow(Nu, 2))) *
                        pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 4) +
                    4 * (201135 + 2 * Nu * (-773107 + Nu * (1214819 + 1157652 * Nu))) *
                        pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 5) +
                    Complex(0, 5) *
                        (333969 +
                         2 * Nu * (-981471 + 4 * Nu * (154039 + 750016 * Nu))) *
                        PhiDOT * r * pow(rDOT, 6) -
                    40 *
                        (13245 +
                         2 * Nu * (-37005 + Nu * (14251 + 130160 * Nu))) *
                        pow(rDOT, 7)) +
               2 * pow(mass, 4) *
                   (4 * rDOT *
                        (269279500 * pow(Nu, 3) +
                         2 * (-174108226 +
                              63063 * S1z *
                                  (Complex(0, 108) +
                                   5 * (24 + 5 * kappa1) * S1z) +
                              63063 * S2z *
                                  (Complex(0, 108) +
                                   5 * (24 + 5 * kappa2) * S2z)) -
                         21 * Nu *
                             (103100846 + 1846845 * pow(M_PI, 2) +
                              Complex(0, 1693692) * S2z -
                              6006 * (S1z * (Complex(0, -282) +
                                             5 * (-180 + 7 * kappa1) * S1z) -
                                      980 * S1z * S2z +
                                      5 * (-180 + 7 * kappa2) * pow(S2z, 2))) -
                         2940 * pow(Nu, 2) *
                             (-122855 +
                              4719 * ((-6 + 7 * kappa1) * pow(S1z, 2) -
                                      26 * S1z * S2z +
                                      (-6 + 7 * kappa2) * pow(S2z, 2)))) +
                    Complex(0, 1) * PhiDOT * r *
                        (-1176172480 * pow(Nu, 3) +
                         8 * (-74084729 +
                              189189 * S1z *
                                  (Complex(0, 99) +
                                   5 * (-8 + 133 * kappa1) * S1z) +
                              189189 * S2z *
                                  (Complex(0, 99) +
                                   5 * (-8 + 133 * kappa2) * S2z)) +
                         35280 * pow(Nu, 2) *
                             (56255 +
                              429 * ((-22 + 65 * kappa1) * pow(S1z, 2) -
                                     174 * S1z * S2z +
                                     (-22 + 65 * kappa2) * pow(S2z, 2))) -
                         147 * Nu *
                             (-65012788 + 4485195 * pow(M_PI, 2) +
                              Complex(0, 3943368) * S2z +
                              10296 *
                                  (S1z * (Complex(0, 383) +
                                          5 * (-96 + 277 * kappa1) * S1z) -
                                   3220 * S1z * S2z +
                                   5 * (-96 + 277 * kappa2) * pow(S2z, 2))))) +
               pow(mass, 3) * r *
                   (Complex(0, -12) * PhiDOT * r * pow(rDOT, 2) *
                        (-1035895280 * pow(Nu, 3) -
                         2 * (-547993687 +
                              63063 * S1z *
                                  (Complex(0, 216) + 545 * kappa1 * S1z) +
                              63063 * S2z *
                                  (Complex(0, 216) + 545 * kappa2 * S2z)) +
                         77 * Nu *
                             (42451610 + 1511055 * pow(M_PI, 2) +
                              Complex(0, 1749384) * S2z +
                              6552 * (S1z * (Complex(0, 267) +
                                             25 * (6 + 11 * kappa1) * S1z) -
                                      5 * S1z * S2z +
                                      25 * (6 + 11 * kappa2) * pow(S2z, 2))) +
                         490 * pow(Nu, 2) *
                             (-5802767 +
                              5148 * ((-6 + 23 * kappa1) * pow(S1z, 2) -
                                      58 * S1z * S2z +
                                      (-6 + 23 * kappa2) * pow(S2z, 2)))) +
                    4 * pow(rDOT, 3) *
                        (-1359334480 * pow(Nu, 3) -
                         4 * (-150254558 +
                              63063 * S1z *
                                  (Complex(0, 54) + 145 * kappa1 * S1z) +
                              63063 * S2z *
                                  (Complex(0, 54) + 145 * kappa2 * S2z)) +
                         231 * Nu *
                             (8490448 + 503685 * pow(M_PI, 2) +
                              Complex(0, 242424) * S2z +
                              2184 * (S1z * (Complex(0, 111) +
                                             110 * kappa1 * S1z) +
                                      70 * S1z * S2z +
                                      110 * kappa2 * pow(S2z, 2))) +
                         11760 * pow(Nu, 2) *
                             (-312980 +
                              429 * ((3 + 25 * kappa1) * pow(S1z, 2) -
                                     44 * S1z * S2z +
                                     (3 + 25 * kappa2) * pow(S2z, 2)))) +
                    6 * pow(PhiDOT, 2) * pow(r, 2) * rDOT *
                        (2368900688 * pow(Nu, 3) +
                         8 * (-812986529 +
                              63063 * S1z *
                                  (Complex(0, 297) + 250 * kappa1 * S1z) +
                              63063 * S2z *
                                  (Complex(0, 297) + 250 * kappa2 * S2z)) -
                         1176 * pow(Nu, 2) *
                             (2423171 +
                              4290 * ((-3 + 41 * kappa1) * pow(S1z, 2) -
                                      88 * S1z * S2z +
                                      (-3 + 41 * kappa2) * pow(S2z, 2))) +
                         539 * Nu *
                             (-24139772 + 647595 * pow(M_PI, 2) +
                              Complex(0, 120744) * S2z -
                              936 * (S1z * (Complex(0, -129) + 600 * S1z +
                                            205 * kappa1 * S1z) -
                                     460 * S1z * S2z +
                                     5 * (120 + 41 * kappa2) * pow(S2z, 2)))) +
                    Complex(0, 1) * pow(PhiDOT, 3) * pow(r, 3) *
                        (-4538040136 * pow(Nu, 3) -
                         88 * (259018351 +
                               17199 * S1z *
                                   (Complex(0, -531) + 110 * kappa1 * S1z) +
                               17199 * S2z *
                                   (Complex(0, -531) + 110 * kappa2 * S2z)) +
                         2352 * pow(Nu, 2) *
                             (7332973 +
                              12870 * ((5 + 23 * kappa1) * pow(S1z, 2) -
                                       36 * S1z * S2z +
                                       (5 + 23 * kappa2) * pow(S2z, 2))) +
                         21 * Nu *
                             (49864815 * pow(M_PI, 2) +
                              8 * (-88128538 - Complex(0, 2099097) * S2z +
                                   9009 *
                                       (S1z * (Complex(0, -233) +
                                               40 * (15 + kappa1) * S1z) +
                                        360 * S1z * S2z +
                                        40 * (15 + kappa2) * pow(S2z, 2))))))) +
          74954880 * delta * pow(mass, 3) *
              (Complex(0, 22) * mass * PhiDOT * r +
               Complex(0, 59) * pow(PhiDOT, 3) * pow(r, 4) + 8 * mass * rDOT +
               66 * pow(PhiDOT, 2) * pow(r, 3) * rDOT +
               Complex(0, 24) * PhiDOT * pow(r, 2) * pow(rDOT, 2) -
               4 * r * pow(rDOT, 3)) *
              log(r / r0)) /
         (2.4216192e7 * sqrt(210) *
          pow(r, 4))) + /* Henry et al. QC spinning hereditary terms */
        ((9 * sqrt(0.04285714285714286) * pow(x, 4.5) *
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

static COMPLEX16 hQC_3_m_3(REAL8 Nu, UINT4 vpnorder, REAL8 x) {
  REAL8 delta = sqrt(1 - 4 * Nu);
  double EulerGamma = 0.5772156649015329;

  /* if(vpnorder == 1){
       return(Complex(0,-1.5)*sqrt(1.0714285714285714)*delta*pow(x,1.5));
   } */

  /* else if(vpnorder == 3){
       return((Complex(0,3)*sqrt(4.285714285714286)*delta -
    Complex(0,3)*sqrt(1.0714285714285714)*delta*Nu)*pow(x,2.5));
   } */

  /* else */ if (vpnorder == 4) {
    return ((9 * sqrt(1.0714285714285714) * delta * pow(x, 3) *
             (Complex(0, -1) * M_PI + log(2.25))) /
            2.);
  }

  /* else if(vpnorder == 5){
       return((Complex(0,-8.386363636363637)*sqrt(0.04285714285714286)*delta +
    Complex(0,83.54545454545455)*sqrt(0.04285714285714286)*delta*Nu -
    Complex(0,20.15909090909091)*sqrt(0.04285714285714286)*delta*
     pow(Nu,2))*pow(x,3.5));
   } */

  else if (vpnorder == 6) {
    return (M_PI *
                (Complex(0, 9) * sqrt(4.285714285714286) * delta * pow(x, 4) -
                 Complex(0, 6.75) * sqrt(1.0714285714285714) * delta * Nu *
                     pow(x, 4)) +
            (-18 * sqrt(4.285714285714286) * delta * pow(x, 4) +
             (27 * sqrt(1.0714285714285714) * delta * Nu * pow(x, 4)) / 2.) *
                log(1.5));
  } else if (vpnorder == 7) {
    return (Complex(0, 1.1149564720993292e-6) * sqrt(0.04285714285714286) *
            delta * pow(x, 4.5) *
            (-465315528 + 74954880 * EulerGamma + 13827800 * Nu +
             124985672 * pow(Nu, 2) - 19373424 * pow(Nu, 3) +
             Complex(0, 47279232) * M_PI - 10090080 * pow(M_PI, 2) -
             4309305 * Nu * pow(M_PI, 2) - 94558464 * log(1.5) -
             Complex(0, 121080960) * M_PI * log(1.5) + 37477440 * log(16 * x) +
             121080960 * log(1.5) * log(1.5)));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hl_3_m_3(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z, REAL8 x) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_3_m_3: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_3_m_3(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x) +
            hQC_3_m_3(Nu, vpnorder, x)) *
           cpolar(1, -3 * Phi);
  }
}

static COMPLEX16 hl_3_m_min3(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z, REAL8 x) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_3_m_min3: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) * pow(-1, 3) *
           conj(hGO_3_m_3(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x) +
                hQC_3_m_3(Nu, vpnorder, x)) *
           cpolar(1, 3 * Phi);
  }
}

// H32
static COMPLEX16 hGO_3_m_2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z,
                           REAL8 x) {
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
        (sqrt(0.7142857142857143) * pow(mass, 2) * Nu *
         (4 * PhiDOT * r + Complex(0, 1) * rDOT) * (S1z + S2z)) /
            (3. * pow(r, 2)) -
        ((4 * sqrt(0.7142857142857143) * Nu * (S1z + S2z) * pow(x, 2.5)) / 3.) +
        ((4 * sqrt(0.7142857142857143) * Nu * (S1z + S2z) * pow(x, 2.5)) / 3.));
  }

  else if (vpnorder == 4) {
    return (-(mass * PhiDOT *
              (2 * mass *
                   ((167 - 925 * Nu + 1615 * pow(Nu, 2)) * PhiDOT * r +
                    Complex(0, 5) * (-82 + 239 * Nu + 55 * pow(Nu, 2)) * rDOT) -
               3 * r *
                   (2 * (-13 - 25 * Nu + 355 * pow(Nu, 2)) * pow(PhiDOT, 3) *
                        pow(r, 3) -
                    Complex(0, 60) * (-8 + 25 * Nu + pow(Nu, 2)) *
                        pow(PhiDOT, 2) * pow(r, 2) * rDOT +
                    12 * (-23 + 70 * Nu + 10 * pow(Nu, 2)) * PhiDOT * r *
                        pow(rDOT, 2) +
                    Complex(0, 5) * (-13 + 38 * Nu + 10 * pow(Nu, 2)) *
                        pow(rDOT, 3)))) /
            (108. * sqrt(35) * r));
  }

  else if (vpnorder == 5) {
    return (
        (pow(mass, 2) * Nu * PhiDOT *
         (Complex(0, 7) * mass +
          r * (Complex(0, 49) * pow(PhiDOT, 2) * pow(r, 2) +
               90 * PhiDOT * r * rDOT - Complex(0, 6) * pow(rDOT, 2)))) /
            (4. * sqrt(35) * pow(r, 2))
        /* Henry et al. QC spin terms */
        /* +((sqrt(0.7142857142857143)*(43*delta*Nu*(S1z - S2z) - 3*Nu*(S1z +
  S2z)
  - 26*pow(Nu,2)*(S1z + S2z) - 8*(S1z + delta*S1z + S2z -
  delta*S2z))*pow(x,3.5))/9.) */
        +
        /* Henry et al. ecc spin terms */ (
            (sqrt(0.7142857142857143) * pow(mass, 2) *
             (Complex(0, 2) * mass * rDOT *
                  ((-12 + Nu * (97 + 4 * Nu) + delta * (-12 + 5 * Nu)) * S1z +
                   (-12 + delta * (12 - 5 * Nu) + Nu * (97 + 4 * Nu)) * S2z) +
              4 * mass * PhiDOT * r *
                  (-((12 + delta * (12 - 23 * Nu) + Nu * (53 + 8 * Nu)) * S1z) -
                   (12 + Nu * (53 + 8 * Nu) + delta * (-12 + 23 * Nu)) * S2z) -
              3 * r *
                  (16 * pow(Nu, 2) * PhiDOT * r *
                       (4 * pow(PhiDOT, 2) * pow(r, 2) -
                        Complex(0, 2) * PhiDOT * r * rDOT + pow(rDOT, 2)) *
                       (S1z + S2z) +
                   Complex(0, 30) * pow(PhiDOT, 2) * pow(r, 2) * rDOT *
                       (S1z + delta * S1z + S2z - delta * S2z) +
                   Nu * (-4 * pow(PhiDOT, 3) * pow(r, 3) *
                             ((5 + 17 * delta) * S1z + (5 - 17 * delta) * S2z) -
                         Complex(0, 1) * pow(PhiDOT, 2) * pow(r, 2) * rDOT *
                             ((189 + 17 * delta) * S1z +
                              (189 - 17 * delta) * S2z) +
                         20 * PhiDOT * r * pow(rDOT, 2) *
                             (-((-3 + delta) * S1z) + (3 + delta) * S2z) -
                         Complex(0, 4) * pow(rDOT, 3) *
                             ((-4 + delta) * S1z - (4 + delta) * S2z))))) /
            (72. * pow(r, 3))));
  }

  else if (vpnorder == 6) {
    return (
        -(mass * PhiDOT *
          (4 * pow(mass, 2) *
               (2 *
                    (5377 + 6438 * Nu - 79866 * pow(Nu, 2) +
                     37348 * pow(Nu, 3)) *
                    PhiDOT * r -
                Complex(0, 5) *
                    (-4115 + 18399 * Nu - 20276 * pow(Nu, 2) + 7 * pow(Nu, 3)) *
                    rDOT) -
           4 * mass * r *
               ((4599 - 15737 * Nu + 36259 * pow(Nu, 2) + 108563 * pow(Nu, 3)) *
                    pow(PhiDOT, 3) * pow(r, 3) -
                Complex(0, 1) *
                    (-34053 + 59698 * Nu + 192949 * pow(Nu, 2) +
                     16193 * pow(Nu, 3)) *
                    pow(PhiDOT, 2) * pow(r, 2) * rDOT +
                (-59058 + 77983 * Nu + 322468 * pow(Nu, 2) -
                 4264 * pow(Nu, 3)) *
                    PhiDOT * r * pow(rDOT, 2) +
                Complex(0, 5) *
                    (-3387 + 8518 * Nu + 8968 * pow(Nu, 2) + 884 * pow(Nu, 3)) *
                    pow(rDOT, 3)) +
           3 * pow(r, 2) *
               (4 *
                    (-710 + 3892 * Nu - 10655 * pow(Nu, 2) +
                     24000 * pow(Nu, 3)) *
                    pow(PhiDOT, 5) * pow(r, 5) +
                Complex(0, 11) *
                    (-1484 + 11693 * Nu - 25006 * pow(Nu, 2) +
                     428 * pow(Nu, 3)) *
                    pow(PhiDOT, 4) * pow(r, 4) * rDOT +
                4 *
                    (4161 - 25618 * Nu + 29489 * pow(Nu, 2) +
                     22078 * pow(Nu, 3)) *
                    pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) +
                Complex(0, 44) *
                    (-151 + 1067 * Nu - 2419 * pow(Nu, 2) + 57 * pow(Nu, 3)) *
                    pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 3) +
                4 *
                    (2041 - 11680 * Nu + 19334 * pow(Nu, 2) +
                     3368 * pow(Nu, 3)) *
                    PhiDOT * r * pow(rDOT, 4) +
                Complex(0, 5) *
                    (477 - 2624 * Nu + 3862 * pow(Nu, 2) + 1160 * pow(Nu, 3)) *
                    pow(rDOT, 5)))) /
            (4752. * sqrt(35) * pow(r, 2))
        /* Henry et al. QC spin terms */
        /* +((2*sqrt(0.7142857142857143)*((2*Nu*(-5 - 5*delta + 4*Nu) +
  3*kappa1*(1 + delta
  - 2*(2 + delta)*Nu + 6*pow(Nu,2)))*pow(S1z,2) + S2z*(2*(4 +
  9*kappa2)*pow(Nu,2)*S2z
  - 3*(-1 + delta)*(Complex(0,2) + kappa2*S2z) + 2*Nu*(Complex(0,-15) + 6*M_PI
  + 5*(-1 + delta)*S2z + 3*(-2 + delta)*kappa2*S2z)) + 2*S1z*(Complex(0,3)
  + Complex(0,3)*delta + Nu*(Complex(0,-15) + 6*M_PI + 2*S2z -
  10*Nu*S2z)))*pow(x,4))/9.) */
        +
        /* Henry et al. ecc spin terms */ (
            (sqrt(0.7142857142857143) * pow(mass, 3) *
             (2 * mass *
                  (Complex(0, 2) * (1 + delta - 2 * Nu) * S1z +
                   Nu * ((1 + delta) * (-6 + kappa1) + 12 * Nu) * pow(S1z, 2) +
                   8 * Nu * (-1 + 3 * Nu) * S1z * S2z +
                   S2z * (Complex(0, -2) * (-1 + delta + 2 * Nu) +
                          Nu * (-6 - delta * (-6 + kappa2) + kappa2 + 12 * Nu) *
                              S2z)) +
              r * (4 * pow(rDOT, 2) *
                       (Complex(0, 2) * (1 + delta - 2 * Nu) * S1z +
                        Nu * ((1 + delta) * (-2 + kappa1) + 4 * Nu) *
                            pow(S1z, 2) +
                        8 * pow(Nu, 2) * S1z * S2z +
                        S2z * (Complex(0, -2) * (-1 + delta + 2 * Nu) +
                               Nu *
                                   (-2 - delta * (-2 + kappa2) + kappa2 +
                                    4 * Nu) *
                                   S2z)) +
                   2 * pow(PhiDOT, 2) * pow(r, 2) *
                       (Complex(0, 14) * (1 + delta - 2 * Nu) * S1z +
                        (6 * (1 + delta) * kappa1 -
                         (26 + 23 * kappa1 + delta * (26 + 11 * kappa1)) * Nu +
                         4 * (1 + 9 * kappa1) * pow(Nu, 2)) *
                            pow(S1z, 2) -
                        64 * pow(Nu, 2) * S1z * S2z +
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
                            pow(S1z, 2) +
                        Complex(0, 20) * (-3 + Nu) * Nu * S1z * S2z +
                        S2z * (-40 * (-1 + delta + 2 * Nu) -
                               Complex(0, 1) *
                                   (8 * Nu * (-2 + 2 * delta + Nu) +
                                    kappa2 * (3 - delta * (3 + 11 * Nu) +
                                              Nu * (5 + 18 * Nu))) *
                                   S2z))))) /
            (24. * pow(r, 4))) + /* Henry et al. QC spinning hereditary terms */
        ((8 * sqrt(0.7142857142857143) * Nu * M_PI * (S1z + S2z) * pow(x, 4)) /
         3.));
  }

  else if (vpnorder == 7) {

    return (
        /* Henry et al. QC spin terms */ /* ((-26902*pow(Nu,3)*(S1z + S2z) -
 4664*(S1z + delta*S1z + S2z - delta*S2z) + Nu*(3960*(1 +
 delta)*kappa1*pow(S1z,3) + 3960*(1 + delta)*kappa1*pow(S1z,2)*S2z + S1z*(28921
 - 18889*delta - 3960*(-1 + delta)*kappa2*pow(S2z,2)) + S2z*(28921 + 18889*delta
 - 3960*(-1 + delta)*kappa2*pow(S2z,2))) - 2*pow(Nu,2)*(1351*delta*(S1z - S2z) +
 6*(S1z + S2z)*(6773 + 660*(kappa1*pow(S1z,2)
 + S2z*(-2*S1z + kappa2*S2z)))))*pow(x,4.5))/(1188.*sqrt(35))
 + */
                                         /* Henry et al. ecc + spin terms */
        (-0.000014029180695847363 *
         (pow(mass, 2) *
          (3 * pow(r, 2) *
               (-120 * pow(Nu, 3) *
                    (3565 * pow(PhiDOT, 5) * pow(r, 5) +
                     Complex(0, 2321) * pow(PhiDOT, 4) * pow(r, 4) * rDOT +
                     8244 * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) -
                     Complex(0, 869) * pow(PhiDOT, 2) * pow(r, 2) *
                         pow(rDOT, 3) -
                     56 * PhiDOT * r * pow(rDOT, 4) -
                     Complex(0, 120) * pow(rDOT, 5)) *
                    (S1z + S2z) +
                2475 * pow(PhiDOT, 2) * pow(r, 2) *
                    (6 * pow(PhiDOT, 3) * pow(r, 3) +
                     Complex(0, 77) * pow(PhiDOT, 2) * pow(r, 2) * rDOT -
                     72 * PhiDOT * r * pow(rDOT, 2) +
                     Complex(0, 6) * pow(rDOT, 3)) *
                    (S1z + delta * S1z + S2z - delta * S2z) -
                3 * Nu *
                    (Complex(0, 22) * pow(PhiDOT, 4) * pow(r, 4) * rDOT *
                         (Complex(0, 36322) + 5 * (2993 + 3893 * delta) * S1z +
                          5 * (2993 - 3893 * delta) * S2z) -
                     Complex(0, 25) * pow(rDOT, 5) *
                         ((1053 + 443 * delta) * S1z +
                          (1053 - 443 * delta) * S2z) +
                     Complex(0, 44) * pow(PhiDOT, 2) * pow(r, 2) *
                         pow(rDOT, 3) *
                         (Complex(0, -5444) + 5 * (1424 + 849 * delta) * S1z +
                          7120 * S2z - 4245 * delta * S2z) -
                     20 * PhiDOT * r * pow(rDOT, 4) *
                         (Complex(0, -1782) + (5963 + 2969 * delta) * S1z +
                          5963 * S2z - 2969 * delta * S2z) +
                     4 * pow(PhiDOT, 5) * pow(r, 5) *
                         (Complex(0, -86889) + 10 * (2063 + 225 * delta) * S1z +
                          20630 * S2z - 2250 * delta * S2z) +
                     4 * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) *
                         (Complex(0, 234861) + 40 * (-1824 + 97 * delta) * S1z -
                          40 * (1824 + 97 * delta) * S2z)) +
                2 * pow(Nu, 2) *
                    (pow(PhiDOT, 5) * pow(r, 5) *
                         (Complex(0, -1549757) +
                          300 * (1448 + 1311 * delta) * S1z +
                          300 * (1448 - 1311 * delta) * S2z) +
                     Complex(0, 11) * pow(PhiDOT, 4) * pow(r, 4) * rDOT *
                         (Complex(0, 329548) +
                          15 * (4113 + 1411 * delta) * S1z + 61695 * S2z -
                          21165 * delta * S2z) +
                     Complex(0, 22) * pow(PhiDOT, 2) * pow(r, 2) *
                         pow(rDOT, 3) *
                         (Complex(0, -23971) + 15 * (3829 + 243 * delta) * S1z +
                          57435 * S2z - 3645 * delta * S2z) +
                     Complex(0, 150) * pow(rDOT, 5) *
                         ((-503 + 92 * delta) * S1z -
                          (503 + 92 * delta) * S2z) +
                     10 * PhiDOT * r * pow(rDOT, 4) *
                         (Complex(0, 4565) + 6 * (-6327 + 991 * delta) * S1z -
                          6 * (6327 + 991 * delta) * S2z) +
                     21 * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) *
                         (Complex(0, 161403) +
                          10 * (-1897 + 1471 * delta) * S1z -
                          10 * (1897 + 1471 * delta) * S2z))) -
           6 * mass * r *
               (60 * pow(Nu, 3) *
                    (2417 * pow(PhiDOT, 3) * pow(r, 3) +
                     Complex(0, 7258) * pow(PhiDOT, 2) * pow(r, 2) * rDOT -
                     4381 * PhiDOT * r * pow(rDOT, 2) +
                     Complex(0, 480) * pow(rDOT, 3)) *
                    (S1z + S2z) -
                165 *
                    (1161 * pow(PhiDOT, 3) * pow(r, 3) -
                     Complex(0, 536) * pow(PhiDOT, 2) * pow(r, 2) * rDOT -
                     2412 * PhiDOT * r * pow(rDOT, 2) -
                     Complex(0, 270) * pow(rDOT, 3)) *
                    (S1z + delta * S1z + S2z - delta * S2z) +
                2 * Nu *
                    (Complex(0, 1) * pow(PhiDOT, 2) * pow(r, 2) * rDOT *
                         (Complex(0, 1015784) +
                          5 * (30849 + 88721 * delta) * S1z +
                          5 * (30849 - 88721 * delta) * S2z) +
                     pow(PhiDOT, 3) * pow(r, 3) *
                         (Complex(0, -173371) +
                          5 * (61569 + 10789 * delta) * S1z + 307845 * S2z -
                          53945 * delta * S2z) -
                     5 * PhiDOT * r * pow(rDOT, 2) *
                         (Complex(0, -115368) + (177417 + 52307 * delta) * S1z +
                          177417 * S2z - 52307 * delta * S2z) +
                     Complex(0, 100) * pow(rDOT, 3) *
                         ((-1545 + 181 * delta) * S1z -
                          (1545 + 181 * delta) * S2z)) +
                pow(Nu, 2) *
                    (20 * PhiDOT * r * pow(rDOT, 2) *
                         (Complex(0, -11187) - 48074 * S1z +
                          11057 * delta * S1z - 48074 * S2z -
                          11057 * delta * S2z) +
                     Complex(0, 725) * pow(rDOT, 3) *
                         (-73 * S1z + 31 * delta * S1z - 73 * S2z -
                          31 * delta * S2z) +
                     pow(PhiDOT, 3) * pow(r, 3) *
                         (Complex(0, 603141) - 543040 * S1z +
                          404620 * delta * S1z -
                          20 * (27152 + 20231 * delta) * S2z) +
                     Complex(0, 1) * pow(PhiDOT, 2) * pow(r, 2) * rDOT *
                         (Complex(0, -1798104) - 648485 * S1z +
                          105755 * delta * S1z -
                          5 * (129697 + 21151 * delta) * S2z))) +
           10 * pow(mass, 2) *
               (24 * pow(Nu, 3) *
                    (6981 * PhiDOT * r + Complex(0, 1600) * rDOT) *
                    (S1z + S2z) -
                66 * (2027 * PhiDOT * r + Complex(0, 380) * rDOT) *
                    (S1z + delta * S1z + S2z - delta * S2z) +
                Complex(0, 30) * Nu * rDOT *
                    (297 * (1 + delta) * kappa1 * pow(S1z, 3) +
                     297 * (1 + delta) * kappa1 * pow(S1z, 2) * S2z +
                     S1z * (17261 - 1641 * delta -
                            297 * (-1 + delta) * kappa2 * pow(S2z, 2)) +
                     S2z * (17261 + 1641 * delta -
                            297 * (-1 + delta) * kappa2 * pow(S2z, 2))) +
                8 * Nu * PhiDOT * r *
                    (Complex(0, -7315) -
                     4455 * (1 + delta) * kappa1 * pow(S1z, 3) +
                     3 * (3881 - 7757 * delta) * S2z -
                     4455 * (1 + delta) * kappa1 * pow(S1z, 2) * S2z +
                     4455 * (-1 + delta) * kappa2 * pow(S2z, 3) +
                     3 * S1z *
                         (3881 + 7757 * delta +
                          1485 * (-1 + delta) * kappa2 * pow(S2z, 2))) +
                3 * pow(Nu, 2) *
                    (Complex(0, -5) * rDOT *
                         (S1z * (18793 + 223 * delta +
                                 1188 * kappa1 * pow(S1z, 2)) +
                          (18793 - 223 * delta +
                           1188 * (-2 + kappa1) * pow(S1z, 2)) *
                              S2z +
                          1188 * (-2 + kappa2) * S1z * pow(S2z, 2) +
                          1188 * kappa2 * pow(S2z, 3)) +
                     4 * PhiDOT * r *
                         (Complex(0, -4939) - 23359 * S1z - 5563 * delta * S1z +
                          5940 * kappa1 * pow(S1z, 3) +
                          (-23359 + 5563 * delta +
                           5940 * (-2 + kappa1) * pow(S1z, 2)) *
                              S2z +
                          5940 * (-2 + kappa2) * S1z * pow(S2z, 2) +
                          5940 * kappa2 * pow(S2z, 3)))))) /
         (sqrt(35) * pow(r, 4))));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hQC_3_m_2(REAL8 Nu, UINT4 vpnorder, REAL8 x) {

  /* if(vpnorder == 2){
       return(((2*sqrt(0.7142857142857143))/3. - 2*sqrt(0.7142857142857143)*Nu)*
       pow(x,2));
   }

   else if(vpnorder == 4){
       return((-193/(27.*sqrt(35)) + (145*sqrt(0.7142857142857143)*Nu)/27. -
       (73*sqrt(0.7142857142857143)*pow(Nu,2))/27.)*pow(x,3));
   } */

  /* else */ if (vpnorder == 5) {
    return ((4 * sqrt(0.7142857142857143) * (1 - 3 * Nu) * M_PI * pow(x, 3.5)) /
            3.);
  }

  /* else if(vpnorder == 6){
       return((-1451/(1188.*sqrt(35)) - (17387*Nu)/(1188.*sqrt(35)) +
       (5557*pow(Nu,2))/(66.*sqrt(35)) -
       (763*sqrt(1.4)*pow(Nu,3))/396.)*pow(x,4));
   } */

  else {
    return 0;
  }
}

static COMPLEX16 hl_3_m_2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z, REAL8 x) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_3_m_2: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_3_m_2(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x) +
            hQC_3_m_2(Nu, vpnorder, x)) *
           cpolar(1, -2 * Phi);
  }
}

static COMPLEX16 hl_3_m_min2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z, REAL8 x) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_3_m_min2: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) * pow(-1, 3) *
           conj(hGO_3_m_2(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x) +
                hQC_3_m_2(Nu, vpnorder, x)) *
           cpolar(1, 2 * Phi);
  }
}

// H31
static COMPLEX16 hGO_3_m_1(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z,
                           REAL8 x) {
  REAL8 delta = sqrt(1 - 4 * Nu);
  REAL8 kappa1 = 1.0;
  REAL8 kappa2 = 1.0;
  REAL8 r0 = 1.0;

  if (vpnorder == 1) {
    return (delta * (mass * (Complex(0, 7) * PhiDOT * r - 12 * rDOT) -
                     Complex(0, 6) * r * (PhiDOT * r - Complex(0, 1) * rDOT) *
                         pow(PhiDOT * r + Complex(0, 1) * rDOT, 2))) /
           (6. * sqrt(14) * r);
  }

  else if (vpnorder == 3) {
    return (delta *
            (Complex(0, 6) * (-5 + 19 * Nu) * pow(r, 2) *
                 pow(PhiDOT * r - Complex(0, 1) * rDOT, 2) *
                 pow(PhiDOT * r + Complex(0, 1) * rDOT, 3) +
             2 * pow(mass, 2) *
                 (Complex(0, 1) * (-101 + 43 * Nu) * PhiDOT * r +
                  (109 - 86 * Nu) * rDOT) +
             3 * mass * r *
                 (Complex(0, -4) * (-9 + 14 * Nu) * pow(PhiDOT, 3) * pow(r, 3) +
                  6 * (2 + 9 * Nu) * pow(PhiDOT, 2) * pow(r, 2) * rDOT -
                  Complex(0, 1) * (33 + 62 * Nu) * PhiDOT * r * pow(rDOT, 2) +
                  4 * (8 + 17 * Nu) * pow(rDOT, 3)))) /
           (36. * sqrt(14) * pow(r, 2));
  }

  else if (vpnorder == 4) {
    return (
        (Complex(0, 0.041666666666666664) * pow(mass, 2) *
         (4 * mass * (-1 + 5 * Nu) * ((1 + delta) * S1z + (-1 + delta) * S2z) -
          r * (2 * pow(rDOT, 2) *
                   (-6 * (1 + delta) * S1z + 5 * (5 + 3 * delta) * Nu * S1z +
                    6 * S2z - 6 * delta * S2z +
                    5 * (-5 + 3 * delta) * Nu * S2z) +
               pow(PhiDOT, 2) * pow(r, 2) *
                   ((24 + 24 * delta - 87 * Nu + 31 * delta * Nu) * S1z +
                    (-24 + 24 * delta + 87 * Nu + 31 * delta * Nu) * S2z) +
               Complex(0, 2) * PhiDOT * r * rDOT *
                   ((6 + 6 * delta - 31 * Nu + 35 * delta * Nu) * S1z +
                    (-6 + 6 * delta + 31 * Nu + 35 * delta * Nu) * S2z)))) /
            (sqrt(14) * pow(r, 3)) -
        ((Complex(0, 0.041666666666666664) *
          (-4 * S1z + 11 * Nu * (S1z - S2z) + 4 * S2z -
           4 * delta * (S1z + S2z) + 13 * delta * Nu * (S1z + S2z)) *
          pow(x, 3)) /
         sqrt(14)) +
        ((Complex(0, 0.041666666666666664) *
          (-4 * S1z + 11 * Nu * (S1z - S2z) + 4 * S2z -
           4 * delta * (S1z + S2z) + 13 * delta * Nu * (S1z + S2z)) *
          pow(x, 3)) /
         sqrt(14)));
  }

  else if (vpnorder == 5) {
    return (
        (delta *
         (Complex(0, -18) * (183 - 1579 * Nu + 3387 * pow(Nu, 2)) * pow(r, 3) *
              pow(PhiDOT * r - Complex(0, 1) * rDOT, 3) *
              pow(PhiDOT * r + Complex(0, 1) * rDOT, 4) +
          2 * pow(mass, 3) *
              (Complex(0, 1) * (26473 - 27451 * Nu + 9921 * pow(Nu, 2)) *
                   PhiDOT * r -
               12 * (623 - 732 * Nu + 1913 * pow(Nu, 2)) * rDOT) +
          2 * pow(mass, 2) * r *
              (Complex(0, -1) * (-8641 - 59189 * Nu + 31959 * pow(Nu, 2)) *
                   pow(PhiDOT, 3) * pow(r, 3) +
               (-32635 - 29345 * Nu + 29541 * pow(Nu, 2)) * pow(PhiDOT, 2) *
                   pow(r, 2) * rDOT -
               Complex(0, 44) * (-256 + 781 * Nu + 840 * pow(Nu, 2)) * PhiDOT *
                   r * pow(rDOT, 2) +
               6 * (-756 + 8238 * Nu + 7357 * pow(Nu, 2)) * pow(rDOT, 3)) +
          3 * mass * pow(r, 2) *
              (Complex(0, 2) * (-2479 - 4505 * Nu + 16785 * pow(Nu, 2)) *
                   pow(PhiDOT, 5) * pow(r, 5) +
               4 * (817 + 1220 * Nu - 7449 * pow(Nu, 2)) * pow(PhiDOT, 4) *
                   pow(r, 4) * rDOT +
               Complex(0, 6) * (-1679 + 1469 * Nu + 12233 * pow(Nu, 2)) *
                   pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) -
               32 * (-460 + 421 * Nu + 2514 * pow(Nu, 2)) * pow(PhiDOT, 2) *
                   pow(r, 2) * pow(rDOT, 3) +
               Complex(0, 1) * (-9851 + 17954 * Nu + 40968 * pow(Nu, 2)) *
                   PhiDOT * r * pow(rDOT, 4) -
               12 * (-771 + 1126 * Nu + 3616 * pow(Nu, 2)) * pow(rDOT, 5)))) /
            (4752. * sqrt(14) *
             pow(r,
                 3)) /* Henry et al. QC spin terms */ /* +((Complex(0,-0.041666666666666664)*(kappa1*(5
- 4*Nu + delta*(5 + 6*Nu))*pow(S1z,2) + S2z*(-12*delta*Nu*S1z + kappa2*(-5 +
5*delta + 4*Nu + 6*delta*Nu)*S2z))* pow(x,3.5))/sqrt(14)) */
        +
        /* Henry et al. ecc spin terms */ (
            (pow(mass, 3) *
             (kappa1 *
                  (Complex(0, -1) *
                       (-13 - 13 * delta + 68 * Nu + 42 * delta * Nu) * PhiDOT *
                       r -
                   34 * (1 + delta) * rDOT + 4 * (26 + 9 * delta) * Nu * rDOT) *
                  pow(S1z, 2) +
              S2z * (12 * delta * Nu * (Complex(0, 7) * PhiDOT * r - 6 * rDOT) *
                         S1z +
                     kappa2 *
                         (Complex(0, -1) *
                              (13 - 13 * delta - 68 * Nu + 42 * delta * Nu) *
                              PhiDOT * r +
                          2 * (17 - 17 * delta - 52 * Nu + 18 * delta * Nu) *
                              rDOT) *
                         S2z))) /
            (24. * sqrt(14) * pow(r, 3))));
  }

  else if (vpnorder == 6) {
    return (
        (delta * pow(mass, 2) * Nu *
         (668 * pow(mass, 2) -
          2 * mass * r *
              (727 * pow(PhiDOT, 2) * pow(r, 2) -
               Complex(0, 99) * PhiDOT * r * rDOT + 452 * pow(rDOT, 2)) +
          pow(r, 2) * (-499 * pow(PhiDOT, 4) * pow(r, 4) +
                       Complex(0, 1534) * pow(PhiDOT, 3) * pow(r, 3) * rDOT +
                       3072 * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) -
                       Complex(0, 680) * PhiDOT * r * pow(rDOT, 3) +
                       1000 * pow(rDOT, 4)))) /
            (180. * sqrt(14) *
             pow(r,
                 4)) /*Henry et al QC spin terms */ /* +((Complex(0,0.004629629629629629)*(70*(S1z
- S2z
+ delta*(S1z + S2z)) - pow(Nu,2)*(931*(S1z - S2z) + 45*delta*(S1z + S2z)) -
Nu*(59*(-S1z + S2z) + 99*delta*(S1z + S2z)))*pow(x,4))/sqrt(14)) */
        + /* Henry et al. ecc spin terms */ (
              (Complex(0, 0.0023148148148148147) * pow(mass, 2) *
               (2 * pow(mass, 2) *
                    ((252 * (1 + delta) - (1277 + 1279 * delta) * Nu +
                      8 * (12 + 47 * delta) * pow(Nu, 2)) *
                         S1z +
                     (252 * (-1 + delta) + (1277 - 1279 * delta) * Nu +
                      8 * (-12 + 47 * delta) * pow(Nu, 2)) *
                         S2z) +
                3 * pow(r, 2) *
                    (2 * pow(rDOT, 4) *
                         ((30 + Nu * (-187 + 318 * Nu) +
                           delta * (30 + Nu * (-101 + 122 * Nu))) *
                              S1z +
                          (-30 + (187 - 318 * Nu) * Nu +
                           delta * (30 + Nu * (-101 + 122 * Nu))) *
                              S2z) +
                     2 * pow(PhiDOT, 4) * pow(r, 4) *
                         ((90 - Nu * (28 + 579 * Nu) +
                           delta * (90 + Nu * (-800 + 551 * Nu))) *
                              S1z +
                          (-90 + Nu * (28 + 579 * Nu) +
                           delta * (90 + Nu * (-800 + 551 * Nu))) *
                              S2z) +
                     Complex(0, 2) * PhiDOT * r * pow(rDOT, 3) *
                         ((186 - Nu * (745 + 354 * Nu) +
                           delta * (186 + Nu * (-191 + 554 * Nu))) *
                              S1z +
                          (-186 + Nu * (745 + 354 * Nu) +
                           delta * (186 + Nu * (-191 + 554 * Nu))) *
                              S2z) +
                     3 * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) *
                         ((32 + Nu * (-451 + 230 * Nu) +
                           delta * (32 + Nu * (691 + 626 * Nu))) *
                              S1z +
                          (-32 + (451 - 230 * Nu) * Nu +
                           delta * (32 + Nu * (691 + 626 * Nu))) *
                              S2z) +
                     Complex(0, 2) * pow(PhiDOT, 3) * pow(r, 3) * rDOT *
                         ((-12 + Nu * (-341 + 315 * Nu) +
                           delta * (-12 + Nu * (-91 + 1213 * Nu))) *
                              S1z +
                          (12 + (341 - 315 * Nu) * Nu +
                           delta * (-12 + Nu * (-91 + 1213 * Nu))) *
                              S2z)) -
                2 * mass * r *
                    (2 * pow(PhiDOT, 2) * pow(r, 2) *
                         ((-312 * (1 + delta) + 2 * (827 - 923 * delta) * Nu +
                           5 * (-201 + 131 * delta) * pow(Nu, 2)) *
                              S1z +
                          (-312 * (-1 + delta) - 2 * (827 + 923 * delta) * Nu +
                           5 * (201 + 131 * delta) * pow(Nu, 2)) *
                              S2z) +
                     2 * pow(rDOT, 2) *
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
              (sqrt(14) * pow(r, 4))));
  }

  else if (vpnorder == 7) {

    return (/* (Complex(0,0.001388888888888889)*pow(x,4.5)*(10*(32*(1 +
      delta)*kappa1
      - (36*(1 + delta) + (65 + delta)*kappa1)*Nu + 2*(150 + 33*delta + (76 +
      22*delta)*kappa1)*pow(Nu,2))*pow(S1z,2) + S2z*(Complex(0,-315)*delta*Nu +
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
                  (1513512 * pow(mass, 3) *
                       (2 * mass *
                            (PhiDOT * r *
                                 (Complex(0, 1) * (-97 + 631 * Nu) * S1z +
                                  5 *
                                      (8 + 16 * Nu * (-7 + 15 * Nu) +
                                       3 * kappa1 *
                                           (-39 + Nu * (149 + 4 * Nu))) *
                                      pow(S1z, 2) +
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
                                      pow(S1z, 2) +
                                  S2z *
                                      (36 - 166 * Nu +
                                       Complex(0, 5) *
                                           (-4 * (6 + Nu * (-25 + 7 * Nu)) +
                                            kappa2 * (155 +
                                                      Nu * (-373 + 164 * Nu))) *
                                           S2z))) +
                        r * (2 * pow(rDOT, 3) *
                                 (S1z *
                                      (18 - 110 * Nu -
                                       Complex(0, 5) *
                                           (69 * kappa1 - 214 * kappa1 * Nu +
                                            4 * (5 + 4 * kappa1) * pow(Nu, 2)) *
                                           S1z) +
                                  2 * (-9 + 55 * Nu) * S2z +
                                  Complex(0, 5) *
                                      (69 * kappa2 - 214 * kappa2 * Nu +
                                       4 * (5 + 4 * kappa2) * pow(Nu, 2)) *
                                      pow(S2z, 2)) +
                             pow(PhiDOT, 3) * pow(r, 3) *
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
                                      pow(S2z, 2)) +
                             2 * pow(PhiDOT, 2) * pow(r, 2) * rDOT *
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
                                      pow(S2z, 2)) +
                             PhiDOT * r * pow(rDOT, 2) *
                                 (Complex(0, 4) * (-114 + 781 * Nu) * S1z +
                                  5 *
                                      (-213 * kappa1 - 72 * Nu +
                                       278 * kappa1 * Nu +
                                       8 * (7 + 44 * kappa1) * pow(Nu, 2)) *
                                      pow(S1z, 2) +
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
                            pow(r, 4) *
                            pow(PhiDOT * r - Complex(0, 1) * rDOT, 4) *
                            pow(PhiDOT * r + Complex(0, 1) * rDOT, 5) +
                        840 * pow(mass, 2) * pow(r, 2) *
                            ((-2555489 +
                              7 * Nu *
                                  (820078 + Nu * (-6623390 + 4948497 * Nu))) *
                                 pow(PhiDOT, 5) * pow(r, 5) +
                             Complex(0, 1) *
                                 (3537631 +
                                  7 * Nu *
                                      (-2817653 +
                                       Nu * (-7052042 + 4017147 * Nu))) *
                                 pow(PhiDOT, 4) * pow(r, 4) * rDOT +
                             3 * (-1428997 + 7 * Nu * (-1230747 + Nu * (-237418 + 4061717 * Nu))) *
                                 pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) +
                             Complex(0, 1) *
                                 (-5153011 +
                                  7 * Nu *
                                      (-2375327 +
                                       9 * Nu * (218846 + 1640185 * Nu))) *
                                 pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 3) +
                             (-7761899 + 7 * Nu *
                                             (2892563 + 5998602 * Nu +
                                              7493619 * pow(Nu, 2))) *
                                 PhiDOT * r * pow(rDOT, 4) +
                             Complex(0, 3) *
                                 (-2422057 +
                                  7 * Nu *
                                      (501045 +
                                       Nu * (2033141 + 2771816 * Nu))) *
                                 pow(rDOT, 5)) -
                        8820 * mass * pow(r, 3) *
                            (2 *
                                 (111737 +
                                  Nu * (-366573 +
                                        Nu * (-618923 + 2278593 * Nu))) *
                                 pow(PhiDOT, 7) * pow(r, 7) +
                             Complex(0, 2) *
                                 (101844 + Nu * (-273675 - 871630 * Nu +
                                                 2069774 * pow(Nu, 2))) *
                                 pow(PhiDOT, 6) * pow(r, 6) * rDOT +
                             2 * (341322 + Nu * (-1429938 + Nu * (-1206083 + 7690681 * Nu))) *
                                 pow(PhiDOT, 5) * pow(r, 5) * pow(rDOT, 2) +
                             Complex(0, 8) *
                                 (90241 + 2 * Nu *
                                              (-206022 +
                                               Nu * (-62113 + 1003558 * Nu))) *
                                 pow(PhiDOT, 4) * pow(r, 4) * pow(rDOT, 3) +
                             2 * (410547 + Nu * (-2269686 + Nu * (762091 + 8400052 * Nu))) *
                                 pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 4) +
                             Complex(0, 4) *
                                 (217935 +
                                  2 * Nu *
                                      (-573699 +
                                       5 * Nu * (18671 + 445748 * Nu))) *
                                 pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 5) +
                             (333969 +
                              2 * Nu *
                                  (-981471 + 4 * Nu * (154039 + 750016 * Nu))) *
                                 PhiDOT * r * pow(rDOT, 6) +
                             Complex(0, 24) *
                                 (13245 +
                                  2 * Nu *
                                      (-37005 + Nu * (14251 + 130160 * Nu))) *
                                 pow(rDOT, 7)) +
                        2 * pow(mass, 4) *
                            (Complex(0, -4178597424) * rDOT +
                             Complex(0, 84) * rDOT *
                                 (38468500 * pow(Nu, 3) +
                                  Complex(0, 648648) * (S1z + S2z) -
                                  90090 * ((-24 + 155 * kappa1) * pow(S1z, 2) +
                                           (-24 + 155 * kappa2) * pow(S2z, 2)) -
                                  420 * pow(Nu, 2) *
                                      (-122855 +
                                       3003 *
                                           ((2 + 11 * kappa1) * pow(S1z, 2) -
                                            18 * S1z * S2z +
                                            (2 + 11 * kappa2) * pow(S2z, 2))) +
                                  3 * Nu *
                                      (-103100846 - 1846845 * pow(M_PI, 2) -
                                       Complex(0, 564564) * S2z +
                                       6006 * (S1z * (Complex(0, -94) +
                                                      5 * (-52 + 63 * kappa1) *
                                                          S1z) -
                                               20 * S1z * S2z +
                                               5 * (-52 + 63 * kappa2) *
                                                   pow(S2z, 2)))) +
                             PhiDOT * r *
                                 (1176172480 * pow(Nu, 3) +
                                  8 * (74084729 -
                                       189189 * S1z *
                                           (Complex(0, 97) +
                                            5 * (-8 + 117 * kappa1) * S1z) -
                                       189189 * S2z *
                                           (Complex(0, 97) +
                                            5 * (-8 + 117 * kappa2) * S2z)) -
                                  176400 * pow(Nu, 2) *
                                      (11251 +
                                       429 *
                                           ((2 + 13 * kappa1) * pow(S1z, 2) -
                                            22 * S1z * S2z +
                                            (2 + 13 * kappa2) * pow(S2z, 2))) +
                                  147 * Nu *
                                      (-65012788 + 4485195 * pow(M_PI, 2) +
                                       Complex(0, 4499352) * S2z +
                                       10296 *
                                           (S1z * (Complex(0, 437) +
                                                   15 * (-32 + 71 * kappa1) *
                                                       S1z) -
                                            3860 * S1z * S2z +
                                            15 * (-32 + 71 * kappa2) *
                                                pow(S2z, 2))))) -
                        3 * pow(mass, 3) * r *
                            (Complex(0, -4) * pow(rDOT, 3) *
                                 (601018232 - 1359334480 * pow(Nu, 3) -
                                  756756 * S1z *
                                      (Complex(0, 6) + 115 * kappa1 * S1z) -
                                  756756 * S2z *
                                      (Complex(0, 6) + 115 * kappa2 * S2z) +
                                  231 * Nu *
                                      (8490448 + 503685 * pow(M_PI, 2) +
                                       Complex(0, 80808) * S2z +
                                       2184 * (S1z * (Complex(0, 37) +
                                                      190 * kappa1 * S1z) +
                                               70 * S1z * S2z +
                                               190 * kappa2 * pow(S2z, 2))) +
                                  58800 * pow(Nu, 2) *
                                      (-62596 +
                                       429 *
                                           ((-1 + 5 * kappa1) * pow(S1z, 2) -
                                            12 * S1z * S2z +
                                            (-1 + 5 * kappa2) * pow(S2z, 2)))) -
                             Complex(0, 14) * pow(PhiDOT, 2) * pow(r, 2) *
                                 rDOT *
                                 (-229522160 * pow(Nu, 3) +
                                  8 * (48303859 +
                                       135135 * S1z *
                                           (Complex(0, -17) +
                                            2 * kappa1 * S1z) +
                                       135135 * S2z *
                                           (Complex(0, -17) +
                                            2 * kappa2 * S2z)) +
                                  2520 * pow(Nu, 2) *
                                      (100913 +
                                       286 *
                                           ((-31 + 5 * kappa1) * pow(S1z, 2) -
                                            72 * S1z * S2z +
                                            (-31 + 5 * kappa2) * pow(S2z, 2))) +
                                  7 * Nu *
                                      (125038052 + 2374515 * pow(M_PI, 2) +
                                       Complex(0, 5858424) * S2z -
                                       10296 * (S1z * (Complex(0, -569) +
                                                       25 * (-24 + 7 * kappa1) *
                                                           S1z) +
                                                700 * S1z * S2z +
                                                25 * (-24 + 7 * kappa2) *
                                                    pow(S2z, 2)))) +
                             4 * PhiDOT * r * pow(rDOT, 2) *
                                 (-1095987374 + 1035895280 * pow(Nu, 3) +
                                  378378 * S1z *
                                      (Complex(0, 152) + 355 * kappa1 * S1z) +
                                  378378 * S2z *
                                      (Complex(0, 152) + 355 * kappa2 * S2z) -
                                  490 * pow(Nu, 2) *
                                      (-5802767 +
                                       5148 *
                                           ((2 + 23 * kappa1) * pow(S1z, 2) -
                                            42 * S1z * S2z +
                                            (2 + 23 * kappa2) * pow(S2z, 2))) -
                                  77 * Nu *
                                      (42451610 + 1511055 * pow(M_PI, 2) +
                                       Complex(0, 3623256) * S2z -
                                       6552 * (S1z * (Complex(0, -553) +
                                                      5 * (18 + 37 * kappa1) *
                                                          S1z) +
                                               965 * S1z * S2z +
                                               5 * (18 + 37 * kappa2) *
                                                   pow(S2z, 2)))) +
                             7 * pow(PhiDOT, 3) * pow(r, 3) *
                                 (512893080 * pow(Nu, 3) -
                                  136 *
                                      (-2089567 +
                                       135135 * S1z *
                                           (Complex(0, 1) + 2 * kappa1 * S1z) +
                                       135135 * S2z *
                                           (Complex(0, 1) + 2 * kappa2 * S2z)) -
                                  560 * pow(Nu, 2) *
                                      (2457671 +
                                       2574 *
                                           ((11 + 53 * kappa1) * pow(S1z, 2) -
                                            84 * S1z * S2z +
                                            (11 + 53 * kappa2) * pow(S2z, 2))) +
                                  3 * Nu *
                                      (16621605 * pow(M_PI, 2) +
                                       8 * (27468722 +
                                            Complex(0, 2681679) * S2z +
                                            3003 * (S1z * (Complex(0, 893) -
                                                           840 * S1z +
                                                           800 * kappa1 * S1z) -
                                                    3160 * S1z * S2z +
                                                    40 * (-21 + 20 * kappa2) *
                                                        pow(S2z, 2)))))))) +
              74954880 * delta * pow(mass, 3) *
                  (mass * (Complex(0, -22) * PhiDOT * r - 24 * rDOT) +
                   3 * r *
                       (Complex(0, 7) * pow(PhiDOT, 3) * pow(r, 3) +
                        14 * pow(PhiDOT, 2) * pow(r, 2) * rDOT -
                        Complex(0, 8) * PhiDOT * r * pow(rDOT, 2) +
                        4 * pow(rDOT, 3))) *
                  log(r / r0)) /
             (3.6324288e8 * sqrt(14) *
              pow(r, 4))) + /* Henry et al. QC spinning hereditary terms */
            ((((-4 + 11 * Nu + delta * (-4 + 13 * Nu)) * S1z +
               (4 - 11 * Nu + delta * (-4 + 13 * Nu)) * S2z) *
              pow(x, 4.5) * (Complex(0, 1) * M_PI + log(4))) /
             (24. * sqrt(14))));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hQC_3_m_1(REAL8 Nu, UINT4 vpnorder, REAL8 x) {
  REAL8 delta = sqrt(1 - 4 * Nu);
  double EulerGamma = 0.5772156649015329;

  /* if(vpnorder == 1){
        return((Complex(0,0.16666666666666666)*delta*pow(x,1.5))/sqrt(14));
    } */

  /* else if(vpnorder == 3){
        return((Complex(0,-0.2222222222222222)*sqrt(0.2857142857142857)*delta -
        (Complex(0,0.1111111111111111)*delta*Nu)/sqrt(14))*pow(x,2.5));
    } */

  /* else */ if (vpnorder == 4) {
    return ((Complex(0, 0.16666666666666666) * delta * M_PI * pow(x, 3)) /
                sqrt(14) +
            (delta * pow(x, 3) * log(2)) / (3. * sqrt(14)));
  }

  /* else if(vpnorder == 5){
        return(((Complex(0,0.5109427609427609)*delta)/sqrt(14) -
        Complex(0,0.11447811447811448)*sqrt(0.2857142857142857)*delta*Nu -
        (Complex(0,0.2079124579124579)*delta*pow(Nu,2))/sqrt(14))*
        pow(x,3.5));
    } */

  else if (vpnorder == 6) {
    return ((Complex(0, -0.027777777777777776) * delta * (16 + 7 * Nu) *
             pow(x, 4) * (M_PI - Complex(0, 2) * log(2))) /
            sqrt(14));
  }

  else if (vpnorder == 7) {
    return ((Complex(0, -2.752978943455134e-9) * delta * pow(x, 4.5) *
             (-430135880 + 74954880 * EulerGamma + 681626456 * Nu -
              641035640 * pow(Nu, 2) + 68698000 * pow(Nu, 3) +
              Complex(0, 47279232) * M_PI - 10090080 * pow(M_PI, 2) -
              38783745 * Nu * pow(M_PI, 2) + 244468224 * log(2) +
              Complex(0, 121080960) * M_PI * log(2) +
              121080960 * pow(log(2), 2) + 37477440 * log(x))) /
            sqrt(14));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hl_3_m_1(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z, REAL8 x) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_4_m_1: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_3_m_1(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x) +
            hQC_3_m_1(Nu, vpnorder, x)) *
           cpolar(1, -1 * Phi);
  }
}

static COMPLEX16 hl_3_m_min1(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z, REAL8 x) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_4_m_min1: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) * pow(-1, 3) *
           conj(hGO_3_m_1(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x) +
                hQC_3_m_1(Nu, vpnorder, x)) *
           cpolar(1, 1 * Phi);
  }
}

// H44
static COMPLEX16 hGO_4_m_4(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z) {
  REAL8 delta = sqrt(1 - 4 * Nu);
  REAL8 kappa1 = 1.0;
  REAL8 kappa2 = 1.0;

  if (vpnorder == 2) {
    return ((sqrt(0.7142857142857143) * (-1 + 3 * Nu) *
             (7 * pow(mass, 2) +
              6 * pow(r, 2) * pow(PhiDOT * r + Complex(0, 1) * rDOT, 4) +
              3 * mass * r *
                  (17 * pow(PhiDOT, 2) * pow(r, 2) +
                   Complex(0, 18) * PhiDOT * r * rDOT - 6 * pow(rDOT, 2)))) /
            (36. * pow(r, 2)));
  }

  else if (vpnorder == 4) {
    return ((40 * pow(mass, 3) * (314 - 987 * Nu + 195 * pow(Nu, 2)) -
             60 * (23 - 159 * Nu + 291 * pow(Nu, 2)) * pow(r, 3) *
                 (PhiDOT * r - Complex(0, 1) * rDOT) *
                 pow(PhiDOT * r + Complex(0, 1) * rDOT, 5) +
             pow(mass, 2) * r *
                 ((53143 - 199660 * Nu + 127500 * pow(Nu, 2)) * pow(PhiDOT, 2) *
                      pow(r, 2) +
                  Complex(0, 24) * (967 - 4615 * Nu + 5935 * pow(Nu, 2)) *
                      PhiDOT * r * rDOT -
                  10 * (290 - 2033 * Nu + 4365 * pow(Nu, 2)) * pow(rDOT, 2)) -
             3 * mass * pow(r, 2) *
                 ((613 - 920 * Nu + 6420 * pow(Nu, 2)) * pow(PhiDOT, 4) *
                      pow(r, 4) -
                  Complex(0, 8) * (-976 + 1745 * Nu + 3150 * pow(Nu, 2)) *
                      pow(PhiDOT, 3) * pow(r, 3) * rDOT +
                  2 * (-6141 + 8980 * Nu + 31500 * pow(Nu, 2)) *
                      pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) +
                  Complex(0, 4) * (-1853 + 1730 * Nu + 13230 * pow(Nu, 2)) *
                      PhiDOT * r * pow(rDOT, 3) -
                  20 * (-83 + 30 * Nu + 762 * pow(Nu, 2)) * pow(rDOT, 4))) /
            (1584. * sqrt(35) * pow(r, 3)));
  }

  else if (vpnorder == 5) {
    return ((pow(mass, 2) * Nu *
             (6 * mass * (Complex(0, -43) * PhiDOT * r + 9 * rDOT) +
              r * (Complex(0, -734) * pow(PhiDOT, 3) * pow(r, 3) +
                   129 * pow(PhiDOT, 2) * pow(r, 2) * rDOT +
                   Complex(0, 156) * PhiDOT * r * pow(rDOT, 2) -
                   26 * pow(rDOT, 3)))) /
                (24. * sqrt(35) * pow(r, 3))
            /* Henry et al. QC spin terms */
            /* +((-32*(39*delta*Nu*(S1z - S2z) + 41*Nu*(S1z + S2z) -
      42*pow(Nu,2)*(S1z + S2z)
      - 10*(S1z + delta*S1z + S2z - delta*S2z))*pow(x,3.5))/(27.*sqrt(35))) */
            + /* Henry et al. ecc spin terms */ (
                  (pow(mass, 2) *
                   (Complex(0, -3) * pow(PhiDOT, 2) * pow(r, 3) * rDOT *
                        ((-250 + 1221 * Nu - 1512 * pow(Nu, 2) +
                          delta * (-250 + 849 * Nu)) *
                             S1z +
                         (-250 + delta * (250 - 849 * Nu) + 1221 * Nu -
                          1512 * pow(Nu, 2)) *
                             S2z) -
                    2 * mass * PhiDOT * r *
                        ((-130 + 757 * Nu - 1224 * pow(Nu, 2) +
                          delta * (-130 + 513 * Nu)) *
                             S1z +
                         (-130 + delta * (130 - 513 * Nu) + 757 * Nu -
                          1224 * pow(Nu, 2)) *
                             S2z) -
                    Complex(0, 2) * mass * rDOT *
                        ((-100 + 577 * Nu - 864 * pow(Nu, 2) +
                          delta * (-100 + 333 * Nu)) *
                             S1z +
                         (-100 + delta * (100 - 333 * Nu) + 577 * Nu -
                          864 * pow(Nu, 2)) *
                             S2z) -
                    6 * pow(PhiDOT, 3) * pow(r, 4) *
                        ((-65 + 263 * Nu - 291 * pow(Nu, 2) +
                          delta * (-65 + 282 * Nu)) *
                             S1z +
                         (-65 + delta * (65 - 282 * Nu) + 263 * Nu -
                          291 * pow(Nu, 2)) *
                             S2z) +
                    12 * PhiDOT * pow(r, 2) * pow(rDOT, 2) *
                        ((-40 + 201 * Nu - 252 * pow(Nu, 2) +
                          delta * (-40 + 129 * Nu)) *
                             S1z +
                         (-40 + delta * (40 - 129 * Nu) + 201 * Nu -
                          252 * pow(Nu, 2)) *
                             S2z) +
                    Complex(0, 6) * r * pow(rDOT, 3) *
                        ((-20 + 107 * Nu - 144 * pow(Nu, 2) +
                          delta * (-20 + 63 * Nu)) *
                             S1z +
                         (-20 + delta * (20 - 63 * Nu) + 107 * Nu -
                          144 * pow(Nu, 2)) *
                             S2z))) /
                  (72. * sqrt(35) * pow(r, 3))));
  }

  else if (vpnorder == 6) {
    return (
        (10 * pow(mass, 4) *
             (-4477296 + 12734393 * Nu - 6895 * pow(Nu, 2) +
              1043805 * pow(Nu, 3)) +
         3150 * (-367 + 4337 * Nu - 17462 * pow(Nu, 2) + 23577 * pow(Nu, 3)) *
             pow(r, 4) * pow(PhiDOT * r - Complex(0, 1) * rDOT, 2) *
             pow(PhiDOT * r + Complex(0, 1) * rDOT, 6) +
         2 * pow(mass, 3) * r *
             ((-36967579 + 245501977 * Nu - 459916170 * pow(Nu, 2) +
               150200680 * pow(Nu, 3)) *
                  pow(PhiDOT, 2) * pow(r, 2) +
              Complex(0, 4) *
                  (7571073 - 10780154 * Nu - 56898800 * pow(Nu, 2) +
                   43665510 * pow(Nu, 3)) *
                  PhiDOT * r * rDOT -
              10 *
                  (1283609 - 5800627 * Nu + 3725295 * pow(Nu, 2) +
                   4771935 * pow(Nu, 3)) *
                  pow(rDOT, 2)) -
         pow(mass, 2) * pow(r, 2) *
             ((-28258134 + 3245207 * Nu + 144051250 * pow(Nu, 2) +
               136991820 * pow(Nu, 3)) *
                  pow(PhiDOT, 4) * pow(r, 4) -
              Complex(0, 24) *
                  (2371982 - 7733376 * Nu - 7948185 * pow(Nu, 2) +
                   9074870 * pow(Nu, 3)) *
                  pow(PhiDOT, 3) * pow(r, 3) * rDOT +
              7 *
                  (6557973 - 50558069 * Nu + 59901380 * pow(Nu, 2) +
                   104752320 * pow(Nu, 3)) *
                  pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) +
              Complex(0, 168) *
                  (52044 - 1084807 * Nu + 1849450 * pow(Nu, 2) +
                   4171730 * pow(Nu, 3)) *
                  PhiDOT * r * pow(rDOT, 3) -
              35 *
                  (1083 - 1246819 * Nu + 2524240 * pow(Nu, 2) +
                   5995845 * pow(Nu, 3)) *
                  pow(rDOT, 4)) -
         105 * mass * pow(r, 3) *
             ((116396 - 551405 * Nu + 560658 * pow(Nu, 2) +
               293036 * pow(Nu, 3)) *
                  pow(PhiDOT, 6) * pow(r, 6) +
              Complex(0, 2) *
                  (158192 - 670661 * Nu + 177718 * pow(Nu, 2) +
                   2163976 * pow(Nu, 3)) *
                  pow(PhiDOT, 5) * pow(r, 5) * rDOT +
              (-393665 + 1322392 * Nu + 1589680 * pow(Nu, 2) -
               8622660 * pow(Nu, 3)) *
                  pow(PhiDOT, 4) * pow(r, 4) * pow(rDOT, 2) -
              Complex(0, 8) *
                  (-23048 + 209397 * Nu - 487057 * pow(Nu, 2) +
                   260396 * pow(Nu, 3)) *
                  pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 3) -
              (630647 - 3391000 * Nu + 2501958 * pow(Nu, 2) +
               7664096 * pow(Nu, 3)) *
                  pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 4) -
              Complex(0, 2) *
                  (218975 - 1037408 * Nu + 148970 * pow(Nu, 2) +
                   3699480 * pow(Nu, 3)) *
                  PhiDOT * r * pow(rDOT, 5) +
              10 *
                  (10233 - 44864 * Nu - 13050 * pow(Nu, 2) +
                   203280 * pow(Nu, 3)) *
                  pow(rDOT, 6))) /
            (1.44144e6 * sqrt(35) * pow(r, 4))
        /* Henry et al. QC spin terms */ /* +((16*sqrt(0.7142857142857143)*(-1 +
  3*Nu)*(kappa1*(1 + delta - 2*Nu)*pow(S1z,2)
  + S2z*(4*Nu*S1z - kappa2*(-1 + delta + 2*Nu)*S2z))*pow(x,4))/9.) */
        + /* Henry et al. ecc spin terms */ (
              (sqrt(0.7142857142857143) * pow(mass, 3) * (-1 + 3 * Nu) *
               (12 * mass +
                r * (53 * pow(PhiDOT, 2) * pow(r, 2) +
                     Complex(0, 26) * PhiDOT * r * rDOT - 8 * pow(rDOT, 2))) *
               (kappa1 * (1 + delta - 2 * Nu) * pow(S1z, 2) +
                S2z * (4 * Nu * S1z + kappa2 * S2z - delta * kappa2 * S2z -
                       2 * kappa2 * Nu * S2z))) /
              (48. * pow(r, 4))));
  }

  else if (vpnorder == 7) {

    return (
        /* Henry et al. QC spin terms */ /* (8*(6630*pow(Nu,3)*(S1z + S2z) -
  786*(S1z + delta*S1z + S2z - delta*S2z)
  + 21*Nu*(273*delta*(S1z - S2z) + 173*(S1z + S2z)) -
  2*pow(Nu,2)*(1557*delta*(S1z - S2z)
   + 2960*(S1z + S2z)))*pow(x,4.5))/(297.*sqrt(35))
   + */
        /* Henry et al. ecc+spin terms */ (
            (pow(mass, 2) *
             (14 * pow(mass, 2) *
                  (120 * pow(Nu, 3) *
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
                   6 * pow(Nu, 2) *
                       (2 * PhiDOT * r *
                            (Complex(0, 217374) + 549175 * S1z +
                             61075 * delta * S1z + 549175 * S2z -
                             61075 * delta * S2z) +
                        Complex(0, 5) * rDOT *
                            (Complex(0, 9861) + 132241 * S1z +
                             18825 * delta * S1z + 132241 * S2z -
                             18825 * delta * S2z))) -
              3 * pow(r, 2) *
                  (1680 * pow(Nu, 3) *
                       (2833 * pow(PhiDOT, 5) * pow(r, 5) +
                        Complex(0, 18796) * pow(PhiDOT, 4) * pow(r, 4) * rDOT -
                        13185 * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) +
                        Complex(0, 1186) * pow(PhiDOT, 2) * pow(r, 2) *
                            pow(rDOT, 3) -
                        5863 * PhiDOT * r * pow(rDOT, 4) -
                        Complex(0, 2100) * pow(rDOT, 5)) *
                       (S1z + S2z) -
                   1050 *
                       (454 * pow(PhiDOT, 5) * pow(r, 5) +
                        Complex(0, 1195) * pow(PhiDOT, 4) * pow(r, 4) * rDOT -
                        1950 * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) -
                        Complex(0, 442) * pow(PhiDOT, 2) * pow(r, 2) *
                            pow(rDOT, 3) -
                        384 * PhiDOT * r * pow(rDOT, 4) -
                        Complex(0, 184) * pow(rDOT, 5)) *
                       (S1z + delta * S1z + S2z - delta * S2z) -
                   6 * pow(Nu, 2) *
                       (2 * pow(PhiDOT, 5) * pow(r, 5) *
                            (Complex(0, -2459811) +
                             35 * (26517 + 10223 * delta) * S1z +
                             (928095 - 357805 * delta) * S2z) -
                        Complex(0, 10) * pow(rDOT, 5) *
                            (Complex(0, 35291) +
                             28 * (2183 + 1155 * delta) * S1z +
                             (61124 - 32340 * delta) * S2z) +
                        Complex(0, 1) * pow(PhiDOT, 2) * pow(r, 2) *
                            pow(rDOT, 3) *
                            (Complex(0, 917901) +
                             1120 * (-191 + 1616 * delta) * S1z -
                             1120 * (191 + 1616 * delta) * S2z) +
                        5 * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) *
                            (Complex(0, 85426) +
                             7 * (-148363 + 4365 * delta) * S1z -
                             7 * (148363 + 4365 * delta) * S2z) -
                        4 * PhiDOT * r * pow(rDOT, 4) *
                            (Complex(0, 280067) +
                             70 * (5844 + 4817 * delta) * S1z -
                             70 * (-5844 + 4817 * delta) * S2z) +
                        Complex(0, 1) * pow(PhiDOT, 4) * pow(r, 4) * rDOT *
                            (Complex(0, 10375501) +
                             70 * (65831 + 22871 * delta) * S1z -
                             70 * (-65831 + 22871 * delta) * S2z)) +
                   Nu * (Complex(0, -40) * pow(rDOT, 5) *
                             (Complex(0, 12203) +
                              42 * (843 + 656 * delta) * S1z -
                              42 * (-843 + 656 * delta) * S2z) +
                         Complex(0, 4) * pow(PhiDOT, 2) * pow(r, 2) *
                             pow(rDOT, 3) *
                             (Complex(0, -266071) +
                              210 * (-2279 + 1010 * delta) * S1z -
                              210 * (2279 + 1010 * delta) * S2z) -
                         16 * PhiDOT * r * pow(rDOT, 4) *
                             (Complex(0, 58753) +
                              105 * (2030 + 1947 * delta) * S1z -
                              105 * (-2030 + 1947 * delta) * S2z) +
                         pow(PhiDOT, 5) * pow(r, 5) *
                             (Complex(0, -8997592) +
                              105 * (39959 + 15835 * delta) * S1z -
                              105 * (-39959 + 15835 * delta) * S2z) -
                         10 * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) *
                             (Complex(0, 254228) +
                              21 * (65293 + 34551 * delta) * S1z -
                              21 * (-65293 + 34551 * delta) * S2z) +
                         Complex(0, 2) * pow(PhiDOT, 4) * pow(r, 4) * rDOT *
                             (Complex(0, 12351083) +
                              105 * (43193 + 37913 * delta) * S1z -
                              105 * (-43193 + 37913 * delta) * S2z))) +
              mass * r *
                  (2520 * pow(Nu, 3) *
                       (8036 * pow(PhiDOT, 3) * pow(r, 3) +
                        Complex(0, 41814) * pow(PhiDOT, 2) * pow(r, 2) * rDOT -
                        30537 * PhiDOT * r * pow(rDOT, 2) -
                        Complex(0, 9064) * pow(rDOT, 3)) *
                       (S1z + S2z) -
                   210 *
                       (11849 * pow(PhiDOT, 3) * pow(r, 3) +
                        Complex(0, 31868) * pow(PhiDOT, 2) * pow(r, 2) * rDOT +
                        3572 * PhiDOT * r * pow(rDOT, 2) +
                        Complex(0, 1508) * pow(rDOT, 3)) *
                       (S1z + delta * S1z + S2z - delta * S2z) -
                   12 * pow(Nu, 2) *
                       (Complex(0, 1) * pow(PhiDOT, 2) * pow(r, 2) * rDOT *
                            (Complex(0, -1150397) +
                             35 * (154765 + 157449 * delta) * S1z +
                             5416775 * S2z - 5510715 * delta * S2z) -
                        PhiDOT * r * pow(rDOT, 2) *
                            (Complex(0, 2306552) +
                             35 * (4931 + 110733 * delta) * S1z + 172585 * S2z -
                             3875655 * delta * S2z) +
                        pow(PhiDOT, 3) * pow(r, 3) *
                            (Complex(0, 12461121) +
                             35 * (18331 + 87381 * delta) * S1z + 641585 * S2z -
                             3058335 * delta * S2z) -
                        Complex(0, 25) * pow(rDOT, 3) *
                            (Complex(0, 36676) +
                             35 * (137 + 1053 * delta) * S1z + 4795 * S2z -
                             36855 * delta * S2z)) +
                   Nu * (Complex(0, -200) * pow(rDOT, 3) *
                             (Complex(0, -2501) +
                              42 * (-308 + 51 * delta) * S1z -
                              42 * (308 + 51 * delta) * S2z) +
                         2 * PhiDOT * r * pow(rDOT, 2) *
                             (Complex(0, -1951984) -
                              105 * (-37907 + 14661 * delta) * S1z +
                              105 * (37907 + 14661 * delta) * S2z) +
                         Complex(0, 8) * pow(PhiDOT, 2) * pow(r, 2) * rDOT *
                             (Complex(0, -10005028) +
                              105 * (40991 + 31689 * delta) * S1z -
                              105 * (-40991 + 31689 * delta) * S2z) +
                         pow(PhiDOT, 3) * pow(r, 3) *
                             (Complex(0, 88418488) +
                              105 * (57793 + 266391 * delta) * S1z -
                              105 * (-57793 + 266391 * delta) * S2z))))) /
            (332640. * sqrt(35) * pow(r, 4))));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hQC_4_m_4(REAL8 Nu, UINT4 vpnorder, REAL8 x) {

  /* if(vpnorder == 2){
       return(((-16*sqrt(0.7142857142857143))/9. +
       (16*sqrt(0.7142857142857143)*Nu)/3.)*pow(x,2));
   } */

  /* else if(vpnorder == 4){
       return((4744/(99.*sqrt(35)) - (10184*sqrt(0.7142857142857143)*Nu)/297. +
       (200*sqrt(35)*pow(Nu,2))/99.)*pow(x,3));
   } */

  /* else */ if (vpnorder == 5) {
    return ((64 * sqrt(0.7142857142857143) * (-1 + 3 * Nu) * pow(x, 3.5) *
             (M_PI + Complex(0, 2) * log(2))) /
            9.);
  }

  /* else if(vpnorder == 6){
       return((-2137342/(45045.*sqrt(35)) + (2176238*Nu)/(6435.*sqrt(35)) -
       (587516*pow(Nu,2))/(1053.*sqrt(35)) +
       (452194*pow(Nu,3))/(3861.*sqrt(35)))*pow(x,4));
   } */

  else {
    return 0;
  }
}

static COMPLEX16 hl_4_m_4(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z, REAL8 x) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_4_m_4: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_4_m_4(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z) +
            hQC_4_m_4(Nu, vpnorder, x)) *
           cpolar(1, -4 * Phi);
  }
}

static COMPLEX16 hl_4_m_min4(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z, REAL8 x) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_4_m_min4: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) * pow(-1, 4) *
           conj(hGO_4_m_4(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z) +
                hQC_4_m_4(Nu, vpnorder, x)) *
           cpolar(1, 4 * Phi);
  }
}

// H43
static COMPLEX16 hGO_4_m_3(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z,
                           REAL8 x) {
  REAL8 delta = sqrt(1 - 4 * Nu);
  REAL8 kappa1 = 1.0;
  REAL8 kappa2 = 1.0;

  if (vpnorder == 3) {
    return (Complex(0, 0.16666666666666666) * delta * mass * (-1 + 2 * Nu) *
            PhiDOT *
            (4 * mass +
             r * (23 * pow(PhiDOT, 2) * pow(r, 2) +
                  Complex(0, 10) * PhiDOT * r * rDOT - 2 * pow(rDOT, 2)))) /
           (sqrt(70) * r);
  }

  else if (vpnorder == 4) {
    return ((Complex(0, -0.041666666666666664) * sqrt(0.35714285714285715) *
             pow(mass, 2) * Nu *
             (4 * mass +
              r * (23 * pow(PhiDOT, 2) * pow(r, 2) +
                   Complex(0, 10) * PhiDOT * r * rDOT - 2 * pow(rDOT, 2))) *
             ((-1 + delta) * S1z + S2z + delta * S2z)) /
                pow(r, 3) -
            (Complex(0, -1.125) * sqrt(0.35714285714285715) * Nu *
             ((-1 + delta) * S1z + S2z + delta * S2z) * pow(x, 3)) +
            (Complex(0, -1.125) * sqrt(0.35714285714285715) * Nu *
             (-S1z + S2z + delta * (S1z + S2z)) * pow(x, 3)));
  }

  else if (vpnorder == 5) {
    return (Complex(0, 0.0012626262626262627) * delta * mass * PhiDOT *
            (2 * pow(mass, 2) * (972 - 2293 * Nu + 1398 * pow(Nu, 2)) +
             2 * mass * r *
                 ((1788 - 9077 * Nu + 13416 * pow(Nu, 2)) * pow(PhiDOT, 2) *
                      pow(r, 2) +
                  Complex(0, 3) * (-2796 + 5299 * Nu + 1622 * pow(Nu, 2)) *
                      PhiDOT * r * rDOT -
                  2 * (-1200 + 2545 * Nu + 162 * pow(Nu, 2)) * pow(rDOT, 2)) -
             3 * pow(r, 2) *
                 ((-524 - 489 * Nu + 6392 * pow(Nu, 2)) * pow(PhiDOT, 4) *
                      pow(r, 4) +
                  Complex(0, 4) * (796 - 1864 * Nu + 133 * pow(Nu, 2)) *
                      pow(PhiDOT, 3) * pow(r, 3) * rDOT +
                  42 * (-51 + 94 * Nu + 56 * pow(Nu, 2)) * pow(PhiDOT, 2) *
                      pow(r, 2) * pow(rDOT, 2) +
                  Complex(0, 4) * (-229 + 366 * Nu + 358 * pow(Nu, 2)) *
                      PhiDOT * r * pow(rDOT, 3) -
                  4 * (-43 + 62 * Nu + 80 * pow(Nu, 2)) * pow(rDOT, 4)))) /
           (sqrt(70) * pow(r, 2));
  }

  else if (vpnorder == 6) {
    return (
        (delta * pow(mass, 2) * Nu * PhiDOT *
         (6 * mass * (181 * PhiDOT * r - Complex(0, 89) * rDOT) +
          r * (4847 * pow(PhiDOT, 3) * pow(r, 3) -
               Complex(0, 7338) * pow(PhiDOT, 2) * pow(r, 2) * rDOT -
               408 * PhiDOT * r * pow(rDOT, 2) +
               Complex(0, 112) * pow(rDOT, 3)))) /
            (180. * sqrt(70) *
             pow(r,
                 2)) /* Henry et al. QC spin terms */ /* +((Complex(0,0.03409090909090909)
*(220*(S1z - S2z + delta*(S1z + S2z)) + Nu*(-2403*(S1z - S2z) + 203*delta*(S1z +
S2z)) + pow(Nu,2)*(3359*(S1z - S2z) + 457*delta*(S1z +
S2z)))*pow(x,4))/sqrt(70)) */
        + /* Henry et al. ecc spin terms */ (
              (Complex(0, -0.0006313131313131314) * pow(mass, 2) *
               (2 * pow(mass, 2) *
                    ((-440 + 6801 * Nu - 1428 * pow(Nu, 2) +
                      delta * (-440 - 3193 * Nu + 300 * pow(Nu, 2))) *
                         S1z +
                     (440 - 6801 * Nu + 1428 * pow(Nu, 2) +
                      delta * (-440 - 3193 * Nu + 300 * pow(Nu, 2))) *
                         S2z) -
                2 * mass * r *
                    (Complex(0, -3) * PhiDOT * r * rDOT *
                         (delta * (-1320 + 9093 * Nu + 59 * pow(Nu, 2)) * S1z -
                          5 * (264 - 311 * Nu + 823 * pow(Nu, 2)) * S1z +
                          delta * (-1320 + 9093 * Nu + 59 * pow(Nu, 2)) * S2z +
                          5 * (264 - 311 * Nu + 823 * pow(Nu, 2)) * S2z) -
                     2 * pow(rDOT, 2) *
                         ((220 + 1659 * Nu - 1512 * pow(Nu, 2) +
                           delta * (220 - 3067 * Nu + 240 * pow(Nu, 2))) *
                              S1z +
                          (-220 - 1659 * Nu + 1512 * pow(Nu, 2) +
                           delta * (220 - 3067 * Nu + 240 * pow(Nu, 2))) *
                              S2z) +
                     2 * pow(PhiDOT, 2) * pow(r, 2) *
                         ((1826 - 19530 * Nu + 20145 * pow(Nu, 2) +
                           delta * (1826 + 1534 * Nu + 567 * pow(Nu, 2))) *
                              S1z +
                          (-1826 + 19530 * Nu - 20145 * pow(Nu, 2) +
                           delta * (1826 + 1534 * Nu + 567 * pow(Nu, 2))) *
                              S2z)) -
                3 * pow(r, 2) *
                    (Complex(0, 3080) * pow(PhiDOT, 3) * pow(r, 3) * rDOT *
                         ((1 + delta) * S1z + (-1 + delta) * S2z) +
                     2 * pow(Nu, 2) *
                         (129 * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) *
                              (41 * S1z + 23 * delta * S1z - 41 * S2z +
                               23 * delta * S2z) -
                          2 * pow(rDOT, 4) *
                              (149 * S1z + 75 * delta * S1z - 149 * S2z +
                               75 * delta * S2z) +
                          Complex(0, 2) * PhiDOT * r * pow(rDOT, 3) *
                              (925 * S1z + 491 * delta * S1z - 925 * S2z +
                               491 * delta * S2z) -
                          Complex(0, 1) * pow(PhiDOT, 3) * pow(r, 3) * rDOT *
                              (-2105 * S1z + 753 * delta * S1z + 2105 * S2z +
                               753 * delta * S2z) +
                          pow(PhiDOT, 4) * pow(r, 4) *
                              (11617 * S1z + 4847 * delta * S1z - 11617 * S2z +
                               4847 * delta * S2z)) +
                     Nu * (16 * pow(PhiDOT, 4) * pow(r, 4) *
                               (-413 * S1z + 127 * delta * S1z + 413 * S2z +
                                127 * delta * S2z) -
                           2 * pow(rDOT, 4) *
                               (-301 * S1z + 213 * delta * S1z + 301 * S2z +
                                213 * delta * S2z) +
                           Complex(0, 2) * PhiDOT * r * pow(rDOT, 3) *
                               (-1625 * S1z + 1009 * delta * S1z + 1625 * S2z +
                                1009 * delta * S2z) +
                           3 * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) *
                               (-2587 * S1z + 1267 * delta * S1z + 2587 * S2z +
                                1267 * delta * S2z) -
                           Complex(0, 2) * pow(PhiDOT, 3) * pow(r, 3) * rDOT *
                               (3635 * S1z + 7981 * delta * S1z - 3635 * S2z +
                                7981 * delta * S2z))))) /
              (sqrt(70) * pow(r, 4))));
  }

  else if (vpnorder == 7) {

    return (
        /*Henry et al. QC spin terms */ /* (Complex(0,-0.020833333333333332)*pow(x,4.5)*(54*(10*(4
  + delta)*pow(Nu,2)
  + kappa1*(5 + 5*delta*pow(1 - 2*Nu,2) - 30*Nu + 28*pow(Nu,2)))*pow(S1z,2) +
  S1z*(Complex(0,297) - Nu*(Complex(0,283) + 810*M_PI) + delta*(Complex(0,297) -
  5*Nu*(Complex(0,629) - 162*M_PI + 216*Nu*S2z)) + Complex(0,3240)*(-1 +
  delta)*Nu*atanh(1/5)) + S2z*(Complex(0,297)*(-1 + delta) + Nu*(Complex(0,283)
  + 810*M_PI) + 54*kappa2*(-5 + 5*delta*pow(1 - 2*Nu,2) + 30*Nu -
  28*pow(Nu,2))*S2z + 5*Nu*(Complex(0,-629)*delta + 162*delta*M_PI - 432*Nu*S2z
  + 108*delta*Nu*S2z + Complex(0,648)*(1 + delta)*atanh(1/5)))))/sqrt(70)
  + */
                                        /* Henry et al. ecc+spin terms */
        ((Complex(0, -1.3875013875013875e-6) * mass *
          (5005 * pow(mass, 2) *
               (-24 * PhiDOT * pow(r, 2) * pow(rDOT, 2) *
                    (Complex(0, 1) * (-11 + 48 * Nu) * S1z +
                     6 * (-5 + 4 * kappa1) * pow(Nu, 2) * pow(S1z, 2) +
                     S2z * (Complex(0, 11) - Complex(0, 48) * Nu -
                            6 * (-5 + 4 * kappa2) * pow(Nu, 2) * S2z)) +
                4 * r * pow(rDOT, 3) *
                    ((-11 + 93 * Nu) * S1z -
                     Complex(0, 6) * (-5 + 4 * kappa1) * pow(Nu, 2) *
                         pow(S1z, 2) +
                     S2z * (11 - 93 * Nu +
                            Complex(0, 6) * (-5 + 4 * kappa2) * pow(Nu, 2) *
                                S2z)) +
                4 * mass * rDOT *
                    ((22 - 111 * Nu) * S1z +
                     Complex(0, 6) * (15 + 8 * kappa1) * pow(Nu, 2) *
                         pow(S1z, 2) +
                     S2z * (-22 + 111 * Nu -
                            Complex(0, 6) * (15 + 8 * kappa2) * pow(Nu, 2) *
                                S2z)) +
                Complex(0, 6) * pow(PhiDOT, 2) * pow(r, 3) * rDOT *
                    (Complex(0, 1) * (-121 + 963 * Nu) * S1z +
                     6 *
                         (-55 * pow(Nu, 2) +
                          kappa1 * (-5 + 30 * Nu + 4 * pow(Nu, 2))) *
                         pow(S1z, 2) +
                     S2z * (Complex(0, 121) + 30 * kappa2 * S2z -
                            6 * (-55 + 4 * kappa2) * pow(Nu, 2) * S2z -
                            9 * Nu * (Complex(0, 107) + 20 * kappa2 * S2z))) +
                2 * mass * PhiDOT * r *
                    (Complex(0, -1) * (-121 + 633 * Nu) * S1z +
                     6 *
                         (180 * pow(Nu, 2) +
                          kappa1 * (9 - 54 * Nu + 28 * pow(Nu, 2))) *
                         pow(S1z, 2) -
                     S2z * (Complex(0, 121) + 54 * kappa2 * S2z +
                            24 * (45 + 7 * kappa2) * pow(Nu, 2) * S2z -
                            3 * Nu * (Complex(0, 211) + 108 * kappa2 * S2z))) +
                pow(PhiDOT, 3) * pow(r, 4) *
                    (Complex(0, -1) * (-649 + 4767 * Nu) * S1z +
                     6 *
                         (820 * pow(Nu, 2) +
                          kappa1 * (75 - 450 * Nu + 364 * pow(Nu, 2))) *
                         pow(S1z, 2) -
                     S2z *
                         (Complex(0, 649) + 450 * kappa2 * S2z +
                          24 * (205 + 91 * kappa2) * pow(Nu, 2) * S2z -
                          3 * Nu * (Complex(0, 1589) + 900 * kappa2 * S2z)))) +
           delta *
               (2 * mass * PhiDOT * pow(r, 3) *
                    (10 *
                         (234744 - 1010534 * Nu + 1024443 * pow(Nu, 2) +
                          5451096 * pow(Nu, 3)) *
                         pow(PhiDOT, 4) * pow(r, 4) -
                     Complex(0, 3) *
                         (-2426804 + 1512854 * Nu + 4994115 * pow(Nu, 2) +
                          610960 * pow(Nu, 3)) *
                         pow(PhiDOT, 3) * pow(r, 3) * rDOT +
                     (-30341028 + 23936528 * Nu + 89326545 * pow(Nu, 2) +
                      19329660 * pow(Nu, 3)) *
                         pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) +
                     Complex(0, 21) *
                         (-668008 + 803028 * Nu + 1908955 * pow(Nu, 2) +
                          540370 * pow(Nu, 3)) *
                         PhiDOT * r * pow(rDOT, 3) -
                     14 *
                         (-172143 + 155683 * Nu + 680580 * pow(Nu, 2) +
                          111840 * pow(Nu, 3)) *
                         pow(rDOT, 4)) -
                105 * PhiDOT * pow(r, 4) *
                    ((-8280 + 24681 * Nu - 151973 * pow(Nu, 2) +
                      624074 * pow(Nu, 3)) *
                         pow(PhiDOT, 6) * pow(r, 6) +
                     Complex(0, 2) *
                         (-32208 + 248485 * Nu - 524074 * pow(Nu, 2) +
                          24546 * pow(Nu, 3)) *
                         pow(PhiDOT, 5) * pow(r, 5) * rDOT +
                     2 *
                         (48924 - 239802 * Nu + 137447 * pow(Nu, 2) +
                          358156 * pow(Nu, 3)) *
                         pow(PhiDOT, 4) * pow(r, 4) * pow(rDOT, 2) +
                     Complex(0, 4) *
                         (174 + 24488 * Nu - 102039 * pow(Nu, 2) +
                          44882 * pow(Nu, 3)) *
                         pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 3) +
                     3 *
                         (10455 - 56490 * Nu + 84504 * pow(Nu, 2) +
                          54016 * pow(Nu, 3)) *
                         pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 4) +
                     Complex(0, 2) *
                         (11175 - 52698 * Nu + 57436 * pow(Nu, 2) +
                          60808 * pow(Nu, 3)) *
                         PhiDOT * r * pow(rDOT, 5) -
                     6 *
                         (829 - 3726 * Nu + 3480 * pow(Nu, 2) +
                          4640 * pow(Nu, 3)) *
                         pow(rDOT, 6)) +
                pow(mass, 2) * r *
                    (20020 * pow(rDOT, 3) * (S1z + S2z) *
                         (-11 + 71 * Nu +
                          Complex(0, 30) * pow(Nu, 2) * (S1z + S2z)) -
                     52 * PhiDOT * r * pow(rDOT, 2) *
                         (129150 * pow(Nu, 3) +
                          7 * Nu *
                              (31961 + Complex(0, 8580) * S1z +
                               Complex(0, 8580) * S2z) -
                          Complex(0, 3) *
                              (Complex(0, 32671) + 8470 * S1z + 8470 * S2z) -
                          35 * pow(Nu, 2) *
                              (33313 + 1980 * pow(S1z, 2) + 3960 * S1z * S2z +
                               1980 * pow(S2z, 2))) +
                     pow(PhiDOT, 3) * pow(r, 3) *
                         (-10566168 - 70869960 * pow(Nu, 3) +
                          Complex(0, 3248245) * S1z +
                          2252250 * kappa1 * pow(S1z, 2) +
                          Complex(0, 3248245) * S2z +
                          2252250 * kappa2 * pow(S2z, 2) -
                          7 * Nu *
                              (4818166 + Complex(0, 2480335) * S1z +
                               1287000 * kappa1 * pow(S1z, 2) +
                               Complex(0, 2480335) * S2z +
                               1287000 * kappa2 * pow(S2z, 2)) +
                          3850 * pow(Nu, 2) *
                              (38873 + 78 * (7 + 30 * kappa1) * pow(S1z, 2) -
                               3588 * S1z * S2z +
                               78 * (7 + 30 * kappa2) * pow(S2z, 2))) -
                     Complex(0, 2) * pow(PhiDOT, 2) * pow(r, 2) * rDOT *
                         (6531280 * pow(Nu, 3) +
                          3 * (5382288 + Complex(0, 605605) * S1z +
                               150150 * kappa1 * pow(S1z, 2) +
                               Complex(0, 605605) * S2z +
                               150150 * kappa2 * pow(S2z, 2)) +
                          210 * pow(Nu, 2) *
                              (322144 + 2145 * (1 + 4 * kappa1) * pow(S1z, 2) -
                               12870 * S1z * S2z +
                               2145 * (1 + 4 * kappa2) * pow(S2z, 2)) -
                          21 * Nu *
                              (2826484 + 85800 * kappa1 * pow(S1z, 2) +
                               Complex(0, 515515) * S2z +
                               85800 * kappa2 * pow(S2z, 2) -
                               5005 * S1z * (Complex(0, -103) + 60 * S2z)))) -
                26 * pow(mass, 3) *
                    (770 * rDOT *
                         (Complex(0, -90) * pow(Nu, 2) * pow(S1z, 2) +
                          S2z * (-22 + 67 * Nu -
                                 Complex(0, 90) * pow(Nu, 2) * S2z) +
                          S1z * (-22 + Nu * (67 + Complex(0, 300) * S2z) -
                                 Complex(0, 180) * pow(Nu, 2) * S2z)) +
                     PhiDOT * r *
                         (-38076 + 174720 * pow(Nu, 3) -
                          Complex(0, 46585) * S1z -
                          20790 * kappa1 * pow(S1z, 2) -
                          Complex(0, 46585) * S2z -
                          20790 * kappa2 * pow(S2z, 2) -
                          1260 * pow(Nu, 2) *
                              (158 + 33 * (5 + 2 * kappa1) * pow(S1z, 2) +
                               198 * S1z * S2z +
                               33 * (5 + 2 * kappa2) * pow(S2z, 2)) +
                          7 * Nu *
                              (6188 + 11880 * kappa1 * pow(S1z, 2) +
                               Complex(0, 21505) * S2z +
                               11880 * kappa2 * pow(S2z, 2) +
                               55 * S1z * (Complex(0, 391) + 744 * S2z))))))) /
         (sqrt(70) *
          pow(r, 4))) + /* Henry et al. QC spinning hereditary terms */
        ((27 * sqrt(0.35714285714285715) * Nu *
          ((-1 + delta) * S1z + S2z + delta * S2z) * pow(x, 4.5) *
          (Complex(0, -1) * M_PI + 4 * atanh(1 / 5))) /
         8.));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hQC_4_m_3(REAL8 Nu, UINT4 vpnorder, REAL8 x) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  /* if(vpnorder == 3){
       return(((Complex(0,-4.5)*delta)/sqrt(70) +
   (Complex(0,9)*delta*Nu)/sqrt(70))* pow(x,2.5));
   } */

  /*  else if(vpnorder == 5){
       return(((Complex(0,15.954545454545455)*delta)/sqrt(70) -
        Complex(0,6.170454545454546)*sqrt(0.7)*delta*Nu +
        (Complex(0,17.863636363636363)*delta*pow(Nu,2))/sqrt(70))*
        pow(x,3.5));
   } */

  /* else */ if (vpnorder == 6) {
    return ((Complex(0, 13.5) * delta * (-1 + 2 * Nu) * pow(x, 4) *
             (M_PI + Complex(0, 2) * log(1.5))) /
            sqrt(70));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hl_4_m_3(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z, REAL8 x) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_4_m_3: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_4_m_3(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x) +
            hQC_4_m_3(Nu, vpnorder, x)) *
           cpolar(1, -3 * Phi);
  }
}

static COMPLEX16 hl_4_m_min3(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z, REAL8 x) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_4_m_min3: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) * pow(-1, 4) *
           conj(hGO_4_m_3(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x) +
                hQC_4_m_3(Nu, vpnorder, x)) *
           cpolar(1, 3 * Phi);
  }
}

// H42
static COMPLEX16 hGO_4_m_2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z) {
  REAL8 delta = sqrt(1 - 4 * Nu);
  REAL8 kappa1 = 1.0;
  REAL8 kappa2 = 1.0;

  if (vpnorder == 2) {
    return (-(sqrt(5) * (-1 + 3 * Nu) *
              (7 * pow(mass, 2) -
               6 * pow(r, 2) * (PhiDOT * r - Complex(0, 1) * rDOT) *
                   pow(PhiDOT * r + Complex(0, 1) * rDOT, 3) +
               3 * mass * r *
                   (pow(PhiDOT, 2) * pow(r, 2) +
                    Complex(0, 9) * PhiDOT * r * rDOT - 6 * pow(rDOT, 2)))) /
            (126. * pow(r, 2)));
  }

  else if (vpnorder == 4) {
    return (-(40 * pow(mass, 3) * (314 - 987 * Nu + 195 * pow(Nu, 2)) +
              60 * (23 - 159 * Nu + 291 * pow(Nu, 2)) * pow(r, 3) *
                  pow(PhiDOT * r - Complex(0, 1) * rDOT, 2) *
                  pow(PhiDOT * r + Complex(0, 1) * rDOT, 4) +
              pow(mass, 2) * r *
                  ((1987 - 11200 * Nu + 12960 * pow(Nu, 2)) * pow(PhiDOT, 2) *
                       pow(r, 2) +
                   Complex(0, 12) * (967 - 4615 * Nu + 5935 * pow(Nu, 2)) *
                       PhiDOT * r * rDOT -
                   10 * (290 - 2033 * Nu + 4365 * pow(Nu, 2)) * pow(rDOT, 2)) -
              3 * mass * pow(r, 2) *
                  ((1577 - 7940 * Nu + 9920 * pow(Nu, 2)) * pow(PhiDOT, 4) *
                       pow(r, 4) +
                   Complex(0, 4) * (-454 - 315 * Nu + 5980 * pow(Nu, 2)) *
                       pow(PhiDOT, 3) * pow(r, 3) * rDOT -
                   2 * (549 - 2140 * Nu + 2140 * pow(Nu, 2)) * pow(PhiDOT, 2) *
                       pow(r, 2) * pow(rDOT, 2) +
                   Complex(0, 2) * (-1853 + 1730 * Nu + 13230 * pow(Nu, 2)) *
                       PhiDOT * r * pow(rDOT, 3) -
                   20 * (-83 + 30 * Nu + 762 * pow(Nu, 2)) * pow(rDOT, 4))) /
            (5544. * sqrt(5) * pow(r, 3)));
  }

  else if (vpnorder == 5) {
    return ((pow(mass, 2) * Nu *
             (mass * (Complex(0, 129) * PhiDOT * r - 54 * rDOT) +
              r * (Complex(0, -73) * pow(PhiDOT, 3) * pow(r, 3) +
                   21 * pow(PhiDOT, 2) * pow(r, 2) * rDOT -
                   Complex(0, 78) * PhiDOT * r * pow(rDOT, 2) +
                   26 * pow(rDOT, 3)))) /
                (84. * sqrt(5) * pow(r, 3))
            /* Henry et al. QC spin terms */
            /* +((4*(21*delta*Nu*(S1z - S2z) + 59*Nu*(S1z + S2z) -
      78*pow(Nu,2)*(S1z + S2z)
      - 10*(S1z + delta*S1z + S2z - delta*S2z))*pow(x,3.5))/(189.*sqrt(5))) */
            + /* Henry et al. ecc spin terms */ (
                  (pow(mass, 2) *
                   (Complex(0, 2) * mass * rDOT *
                        ((-100 - 100 * delta + 577 * Nu + 333 * delta * Nu -
                          864 * pow(Nu, 2)) *
                             S1z +
                         (-100 + 100 * delta + 577 * Nu - 333 * delta * Nu -
                          864 * pow(Nu, 2)) *
                             S2z) +
                    4 * mass * PhiDOT * r *
                        ((-70 - 70 * delta + 253 * Nu + 237 * delta * Nu -
                          156 * pow(Nu, 2)) *
                             S1z +
                         (-70 + 70 * delta + 253 * Nu - 237 * delta * Nu -
                          156 * pow(Nu, 2)) *
                             S2z) -
                    3 * r *
                        (Complex(0, 2) * pow(rDOT, 3) *
                             ((-20 - 20 * delta + 107 * Nu + 63 * delta * Nu -
                               144 * pow(Nu, 2)) *
                                  S1z +
                              (-20 + 20 * delta + 107 * Nu - 63 * delta * Nu -
                               144 * pow(Nu, 2)) *
                                  S2z) +
                         4 * PhiDOT * r * pow(rDOT, 2) *
                             ((-20 - 20 * delta + 63 * Nu + 67 * delta * Nu -
                               16 * pow(Nu, 2)) *
                                  S1z +
                              (-20 + 20 * delta + 63 * Nu - 67 * delta * Nu -
                               16 * pow(Nu, 2)) *
                                  S2z) +
                         4 * pow(PhiDOT, 3) * pow(r, 3) *
                             ((5 + 5 * delta - 7 * Nu + 2 * delta * Nu -
                               41 * pow(Nu, 2)) *
                                  S1z -
                              (-5 + 5 * delta + 7 * Nu + 2 * delta * Nu +
                               41 * pow(Nu, 2)) *
                                  S2z) -
                         Complex(0, 1) * pow(PhiDOT, 2) * pow(r, 2) * rDOT *
                             ((-10 - 10 * delta - 339 * Nu + 89 * delta * Nu +
                               1048 * pow(Nu, 2)) *
                                  S1z +
                              (-10 + 10 * delta - 339 * Nu - 89 * delta * Nu +
                               1048 * pow(Nu, 2)) *
                                  S2z)))) /
                  (504. * sqrt(5) * pow(r, 3))));
  }

  else if (vpnorder == 6) {
    return (
        (-10 * pow(mass, 4) *
             (-4477296 + 12734393 * Nu - 6895 * pow(Nu, 2) +
              1043805 * pow(Nu, 3)) +
         3150 * (-367 + 4337 * Nu - 17462 * pow(Nu, 2) + 23577 * pow(Nu, 3)) *
             pow(r, 4) * pow(PhiDOT * r - Complex(0, 1) * rDOT, 3) *
             pow(PhiDOT * r + Complex(0, 1) * rDOT, 5) -
         2 * pow(mass, 3) * r *
             (7 *
                  (-100473 + 3430399 * Nu - 9132990 * pow(Nu, 2) +
                   2885660 * pow(Nu, 3)) *
                  pow(PhiDOT, 2) * pow(r, 2) +
              Complex(0, 2) *
                  (7571073 - 10780154 * Nu - 56898800 * pow(Nu, 2) +
                   43665510 * pow(Nu, 3)) *
                  PhiDOT * r * rDOT -
              10 *
                  (1283609 - 5800627 * Nu + 3725295 * pow(Nu, 2) +
                   4771935 * pow(Nu, 3)) *
                  pow(rDOT, 2)) +
         7 * pow(mass, 2) * pow(r, 2) *
             ((1071402 + 3846989 * Nu - 27339110 * pow(Nu, 2) +
               17538420 * pow(Nu, 3)) *
                  pow(PhiDOT, 4) * pow(r, 4) +
              Complex(0, 12) *
                  (671714 - 1645932 * Nu - 1903365 * pow(Nu, 2) +
                   3346250 * pow(Nu, 3)) *
                  pow(PhiDOT, 3) * pow(r, 3) * rDOT -
              (481563 + 4291961 * Nu - 17137220 * pow(Nu, 2) +
               9315720 * pow(Nu, 3)) *
                  pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) +
              Complex(0, 12) *
                  (52044 - 1084807 * Nu + 1849450 * pow(Nu, 2) +
                   4171730 * pow(Nu, 3)) *
                  PhiDOT * r * pow(rDOT, 3) -
              5 *
                  (1083 - 1246819 * Nu + 2524240 * pow(Nu, 2) +
                   5995845 * pow(Nu, 3)) *
                  pow(rDOT, 4)) -
         105 * mass * pow(r, 3) *
             ((54272 - 58271 * Nu - 815454 * pow(Nu, 2) +
               1435572 * pow(Nu, 3)) *
                  pow(PhiDOT, 6) * pow(r, 6) +
              Complex(0, 1) *
                  (73976 - 157355 * Nu - 811766 * pow(Nu, 2) +
                   2935488 * pow(Nu, 3)) *
                  pow(PhiDOT, 5) * pow(r, 5) * rDOT +
              (72637 - 400832 * Nu + 282028 * pow(Nu, 2) +
               1063956 * pow(Nu, 3)) *
                  pow(PhiDOT, 4) * pow(r, 4) * pow(rDOT, 2) +
              Complex(0, 4) *
                  (90174 - 385167 * Nu - 126419 * pow(Nu, 2) +
                   1739072 * pow(Nu, 3)) *
                  pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 3) +
              (-55817 + 58920 * Nu + 989942 * pow(Nu, 2) -
               2334016 * pow(Nu, 3)) *
                  pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 4) +
              Complex(0, 1) *
                  (218975 - 1037408 * Nu + 148970 * pow(Nu, 2) +
                   3699480 * pow(Nu, 3)) *
                  PhiDOT * r * pow(rDOT, 5) +
              10 *
                  (-10233 + 44864 * Nu + 13050 * pow(Nu, 2) -
                   203280 * pow(Nu, 3)) *
                  pow(rDOT, 6))) /
            (5.04504e6 * sqrt(5) *
             pow(r,
                 4)) /* Henry et al. QC spin terms */ /* +((2*sqrt(5)*(kappa1*(1
+ delta - 2*Nu
+ 6*pow(Nu,2))*pow(S1z,2) + S2z*(4*(1 - 3*Nu)*Nu*S1z - kappa2*(-1 + delta + 2*Nu
- 6*pow(Nu,2))*S2z))*pow(x,4))/63.) */
        +
        /* Henry et al. ecc spin terms */ (
            (sqrt(5) * pow(mass, 3) *
             (kappa1 *
                  (mass * (12 + delta * (12 - 34 * Nu) - 58 * Nu +
                           72 * pow(Nu, 2)) +
                   r * ((5 - delta * (-5 + Nu) - 11 * Nu + 30 * pow(Nu, 2)) *
                            pow(PhiDOT, 2) * pow(r, 2) +
                        Complex(0, 1) *
                            (13 + delta * (13 - 59 * Nu) - 85 * Nu +
                             78 * pow(Nu, 2)) *
                            PhiDOT * r * rDOT +
                        4 *
                            (-2 + 11 * Nu - 12 * pow(Nu, 2) +
                             delta * (-2 + 7 * Nu)) *
                            pow(rDOT, 2))) *
                  pow(S1z, 2) +
              S2z *
                  (-2 * mass *
                       (-24 * Nu * S1z + 6 * (-1 + delta) * kappa2 * S2z +
                        (29 - 17 * delta) * kappa2 * Nu * S2z +
                        36 * pow(Nu, 2) * (2 * S1z - kappa2 * S2z)) +
                   r * (-((-1 + delta) * kappa2 *
                          (5 * pow(PhiDOT, 2) * pow(r, 2) +
                           Complex(0, 13) * PhiDOT * r * rDOT -
                           8 * pow(rDOT, 2)) *
                          S2z) -
                        6 * pow(Nu, 2) *
                            (5 * pow(PhiDOT, 2) * pow(r, 2) +
                             Complex(0, 13) * PhiDOT * r * rDOT -
                             8 * pow(rDOT, 2)) *
                            (2 * S1z - kappa2 * S2z) +
                        Nu * (pow(PhiDOT, 2) * pow(r, 2) *
                                  (20 * S1z + (-11 + delta) * kappa2 * S2z) -
                              4 * pow(rDOT, 2) *
                                  (8 * S1z + (-11 + 7 * delta) * kappa2 * S2z) +
                              Complex(0, 1) * PhiDOT * r * rDOT *
                                  (52 * S1z +
                                   (-85 + 59 * delta) * kappa2 * S2z)))))) /
            (168. * pow(r, 4))));
  }

  else if (vpnorder == 7) {

    return (
        /* Henry et al. QC spin terms */ /* ((-5958*delta*pow(Nu,2)*(S1z - S2z)
 + 45*delta*Nu*(-S1z + S2z)
  - 6009*Nu*(S1z + S2z) + 7528*pow(Nu,2)*(S1z + S2z) - 1134*pow(Nu,3)*(S1z +
 S2z) + 1104*(S1z + delta*S1z + S2z - delta*S2z))*pow(x,4.5))/(2079.*sqrt(5))
 + */
        /* Henry et al. ecc+spin terms */ (
            (pow(mass, 2) *
             (-14 * pow(mass, 2) *
                  (60 * pow(Nu, 3) *
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
                   3 * pow(Nu, 2) *
                       (4 * PhiDOT * r *
                            (Complex(0, 108687) + 41225 * S1z +
                             29825 * delta * S1z + 41225 * S2z -
                             29825 * delta * S2z) +
                        Complex(0, 5) * rDOT *
                            (Complex(0, 19722) + 132241 * S1z +
                             18825 * delta * S1z + 132241 * S2z -
                             18825 * delta * S2z))) -
              3 * pow(r, 2) *
                  (840 * pow(Nu, 3) *
                       (3159 * pow(PhiDOT, 5) * pow(r, 5) +
                        Complex(0, 18450) * pow(PhiDOT, 4) * pow(r, 4) * rDOT +
                        3678 * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) +
                        Complex(0, 13970) * pow(PhiDOT, 2) * pow(r, 2) *
                            pow(rDOT, 3) -
                        1546 * PhiDOT * r * pow(rDOT, 4) +
                        Complex(0, 2100) * pow(rDOT, 5)) *
                       (S1z + S2z) -
                   525 *
                       (354 * pow(PhiDOT, 5) * pow(r, 5) +
                        Complex(0, 453) * pow(PhiDOT, 4) * pow(r, 4) * rDOT -
                        60 * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) +
                        Complex(0, 1546) * pow(PhiDOT, 2) * pow(r, 2) *
                            pow(rDOT, 3) -
                        336 * PhiDOT * r * pow(rDOT, 4) +
                        Complex(0, 184) * pow(rDOT, 5)) *
                       (S1z + delta * S1z + S2z - delta * S2z) +
                   Nu * (Complex(0, 7) * pow(PhiDOT, 4) * pow(r, 4) * rDOT *
                             (Complex(0, 1062898) +
                              15 * (11275 + 23003 * delta) * S1z +
                              (169125 - 345045 * delta) * S2z) -
                         8 * PhiDOT * r * pow(rDOT, 4) *
                             (Complex(0, -58753) +
                              420 * (145 + 6 * delta) * S1z -
                              420 * (-145 + 6 * delta) * S2z) +
                         Complex(0, 40) * pow(rDOT, 5) *
                             (Complex(0, 12203) +
                              21 * (843 + 656 * delta) * S1z -
                              21 * (-843 + 656 * delta) * S2z) -
                         5 * pow(PhiDOT, 5) * pow(r, 5) *
                             (Complex(0, 240964) +
                              21 * (-20361 + 2747 * delta) * S1z -
                              21 * (20361 + 2747 * delta) * S2z) +
                         2 * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) *
                             (Complex(0, 1832842) +
                              105 * (-13567 + 6603 * delta) * S1z -
                              105 * (13567 + 6603 * delta) * S2z) +
                         Complex(0, 4) * pow(PhiDOT, 2) * pow(r, 2) *
                             pow(rDOT, 3) *
                             (Complex(0, 97999) +
                              105 * (9121 + 7144 * delta) * S1z -
                              105 * (-9121 + 7144 * delta) * S2z)) +
                   6 * pow(Nu, 2) *
                       (Complex(0, -10) * pow(rDOT, 5) *
                            (Complex(0, 35291) +
                             14 * (2183 + 1155 * delta) * S1z +
                             (30562 - 16170 * delta) * S2z) -
                        Complex(0, 3) * pow(PhiDOT, 2) * pow(r, 2) *
                            pow(rDOT, 3) *
                            (Complex(0, 120747) +
                             280 * (1571 + 165 * delta) * S1z -
                             280 * (-1571 + 165 * delta) * S2z) +
                        2 * PhiDOT * r * pow(rDOT, 4) *
                            (Complex(0, -280067) -
                             140 * (-124 + 1293 * delta) * S1z +
                             140 * (124 + 1293 * delta) * S2z) +
                        5 * pow(PhiDOT, 5) * pow(r, 5) *
                            (Complex(0, 45393) +
                             14 * (-12899 + 1573 * delta) * S1z -
                             14 * (12899 + 1573 * delta) * S2z) -
                        Complex(0, 7) * pow(PhiDOT, 4) * pow(r, 4) * rDOT *
                            (Complex(0, 213613) +
                             5 * (31897 + 4025 * delta) * S1z -
                             5 * (-31897 + 4025 * delta) * S2z) +
                        pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) *
                            (Complex(0, -2032699) -
                             35 * (-33139 + 11613 * delta) * S1z +
                             35 * (33139 + 11613 * delta) * S2z))) +
              mass * r *
                  (2520 * pow(Nu, 3) *
                       (3772 * pow(PhiDOT, 3) * pow(r, 3) +
                        Complex(0, 22365) * pow(PhiDOT, 2) * pow(r, 2) * rDOT -
                        1427 * PhiDOT * r * pow(rDOT, 2) +
                        Complex(0, 4532) * pow(rDOT, 3)) *
                       (S1z + S2z) +
                   420 *
                       (643 * pow(PhiDOT, 3) * pow(r, 3) -
                        Complex(0, 17143) * pow(PhiDOT, 2) * pow(r, 2) * rDOT +
                        4684 * PhiDOT * r * pow(rDOT, 2) +
                        Complex(0, 377) * pow(rDOT, 3)) *
                       (S1z + delta * S1z + S2z - delta * S2z) -
                   6 * pow(Nu, 2) *
                       (Complex(0, 25) * pow(rDOT, 3) *
                            (Complex(0, 73352) +
                             35 * (137 + 1053 * delta) * S1z +
                             (4795 - 36855 * delta) * S2z) +
                        pow(PhiDOT, 3) * pow(r, 3) *
                            (Complex(0, 3625463) +
                             70 * (45113 + 1515 * delta) * S1z -
                             70 * (-45113 + 1515 * delta) * S2z) +
                        4 * PhiDOT * r * pow(rDOT, 2) *
                            (Complex(0, 576638) +
                             245 * (-884 + 2103 * delta) * S1z -
                             245 * (884 + 2103 * delta) * S2z) -
                        Complex(0, 1) * pow(PhiDOT, 2) * pow(r, 2) * rDOT *
                            (Complex(0, -3147986) +
                             35 * (-290963 + 34473 * delta) * S1z -
                             35 * (290963 + 34473 * delta) * S2z)) +
                   Nu * (Complex(0, 200) * pow(rDOT, 3) *
                             (Complex(0, -2501) +
                              21 * (-308 + 51 * delta) * S1z -
                              21 * (308 + 51 * delta) * S2z) +
                         5 * pow(PhiDOT, 3) * pow(r, 3) *
                             (Complex(0, 1501652) -
                              21 * (-37535 + 8031 * delta) * S1z +
                              21 * (37535 + 8031 * delta) * S2z) -
                         2 * PhiDOT * r * pow(rDOT, 2) *
                             (Complex(0, -975992) +
                              105 * (33593 + 24921 * delta) * S1z -
                              105 * (-33593 + 24921 * delta) * S2z) +
                         Complex(0, 4) * pow(PhiDOT, 2) * pow(r, 2) * rDOT *
                             (Complex(0, 4843364) +
                              105 * (88681 + 40215 * delta) * S1z -
                              105 * (-88681 + 40215 * delta) * S2z))))) /
            (1.16424e6 * sqrt(5) * pow(r, 4))));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hQC_4_m_2(REAL8 Nu, UINT4 vpnorder, REAL8 x) {

  /* if(vpnorder == 2){
       return(((2*sqrt(5))/63. - (2*sqrt(5)*Nu)/21.)*pow(x,2));
   }

   else if(vpnorder == 4){
       return((-437/(693.*sqrt(5)) + (115*sqrt(5)*Nu)/297. -
       (19*sqrt(5)*pow(Nu,2))/693.)*pow(x,3));
   } */

  /* else */ if (vpnorder == 5) {
    return (M_PI * ((4 * sqrt(5) * pow(x, 3.5)) / 63. -
                    (4 * sqrt(5) * Nu * pow(x, 3.5)) / 21.));
  }
  /* else if(vpnorder == 6){
       return((346013/(420420.*sqrt(5)) - (606751*Nu)/(180180.*sqrt(5)) +
       (400453*pow(Nu,2))/(162162.*sqrt(5)) +
       (25783*pow(Nu,3))/(108108.*sqrt(5)))*pow(x,4));
   } */
  else {
    return 0;
  }
}

static COMPLEX16 hl_4_m_2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z, REAL8 x) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_5_m_2: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_4_m_2(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z) +
            hQC_4_m_2(Nu, vpnorder, x)) *
           cpolar(1, -2 * Phi);
  }
}

static COMPLEX16 hl_4_m_min2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z, REAL8 x) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_5_m_min2: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) * pow(-1, 4) *
           conj(hGO_4_m_2(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z) +
                hQC_4_m_2(Nu, vpnorder, x)) *
           cpolar(1, 2 * Phi);
  }
}

// H41
static COMPLEX16 hGO_4_m_1(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z,
                           REAL8 x) {
  REAL8 delta = sqrt(1 - 4 * Nu);
  REAL8 kappa1 = 1.0;
  REAL8 kappa2 = 1.0;

  if (vpnorder == 3) {
    return (delta * mass * (-1 + 2 * Nu) * PhiDOT *
            (Complex(0, -12) * mass +
             r * (Complex(0, 11) * pow(PhiDOT, 2) * pow(r, 2) +
                  10 * PhiDOT * r * rDOT + Complex(0, 6) * pow(rDOT, 2)))) /
           (42. * sqrt(10) * r);
  }

  else if (vpnorder == 4) {
    return ((Complex(0, 0.005952380952380952) * sqrt(2.5) * pow(mass, 2) * Nu *
             (12 * mass +
              r * (-11 * pow(PhiDOT, 2) * pow(r, 2) +
                   Complex(0, 10) * PhiDOT * r * rDOT - 6 * pow(rDOT, 2))) *
             ((-1 + delta) * S1z + S2z + delta * S2z)) /
                pow(r, 3) -
            (Complex(0, 0.005952380952380952) * sqrt(2.5) * Nu *
             (-S1z + S2z + delta * (S1z + S2z)) * pow(x, 3)) +
            (Complex(0, 0.005952380952380952) * sqrt(2.5) * Nu *
             (-S1z + S2z + delta * (S1z + S2z)) * pow(x, 3)));
  }

  else if (vpnorder == 5) {
    return (Complex(0, -0.00018037518037518038) * delta * mass * PhiDOT *
            (6 * pow(mass, 2) * (972 - 2293 * Nu + 1398 * pow(Nu, 2)) -
             2 * mass * r *
                 ((-340 - 785 * Nu + 6504 * pow(Nu, 2)) * pow(PhiDOT, 2) *
                      pow(r, 2) -
                  Complex(0, 3) * (-2796 + 5299 * Nu + 1622 * pow(Nu, 2)) *
                      PhiDOT * r * rDOT +
                  6 * (-1200 + 2545 * Nu + 162 * pow(Nu, 2)) * pow(rDOT, 2)) +
             3 * pow(r, 2) *
                 ((-540 + 235 * Nu + 2648 * pow(Nu, 2)) * pow(PhiDOT, 4) *
                      pow(r, 4) -
                  Complex(0, 4) * (-1764 + 3536 * Nu + 373 * pow(Nu, 2)) *
                      pow(PhiDOT, 3) * pow(r, 3) * rDOT +
                  2 * (-723 + 1022 * Nu + 1384 * pow(Nu, 2)) * pow(PhiDOT, 2) *
                      pow(r, 2) * pow(rDOT, 2) -
                  Complex(0, 4) * (-229 + 366 * Nu + 358 * pow(Nu, 2)) *
                      PhiDOT * r * pow(rDOT, 3) +
                  12 * (-43 + 62 * Nu + 80 * pow(Nu, 2)) * pow(rDOT, 4)))) /
           (sqrt(10) * pow(r, 2));
  }

  else if (vpnorder == 6) {
    return (
        -(delta * pow(mass, 2) * Nu * PhiDOT *
          (mass * (362 * PhiDOT * r - Complex(0, 534) * rDOT) +
           r * (149 * pow(PhiDOT, 3) * pow(r, 3) +
                Complex(0, 182) * pow(PhiDOT, 2) * pow(r, 2) * rDOT -
                136 * PhiDOT * r * pow(rDOT, 2) +
                Complex(0, 112) * pow(rDOT, 3)))) /
            (420. * sqrt(10) * pow(r, 2))
        /* Henry et al. QC spin terms */
        /* +((Complex(0,-0.00018037518037518038)*(220*(S1z - S2z + delta*(S1z +
  S2z))
  + Nu*(-2247*(S1z - S2z) + 47*delta*(S1z + S2z)) +
  pow(Nu,2)*(2891*(S1z - S2z) + 613*delta*(S1z + S2z)))*pow(x,4))/sqrt(10)) */
        + /*Henry et al. ecc spin terms*/ (
              (Complex(0, 0.00027056277056277056) * pow(mass, 2) *
               (2 * pow(mass, 2) *
                    ((-440 + 6801 * Nu - 1428 * pow(Nu, 2) +
                      delta * (-440 - 3193 * Nu + 300 * pow(Nu, 2))) *
                         S1z +
                     (440 - 6801 * Nu + 1428 * pow(Nu, 2) +
                      delta * (-440 - 3193 * Nu + 300 * pow(Nu, 2))) *
                         S2z) +
                2 * mass * r *
                    (Complex(0, 1) * PhiDOT * r * rDOT *
                         (delta * (-1320 + 9093 * Nu + 59 * pow(Nu, 2)) * S1z -
                          5 * (264 - 311 * Nu + 823 * pow(Nu, 2)) * S1z +
                          delta * (-1320 + 9093 * Nu + 59 * pow(Nu, 2)) * S2z +
                          5 * (264 - 311 * Nu + 823 * pow(Nu, 2)) * S2z) +
                     2 * pow(rDOT, 2) *
                         ((220 + 1659 * Nu - 1512 * pow(Nu, 2) +
                           delta * (220 - 3067 * Nu + 240 * pow(Nu, 2))) *
                              S1z +
                          (-220 - 1659 * Nu + 1512 * pow(Nu, 2) +
                           delta * (220 - 3067 * Nu + 240 * pow(Nu, 2))) *
                              S2z) +
                     2 * pow(PhiDOT, 2) * pow(r, 2) *
                         ((-66 + 158 * Nu - 3389 * pow(Nu, 2) +
                           delta * (-66 + 238 * Nu + 301 * pow(Nu, 2))) *
                              S1z +
                          (66 - 158 * Nu + 3389 * pow(Nu, 2) +
                           delta * (-66 + 238 * Nu + 301 * pow(Nu, 2))) *
                              S2z)) +
                pow(r, 2) *
                    (Complex(0, 3960) * pow(PhiDOT, 3) * pow(r, 3) * rDOT *
                         ((1 + delta) * S1z + (-1 + delta) * S2z) +
                     2 * pow(Nu, 2) *
                         (6 * pow(rDOT, 4) *
                              (149 * S1z + 75 * delta * S1z - 149 * S2z +
                               75 * delta * S2z) +
                          7 * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) *
                              (331 * S1z + 293 * delta * S1z - 331 * S2z +
                               293 * delta * S2z) -
                          Complex(0, 2) * PhiDOT * r * pow(rDOT, 3) *
                              (925 * S1z + 491 * delta * S1z - 925 * S2z +
                               491 * delta * S2z) +
                          pow(PhiDOT, 4) * pow(r, 4) *
                              (2593 * S1z + 551 * delta * S1z - 2593 * S2z +
                               551 * delta * S2z) +
                          Complex(0, 1) * pow(PhiDOT, 3) * pow(r, 3) * rDOT *
                              (8815 * S1z + 2393 * delta * S1z - 8815 * S2z +
                               2393 * delta * S2z)) +
                     Nu * (6 * pow(rDOT, 4) *
                               (-301 * S1z + 213 * delta * S1z + 301 * S2z +
                                213 * delta * S2z) +
                           8 * pow(PhiDOT, 4) * pow(r, 4) *
                               (-349 * S1z + 305 * delta * S1z + 349 * S2z +
                                305 * delta * S2z) -
                           Complex(0, 2) * PhiDOT * r * pow(rDOT, 3) *
                               (-1625 * S1z + 1009 * delta * S1z + 1625 * S2z +
                                1009 * delta * S2z) +
                           pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) *
                               (-4973 * S1z + 2773 * delta * S1z + 4973 * S2z +
                                2773 * delta * S2z) -
                           Complex(0, 2) * pow(PhiDOT, 3) * pow(r, 3) * rDOT *
                               (1765 * S1z + 14779 * delta * S1z - 1765 * S2z +
                                14779 * delta * S2z))))) /
              (sqrt(10) * pow(r, 4))));
  }

  else if (vpnorder == 7) {

    return (
        /* Henry et al. QC spiun terms */ /* (Complex(0,0.000992063492063492)*pow(x,4.5)*(6*(10*(4
 + delta)*pow(Nu,2)
   + kappa1*(5 + 5*delta*pow(1 - 2*Nu,2) + 6*Nu*(-5 + 2*Nu)))*pow(S1z,2) +
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
             (14 * mass * PhiDOT * pow(r, 3) *
                  (2 *
                       (243784 - 839354 * Nu + 840205 * pow(Nu, 2) +
                        1464600 * pow(Nu, 3)) *
                       pow(PhiDOT, 4) * pow(r, 4) -
                   Complex(0, 3) *
                       (-14468 - 3185682 * Nu + 8381955 * pow(Nu, 2) +
                        30000 * pow(Nu, 3)) *
                       pow(PhiDOT, 3) * pow(r, 3) * rDOT +
                   5 *
                       (-657932 + 496144 * Nu + 2150479 * pow(Nu, 2) +
                        542628 * pow(Nu, 3)) *
                       pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) -
                   Complex(0, 3) *
                       (-668008 + 803028 * Nu + 1908955 * pow(Nu, 2) +
                        540370 * pow(Nu, 3)) *
                       PhiDOT * r * pow(rDOT, 3) +
                   6 *
                       (-172143 + 155683 * Nu + 680580 * pow(Nu, 2) +
                        111840 * pow(Nu, 3)) *
                       pow(rDOT, 4)) -
              105 * PhiDOT * pow(r, 4) *
                  ((1320 - 35155 * Nu + 36743 * pow(Nu, 2) +
                    241538 * pow(Nu, 3)) *
                       pow(PhiDOT, 6) * pow(r, 6) +
                   Complex(0, 26) *
                       (-7536 + 50167 * Nu - 84782 * pow(Nu, 2) +
                        22 * pow(Nu, 3)) *
                       pow(PhiDOT, 5) * pow(r, 5) * rDOT +
                   2 *
                       (46908 - 203266 * Nu + 59723 * pow(Nu, 2) +
                        272204 * pow(Nu, 3)) *
                       pow(PhiDOT, 4) * pow(r, 4) * pow(rDOT, 2) -
                   Complex(0, 4) *
                       (38798 - 217656 * Nu + 324761 * pow(Nu, 2) +
                        67186 * pow(Nu, 3)) *
                       pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 3) +
                   (58177 - 260454 * Nu + 239688 * pow(Nu, 2) +
                    324160 * pow(Nu, 3)) *
                       pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 4) -
                   Complex(0, 2) *
                       (11175 - 52698 * Nu + 57436 * pow(Nu, 2) +
                        60808 * pow(Nu, 3)) *
                       PhiDOT * r * pow(rDOT, 5) +
                   18 *
                       (829 - 3726 * Nu + 3480 * pow(Nu, 2) +
                        4640 * pow(Nu, 3)) *
                       pow(rDOT, 6)) +
              26 * pow(mass, 3) *
                  (2310 * rDOT *
                       (Complex(0, -30) * pow(Nu, 2) * pow(S1z, 2) +
                        S2z * (-22 + 67 * Nu -
                               Complex(0, 30) * pow(Nu, 2) * S2z) +
                        S1z * (-22 + Nu * (67 + Complex(0, 100) * S2z) -
                               Complex(0, 60) * pow(Nu, 2) * S2z)) +
                   PhiDOT * r *
                       (-114228 + 524160 * pow(Nu, 3) -
                        Complex(0, 46585) * S1z - 62370 * kappa1 * pow(S1z, 2) -
                        Complex(0, 46585) * S2z - 62370 * kappa2 * pow(S2z, 2) -
                        420 * pow(Nu, 2) *
                            (1422 + 11 * (-25 + 54 * kappa1) * pow(S1z, 2) -
                             1738 * S1z * S2z +
                             11 * (-25 + 54 * kappa2) * pow(S2z, 2)) +
                        7 * Nu *
                            (18564 + 35640 * kappa1 * pow(S1z, 2) +
                             Complex(0, 21505) * S2z +
                             35640 * kappa2 * pow(S2z, 2) +
                             55 * S1z * (Complex(0, 391) + 312 * S2z)))) +
              pow(mass, 2) * r *
                  (60060 * pow(rDOT, 3) * (S1z + S2z) *
                       (11 - 71 * Nu -
                        Complex(0, 10) * pow(Nu, 2) * (S1z + S2z)) +
                   156 * PhiDOT * r * pow(rDOT, 2) *
                       (98013 + 129150 * pow(Nu, 3) - Complex(0, 8470) * S1z +
                        7 * Nu *
                            (31961 + Complex(0, 2860) * S1z +
                             Complex(0, 2860) * S2z) -
                        Complex(0, 8470) * S2z +
                        35 * pow(Nu, 2) *
                            (-33313 + 220 * pow(S1z, 2) + 440 * S1z * S2z +
                             220 * pow(S2z, 2))) -
                   5 * pow(PhiDOT, 3) * pow(r, 3) *
                       (3966680 + 5749576 * pow(Nu, 3) -
                        Complex(0, 231231) * S1z -
                        90090 * kappa1 * pow(S1z, 2) -
                        Complex(0, 231231) * S2z -
                        90090 * kappa2 * pow(S2z, 2) +
                        7 * Nu *
                            (-921242 + Complex(0, 110253) * S1z +
                             51480 * kappa1 * pow(S1z, 2) +
                             Complex(0, 110253) * S2z +
                             51480 * kappa2 * pow(S2z, 2)) -
                        462 * pow(Nu, 2) *
                            (10483 + 130 * (-5 + 6 * kappa1) * pow(S1z, 2) -
                             2860 * S1z * S2z +
                             130 * (-5 + 6 * kappa2) * pow(S2z, 2))) +
                   Complex(0, 2) * pow(PhiDOT, 2) * pow(r, 2) * rDOT *
                       (6531280 * pow(Nu, 3) +
                        3 * (5382288 - Complex(0, 385385) * S1z +
                             150150 * kappa1 * pow(S1z, 2) -
                             Complex(0, 385385) * S2z +
                             150150 * kappa2 * pow(S2z, 2)) +
                        210 * pow(Nu, 2) *
                            (322144 + 715 * (19 + 12 * kappa1) * pow(S1z, 2) +
                             10010 * S1z * S2z +
                             715 * (19 + 12 * kappa2) * pow(S2z, 2)) -
                        21 * Nu *
                            (2826484 + 85800 * kappa1 * pow(S1z, 2) -
                             Complex(0, 140855) * S2z +
                             85800 * kappa2 * pow(S2z, 2) +
                             715 * S1z * (Complex(0, -197) + 540 * S2z)))))) +
           5005 * pow(mass, 2) *
               (-3 * r *
                    (-8 * PhiDOT * r * pow(rDOT, 2) *
                         (Complex(0, -1) * (-11 + 48 * Nu) * S1z +
                          2 * (-5 + 16 * kappa1) * pow(Nu, 2) * pow(S1z, 2) +
                          S2z * (Complex(0, -11) + Complex(0, 48) * Nu +
                                 2 * (5 - 16 * kappa2) * pow(Nu, 2) * S2z)) +
                     4 * pow(rDOT, 3) *
                         ((11 - 93 * Nu) * S1z +
                          Complex(0, 2) * (-5 + 4 * kappa1) * pow(Nu, 2) *
                              pow(S1z, 2) +
                          S2z * (-11 + 93 * Nu -
                                 Complex(0, 2) * (-5 + 4 * kappa2) *
                                     pow(Nu, 2) * S2z)) +
                     Complex(0, 2) * pow(PhiDOT, 2) * pow(r, 2) * rDOT *
                         (Complex(0, 1) * (-77 + 351 * Nu) * S1z +
                          2 * (245 * pow(Nu, 2) + kappa1 * (15 - 90 * Nu + 68 * pow(Nu, 2))) *
                              pow(S1z, 2) -
                          S2z * (Complex(0, -77) + 30 * kappa2 * S2z +
                                 2 * (245 + 68 * kappa2) * pow(Nu, 2) * S2z -
                                 9 * Nu *
                                     (Complex(0, -39) + 20 * kappa2 * S2z))) +
                     pow(PhiDOT, 3) * pow(r, 3) *
                         (Complex(0, -1) * (-77 + 411 * Nu) * S1z +
                          2 * (100 * pow(Nu, 2) + kappa1 * (15 - 90 * Nu + 172 * pow(Nu, 2))) *
                              pow(S1z, 2) -
                          S2z * (Complex(0, 77) + 30 * kappa2 * S2z +
                                 8 * (25 + 43 * kappa2) * pow(Nu, 2) * S2z -
                                 3 * Nu *
                                     (Complex(0, 137) + 60 * kappa2 * S2z)))) +
                2 * mass *
                    (6 * rDOT *
                         ((22 - 111 * Nu) * S1z +
                          Complex(0, 2) * (15 + 8 * kappa1) * pow(Nu, 2) *
                              pow(S1z, 2) +
                          S2z * (-22 + 111 * Nu -
                                 Complex(0, 2) * (15 + 8 * kappa2) *
                                     pow(Nu, 2) * S2z)) +
                     PhiDOT * r *
                         (Complex(0, -1) * (-121 + 633 * Nu) * S1z +
                          6 * (220 * pow(Nu, 2) + 3 * kappa1 * (9 - 54 * Nu + 76 * pow(Nu, 2))) *
                              pow(S1z, 2) -
                          S2z * (Complex(0, 121) + 162 * kappa2 * S2z +
                                 24 * (55 + 57 * kappa2) * pow(Nu, 2) * S2z -
                                 3 * Nu *
                                     (Complex(0, 211) +
                                      324 * kappa2 * S2z))))))) /
         (sqrt(10) * pow(r, 4))) +
        /* Henry et al. QC spinning hereditary terms */ (
            (sqrt(2.5) * Nu * ((-1 + delta) * S1z + S2z + delta * S2z) *
             pow(x, 4.5) * (Complex(0, 1) * M_PI + log(4))) /
            168.));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hQC_4_m_1(REAL8 Nu, UINT4 vpnorder, REAL8 x) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  /* if(vpnorder == 3){
       return(((Complex(0,0.023809523809523808)*delta)/sqrt(10) -
       (Complex(0,0.047619047619047616)*delta*Nu)/sqrt(10))*pow(x,2.5));
   }
   else if(vpnorder == 5){
       return(((Complex(0,-0.07287157287157287)*delta)/sqrt(10) +
       (Complex(0,0.18235930735930736)*delta*Nu)/sqrt(10) -
       (Complex(0,0.05988455988455989)*delta*pow(Nu,2))/sqrt(10))*
       pow(x,3.5));
   } */
  /* else */ if (vpnorder == 6) {
    return (
        M_PI *
            ((Complex(0, 0.023809523809523808) * delta * pow(x, 4)) / sqrt(10) -
             (Complex(0, 0.047619047619047616) * delta * Nu * pow(x, 4)) /
                 sqrt(10)) +
        ((delta * pow(x, 4)) / (21. * sqrt(10)) -
         (sqrt(0.4) * delta * Nu * pow(x, 4)) / 21.) *
            log(2));
  } else {
    return 0;
  }
}

static COMPLEX16 hl_4_m_1(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z, REAL8 x) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_4_m_1: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_4_m_1(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x) +
            hQC_4_m_1(Nu, vpnorder, x)) *
           cpolar(1, -1 * Phi);
  }
}

static COMPLEX16 hl_4_m_min1(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z, REAL8 x) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_4_m_min1: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) * pow(-1, 4) *
           conj(hGO_4_m_1(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x) +
                hQC_4_m_1(Nu, vpnorder, x)) *
           cpolar(1, 1 * Phi);
  }
}

// H55
static COMPLEX16 hGO_5_m_5(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z) {
  REAL8 delta = sqrt(1 - 4 * Nu);
  REAL8 kappa1 = 1.0;
  REAL8 kappa2 = 1.0;

  if (vpnorder == 3) {
    return ((delta * (-1 + 2 * Nu) *
             (24 * pow(r, 2) * pow(Complex(0, -1) * PhiDOT * r + rDOT, 5) +
              2 * pow(mass, 2) * (Complex(0, -86) * PhiDOT * r + 41 * rDOT) +
              3 * mass * r *
                  (Complex(0, -143) * pow(PhiDOT, 3) * pow(r, 3) +
                   208 * pow(PhiDOT, 2) * pow(r, 2) * rDOT +
                   Complex(0, 132) * PhiDOT * r * pow(rDOT, 2) -
                   32 * pow(rDOT, 3)))) /
            (48. * sqrt(66) * pow(r, 2)));
  }

  else if (vpnorder == 5) {
    return (
        (delta *
         (360 * (33 - 197 * Nu + 294 * pow(Nu, 2)) * pow(r, 3) *
              pow(PhiDOT * r + Complex(0, 1) * rDOT, 6) *
              (Complex(0, 1) * PhiDOT * r + rDOT) +
          2 * pow(mass, 3) *
              (Complex(0, -5) * (53311 - 121906 * Nu + 42816 * pow(Nu, 2)) *
                   PhiDOT * r +
               78 * (1141 - 2760 * Nu + 1420 * pow(Nu, 2)) * rDOT) +
          2 * pow(mass, 2) * r *
              (Complex(0, -10) * (40826 - 125981 * Nu + 87534 * pow(Nu, 2)) *
                   pow(PhiDOT, 3) * pow(r, 3) +
               (112966 - 818425 * Nu + 1385970 * pow(Nu, 2)) * pow(PhiDOT, 2) *
                   pow(r, 2) * rDOT +
               Complex(0, 20) * (-2636 - 11335 * Nu + 43962 * pow(Nu, 2)) *
                   PhiDOT * r * pow(rDOT, 2) -
               39 * (-639 - 735 * Nu + 5690 * pow(Nu, 2)) * pow(rDOT, 3)) +
          15 * mass * pow(r, 2) *
              (Complex(0, 4) * (139 - 1687 * Nu + 11376 * pow(Nu, 2)) *
                   pow(PhiDOT, 5) * pow(r, 5) -
               2 * (11276 - 19559 * Nu + 5982 * pow(Nu, 2)) * pow(PhiDOT, 4) *
                   pow(r, 4) * rDOT +
               Complex(0, 3) * (-14615 + 12440 * Nu + 37132 * pow(Nu, 2)) *
                   pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) -
               8 * (-4666 + 139 * Nu + 22194 * pow(Nu, 2)) * pow(PhiDOT, 2) *
                   pow(r, 2) * pow(rDOT, 3) -
               Complex(0, 4) * (-3971 - 3226 * Nu + 27360 * pow(Nu, 2)) *
                   PhiDOT * r * pow(rDOT, 4) +
               48 * (-57 - 97 * Nu + 518 * pow(Nu, 2)) * pow(rDOT, 5)))) /
        (18720. * sqrt(66) * pow(r, 3)));
  }

  else if (vpnorder == 6) {
    return (
        -(delta * pow(mass, 2) * Nu *
          (3566 * pow(mass, 2) +
           6 * mass * r *
               (11305 * pow(PhiDOT, 2) * pow(r, 2) +
                Complex(0, 3921) * PhiDOT * r * rDOT - 906 * pow(rDOT, 2)) +
           pow(r, 2) * (104681 * pow(PhiDOT, 4) * pow(r, 4) +
                        Complex(0, 17192) * pow(PhiDOT, 3) * pow(r, 3) * rDOT -
                        27840 * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) -
                        Complex(0, 9968) * PhiDOT * r * pow(rDOT, 3) +
                        1424 * pow(rDOT, 4)))) /
            (1680. * sqrt(66) * pow(r, 4))
        /* +((Complex(0,-21.70138888888889)*((1 + delta + 3*delta*(-1 + Nu)*Nu +
  Nu*(-7 + 11*Nu))*S1z
  + (-1 + delta + (7 - 11*Nu)*Nu + 3*delta*(-1 +
  Nu)*Nu)*S2z)*pow(x,4))/sqrt(66)) */
        +
        /* Henry et al. ecc spin terms */ (
            (Complex(0, -0.006944444444444444) * pow(mass, 2) *
             (pow(mass, 2) * ((90 - 622 * Nu + 966 * pow(Nu, 2) +
                               delta * (90 - 524 * Nu + 770 * pow(Nu, 2))) *
                                  S1z +
                              2 *
                                  (-45 + 311 * Nu - 483 * pow(Nu, 2) +
                                   delta * (45 - 262 * Nu + 385 * pow(Nu, 2))) *
                                  S2z) +
              3 * pow(r, 2) *
                  (4 * pow(rDOT, 4) *
                       ((15 - 92 * Nu + 126 * pow(Nu, 2) +
                         delta * (15 - 64 * Nu + 70 * pow(Nu, 2))) *
                            S1z +
                        (-15 + 92 * Nu - 126 * pow(Nu, 2) +
                         delta * (15 - 64 * Nu + 70 * pow(Nu, 2))) *
                            S2z) -
                   15 * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) *
                       ((41 - 258 * Nu + 364 * pow(Nu, 2) +
                         delta * (41 - 154 * Nu + 156 * pow(Nu, 2))) *
                            S1z +
                        (-41 + 258 * Nu - 364 * pow(Nu, 2) +
                         delta * (41 - 154 * Nu + 156 * pow(Nu, 2))) *
                            S2z) -
                   Complex(0, 4) * PhiDOT * r * pow(rDOT, 3) *
                       ((75 - 464 * Nu + 642 * pow(Nu, 2) +
                         delta * (75 - 298 * Nu + 310 * pow(Nu, 2))) *
                            S1z +
                        (-75 + 464 * Nu - 642 * pow(Nu, 2) +
                         delta * (75 - 298 * Nu + 310 * pow(Nu, 2))) *
                            S2z) +
                   pow(PhiDOT, 4) * pow(r, 4) *
                       ((300 - 2347 * Nu + 4041 * pow(Nu, 2) +
                         delta * (300 - 869 * Nu + 1085 * pow(Nu, 2))) *
                            S1z +
                        (-300 + 2347 * Nu - 4041 * pow(Nu, 2) +
                         delta * (300 - 869 * Nu + 1085 * pow(Nu, 2))) *
                            S2z) +
                   Complex(0, 1) * pow(PhiDOT, 3) * pow(r, 3) * rDOT *
                       ((675 - 4414 * Nu + 6492 * pow(Nu, 2) +
                         delta * (675 - 2558 * Nu + 2780 * pow(Nu, 2))) *
                            S1z +
                        (-675 + 4414 * Nu - 6492 * pow(Nu, 2) +
                         delta * (675 - 2558 * Nu + 2780 * pow(Nu, 2))) *
                            S2z)) +
              mass * r *
                  (-2 * pow(rDOT, 2) *
                       ((255 - 1592 * Nu + 2226 * pow(Nu, 2) +
                         delta * (255 - 1144 * Nu + 1330 * pow(Nu, 2))) *
                            S1z +
                        (-255 + 1592 * Nu - 2226 * pow(Nu, 2) +
                         delta * (255 - 1144 * Nu + 1330 * pow(Nu, 2))) *
                            S2z) +
                   Complex(0, 2) * PhiDOT * r * rDOT *
                       ((870 - 5563 * Nu + 7989 * pow(Nu, 2) +
                         delta * (870 - 3743 * Nu + 4349 * pow(Nu, 2))) *
                            S1z +
                        (-870 + 5563 * Nu - 7989 * pow(Nu, 2) +
                         delta * (870 - 3743 * Nu + 4349 * pow(Nu, 2))) *
                            S2z) +
                   pow(PhiDOT, 2) * pow(r, 2) *
                       ((1329 - 9376 * Nu + 14838 * pow(Nu, 2) +
                         delta * (1329 - 5438 * Nu + 6962 * pow(Nu, 2))) *
                            S1z +
                        (-1329 + 9376 * Nu - 14838 * pow(Nu, 2) +
                         delta * (1329 - 5438 * Nu + 6962 * pow(Nu, 2))) *
                            S2z)))) /
            (sqrt(66) * pow(r, 4))));
  }

  else if (vpnorder == 7) {

    return (/* (Complex(0,16.276041666666668)*(-1 + 2*Nu)*(kappa1*(-1 - delta +
      2*(2 + delta)*Nu)*pow(S1z,2)
      + S2z*(-4*delta*Nu*S1z + kappa2*(1 - delta + 2*(-2 +
      delta)*Nu)*S2z))*pow(x,4.5))/sqrt(66) */
            /*  + */ /* Henry et al. ecc+spin terms */ (
                (Complex(0, 24570) * pow(mass, 3) *
                     (1 - 6 * Nu + 8 * pow(Nu, 2)) *
                     (mass * (298 * PhiDOT * r + Complex(0, 94) * rDOT) +
                      5 * r *
                          (95 * pow(PhiDOT, 3) * pow(r, 3) +
                           Complex(0, 58) * pow(PhiDOT, 2) * pow(r, 2) * rDOT -
                           34 * PhiDOT * r * pow(rDOT, 2) -
                           Complex(0, 8) * pow(rDOT, 3))) *
                     (kappa1 * pow(S1z, 2) - kappa2 * pow(S2z, 2)) +
                 delta *
                     (7560 *
                          (-135 + 1467 * Nu - 5285 * pow(Nu, 2) +
                           6230 * pow(Nu, 3)) *
                          pow(r, 4) *
                          pow(Complex(0, 1) * PhiDOT * r - rDOT, 7) *
                          pow(PhiDOT * r - Complex(0, 1) * rDOT, 2) +
                      6 * pow(mass, 2) * pow(r, 2) *
                          (Complex(0, 2) *
                               (-2032864 - 4559381 * Nu - 3520590 * pow(Nu, 2) +
                                36950400 * pow(Nu, 3)) *
                               pow(PhiDOT, 5) * pow(r, 5) -
                           8 *
                               (-1581827 + 6211490 * Nu + 607860 * pow(Nu, 2) +
                                6123006 * pow(Nu, 3)) *
                               pow(PhiDOT, 4) * pow(r, 4) * rDOT +
                           Complex(0, 12) *
                               (738669 - 8183419 * Nu + 7682580 * pow(Nu, 2) +
                                9432160 * pow(Nu, 3)) *
                               pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) -
                           2 *
                               (323812 - 40147318 * Nu + 43264365 * pow(Nu, 2) +
                                116543718 * pow(Nu, 3)) *
                               pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 3) -
                           Complex(0, 14) *
                               (-19909 - 2880941 * Nu + 2929365 * pow(Nu, 2) +
                                11172930 * pow(Nu, 3)) *
                               PhiDOT * r * pow(rDOT, 4) +
                           21 *
                               (1285 - 417939 * Nu + 376375 * pow(Nu, 2) +
                                1794742 * pow(Nu, 3)) *
                               pow(rDOT, 5)) +
                      45 * mass * pow(r, 3) *
                          (Complex(0, -4) *
                               (-76218 + 328011 * Nu - 511853 * pow(Nu, 2) +
                                804238 * pow(Nu, 3)) *
                               pow(PhiDOT, 7) * pow(r, 7) -
                           4 *
                               (241824 - 1017903 * Nu + 1022762 * pow(Nu, 2) +
                                918896 * pow(Nu, 3)) *
                               pow(PhiDOT, 6) * pow(r, 6) * rDOT -
                           Complex(0, 4) *
                               (450537 - 1503933 * Nu - 642815 * pow(Nu, 2) +
                                6184642 * pow(Nu, 3)) *
                               pow(PhiDOT, 5) * pow(r, 5) * pow(rDOT, 2) +
                           8 *
                               (53226 + 197748 * Nu - 1673279 * pow(Nu, 2) +
                                3045358 * pow(Nu, 3)) *
                               pow(PhiDOT, 4) * pow(r, 4) * pow(rDOT, 3) -
                           Complex(0, 1) *
                               (1855503 - 9943746 * Nu + 10573448 * pow(Nu, 2) +
                                6270464 * pow(Nu, 3)) *
                               pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 4) +
                           16 *
                               (147339 - 604776 * Nu + 17530 * pow(Nu, 2) +
                                1666996 * pow(Nu, 3)) *
                               pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 5) +
                           Complex(0, 28) *
                               (41703 - 148362 * Nu - 104120 * pow(Nu, 2) +
                                627040 * pow(Nu, 3)) *
                               PhiDOT * r * pow(rDOT, 6) -
                           1344 *
                               (162 - 516 * Nu - 707 * pow(Nu, 2) +
                                2870 * pow(Nu, 3)) *
                               pow(rDOT, 7)) +
                      28 * pow(mass, 4) *
                          (-3 * rDOT *
                               (300238 - 535772 * pow(Nu, 3) +
                                27495 * kappa1 * pow(S1z, 2) +
                                27495 * kappa2 * pow(S2z, 2) +
                                180 * pow(Nu, 2) *
                                    (3335 + 611 * kappa1 * pow(S1z, 2) -
                                     1222 * S1z * S2z +
                                     611 * kappa2 * pow(S2z, 2)) -
                                15 * Nu *
                                    (38797 + 7332 * kappa1 * pow(S1z, 2) -
                                     7332 * S1z * S2z +
                                     7332 * kappa2 * pow(S2z, 2))) -
                           Complex(0, 1) * PhiDOT * r *
                               (-5713898 + 2860180 * pow(Nu, 3) -
                                261495 * kappa1 * pow(S1z, 2) -
                                261495 * kappa2 * pow(S2z, 2) -
                                30 * pow(Nu, 2) *
                                    (278563 + 34866 * kappa1 * pow(S1z, 2) -
                                     69732 * S1z * S2z +
                                     34866 * kappa2 * pow(S2z, 2)) +
                                12 * Nu *
                                    (1137419 + 87165 * kappa1 * pow(S1z, 2) -
                                     87165 * S1z * S2z +
                                     87165 * kappa2 * pow(S2z, 2)))) +
                      6 * pow(mass, 3) * r *
                          (28 * pow(rDOT, 3) *
                               (-175718 - 773906 * pow(Nu, 3) +
                                5850 * kappa1 * pow(S1z, 2) +
                                5850 * kappa2 * pow(S2z, 2) +
                                45 * pow(Nu, 2) *
                                    (-11217 + 520 * kappa1 * pow(S1z, 2) -
                                     1040 * S1z * S2z +
                                     520 * kappa2 * pow(S2z, 2)) -
                                3 * Nu *
                                    (-243517 + 7800 * kappa1 * pow(S1z, 2) -
                                     7800 * S1z * S2z +
                                     7800 * kappa2 * pow(S2z, 2))) +
                           Complex(0, 14) * PhiDOT * r * pow(rDOT, 2) *
                               (1974194 + 5947390 * pow(Nu, 3) -
                                49725 * kappa1 * pow(S1z, 2) -
                                49725 * kappa2 * pow(S2z, 2) -
                                1105 * pow(Nu, 2) *
                                    (-307 + 180 * kappa1 * pow(S1z, 2) -
                                     360 * S1z * S2z +
                                     180 * kappa2 * pow(S2z, 2)) +
                                Nu * (-5423729 + 198900 * kappa1 * pow(S1z, 2) -
                                      198900 * S1z * S2z +
                                      198900 * kappa2 * pow(S2z, 2))) +
                           2 * pow(PhiDOT, 2) * pow(r, 2) * rDOT *
                               (22105852 + 67721336 * pow(Nu, 3) -
                                593775 * kappa1 * pow(S1z, 2) -
                                593775 * kappa2 * pow(S2z, 2) -
                                60 * pow(Nu, 2) *
                                    (1466303 + 39585 * kappa1 * pow(S1z, 2) -
                                     79170 * S1z * S2z +
                                     39585 * kappa2 * pow(S2z, 2)) +
                                4 * Nu *
                                    (-4950311 + 593775 * kappa1 * pow(S1z, 2) -
                                     593775 * S1z * S2z +
                                     593775 * kappa2 * pow(S2z, 2))) -
                           Complex(0, 1) * pow(PhiDOT, 3) * pow(r, 3) *
                               (82632040 * pow(Nu, 3) -
                                3 * (4583428 + 648375 * kappa1 * pow(S1z, 2) +
                                     648375 * kappa2 * pow(S2z, 2)) -
                                10 * pow(Nu, 2) *
                                    (24739703 + 778050 * kappa1 * pow(S1z, 2) -
                                     1556100 * S1z * S2z +
                                     778050 * kappa2 * pow(S2z, 2)) +
                                9 * Nu *
                                    (14885821 + 864500 * kappa1 * pow(S1z, 2) -
                                     864500 * S1z * S2z +
                                     864500 * kappa2 * pow(S2z, 2)))))) /
                (1.57248e6 * sqrt(66) * pow(r, 4))));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hQC_5_m_5(REAL8 Nu, UINT4 vpnorder, REAL8 x) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  /* if(vpnorder == 3){
       return(((Complex(0,13.020833333333334)*delta)/sqrt(66) -
       (Complex(0,26.041666666666668)*delta*Nu)/sqrt(66))*pow(x,2.5));
   }
   else if(vpnorder == 5){
       return(((Complex(0,-87.80715811965813)*delta)/sqrt(66) +
       (Complex(0,229.7008547008547)*delta*Nu)/sqrt(66) -
       Complex(0,42.73504273504273)*sqrt(0.06060606060606061)*delta*
       pow(Nu,2))*pow(x,3.5));
   } */
  /* else */ if (vpnorder == 6) {
    return ((3125 * delta * (-1 + 2 * Nu) * pow(x, 4) *
             (Complex(0, -1) * M_PI + log(6.25))) /
            (48. * sqrt(66)));
  } else {
    return 0;
  }
}

static COMPLEX16 hl_5_m_5(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z, REAL8 x) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_5_m_5: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_5_m_5(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z) +
            hQC_5_m_5(Nu, vpnorder, x)) *
           cpolar(1, -5 * Phi);
  }
}

static COMPLEX16 hl_5_m_min5(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z, REAL8 x) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_5_m_min5: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) * pow(-1, 5) *
           conj(hGO_5_m_5(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z) +
                hQC_5_m_5(Nu, vpnorder, x)) *
           cpolar(1, 5 * Phi);
  }
}

// H54
static COMPLEX16 hGO_5_m_4(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  if (vpnorder == 4) {
    return (-(mass * (1 - 5 * Nu + 5 * pow(Nu, 2)) * PhiDOT *
              (mass * (82 * PhiDOT * r + Complex(0, 22) * rDOT) +
               3 * r *
                   (58 * pow(PhiDOT, 3) * pow(r, 3) +
                    Complex(0, 33) * pow(PhiDOT, 2) * pow(r, 2) * rDOT -
                    12 * PhiDOT * r * pow(rDOT, 2) -
                    Complex(0, 2) * pow(rDOT, 3)))) /
            (36. * sqrt(165) * r));
  }

  else if (vpnorder == 5) {

    return (/* (32*Nu*((-1 + delta + 2*Nu)*S1z - (1 + delta -
      2*Nu)*S2z)*pow(x,3.5))/(3.*sqrt(165))
      + */
            /* Henry et al. ecc spin terms */ (
                (pow(mass, 2) * Nu *
                 (82 * mass * PhiDOT * r + 174 * pow(PhiDOT, 3) * pow(r, 4) +
                  Complex(0, 22) * mass * rDOT +
                  Complex(0, 99) * pow(PhiDOT, 2) * pow(r, 3) * rDOT -
                  36 * PhiDOT * pow(r, 2) * pow(rDOT, 2) -
                  Complex(0, 6) * r * pow(rDOT, 3)) *
                 ((-1 + delta + 2 * Nu) * S1z - (1 + delta - 2 * Nu) * S2z)) /
                (24. * sqrt(165) * pow(r, 3))));
  }

  else if (vpnorder == 6) {
    return (mass * PhiDOT *
            (-4 * pow(mass, 2) *
                 (26 *
                      (-5051 + 28623 * Nu - 46305 * pow(Nu, 2) +
                       29470 * pow(Nu, 3)) *
                      PhiDOT * r +
                  Complex(0, 7) *
                      (5684 - 26697 * Nu + 14225 * pow(Nu, 2) +
                       25355 * pow(Nu, 3)) *
                      rDOT) -
             4 * mass * r *
                 ((-149157 + 1133006 * Nu - 2731750 * pow(Nu, 2) +
                   2085685 * pow(Nu, 3)) *
                      pow(PhiDOT, 3) * pow(r, 3) +
                  Complex(0, 7) *
                      (101118 - 491779 * Nu + 402185 * pow(Nu, 2) +
                       172105 * pow(Nu, 3)) *
                      pow(PhiDOT, 2) * pow(r, 2) * rDOT -
                  (363357 - 1825691 * Nu + 1714720 * pow(Nu, 2) +
                   395255 * pow(Nu, 3)) *
                      PhiDOT * r * pow(rDOT, 2) -
                  Complex(0, 7) *
                      (9717 - 48896 * Nu + 46000 * pow(Nu, 2) +
                       10865 * pow(Nu, 3)) *
                      pow(rDOT, 3)) +
             15 * pow(r, 2) *
                 (4 *
                      (3449 - 6580 * Nu - 56728 * pow(Nu, 2) +
                       115269 * pow(Nu, 3)) *
                      pow(PhiDOT, 5) * pow(r, 5) +
                  Complex(0, 7) *
                      (-8128 + 45859 * Nu - 62702 * pow(Nu, 2) +
                       13996 * pow(Nu, 3)) *
                      pow(PhiDOT, 4) * pow(r, 4) * rDOT +
                  4 *
                      (10125 - 47852 * Nu + 26635 * pow(Nu, 2) +
                       44240 * pow(Nu, 3)) *
                      pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) +
                  Complex(0, 182) *
                      (127 - 548 * Nu + 73 * pow(Nu, 2) + 816 * pow(Nu, 3)) *
                      pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 3) -
                  8 *
                      (1009 - 4060 * Nu - 889 * pow(Nu, 2) +
                       7952 * pow(Nu, 3)) *
                      PhiDOT * r * pow(rDOT, 4) -
                  Complex(0, 28) *
                      (45 - 172 * Nu - 85 * pow(Nu, 2) + 400 * pow(Nu, 3)) *
                      pow(rDOT, 5)))) /
           (65520. * sqrt(165) * pow(r, 2));
  }

  else if (vpnorder == 7) {

    return (/* (16*(-530*pow(Nu,3)*(S1z + S2z) + 104*(S1z + delta*S1z + S2z -
     delta*S2z)
      + 2*pow(Nu,2)*(541*delta*(S1z - S2z) - 120*(S1z + S2z)) -
     Nu*(1139*delta*(S1z - S2z) + 109*(S1z + S2z)))*pow(x,4.5))/(117.*sqrt(165))
     + */
            /* Henry et al. ecc+spin terms */
            ((pow(mass, 2) *
              (4 * pow(mass, 2) *
                   (210 * pow(Nu, 3) *
                        (1950 * PhiDOT * r + Complex(0, 1373) * rDOT) *
                        (S1z + S2z) +
                    10920 * (20 * PhiDOT * r + Complex(0, 7) * rDOT) *
                        (S1z + delta * S1z + S2z - delta * S2z) +
                    21 * pow(Nu, 2) *
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
               3 * pow(r, 2) *
                   (-840 * pow(Nu, 3) *
                        (13072 * pow(PhiDOT, 5) * pow(r, 5) +
                         Complex(0, 679) * pow(PhiDOT, 4) * pow(r, 4) * rDOT +
                         8700 * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) +
                         Complex(0, 4459) * pow(PhiDOT, 2) * pow(r, 2) *
                             pow(rDOT, 3) -
                         1432 * PhiDOT * r * pow(rDOT, 4) -
                         Complex(0, 210) * pow(rDOT, 5)) *
                        (S1z + S2z) +
                    Complex(0, 859950) * pow(PhiDOT, 4) * pow(r, 4) * rDOT *
                        (S1z + delta * S1z + S2z - delta * S2z) +
                    14 * pow(Nu, 2) *
                        (30 * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) *
                             (Complex(0, -3796) + (77 + 4765 * delta) * S1z +
                              (77 - 4765 * delta) * S2z) -
                         Complex(0, 210) * pow(rDOT, 5) *
                             ((-23 + 5 * delta) * S1z -
                              (23 + 5 * delta) * S2z) +
                         Complex(0, 13) * pow(PhiDOT, 2) * pow(r, 2) *
                             pow(rDOT, 3) *
                             (Complex(0, 5296) +
                              105 * (-29 + 39 * delta) * S1z -
                              105 * (29 + 39 * delta) * S2z) +
                         Complex(0, 2) * pow(PhiDOT, 4) * pow(r, 4) * rDOT *
                             (Complex(0, -1058096) +
                              105 * (4821 + 131 * delta) * S1z -
                              105 * (-4821 + 131 * delta) * S2z) -
                         4 * PhiDOT * r * pow(rDOT, 4) *
                             (Complex(0, 3796) +
                              15 * (-399 + 193 * delta) * S1z -
                              15 * (399 + 193 * delta) * S2z) +
                         4 * pow(PhiDOT, 5) * pow(r, 5) *
                             (Complex(0, 387296) +
                              15 * (-623 + 4413 * delta) * S1z -
                              15 * (623 + 4413 * delta) * S2z)) -
                    3 * Nu *
                        (Complex(0, -140) * pow(rDOT, 5) *
                             ((-141 + 193 * delta) * S1z -
                              (141 + 193 * delta) * S2z) +
                         Complex(0, 182) * pow(PhiDOT, 2) * pow(r, 2) *
                             pow(rDOT, 3) *
                             (Complex(0, 632) + 5 * (-319 + 555 * delta) * S1z -
                              5 * (319 + 555 * delta) * S2z) -
                         8 * PhiDOT * r * pow(rDOT, 4) *
                             (Complex(0, 2574) +
                              35 * (-415 + 623 * delta) * S1z -
                              35 * (415 + 623 * delta) * S2z) +
                         2 * pow(PhiDOT, 5) * pow(r, 5) *
                             (Complex(0, 1106001) +
                              140 * (-1069 + 3045 * delta) * S1z -
                              140 * (1069 + 3045 * delta) * S2z) +
                         20 * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) *
                             (Complex(0, -8502) +
                              7 * (-2883 + 6419 * delta) * S1z -
                              7 * (2883 + 6419 * delta) * S2z) +
                         Complex(0, 7) * pow(PhiDOT, 4) * pow(r, 4) * rDOT *
                             (Complex(0, -431288) +
                              5 * (72097 + 7957 * delta) * S1z -
                              5 * (-72097 + 7957 * delta) * S2z))) -
               6 * mass * r *
                   (420 * pow(Nu, 3) *
                        (2698 * pow(PhiDOT, 3) * pow(r, 3) +
                         Complex(0, 307) * pow(PhiDOT, 2) * pow(r, 2) * rDOT +
                         739 * PhiDOT * r * pow(rDOT, 2) +
                         Complex(0, 183) * pow(rDOT, 3)) *
                        (S1z + S2z) -
                    5460 *
                        (86 * pow(PhiDOT, 3) * pow(r, 3) +
                         Complex(0, 137) * pow(PhiDOT, 2) * pow(r, 2) * rDOT -
                         30 * PhiDOT * r * pow(rDOT, 2) -
                         Complex(0, 3) * pow(rDOT, 3)) *
                        (S1z + delta * S1z + S2z - delta * S2z) +
                    Nu * (Complex(0, 700) * pow(rDOT, 3) *
                              ((-651 + 482 * delta) * S1z -
                               (651 + 482 * delta) * S2z) +
                          28 * pow(PhiDOT, 2) * pow(r, 2) * rDOT *
                              (34294 -
                               Complex(0, 5) * (-56523 + 1013 * delta) * S1z +
                               Complex(0, 5) * (56523 + 1013 * delta) * S2z) +
                          pow(PhiDOT, 3) * pow(r, 3) *
                              (Complex(0, 1423669) +
                               280 * (2589 + 18458 * delta) * S1z -
                               280 * (-2589 + 18458 * delta) * S2z) +
                          2 * PhiDOT * r * pow(rDOT, 2) *
                              (Complex(0, 34203) +
                               35 * (-40527 + 19571 * delta) * S1z -
                               35 * (40527 + 19571 * delta) * S2z)) -
                    14 * pow(Nu, 2) *
                        (Complex(0, 50) * pow(rDOT, 3) *
                             ((-1330 + 223 * delta) * S1z -
                              (1330 + 223 * delta) * S2z) +
                         2 * PhiDOT * r * pow(rDOT, 2) *
                             (Complex(0, 7956) +
                              5 * (-39389 + 2528 * delta) * S1z -
                              5 * (39389 + 2528 * delta) * S2z) +
                         Complex(0, 1) * pow(PhiDOT, 2) * pow(r, 2) * rDOT *
                             (Complex(0, -213408) +
                              5 * (202025 + 36547 * delta) * S1z -
                              5 * (-202025 + 36547 * delta) * S2z) +
                         2 * pow(PhiDOT, 3) * pow(r, 3) *
                             (Complex(0, 160212) +
                              5 * (-2831 + 38849 * delta) * S1z -
                              5 * (2831 + 38849 * delta) * S2z))))) /
             (393120. * sqrt(165) * pow(r, 4))));
  }

  else {
    return 0;
  }
}

/* static COMPLEX16 hQC_5_m_4(REAL8 Nu, UINT4 vpnorder, REAL8 x){

    if(vpnorder == 4){
        return((-64/(9.*sqrt(165)) + (64*sqrt(0.15151515151515152)*Nu)/9. -
        (64*sqrt(0.15151515151515152)*pow(Nu,2))/9.)*pow(x,3));
    }
    else if(vpnorder == 6){
        return((142432/(4095.*sqrt(165)) -
(10528*sqrt(0.7333333333333333)*Nu)/585. + (33344*pow(Nu,2))/(117.*sqrt(165)) -
        (3616*pow(Nu,3))/(39.*sqrt(165)))*pow(x,4));
    }
    else{
        return 0;
    }
} */

static COMPLEX16 hl_5_m_4(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_5_m_4: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_5_m_4(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z,
                      S2z) /* +hQC_5_m_4(Nu,vpnorder,x) */) *
           cpolar(1, -4 * Phi);
  }
}

static COMPLEX16 hl_5_m_min4(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_5_m_min4: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) * pow(-1, 5) *
           conj(hGO_5_m_4(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z,
                          S2z) /* +hQC_5_m_4(Nu,vpnorder,x) */) *
           cpolar(1, 4 * Phi);
  }
}

// H53
static COMPLEX16 hGO_5_m_3(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z) {
  REAL8 delta = sqrt(1 - 4 * Nu);
  REAL8 kappa1 = 1.0;
  REAL8 kappa2 = 1.0;

  if (vpnorder == 3) {
    return (delta * (-1 + 2 * Nu) *
            (2 * pow(mass, 2) * (Complex(0, 258) * PhiDOT * r - 205 * rDOT) -
             Complex(0, 120) * pow(r, 2) * (PhiDOT * r - Complex(0, 1) * rDOT) *
                 pow(PhiDOT * r + Complex(0, 1) * rDOT, 4) +
             3 * mass * r *
                 (Complex(0, -51) * pow(PhiDOT, 3) * pow(r, 3) -
                  240 * pow(PhiDOT, 2) * pow(r, 2) * rDOT -
                  Complex(0, 396) * PhiDOT * r * pow(rDOT, 2) +
                  160 * pow(rDOT, 3)))) /
           (144. * sqrt(330) * pow(r, 2));
  }

  else if (vpnorder == 5) {
    return (delta *
            (Complex(0, 120) * (33 - 197 * Nu + 294 * pow(Nu, 2)) * pow(r, 3) *
                 pow(PhiDOT * r - Complex(0, 1) * rDOT, 2) *
                 pow(PhiDOT * r + Complex(0, 1) * rDOT, 5) +
             2 * pow(mass, 3) *
                 (Complex(0, 1) * (53311 - 121906 * Nu + 42816 * pow(Nu, 2)) *
                      PhiDOT * r -
                  26 * (1141 - 2760 * Nu + 1420 * pow(Nu, 2)) * rDOT) +
             2 * pow(mass, 2) * r *
                 (Complex(0, -2) * (6350 - 12803 * Nu + 10314 * pow(Nu, 2)) *
                      pow(PhiDOT, 3) * pow(r, 3) +
                  (-6546 + 64131 * Nu - 109702 * pow(Nu, 2)) * pow(PhiDOT, 2) *
                      pow(r, 2) * rDOT -
                  Complex(0, 4) * (-2636 - 11335 * Nu + 43962 * pow(Nu, 2)) *
                      PhiDOT * r * pow(rDOT, 2) +
                  13 * (-639 - 735 * Nu + 5690 * pow(Nu, 2)) * pow(rDOT, 3)) +
             3 * mass * pow(r, 2) *
                 (Complex(0, -4) * (1223 - 4567 * Nu + 3396 * pow(Nu, 2)) *
                      pow(PhiDOT, 5) * pow(r, 5) +
                  26 * (-412 - 437 * Nu + 3114 * pow(Nu, 2)) * pow(PhiDOT, 4) *
                      pow(r, 4) * rDOT +
                  Complex(0, 1) * (-2331 - 31496 * Nu + 93276 * pow(Nu, 2)) *
                      pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) +
                  8 * (-1994 + 1541 * Nu + 5718 * pow(Nu, 2)) * pow(PhiDOT, 2) *
                      pow(r, 2) * pow(rDOT, 3) +
                  Complex(0, 4) * (-3971 - 3226 * Nu + 27360 * pow(Nu, 2)) *
                      PhiDOT * r * pow(rDOT, 4) -
                  80 * (-57 - 97 * Nu + 518 * pow(Nu, 2)) * pow(rDOT, 5)))) /
           (3744. * sqrt(330) * pow(r, 3));
  }

  else if (vpnorder == 6) {
    return (
        (delta * pow(mass, 2) * Nu *
         (17830 * pow(mass, 2) +
          18 * mass * r *
              (5231 * pow(PhiDOT, 2) * pow(r, 2) +
               Complex(0, 3921) * PhiDOT * r * rDOT - 1510 * pow(rDOT, 2)) -
          pow(r, 2) * (48579 * pow(PhiDOT, 4) * pow(r, 4) +
                       Complex(0, 31304) * pow(PhiDOT, 3) * pow(r, 3) * rDOT +
                       33024 * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) +
                       Complex(0, 29904) * PhiDOT * r * pow(rDOT, 3) -
                       7120 * pow(rDOT, 4)))) /
            (5040. * sqrt(330) * pow(r, 4))
        /* +(Complex(0,0.1875)*sqrt(0.02727272727272727)*((5 + Nu*(-27 + 31*Nu)
   + delta*(5 + 23*(-1 + Nu)*Nu))*S1z
   + (-5 + (27 - 31*Nu)*Nu + delta*(5 + 23*(-1 + Nu)*Nu))*S2z)*pow(x,4)) */
        +
        /* Henry et al. ecc spin terms */ (
            (Complex(0, 0.006944444444444444) * pow(mass, 2) *
             (pow(mass, 2) * ((90 - 622 * Nu + 966 * pow(Nu, 2) +
                               delta * (90 - 524 * Nu + 770 * pow(Nu, 2))) *
                                  S1z +
                              2 *
                                  (-45 + 311 * Nu - 483 * pow(Nu, 2) +
                                   delta * (45 - 262 * Nu + 385 * pow(Nu, 2))) *
                                  S2z) -
              3 * pow(r, 2) *
                  (Complex(0, 4) * PhiDOT * r * pow(rDOT, 3) *
                       (5 * (9 - 56 * Nu + 78 * pow(Nu, 2)) * S1z +
                        delta * (45 - 142 * Nu + 114 * pow(Nu, 2)) * S1z -
                        5 * (9 - 56 * Nu + 78 * pow(Nu, 2)) * S2z +
                        delta * (45 - 142 * Nu + 114 * pow(Nu, 2)) * S2z) +
                   Complex(0, 1) * pow(PhiDOT, 3) * pow(r, 3) * rDOT *
                       (25 * (3 - 14 * Nu + 12 * pow(Nu, 2)) * S1z +
                        delta * (75 - 958 * Nu + 1516 * pow(Nu, 2)) * S1z -
                        25 * (3 - 14 * Nu + 12 * pow(Nu, 2)) * S2z +
                        delta * (75 - 958 * Nu + 1516 * pow(Nu, 2)) * S2z) -
                   4 * pow(rDOT, 4) *
                       ((15 - 92 * Nu + 126 * pow(Nu, 2) +
                         delta * (15 - 64 * Nu + 70 * pow(Nu, 2))) *
                            S1z +
                        (-15 + 92 * Nu - 126 * pow(Nu, 2) +
                         delta * (15 - 64 * Nu + 70 * pow(Nu, 2))) *
                            S2z) +
                   pow(PhiDOT, 4) * pow(r, 4) *
                       ((20 - 317 * Nu + 751 * pow(Nu, 2) +
                         delta * (20 + 13 * Nu + 91 * pow(Nu, 2))) *
                            S1z +
                        (-20 + 317 * Nu - 751 * pow(Nu, 2) +
                         delta * (20 + 13 * Nu + 91 * pow(Nu, 2))) *
                            S2z) -
                   3 * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) *
                       ((-45 + 298 * Nu - 444 * pow(Nu, 2) +
                         delta * (-45 + 2 * Nu + 148 * pow(Nu, 2))) *
                            S1z +
                        (45 - 298 * Nu + 444 * pow(Nu, 2) +
                         delta * (-45 + 2 * Nu + 148 * pow(Nu, 2))) *
                            S2z)) +
              mass * r *
                  (Complex(0, 6) * PhiDOT * r * rDOT *
                       (5 * (38 - 243 * Nu + 349 * pow(Nu, 2)) * S1z +
                        delta * (190 - 643 * Nu + 601 * pow(Nu, 2)) * S1z -
                        5 * (38 - 243 * Nu + 349 * pow(Nu, 2)) * S2z +
                        delta * (190 - 643 * Nu + 601 * pow(Nu, 2)) * S2z) +
                   pow(PhiDOT, 2) * pow(r, 2) *
                       ((849 - 5360 * Nu + 7590 * pow(Nu, 2) +
                         delta * (849 - 1774 * Nu + 418 * pow(Nu, 2))) *
                            S1z +
                        (-849 + 5360 * Nu - 7590 * pow(Nu, 2) +
                         delta * (849 - 1774 * Nu + 418 * pow(Nu, 2))) *
                            S2z) -
                   2 * pow(rDOT, 2) *
                       ((255 - 1592 * Nu + 2226 * pow(Nu, 2) +
                         delta * (255 - 1144 * Nu + 1330 * pow(Nu, 2))) *
                            S1z +
                        (-255 + 1592 * Nu - 2226 * pow(Nu, 2) +
                         delta * (255 - 1144 * Nu + 1330 * pow(Nu, 2))) *
                            S2z)))) /
            (sqrt(330) * pow(r, 4))));
  }

  else if (vpnorder == 7) {

    return (/* Complex(0,-0.140625)*sqrt(0.02727272727272727)*(kappa1*(5 +
      5*delta*pow(1 - 2*Nu,2)
      - 30*Nu + 8*pow(Nu,2))*pow(S1z,2) + S2z*(20*delta*(1 - 2*Nu)*Nu*S1z
      + kappa2*(-5 + 5*delta*pow(1 - 2*Nu,2) + 30*Nu -
      8*pow(Nu,2))*S2z))*pow(x,4.5)
      + */
            /* Henry et al. ecc+spin terms */
            ((8190 * pow(mass, 3) *
                  (Complex(0, -2) * mass *
                       (447 - 2682 * Nu + 3224 * pow(Nu, 2)) * PhiDOT * r +
                   2 * mass * (235 + 2 * Nu * (-705 + 1004 * Nu)) * rDOT +
                   r * (Complex(0, 1) * (15 - 90 * Nu + 2008 * pow(Nu, 2)) *
                            pow(PhiDOT, 3) * pow(r, 3) +
                        66 * (5 + 6 * Nu * (-5 + 12 * Nu)) * pow(PhiDOT, 2) *
                            pow(r, 2) * rDOT +
                        Complex(0, 6) * (85 - 510 * Nu + 808 * pow(Nu, 2)) *
                            PhiDOT * r * pow(rDOT, 2) -
                        8 * (-5 + 12 * Nu) * (-5 + 18 * Nu) * pow(rDOT, 3))) *
                  (kappa1 * pow(S1z, 2) - kappa2 * pow(S2z, 2)) +
              delta *
                  (Complex(0, -12600) *
                       (-135 + Nu * (1467 + 35 * Nu * (-151 + 178 * Nu))) *
                       pow(r, 4) * pow(PhiDOT * r - Complex(0, 1) * rDOT, 3) *
                       pow(PhiDOT * r + Complex(0, 1) * rDOT, 6) +
                   6 * pow(mass, 2) * pow(r, 2) *
                       (Complex(0, -10) *
                            (129184 +
                             Nu *
                                 (917839 + 2 * Nu * (-1601563 + 652752 * Nu))) *
                            pow(PhiDOT, 5) * pow(r, 5) +
                        8 *
                            (3087667 +
                             2 * Nu *
                                 (-3131631 +
                                  5 * Nu * (-385362 + 1212143 * Nu))) *
                            pow(PhiDOT, 4) * pow(r, 4) * rDOT +
                        Complex(0, 4) *
                            (1877161 +
                             Nu * (-1391951 +
                                   20 * Nu * (-594787 + 1654110 * Nu))) *
                            pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) +
                        2 *
                            (-919652 +
                             Nu * (-16365978 +
                                   5 * Nu * (5343599 + 5349122 * Nu))) *
                            pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 3) +
                        Complex(0, 14) *
                            (-19909 + Nu * (-2880941 +
                                            15 * Nu * (195291 + 744862 * Nu))) *
                            PhiDOT * r * pow(rDOT, 4) -
                        35 *
                            (1285 +
                             Nu * (-417939 + Nu * (376375 + 1794742 * Nu))) *
                            pow(rDOT, 5)) +
                   45 * mass * pow(r, 3) *
                       (Complex(0, 4) *
                            (54310 +
                             Nu * (-20821 + Nu * (-636477 + 820286 * Nu))) *
                            pow(PhiDOT, 7) * pow(r, 7) +
                        4 *
                            (-119296 + Nu * (251655 + 893222 * Nu -
                                             3122048 * pow(Nu, 2))) *
                            pow(PhiDOT, 6) * pow(r, 6) * rDOT -
                        Complex(0, 4) *
                            (42783 +
                             Nu * (123661 + 3 * Nu * (-447987 + 841994 * Nu))) *
                            pow(PhiDOT, 5) * pow(r, 5) * pow(rDOT, 2) -
                        8 *
                            (237862 +
                             Nu * (-912196 + 3 * Nu * (-46779 + 933382 * Nu))) *
                            pow(PhiDOT, 4) * pow(r, 4) * pow(rDOT, 3) -
                        Complex(0, 1) *
                            (1550225 +
                             6 * Nu *
                                 (-680405 +
                                  4 * Nu * (-402523 + 1277568 * Nu))) *
                            pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 4) -
                        16 *
                            (39855 +
                             2 * Nu * (-106258 + Nu * (113231 + 75438 * Nu))) *
                            pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 5) -
                        Complex(0, 28) *
                            (41703 +
                             2 * Nu *
                                 (-74181 + 20 * Nu * (-2603 + 15676 * Nu))) *
                            PhiDOT * r * pow(rDOT, 6) +
                        2240 *
                            (162 + Nu * (-516 + 7 * Nu * (-101 + 410 * Nu))) *
                            pow(rDOT, 7)) +
                   28 * pow(mass, 4) *
                       (Complex(0, 1) * PhiDOT * r *
                            (-5713898 - 261495 * kappa1 * pow(S1z, 2) +
                             2 * Nu *
                                 (6824514 - 4178445 * Nu +
                                  1430090 * pow(Nu, 2) -
                                  522990 * kappa1 * (-1 + Nu) * pow(S1z, 2)) +
                             1045980 * Nu * (-1 + 2 * Nu) * S1z * S2z -
                             261495 * kappa2 * pow(1 - 2 * Nu, 2) *
                                 pow(S2z, 2)) +
                        5 * rDOT *
                            (300238 + 27495 * kappa1 * pow(S1z, 2) +
                             Nu * (-581955 + 600300 * Nu - 535772 * pow(Nu, 2) +
                                   109980 * kappa1 * (-1 + Nu) * pow(S1z, 2)) -
                             109980 * Nu * (-1 + 2 * Nu) * S1z * S2z +
                             27495 * kappa2 * pow(1 - 2 * Nu, 2) *
                                 pow(S2z, 2))) +
                   2 * pow(mass, 3) * r *
                       (Complex(0, -1) * pow(PhiDOT, 3) * pow(r, 3) *
                            (-13 * (371212 + 4725 * kappa1 * pow(S1z, 2)) +
                             Nu * (-7716399 - 11038350 * Nu +
                                   21939560 * pow(Nu, 2) -
                                   245700 * kappa1 * (-1 + Nu) * pow(S1z, 2)) +
                             245700 * Nu * (-1 + 2 * Nu) * S1z * S2z -
                             61425 * kappa2 * pow(1 - 2 * Nu, 2) *
                                 pow(S2z, 2)) -
                        Complex(0, 126) * PhiDOT * r * pow(rDOT, 2) *
                            ((1974194 - 5423729 * Nu + 339235 * pow(Nu, 2) +
                              5947390 * pow(Nu, 3) -
                              49725 * kappa1 * pow(1 - 2 * Nu, 2) *
                                  pow(S1z, 2)) /
                                 3. +
                             66300 * Nu * (-1 + 2 * Nu) * S1z * S2z -
                             16575 * kappa2 * pow(1 - 2 * Nu, 2) *
                                 pow(S2z, 2)) +
                        140 * pow(rDOT, 3) *
                            (175718 - 5850 * kappa1 * pow(S1z, 2) +
                             Nu * (-730551 + 504765 * Nu + 773906 * pow(Nu, 2) -
                                   23400 * kappa1 * (-1 + Nu) * pow(S1z, 2)) +
                             23400 * Nu * (-1 + 2 * Nu) * S1z * S2z -
                             5850 * kappa2 * pow(1 - 2 * Nu, 2) * pow(S2z, 2)) +
                        6 * pow(PhiDOT, 2) * pow(r, 2) * rDOT *
                            (-8099876 + 225225 * kappa1 * pow(S1z, 2) +
                             4 * Nu *
                                 (1323649 + 9536535 * Nu -
                                  6899570 * pow(Nu, 2) +
                                  225225 * kappa1 * (-1 + Nu) * pow(S1z, 2)) -
                             900900 * Nu * (-1 + 2 * Nu) * S1z * S2z +
                             225225 * kappa2 * pow(1 - 2 * Nu, 2) *
                                 pow(S2z, 2))))) /
             (1.57248e6 * sqrt(330) * pow(r, 4))));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hQC_5_m_3(REAL8 Nu, UINT4 vpnorder, REAL8 x) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  /* if(vpnorder == 3){
       return((Complex(0,-0.5625)*sqrt(0.02727272727272727)*delta +
       Complex(0,1.125)*sqrt(0.02727272727272727)*delta*Nu)*pow(x,2.5));
   }
   else if(vpnorder == 5){
       return((Complex(0,2.985576923076923)*sqrt(0.02727272727272727)*delta -
       Complex(0,6.6923076923076925)*sqrt(0.02727272727272727)*delta*Nu +
       Complex(0,0.11538461538461539)*sqrt(3.3)*delta*pow(Nu,2))*
       pow(x,3.5));
   } */
  /* else */ if (vpnorder == 6) {
    return (Complex(0, 1.6875) * sqrt(0.02727272727272727) * delta *
            (-1 + 2 * Nu) * pow(x, 4) * (M_PI + Complex(0, 2) * log(1.5)));
  } else {
    return 0;
  }
}

static COMPLEX16 hl_5_m_3(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z, REAL8 x) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_5_m_3: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_5_m_3(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z) +
            hQC_5_m_3(Nu, vpnorder, x)) *
           cpolar(1, -3 * Phi);
  }
}

static COMPLEX16 hl_5_m_min3(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z, REAL8 x) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_5_m_min3: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) * pow(-1, 5) *
           conj(hGO_5_m_3(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z) +
                hQC_5_m_3(Nu, vpnorder, x)) *
           cpolar(1, 3 * Phi);
  }
}

// H52
static COMPLEX16 hGO_5_m_2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  if (vpnorder == 4) {
    return (mass * (1 - 5 * Nu + 5 * pow(Nu, 2)) * PhiDOT *
            (mass * (41 * PhiDOT * r + Complex(0, 22) * rDOT) -
             3 * r *
                 (11 * pow(PhiDOT, 3) * pow(r, 3) -
                  Complex(0, 3) * pow(PhiDOT, 2) * pow(r, 2) * rDOT +
                  6 * PhiDOT * r * pow(rDOT, 2) +
                  Complex(0, 2) * pow(rDOT, 3)))) /
           (54. * sqrt(55) * r);
  }

  else if (vpnorder == 5) {

    return (/* (2*Nu*(S1z - delta*S1z + S2z + delta*S2z - 2*Nu*(S1z +
      S2z))*pow(x,3.5))/(9.*sqrt(55))
      + */
            /* Henry et al. ecc spin terms */ (
                -0.027777777777777776 *
                (pow(mass, 2) * Nu *
                 (mass * (41 * PhiDOT * r + Complex(0, 22) * rDOT) -
                  3 * r *
                      (11 * pow(PhiDOT, 3) * pow(r, 3) -
                       Complex(0, 3) * pow(PhiDOT, 2) * pow(r, 2) * rDOT +
                       6 * PhiDOT * r * pow(rDOT, 2) +
                       Complex(0, 2) * pow(rDOT, 3))) *
                 ((-1 + delta + 2 * Nu) * S1z - (1 + delta - 2 * Nu) * S2z)) /
                (sqrt(55) * pow(r, 3))));
  }

  else if (vpnorder == 6) {
    return (mass * PhiDOT *
            (4 * pow(mass, 2) *
                 (13 *
                      (-5051 + 28623 * Nu - 46305 * pow(Nu, 2) +
                       29470 * pow(Nu, 3)) *
                      PhiDOT * r +
                  Complex(0, 7) *
                      (5684 - 26697 * Nu + 14225 * pow(Nu, 2) +
                       25355 * pow(Nu, 3)) *
                      rDOT) -
             2 * mass * r *
                 ((23157 + 154 * Nu - 648410 * pow(Nu, 2) +
                   1133195 * pow(Nu, 3)) *
                      pow(PhiDOT, 3) * pow(r, 3) -
                  Complex(0, 14) *
                      (6648 - 31729 * Nu + 23795 * pow(Nu, 2) +
                       13225 * pow(Nu, 3)) *
                      pow(PhiDOT, 2) * pow(r, 2) * rDOT +
                  (363357 - 1825691 * Nu + 1714720 * pow(Nu, 2) +
                   395255 * pow(Nu, 3)) *
                      PhiDOT * r * pow(rDOT, 2) +
                  Complex(0, 14) *
                      (9717 - 48896 * Nu + 46000 * pow(Nu, 2) +
                       10865 * pow(Nu, 3)) *
                      pow(rDOT, 3)) +
             15 * pow(r, 2) *
                 (2 *
                      (2207 - 6076 * Nu - 18424 * pow(Nu, 2) +
                       38787 * pow(Nu, 3)) *
                      pow(PhiDOT, 5) * pow(r, 5) -
                  Complex(0, 7) *
                      (4744 - 23965 * Nu + 23906 * pow(Nu, 2) +
                       1892 * pow(Nu, 3)) *
                      pow(PhiDOT, 4) * pow(r, 4) * rDOT +
                  2 *
                      (5667 - 22652 * Nu - 5747 * pow(Nu, 2) +
                       45416 * pow(Nu, 3)) *
                      pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) -
                  Complex(0, 14) *
                      (97 - 536 * Nu + 643 * pow(Nu, 2) + 36 * pow(Nu, 3)) *
                      pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 3) +
                  4 *
                      (1009 - 4060 * Nu - 889 * pow(Nu, 2) +
                       7952 * pow(Nu, 3)) *
                      PhiDOT * r * pow(rDOT, 4) +
                  Complex(0, 28) *
                      (45 - 172 * Nu - 85 * pow(Nu, 2) + 400 * pow(Nu, 3)) *
                      pow(rDOT, 5)))) /
           (98280. * sqrt(55) * pow(r, 2));
  }

  else if (vpnorder == 7) {

    return (
        (/* (698*pow(Nu,3)*(S1z + S2z) - 104*(S1z + delta*S1z + S2z - delta*S2z)
   - 2*pow(Nu,2)*(457*delta*(S1z - S2z) + 48*(S1z + S2z)) +
  Nu*(1055*delta*(S1z - S2z) + 193*(S1z + S2z)))*pow(x,4.5))/(351.*sqrt(55))
  + */
         /* Henry et al. ecc+spin terms */
         (-1.6958350291683625e-6 *
          (pow(mass, 2) *
           (4 * pow(mass, 2) *
                (210 * pow(Nu, 3) *
                     (975 * PhiDOT * r + Complex(0, 1373) * rDOT) *
                     (S1z + S2z) +
                 10920 * (10 * PhiDOT * r + Complex(0, 7) * rDOT) *
                     (S1z + delta * S1z + S2z - delta * S2z) +
                 21 * pow(Nu, 2) *
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
            3 * pow(r, 2) *
                (840 * pow(Nu, 3) *
                     (896 * pow(PhiDOT, 5) * pow(r, 5) +
                      Complex(0, 1390) * pow(PhiDOT, 4) * pow(r, 4) * rDOT +
                      2682 * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) -
                      Complex(0, 235) * pow(PhiDOT, 2) * pow(r, 2) *
                          pow(rDOT, 3) +
                      716 * PhiDOT * r * pow(rDOT, 4) +
                      Complex(0, 210) * pow(rDOT, 5)) *
                     (S1z + S2z) -
                 Complex(0, 286650) * pow(PhiDOT, 4) * pow(r, 4) * rDOT *
                     (S1z + delta * S1z + S2z - delta * S2z) -
                 14 * pow(Nu, 2) *
                     (Complex(0, 2) * pow(PhiDOT, 4) * pow(r, 4) * rDOT *
                          (Complex(0, -41392) +
                           15 * (15219 + 1835 * delta) * S1z +
                           (228285 - 27525 * delta) * S2z) +
                      Complex(0, 210) * pow(rDOT, 5) *
                          ((-23 + 5 * delta) * S1z - (23 + 5 * delta) * S2z) +
                      2 * PhiDOT * r * pow(rDOT, 4) *
                          (Complex(0, 7592) + 15 * (-399 + 193 * delta) * S1z -
                           15 * (399 + 193 * delta) * S2z) +
                      10 * pow(PhiDOT, 5) * pow(r, 5) *
                          (Complex(0, 1508) + 3 * (-1889 + 483 * delta) * S1z -
                           3 * (1889 + 483 * delta) * S2z) +
                      pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 3) *
                          (34424 - Complex(0, 15) * (-161 + 555 * delta) * S1z +
                           Complex(0, 15) * (161 + 555 * delta) * S2z) +
                      3 * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) *
                          (Complex(0, 15288) + 5 * (539 + 907 * delta) * S1z -
                           5 * (-539 + 907 * delta) * S2z)) +
                 3 * Nu *
                     (Complex(0, 140) * pow(rDOT, 5) *
                          ((-141 + 193 * delta) * S1z -
                           (141 + 193 * delta) * S2z) +
                      14 * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 3) *
                          (4108 - Complex(0, 5) * (-193 + 453 * delta) * S1z +
                           Complex(0, 5) * (193 + 453 * delta) * S2z) +
                      4 * PhiDOT * r * pow(rDOT, 4) *
                          (Complex(0, 5148) + 35 * (-415 + 623 * delta) * S1z -
                           35 * (415 + 623 * delta) * S2z) +
                      10 * pow(PhiDOT, 5) * pow(r, 5) *
                          (Complex(0, 2223) + 14 * (-859 + 963 * delta) * S1z -
                           14 * (859 + 963 * delta) * S2z) +
                      7 * pow(PhiDOT, 4) * pow(r, 4) * rDOT *
                          (16276 -
                           Complex(0, 5) * (-31495 + 2453 * delta) * S1z +
                           Complex(0, 5) * (31495 + 2453 * delta) * S2z) +
                      2 * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) *
                          (Complex(0, 39156) +
                           35 * (-1773 + 3437 * delta) * S1z -
                           35 * (1773 + 3437 * delta) * S2z))) -
            6 * mass * r *
                (210 * pow(Nu, 3) *
                     (946 * pow(PhiDOT, 3) * pow(r, 3) -
                      Complex(0, 1150) * pow(PhiDOT, 2) * pow(r, 2) * rDOT +
                      739 * PhiDOT * r * pow(rDOT, 2) +
                      Complex(0, 366) * pow(rDOT, 3)) *
                     (S1z + S2z) -
                 5460 *
                     (3 * pow(PhiDOT, 3) * pow(r, 3) +
                      Complex(0, 17) * pow(PhiDOT, 2) * pow(r, 2) * rDOT -
                      15 * PhiDOT * r * pow(rDOT, 2) -
                      Complex(0, 3) * pow(rDOT, 3)) *
                     (S1z + delta * S1z + S2z - delta * S2z) -
                 14 * pow(Nu, 2) *
                     (Complex(0, 1) * pow(PhiDOT, 2) * pow(r, 2) * rDOT *
                          (Complex(0, -106704) +
                           5 * (13991 + 6349 * delta) * S1z +
                           (69955 - 31745 * delta) * S2z) +
                      Complex(0, 50) * pow(rDOT, 3) *
                          ((-1330 + 223 * delta) * S1z -
                           (1330 + 223 * delta) * S2z) +
                      PhiDOT * r * pow(rDOT, 2) *
                          (Complex(0, 15912) +
                           5 * (-39389 + 2528 * delta) * S1z -
                           5 * (39389 + 2528 * delta) * S2z) +
                      pow(PhiDOT, 3) * pow(r, 3) *
                          (Complex(0, 62868) +
                           5 * (-5543 + 5105 * delta) * S1z -
                           5 * (5543 + 5105 * delta) * S2z)) +
                 Nu * (Complex(0, 28) * pow(PhiDOT, 2) * pow(r, 2) * rDOT *
                           (Complex(0, -17147) +
                            5 * (5007 + 1675 * delta) * S1z +
                            (25035 - 8375 * delta) * S2z) +
                       5 * pow(PhiDOT, 3) * pow(r, 3) *
                           (Complex(0, 57161) + 28 * (495 + 272 * delta) * S1z -
                            28 * (-495 + 272 * delta) * S2z) +
                       Complex(0, 700) * pow(rDOT, 3) *
                           ((-651 + 482 * delta) * S1z -
                            (651 + 482 * delta) * S2z) +
                       PhiDOT * r * pow(rDOT, 2) *
                           (Complex(0, 68406) +
                            35 * (-40527 + 19571 * delta) * S1z -
                            35 * (40527 + 19571 * delta) * S2z))))) /
          (sqrt(55) * pow(r, 4)))));
  } else {
    return 0;
  }
}

/* static COMPLEX16 hQC_5_m_2(REAL8 Nu, UINT4 vpnorder, REAL8 x){

    if(vpnorder == 4){
        return((4/(27.*sqrt(55)) - (4*sqrt(0.45454545454545453)*Nu)/27. +
        (4*sqrt(0.45454545454545453)*pow(Nu,2))/27.)*pow(x,3));
    }
    else if(vpnorder == 6){
        return((-7822/(12285.*sqrt(55)) + (6158*Nu)/(1755.*sqrt(55)) -
        (1652*pow(Nu,2))/(351.*sqrt(55)) +
        (14*sqrt(2.2)*pow(Nu,3))/117.)*pow(x,4));
    }
    else{
        return 0;
    }
} */

static COMPLEX16 hl_5_m_2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_5_m_2: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_5_m_2(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z,
                      S2z) /* +hQC_5_m_2(Nu,vpnorder,x) */) *
           cpolar(1, -2 * Phi);
  }
}

static COMPLEX16 hl_5_m_min2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_5_m_min2: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) * pow(-1, 5) *
           conj(hGO_5_m_2(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z,
                          S2z) /* +hQC_5_m_2(Nu,vpnorder,x) */) *
           cpolar(1, 2 * Phi);
  }
}

// H51
static COMPLEX16 hGO_5_m_1(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z) {
  REAL8 delta = sqrt(1 - 4 * Nu);
  REAL8 kappa1 = 1.0;
  REAL8 kappa2 = 1.0;

  if (vpnorder == 3) {
    return (delta * (-1 + 2 * Nu) *
            (120 * pow(r, 2) * pow(Complex(0, 1) * PhiDOT * r - rDOT, 3) *
                 pow(PhiDOT * r - Complex(0, 1) * rDOT, 2) +
             2 * pow(mass, 2) * (Complex(0, -86) * PhiDOT * r + 205 * rDOT) +
             Complex(0, 3) * mass * r *
                 (97 * pow(PhiDOT, 3) * pow(r, 3) +
                  Complex(0, 160) * pow(PhiDOT, 2) * pow(r, 2) * rDOT +
                  132 * PhiDOT * r * pow(rDOT, 2) +
                  Complex(0, 160) * pow(rDOT, 3)))) /
           (144. * sqrt(385) * pow(r, 2));
  }

  else if (vpnorder == 5) {
    return (delta *
            (-360 * (33 - 197 * Nu + 294 * pow(Nu, 2)) * pow(r, 3) *
                 pow(PhiDOT * r + Complex(0, 1) * rDOT, 4) *
                 pow(Complex(0, 1) * PhiDOT * r + rDOT, 3) +
             2 * pow(mass, 3) *
                 (Complex(0, -1) * (53311 - 121906 * Nu + 42816 * pow(Nu, 2)) *
                      PhiDOT * r +
                  78 * (1141 - 2760 * Nu + 1420 * pow(Nu, 2)) * rDOT) +
             2 * pow(mass, 2) * r *
                 (Complex(0, 2) * (29938 - 82195 * Nu + 59238 * pow(Nu, 2)) *
                      pow(PhiDOT, 3) * pow(r, 3) +
                  (-27026 + 120623 * Nu - 199326 * pow(Nu, 2)) *
                      pow(PhiDOT, 2) * pow(r, 2) * rDOT +
                  Complex(0, 4) * (-2636 - 11335 * Nu + 43962 * pow(Nu, 2)) *
                      PhiDOT * r * pow(rDOT, 2) -
                  39 * (-639 - 735 * Nu + 5690 * pow(Nu, 2)) * pow(rDOT, 3)) +
             3 * mass * pow(r, 2) *
                 (Complex(0, -4) * (3115 - 15385 * Nu + 22098 * pow(Nu, 2)) *
                      pow(PhiDOT, 5) * pow(r, 5) +
                  2 * (-2108 - 17893 * Nu + 56466 * pow(Nu, 2)) *
                      pow(PhiDOT, 4) * pow(r, 4) * rDOT -
                  Complex(0, 3) * (-8473 - 9528 * Nu + 65204 * pow(Nu, 2)) *
                      pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) +
                  8 * (-2692 - 6587 * Nu + 29754 * pow(Nu, 2)) *
                      pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 3) -
                  Complex(0, 4) * (-3971 - 3226 * Nu + 27360 * pow(Nu, 2)) *
                      PhiDOT * r * pow(rDOT, 4) +
                  240 * (-57 - 97 * Nu + 518 * pow(Nu, 2)) * pow(rDOT, 5)))) /
           (11232. * sqrt(385) * pow(r, 3));
  }

  else if (vpnorder == 6) {
    return (
        -(delta * pow(mass, 2) * Nu *
          (17830 * pow(mass, 2) -
           6 * mass * r *
               (4723 * pow(PhiDOT, 2) * pow(r, 2) -
                Complex(0, 3921) * PhiDOT * r * rDOT + 4530 * pow(rDOT, 2)) +
           pow(r, 2) * (12629 * pow(PhiDOT, 4) * pow(r, 4) -
                        Complex(0, 24248) * pow(PhiDOT, 3) * pow(r, 3) * rDOT +
                        20064 * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) -
                        Complex(0, 9968) * PhiDOT * r * pow(rDOT, 3) +
                        7120 * pow(rDOT, 4)))) /
            (5040. * sqrt(385) * pow(r, 4))
        /* +((Complex(0,-0.0023148148148148147)*((5 + Nu*(-23 + 19*Nu)
  + delta*(5 + 27*(-1 + Nu)*Nu))*S1z + (-5 + (23 - 19*Nu)*Nu
  + delta*(5 + 27*(-1 + Nu)*Nu))*S2z)*pow(x,4))/sqrt(385)) */
        +
        /* Henry et al. ecc spin terms */ (
            (Complex(0, -0.0023148148148148147) * pow(mass, 2) *
             (pow(mass, 2) * ((90 - 622 * Nu + 966 * pow(Nu, 2) +
                               delta * (90 - 524 * Nu + 770 * pow(Nu, 2))) *
                                  S1z +
                              2 *
                                  (-45 + 311 * Nu - 483 * pow(Nu, 2) +
                                   delta * (45 - 262 * Nu + 385 * pow(Nu, 2))) *
                                  S2z) +
              3 * pow(r, 2) *
                  (Complex(0, 4) * PhiDOT * r * pow(rDOT, 3) *
                       (-5 * (3 - 20 * Nu + 30 * pow(Nu, 2)) * S1z +
                        delta * (-15 - 106 * Nu + 262 * pow(Nu, 2)) * S1z +
                        5 * (3 - 20 * Nu + 30 * pow(Nu, 2)) * S2z +
                        delta * (-15 - 106 * Nu + 262 * pow(Nu, 2)) * S2z) +
                   Complex(0, 1) * pow(PhiDOT, 3) * pow(r, 3) * rDOT *
                       (-35 * (3 - 22 * Nu + 36 * pow(Nu, 2)) * S1z +
                        delta * (-105 - 494 * Nu + 1268 * pow(Nu, 2)) * S1z +
                        35 * (3 - 22 * Nu + 36 * pow(Nu, 2)) * S2z +
                        delta * (-105 - 494 * Nu + 1268 * pow(Nu, 2)) * S2z) +
                   4 * pow(rDOT, 4) *
                       ((15 - 92 * Nu + 126 * pow(Nu, 2) +
                         delta * (15 - 64 * Nu + 70 * pow(Nu, 2))) *
                            S1z +
                        (-15 + 92 * Nu - 126 * pow(Nu, 2) +
                         delta * (15 - 64 * Nu + 70 * pow(Nu, 2))) *
                            S2z) +
                   pow(PhiDOT, 4) * pow(r, 4) *
                       ((-180 + 1009 * Nu - 1227 * pow(Nu, 2) +
                         delta * (-180 + 95 * Nu + 601 * pow(Nu, 2))) *
                            S1z +
                        (180 - 1009 * Nu + 1227 * pow(Nu, 2) +
                         delta * (-180 + 95 * Nu + 601 * pow(Nu, 2))) *
                            S2z) +
                   3 * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) *
                       ((35 - 198 * Nu + 244 * pow(Nu, 2) +
                         delta * (35 - 382 * Nu + 612 * pow(Nu, 2))) *
                            S1z +
                        (-35 + 198 * Nu - 244 * pow(Nu, 2) +
                         delta * (35 - 382 * Nu + 612 * pow(Nu, 2))) *
                            S2z)) +
              mass * r *
                  (Complex(0, -2) * PhiDOT * r * rDOT *
                       (-5 * (78 - 499 * Nu + 717 * pow(Nu, 2)) * S1z +
                        delta * (-390 - 677 * Nu + 2759 * pow(Nu, 2)) * S1z +
                        5 * (78 - 499 * Nu + 717 * pow(Nu, 2)) * S2z +
                        delta * (-390 - 677 * Nu + 2759 * pow(Nu, 2)) * S2z) +
                   pow(PhiDOT, 2) * pow(r, 2) *
                       ((609 - 3352 * Nu + 3966 * pow(Nu, 2) +
                         delta * (609 + 58 * Nu - 2854 * pow(Nu, 2))) *
                            S1z +
                        (-609 + 3352 * Nu - 3966 * pow(Nu, 2) +
                         delta * (609 + 58 * Nu - 2854 * pow(Nu, 2))) *
                            S2z) -
                   2 * pow(rDOT, 2) *
                       ((255 - 1592 * Nu + 2226 * pow(Nu, 2) +
                         delta * (255 - 1144 * Nu + 1330 * pow(Nu, 2))) *
                            S1z +
                        (-255 + 1592 * Nu - 2226 * pow(Nu, 2) +
                         delta * (255 - 1144 * Nu + 1330 * pow(Nu, 2))) *
                            S2z)))) /
            (sqrt(385) * pow(r, 4))));
  }

  else if (vpnorder == 7) {

    return (/* (Complex(0,0.001736111111111111)*(kappa1*(5 + 5*delta*pow(1 -
      2*Nu,2)
      - 2*Nu*(15 + 4*Nu))*pow(S1z,2) + S2z*(20*delta*(1 - 2*Nu)*Nu*S1z
      + kappa2*(-5 + 5*delta*pow(1 - 2*Nu,2) + 30*Nu +
      8*pow(Nu,2))*S2z))*pow(x,4.5))/sqrt(385)
      + */
            /* Henry et al. ecc+spin terms */
            ((3510 * pow(mass, 3) *
                  (Complex(0, 2) * mass * (149 + 2 * Nu * (-447 + 508 * Nu)) *
                       PhiDOT * r -
                   2 * mass * (235 + 2 * Nu * (-705 + 1036 * Nu)) * rDOT +
                   r * (Complex(0, -7) * (35 - 210 * Nu + 232 * pow(Nu, 2)) *
                            pow(PhiDOT, 3) * pow(r, 3) +
                        2 * (115 + 2 * Nu * (-345 + 628 * Nu)) *
                            pow(PhiDOT, 2) * pow(r, 2) * rDOT -
                        Complex(0, 2) * (85 - 510 * Nu + 872 * pow(Nu, 2)) *
                            PhiDOT * r * pow(rDOT, 2) +
                        8 * (-5 + 14 * Nu) * (-5 + 16 * Nu) * pow(rDOT, 3))) *
                  (kappa1 * pow(S1z, 2) - kappa2 * pow(S2z, 2)) +
              delta *
                  (5400 * (-135 + Nu * (1467 + 35 * Nu * (-151 + 178 * Nu))) *
                       pow(r, 4) * pow(PhiDOT * r - Complex(0, 1) * rDOT, 4) *
                       pow(Complex(0, -1) * PhiDOT * r + rDOT, 5) +
                   6 * pow(mass, 2) * pow(r, 2) *
                       (Complex(0, -2) *
                            (556064 +
                             7 * Nu *
                                 (259253 + 130 * Nu * (-10297 + 8136 * Nu))) *
                            pow(PhiDOT, 5) * pow(r, 5) +
                        8 *
                            (349973 +
                             4 * Nu *
                                 (-94417 + 5 * Nu * (-65429 + 126618 * Nu))) *
                            pow(PhiDOT, 4) * pow(r, 4) * rDOT -
                        Complex(0, 60) *
                            (37369 + Nu * (-136791 - 60188 * Nu +
                                           607348 * pow(Nu, 2))) *
                            pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) +
                        2 *
                            (706852 +
                             Nu * (-3817342 +
                                   15 * Nu * (-114953 + 1628610 * Nu))) *
                            pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 3) -
                        Complex(0, 2) *
                            (-19909 + Nu * (-2880941 +
                                            15 * Nu * (195291 + 744862 * Nu))) *
                            PhiDOT * r * pow(rDOT, 4) +
                        15 *
                            (1285 +
                             Nu * (-417939 + Nu * (376375 + 1794742 * Nu))) *
                            pow(rDOT, 5)) +
                   45 * mass * pow(r, 3) *
                       (Complex(0, 4) *
                            (11270 +
                             Nu * (40035 + Nu * (-338333 + 490662 * Nu))) *
                            pow(PhiDOT, 7) * pow(r, 7) -
                        4 *
                            (12816 +
                             Nu * (1923 - 262738 * Nu + 582488 * pow(Nu, 2))) *
                            pow(PhiDOT, 6) * pow(r, 6) * rDOT +
                        Complex(0, 4) *
                            (76389 +
                             Nu * (-196893 + Nu * (-529811 + 1606402 * Nu))) *
                            pow(PhiDOT, 5) * pow(r, 5) * pow(rDOT, 2) +
                        8 *
                            (-45746 +
                             Nu * (115572 + (313303 - 946222 * Nu) * Nu)) *
                            pow(PhiDOT, 4) * pow(r, 4) * pow(rDOT, 3) +
                        Complex(0, 1) *
                            (464727 +
                             2 * Nu *
                                 (-792537 + 4 * Nu * (-164359 + 877280 * Nu))) *
                            pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 4) -
                        32 *
                            (13500 +
                             Nu * (-39687 + Nu * (-69661 + 249182 * Nu))) *
                            pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 5) +
                        Complex(0, 4) *
                            (41703 +
                             2 * Nu *
                                 (-74181 + 20 * Nu * (-2603 + 15676 * Nu))) *
                            PhiDOT * r * pow(rDOT, 6) -
                        960 * (162 + Nu * (-516 + 7 * Nu * (-101 + 410 * Nu))) *
                            pow(rDOT, 7)) +
                   4 * pow(mass, 4) *
                       (Complex(0, -1) * PhiDOT * r *
                            (-5713898 - 261495 * kappa1 * pow(S1z, 2) +
                             2 * Nu *
                                 (6824514 - 4178445 * Nu +
                                  1430090 * pow(Nu, 2) -
                                  522990 * kappa1 * (-1 + Nu) * pow(S1z, 2)) +
                             1045980 * Nu * (-1 + 2 * Nu) * S1z * S2z -
                             261495 * kappa2 * pow(1 - 2 * Nu, 2) *
                                 pow(S2z, 2)) +
                        15 * rDOT *
                            (-300238 - 27495 * kappa1 * pow(S1z, 2) +
                             Nu * (581955 - 600300 * Nu + 535772 * pow(Nu, 2) -
                                   109980 * kappa1 * (-1 + Nu) * pow(S1z, 2)) +
                             109980 * Nu * (-1 + 2 * Nu) * S1z * S2z -
                             27495 * kappa2 * pow(1 - 2 * Nu, 2) *
                                 pow(S2z, 2))) +
                   6 * pow(mass, 3) * r *
                       (Complex(0, 2) * PhiDOT * r * pow(rDOT, 2) *
                            (1974194 - 49725 * kappa1 * pow(S1z, 2) +
                             Nu * (-5423729 + 339235 * Nu +
                                   5947390 * pow(Nu, 2) -
                                   198900 * kappa1 * (-1 + Nu) * pow(S1z, 2)) +
                             198900 * Nu * (-1 + 2 * Nu) * S1z * S2z -
                             49725 * kappa2 * pow(1 - 2 * Nu, 2) *
                                 pow(S2z, 2)) +
                        Complex(0, 25) * pow(PhiDOT, 3) * pow(r, 3) *
                            ((-265372 +
                              Nu * (1803657 +
                                    2 * Nu * (-1845967 + 746940 * Nu))) /
                                 5. -
                             5733 * kappa1 * pow(1 - 2 * Nu, 2) * pow(S1z, 2) +
                             22932 * Nu * (-1 + 2 * Nu) * S1z * S2z -
                             5733 * kappa2 * pow(1 - 2 * Nu, 2) * pow(S2z, 2)) +
                        20 * pow(rDOT, 3) *
                            (-175718 + 5850 * kappa1 * pow(S1z, 2) +
                             Nu * (730551 - 504765 * Nu - 773906 * pow(Nu, 2) +
                                   23400 * kappa1 * (-1 + Nu) * pow(S1z, 2)) -
                             23400 * Nu * (-1 + 2 * Nu) * S1z * S2z +
                             5850 * kappa2 * pow(1 - 2 * Nu, 2) * pow(S2z, 2)) +
                        2 * pow(PhiDOT, 2) * pow(r, 2) * rDOT *
                            (-2687884 + 67275 * kappa1 * pow(S1z, 2) +
                             4 * Nu *
                                 (917051 + 1724565 * Nu - 1611110 * pow(Nu, 2) +
                                  67275 * kappa1 * (-1 + Nu) * pow(S1z, 2)) -
                             269100 * Nu * (-1 + 2 * Nu) * S1z * S2z +
                             67275 * kappa2 * pow(1 - 2 * Nu, 2) *
                                 pow(S2z, 2))))) /
             (673920. * sqrt(385) * pow(r, 4))));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hQC_5_m_1(REAL8 Nu, UINT4 vpnorder, REAL8 x) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  /* if(vpnorder == 3){
       return(((Complex(0,0.006944444444444444)*delta)/sqrt(385) -
       (Complex(0,0.013888888888888888)*delta*Nu)/sqrt(385))*pow(x,2.5));
   }
   else if(vpnorder == 5){
       return(((Complex(0,-0.03187321937321937)*delta)/sqrt(385) +
        Complex(0,0.005698005698005698)*sqrt(0.3142857142857143)*delta*
        Nu - (Complex(0,0.0007122507122507123)*delta*pow(Nu,2))/
        sqrt(385))*pow(x,3.5));
   } */
  /* else */ if (vpnorder == 6) {
    return ((Complex(0, -0.006944444444444444) * delta * (-1 + 2 * Nu) *
             pow(x, 4) * (M_PI - Complex(0, 2) * log(2))) /
            sqrt(385));
  } else {
    return 0;
  }
}

static COMPLEX16 hl_5_m_1(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z, REAL8 x) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_5_m_1: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_5_m_1(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z) +
            hQC_5_m_1(Nu, vpnorder, x)) *
           cpolar(1, -1 * Phi);
  }
}

static COMPLEX16 hl_5_m_min1(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z, REAL8 x) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_5_m_min1: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) * pow(-1, 5) *
           conj(hGO_5_m_1(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z) +
                hQC_5_m_1(Nu, vpnorder, x)) *
           cpolar(1, 1 * Phi);
  }
}

// 66
static COMPLEX16 hGO_6_m_6(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  if (vpnorder == 4) {
    return (((1 - 5 * Nu + 5 * pow(Nu, 2)) *
             (172 * pow(mass, 3) +
              120 * pow(r, 3) * pow(PhiDOT * r + Complex(0, 1) * rDOT, 6) +
              pow(mass, 2) * r *
                  (3269 * pow(PhiDOT, 2) * pow(r, 2) +
                   Complex(0, 2920) * PhiDOT * r * rDOT - 806 * pow(rDOT, 2)) +
              15 * mass * pow(r, 2) *
                  (281 * pow(PhiDOT, 4) * pow(r, 4) +
                   Complex(0, 494) * pow(PhiDOT, 3) * pow(r, 3) * rDOT -
                   444 * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) -
                   Complex(0, 208) * PhiDOT * r * pow(rDOT, 3) +
                   40 * pow(rDOT, 4)))) /
            (360. * sqrt(143) * pow(r, 3)));
  }

  else if (vpnorder == 6) {
    return (
        (14 * pow(mass, 4) *
             (-5740 + 29361 * Nu - 33348 * pow(Nu, 2) + 7334 * pow(Nu, 3)) -
         7560 * (-1 + 9 * Nu - 26 * pow(Nu, 2) + 23 * pow(Nu, 3)) * pow(r, 4) *
             (PhiDOT * r - Complex(0, 1) * rDOT) *
             pow(PhiDOT * r + Complex(0, 1) * rDOT, 7) +
         2 * pow(mass, 3) * r *
             ((-539645 + 2950311 * Nu - 4086684 * pow(Nu, 2) +
               1644517 * pow(Nu, 3)) *
                  pow(PhiDOT, 2) * pow(r, 2) +
              Complex(0, 6) *
                  (-47984 + 275121 * Nu - 442540 * pow(Nu, 2) +
                   255850 * pow(Nu, 3)) *
                  PhiDOT * r * rDOT +
              14 *
                  (3614 - 21621 * Nu + 39684 * pow(Nu, 2) -
                   29332 * pow(Nu, 3)) *
                  pow(rDOT, 2)) +
         3 * pow(mass, 2) * pow(r, 2) *
             ((-311847 + 1966993 * Nu - 3502751 * pow(Nu, 2) +
               1752968 * pow(Nu, 3)) *
                  pow(PhiDOT, 4) * pow(r, 4) +
              Complex(0, 4) *
                  (629 + 160412 * Nu - 846370 * pow(Nu, 2) +
                   912975 * pow(Nu, 3)) *
                  pow(PhiDOT, 3) * pow(r, 3) * rDOT -
              3 *
                  (65519 - 144403 * Nu - 684796 * pow(Nu, 2) +
                   1205253 * pow(Nu, 3)) *
                  pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) -
              Complex(0, 32) *
                  (3921 - 11389 * Nu - 27265 * pow(Nu, 2) +
                   58450 * pow(Nu, 3)) *
                  PhiDOT * r * pow(rDOT, 3) +
              14 *
                  (1867 - 5501 * Nu - 12824 * pow(Nu, 2) + 28137 * pow(Nu, 3)) *
                  pow(rDOT, 4)) -
         45 * mass * pow(r, 3) *
             ((195 + 3619 * Nu - 36617 * pow(Nu, 2) + 66836 * pow(Nu, 3)) *
                  pow(PhiDOT, 6) * pow(r, 6) +
              Complex(0, 4) *
                  (-1878 + 10969 * Nu - 20741 * pow(Nu, 2) +
                   18263 * pow(Nu, 3)) *
                  pow(PhiDOT, 5) * pow(r, 5) * rDOT +
              (17169 - 75446 * Nu + 35497 * pow(Nu, 2) + 47054 * pow(Nu, 3)) *
                  pow(PhiDOT, 4) * pow(r, 4) * pow(rDOT, 2) +
              Complex(0, 2) *
                  (9183 - 30296 * Nu - 37835 * pow(Nu, 2) +
                   95060 * pow(Nu, 3)) *
                  pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 3) -
              4 * (2781 - 6062 * Nu - 28595 * pow(Nu, 2) + 49070 * pow(Nu, 3)) *
                  pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 4) -
              Complex(0, 16) *
                  (228 - 217 * Nu - 3871 * pow(Nu, 2) + 5803 * pow(Nu, 3)) *
                  PhiDOT * r * pow(rDOT, 5) +
              56 * (9 + 4 * Nu - 221 * pow(Nu, 2) + 308 * pow(Nu, 3)) *
                  pow(rDOT, 6))) /
        (15120. * sqrt(143) * pow(r, 4)));
  }

  else if (vpnorder == 7) {

    return (/* (108*(110*pow(Nu,3)*(S1z + S2z) - 14*(S1z + delta*S1z + S2z -
      delta*S2z)
      - 50*pow(Nu,2)*(2*delta*(S1z - S2z) + 3*(S1z + S2z)) +  Nu*(85*delta*(S1z
      - S2z)
      + 83*(S1z + S2z)))*pow(x,4.5))/(35.*sqrt(143))
      + */
            /* Henry et al. ecc+spin terms */
            ((pow(mass, 2) *
              (8 * pow(mass, 2) *
                   (30 * pow(Nu, 3) *
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
                    3 * pow(Nu, 2) *
                        (Complex(0, 1) * rDOT *
                             (Complex(0, 15107) +
                              40 * (1592 + 485 * delta) * S1z + 63680 * S2z -
                              19400 * delta * S2z) +
                         5 * PhiDOT * r *
                             (Complex(0, 13198) + (21493 + 7303 * delta) * S1z +
                              21493 * S2z - 7303 * delta * S2z))) -
               3 * pow(r, 2) *
                   (-900 * pow(Nu, 3) *
                        (975 * pow(PhiDOT, 5) * pow(r, 5) +
                         Complex(0, 2382) * pow(PhiDOT, 4) * pow(r, 4) * rDOT -
                         2208 * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) -
                         Complex(0, 1592) * pow(PhiDOT, 2) * pow(r, 2) *
                             pow(rDOT, 3) +
                         688 * PhiDOT * r * pow(rDOT, 4) +
                         Complex(0, 128) * pow(rDOT, 5)) *
                        (S1z + S2z) +
                    1260 *
                        (61 * pow(PhiDOT, 5) * pow(r, 5) +
                         Complex(0, 151) * pow(PhiDOT, 4) * pow(r, 4) * rDOT -
                         172 * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) -
                         Complex(0, 122) * pow(PhiDOT, 2) * pow(r, 2) *
                             pow(rDOT, 3) +
                         48 * PhiDOT * r * pow(rDOT, 4) +
                         Complex(0, 8) * pow(rDOT, 5)) *
                        (S1z + delta * S1z + S2z - delta * S2z) -
                    3 * Nu *
                        (Complex(0, 1) * pow(PhiDOT, 4) * pow(r, 4) * rDOT *
                             (Complex(0, 56511) +
                              10 * (43269 + 36265 * delta) * S1z +
                              432690 * S2z - 362650 * delta * S2z) -
                         Complex(0, 4) * pow(PhiDOT, 2) * pow(r, 2) *
                             pow(rDOT, 3) *
                             (Complex(0, 14259) +
                              40 * (2154 + 1675 * delta) * S1z + 86160 * S2z -
                              67000 * delta * S2z) +
                         4 * pow(PhiDOT, 5) * pow(r, 5) *
                             (Complex(0, 94823) +
                              25 * (1494 + 1901 * delta) * S1z + 37350 * S2z -
                              47525 * delta * S2z) -
                         16 * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) *
                             (Complex(0, 6984) +
                              5 * (5967 + 4855 * delta) * S1z + 29835 * S2z -
                              24275 * delta * S2z) +
                         Complex(0, 8) * pow(rDOT, 5) *
                             (Complex(0, 251) + 10 * (303 + 215 * delta) * S1z +
                              3030 * S2z - 2150 * delta * S2z) +
                         64 * PhiDOT * r * pow(rDOT, 4) *
                             (Complex(0, 251) + 5 * (438 + 325 * delta) * S1z +
                              2190 * S2z - 1625 * delta * S2z)) +
                    pow(Nu, 2) *
                        (Complex(0, 8) * pow(rDOT, 5) *
                             (Complex(0, 2449) + 600 * (34 + 11 * delta) * S1z -
                              600 * (-34 + 11 * delta) * S2z) +
                         64 * PhiDOT * r * pow(rDOT, 4) *
                             (Complex(0, 2449) + 75 * (188 + 67 * delta) * S1z -
                              75 * (-188 + 67 * delta) * S2z) -
                         16 * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) *
                             (Complex(0, 66866) +
                              75 * (2423 + 1039 * delta) * S1z -
                              75 * (-2423 + 1039 * delta) * S2z) -
                         Complex(0, 4) * pow(PhiDOT, 2) * pow(r, 2) *
                             pow(rDOT, 3) *
                             (Complex(0, 139541) +
                              150 * (3551 + 1399 * delta) * S1z -
                              150 * (-3551 + 1399 * delta) * S2z) +
                         Complex(0, 1) * pow(PhiDOT, 4) * pow(r, 4) * rDOT *
                             (Complex(0, 515089) +
                              600 * (4703 + 2041 * delta) * S1z -
                              600 * (-4703 + 2041 * delta) * S2z) +
                         4 * pow(PhiDOT, 5) * pow(r, 5) *
                             (Complex(0, 925502) +
                              75 * (2966 + 2521 * delta) * S1z -
                              75 * (-2966 + 2521 * delta) * S2z))) +
               12 * mass * r *
                   (30 * pow(Nu, 3) *
                        (20149 * pow(PhiDOT, 3) * pow(r, 3) +
                         Complex(0, 33176) * pow(PhiDOT, 2) * pow(r, 2) * rDOT -
                         18519 * PhiDOT * r * pow(rDOT, 2) -
                         Complex(0, 4092) * pow(rDOT, 3)) *
                        (S1z + S2z) -
                    42 *
                        (1167 * pow(PhiDOT, 3) * pow(r, 3) +
                         Complex(0, 2042) * pow(PhiDOT, 2) * pow(r, 2) * rDOT -
                         1131 * PhiDOT * r * pow(rDOT, 2) -
                         Complex(0, 236) * pow(rDOT, 3)) *
                        (S1z + delta * S1z + S2z - delta * S2z) +
                    pow(Nu, 2) *
                        (-2 * pow(PhiDOT, 3) * pow(r, 3) *
                             (Complex(0, 452503) +
                              5 * (72902 + 36559 * delta) * S1z + 364510 * S2z -
                              182795 * delta * S2z) +
                         3 * PhiDOT * r * pow(rDOT, 2) *
                             (Complex(0, 53821) +
                              30 * (8510 + 2941 * delta) * S1z + 255300 * S2z -
                              88230 * delta * S2z) -
                         Complex(0, 6) * pow(PhiDOT, 2) * pow(r, 2) * rDOT *
                             (Complex(0, 65461) +
                              30 * (7463 + 2831 * delta) * S1z + 223890 * S2z -
                              84930 * delta * S2z) +
                         Complex(0, 1) * pow(rDOT, 3) *
                             (Complex(0, 26693) +
                              20 * (8419 + 2675 * delta) * S1z + 168380 * S2z -
                              53500 * delta * S2z)) +
                    Nu * (pow(PhiDOT, 3) * pow(r, 3) *
                              (Complex(0, 280937) +
                               (329808 + 305330 * delta) * S1z +
                               (329808 - 305330 * delta) * S2z) -
                          Complex(0, 1) * pow(rDOT, 3) *
                              (Complex(0, 8281) +
                               10 * (7293 + 5153 * delta) * S1z +
                               (72930 - 51530 * delta) * S2z) -
                          3 * PhiDOT * r * pow(rDOT, 2) *
                              (Complex(0, 16877) +
                               90 * (1261 + 930 * delta) * S1z -
                               90 * (-1261 + 930 * delta) * S2z) +
                          Complex(0, 2) * pow(PhiDOT, 2) * pow(r, 2) * rDOT *
                              (Complex(0, 60821) +
                               24 * (12576 + 9775 * delta) * S1z -
                               24 * (-12576 + 9775 * delta) * S2z))))) /
             (30240. * sqrt(143) * pow(r, 4))));
  }

  else {
    return 0;
  }
}

/* static COMPLEX16 hQC_6_m_6(REAL8 Nu, UINT4 vpnorder, REAL8 x){

    if(vpnorder == 4){
        return((108/(5.*sqrt(143)) - (108*Nu)/sqrt(143) +
(108*pow(Nu,2))/sqrt(143))* pow(x,3));
    }
    else if(vpnorder == 6){
        return((-6102/(35.*sqrt(143)) + (378*sqrt(1.1818181818181819)*Nu)/5. -
        (6912*pow(Nu,2))/(5.*sqrt(143)) +
        (162*sqrt(1.1818181818181819)*pow(Nu,3))/5.)*pow(x,4));
    }
    else{
        return 0;
    }
} */

static COMPLEX16 hl_6_m_6(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_6_m_6: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_6_m_6(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z,
                      S2z) /* +hQC_6_m_6(Nu,vpnorder,x) */) *
           cpolar(1, -6 * Phi);
  }
}

static COMPLEX16 hl_6_m_min6(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_6_m_min6: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) * pow(-1, 6) *
           conj(hGO_6_m_6(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z,
                          S2z) /* +hQC_6_m_6(Nu,vpnorder,x) */) *
           cpolar(1, 6 * Phi);
  }
}

// 65
static COMPLEX16 hGO_6_m_5(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  if (vpnorder == 5) {
    return (Complex(0, 0.003968253968253968) * delta * mass *
            (1 - 4 * Nu + 3 * pow(Nu, 2)) * PhiDOT *
            (82 * pow(mass, 2) +
             2 * mass * r *
                 (701 * pow(PhiDOT, 2) * pow(r, 2) +
                  Complex(0, 343) * PhiDOT * r * rDOT - 62 * pow(rDOT, 2)) +
             3 * pow(r, 2) *
                 (547 * pow(PhiDOT, 4) * pow(r, 4) +
                  Complex(0, 364) * pow(PhiDOT, 3) * pow(r, 3) * rDOT -
                  180 * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) -
                  Complex(0, 56) * PhiDOT * r * pow(rDOT, 3) +
                  8 * pow(rDOT, 4)))) /
           (sqrt(429) * pow(r, 2));
  }

  else if (vpnorder == 6) {

    return (/* (Complex(0,-21.70138888888889)*Nu*(S1z + delta*(-1 + Nu)*S1z
      - 3*Nu*S1z - (1 + delta)*S2z + (3 + delta)*Nu*S2z)*pow(x,4))/sqrt(429)
      + */
            /* Henry et al. ecc spin terms */ (
                (Complex(0, -0.006944444444444444) * pow(mass, 2) * Nu *
                 (82 * pow(mass, 2) +
                  2 * mass * r *
                      (701 * pow(PhiDOT, 2) * pow(r, 2) +
                       Complex(0, 343) * PhiDOT * r * rDOT -
                       62 * pow(rDOT, 2)) +
                  3 * pow(r, 2) *
                      (547 * pow(PhiDOT, 4) * pow(r, 4) +
                       Complex(0, 364) * pow(PhiDOT, 3) * pow(r, 3) * rDOT -
                       180 * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) -
                       Complex(0, 56) * PhiDOT * r * pow(rDOT, 3) +
                       8 * pow(rDOT, 4))) *
                 ((1 + delta * (-1 + Nu) - 3 * Nu) * S1z +
                  (-1 + delta * (-1 + Nu) + 3 * Nu) * S2z)) /
                (sqrt(429) * pow(r, 4))));
  }

  else if (vpnorder == 7) {

    return (
        (Complex(0, 0.00003306878306878307) * delta * mass * PhiDOT *
         (16 * pow(mass, 3) *
              (-3893 + 16386 * Nu - 17175 * pow(Nu, 2) + 6922 * pow(Nu, 3)) +
          24 * pow(mass, 2) * r *
              ((-25472 + 127239 * Nu - 198605 * pow(Nu, 2) +
                117623 * pow(Nu, 3)) *
                   pow(PhiDOT, 2) * pow(r, 2) +
               Complex(0, 10) *
                   (1954 - 7090 * Nu + 1431 * pow(Nu, 2) + 5232 * pow(Nu, 3)) *
                   PhiDOT * r * rDOT -
               (6603 - 26446 * Nu + 16325 * pow(Nu, 2) + 7138 * pow(Nu, 3)) *
                   pow(rDOT, 2)) +
          12 * mass * pow(r, 2) *
              (2 *
                   (-23971 + 144172 * Nu - 256320 * pow(Nu, 2) +
                    127374 * pow(Nu, 3)) *
                   pow(PhiDOT, 4) * pow(r, 4) +
               Complex(0, 5) *
                   (40015 - 154144 * Nu + 84408 * pow(Nu, 2) +
                    41694 * pow(Nu, 3)) *
                   pow(PhiDOT, 3) * pow(r, 3) * rDOT -
               3 *
                   (46471 - 183392 * Nu + 109660 * pow(Nu, 2) +
                    47046 * pow(Nu, 3)) *
                   pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) -
               Complex(0, 5) *
                   (9713 - 38048 * Nu + 21348 * pow(Nu, 2) +
                    11562 * pow(Nu, 3)) *
                   PhiDOT * r * pow(rDOT, 3) +
               4 * (1811 - 7007 * Nu + 3585 * pow(Nu, 2) + 2511 * pow(Nu, 3)) *
                   pow(rDOT, 4)) -
          45 * pow(r, 3) *
              ((3177 - 1626 * Nu - 54862 * pow(Nu, 2) + 73376 * pow(Nu, 3)) *
                   pow(PhiDOT, 6) * pow(r, 6) +
               Complex(0, 2) *
                   (-5337 + 26646 * Nu - 34862 * pow(Nu, 2) +
                    11212 * pow(Nu, 3)) *
                   pow(PhiDOT, 5) * pow(r, 5) * rDOT +
               4 * (1917 - 7098 * Nu + 1187 * pow(Nu, 2) + 6278 * pow(Nu, 3)) *
                   pow(PhiDOT, 4) * pow(r, 4) * pow(rDOT, 2) +
               Complex(0, 8) *
                   (639 - 1986 * Nu - 1156 * pow(Nu, 2) + 3296 * pow(Nu, 3)) *
                   pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 3) -
               16 * (153 - 402 * Nu - 577 * pow(Nu, 2) + 1022 * pow(Nu, 3)) *
                   pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 4) -
               Complex(0, 16) *
                   (45 - 102 * Nu - 236 * pow(Nu, 2) + 352 * pow(Nu, 3)) *
                   PhiDOT * r * pow(rDOT, 5) +
               32 * (3 - 6 * Nu - 19 * pow(Nu, 2) + 26 * pow(Nu, 3)) *
                   pow(rDOT, 6)))) /
        (sqrt(429) * pow(r, 4)));
  }

  else {
    return 0;
  }
}

/* static COMPLEX16 hQC_6_m_5(REAL8 Nu, UINT4 vpnorder, REAL8 x){
    REAL8 delta = sqrt(1-4*Nu);

    if(vpnorder == 5){
        return(((Complex(0,12.40079365079365)*delta)/sqrt(429) -
        (Complex(0,49.6031746031746)*delta*Nu)/sqrt(429) +
        (Complex(0,37.20238095238095)*delta*pow(Nu,2))/sqrt(429))*
        pow(x,3.5));
    }
    else{
        return 0;
    }
} */

static COMPLEX16 hl_6_m_5(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_6_m_m5: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_6_m_5(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z,
                      S2z) /* +hQC_6_m_5(Nu,vpnorder,x) */) *
           cpolar(1, -5 * Phi);
  }
}

static COMPLEX16 hl_6_m_min5(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_6_m_min5: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) * pow(-1, 6) *
           conj(hGO_6_m_5(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z,
                          S2z) /* +hQC_6_m_5(Nu,vpnorder,x) */) *
           cpolar(1, 5 * Phi);
  }
}

// 64
static COMPLEX16 hGO_6_m_4(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  if (vpnorder == 4) {
    return (((1 - 5 * Nu + 5 * pow(Nu, 2)) *
             (-516 * pow(mass, 3) +
              360 * pow(r, 3) * (PhiDOT * r - Complex(0, 1) * rDOT) *
                  pow(PhiDOT * r + Complex(0, 1) * rDOT, 5) +
              pow(mass, 2) * r *
                  (-3587 * pow(PhiDOT, 2) * pow(r, 2) -
                   Complex(0, 5840) * PhiDOT * r * rDOT + 2418 * pow(rDOT, 2)) +
              15 * mass * pow(r, 2) *
                  (113 * pow(PhiDOT, 4) * pow(r, 4) -
                   Complex(0, 108) * pow(PhiDOT, 3) * pow(r, 3) * rDOT +
                   468 * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) +
                   Complex(0, 416) * PhiDOT * r * pow(rDOT, 3) -
                   120 * pow(rDOT, 4)))) /
            (1980. * sqrt(78) * pow(r, 3)));
  }

  else if (vpnorder == 6) {
    return (
        -(14 * pow(mass, 4) *
              (-5740 + 29361 * Nu - 33348 * pow(Nu, 2) + 7334 * pow(Nu, 3)) +
          7560 * (-1 + 9 * Nu - 26 * pow(Nu, 2) + 23 * pow(Nu, 3)) * pow(r, 4) *
              pow(PhiDOT * r - Complex(0, 1) * rDOT, 2) *
              pow(PhiDOT * r + Complex(0, 1) * rDOT, 6) +
          2 * pow(mass, 3) * r *
              ((-196625 + 1082991 * Nu - 1522164 * pow(Nu, 2) +
                618457 * pow(Nu, 3)) *
                   pow(PhiDOT, 2) * pow(r, 2) +
               Complex(0, 4) *
                   (-47984 + 275121 * Nu - 442540 * pow(Nu, 2) +
                    255850 * pow(Nu, 3)) *
                   PhiDOT * r * rDOT +
               14 *
                   (3614 - 21621 * Nu + 39684 * pow(Nu, 2) -
                    29332 * pow(Nu, 3)) *
                   pow(rDOT, 2)) +
          pow(mass, 2) * pow(r, 2) *
              ((133599 - 779681 * Nu + 1417087 * pow(Nu, 2) -
                1130416 * pow(Nu, 3)) *
                   pow(PhiDOT, 4) * pow(r, 4) +
               Complex(0, 8) *
                   (3849 + 4172 * Nu - 80290 * pow(Nu, 2) +
                    64435 * pow(Nu, 3)) *
                   pow(PhiDOT, 3) * pow(r, 3) * rDOT +
               (-226971 + 551047 * Nu + 2049124 * pow(Nu, 2) -
                3713857 * pow(Nu, 3)) *
                   pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) -
               Complex(0, 64) *
                   (3921 - 11389 * Nu - 27265 * pow(Nu, 2) +
                    58450 * pow(Nu, 3)) *
                   PhiDOT * r * pow(rDOT, 3) +
               42 *
                   (1867 - 5501 * Nu - 12824 * pow(Nu, 2) +
                    28137 * pow(Nu, 3)) *
                   pow(rDOT, 4)) +
          15 * mass * pow(r, 3) *
              ((2267 - 12733 * Nu + 13895 * pow(Nu, 2) + 6300 * pow(Nu, 3)) *
                   pow(PhiDOT, 6) * pow(r, 6) -
               Complex(0, 8) *
                   (908 - 2597 * Nu - 5873 * pow(Nu, 2) + 11809 * pow(Nu, 3)) *
                   pow(PhiDOT, 5) * pow(r, 5) * rDOT +
               (5241 + 10066 * Nu - 173159 * pow(Nu, 2) + 235382 * pow(Nu, 3)) *
                   pow(PhiDOT, 4) * pow(r, 4) * pow(rDOT, 2) +
               Complex(0, 4) *
                   (-1651 + 11312 * Nu - 25417 * pow(Nu, 2) +
                    20916 * pow(Nu, 3)) *
                   pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 3) +
               4 *
                   (3127 - 8386 * Nu - 23569 * pow(Nu, 2) +
                    45122 * pow(Nu, 3)) *
                   pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 4) +
               Complex(0, 32) *
                   (228 - 217 * Nu - 3871 * pow(Nu, 2) + 5803 * pow(Nu, 3)) *
                   PhiDOT * r * pow(rDOT, 5) -
               168 * (9 + 4 * Nu - 221 * pow(Nu, 2) + 308 * pow(Nu, 3)) *
                   pow(rDOT, 6))) /
        (27720. * sqrt(78) * pow(r, 4)));
  }

  else if (vpnorder == 7) {

    return (/* (256*sqrt(0.05128205128205128)*(-150*pow(Nu,3)*(S1z + S2z)
      + 14*(S1z + delta*S1z + S2z - delta*S2z) + 10*pow(Nu,2)*(6*delta*(S1z -
      S2z) + 23*(S1z + S2z)) - Nu*(65*delta*(S1z - S2z) + 103*(S1z +
      S2z)))*pow(x,4.5))/3465.
      + */
            /* Henry et al. ecc+spin terms */
            ((pow(mass, 2) *
              (-8 * pow(mass, 2) *
                   (30 * pow(Nu, 3) *
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
                    3 * pow(Nu, 2) *
                        (5 * PhiDOT * r *
                             (Complex(0, 26396) +
                              (27731 + 11513 * delta) * S1z +
                              (27731 - 11513 * delta) * S2z) +
                         Complex(0, 1) * rDOT *
                             (Complex(0, 45321) +
                              80 * (1592 + 485 * delta) * S1z -
                              80 * (-1592 + 485 * delta) * S2z))) -
               3 * pow(r, 2) *
                   (-60 * pow(Nu, 3) *
                        (401 * pow(PhiDOT, 5) * pow(r, 5) +
                         Complex(0, 30572) * pow(PhiDOT, 4) * pow(r, 4) * rDOT -
                         18000 * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) +
                         Complex(0, 3632) * pow(PhiDOT, 2) * pow(r, 2) *
                             pow(rDOT, 3) -
                         10736 * PhiDOT * r * pow(rDOT, 4) -
                         Complex(0, 3840) * pow(rDOT, 5)) *
                        (S1z + S2z) +
                    1260 *
                        (17 * pow(PhiDOT, 5) * pow(r, 5) +
                         Complex(0, 58) * pow(PhiDOT, 4) * pow(r, 4) * rDOT +
                         16 * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) +
                         Complex(0, 84) * pow(PhiDOT, 2) * pow(r, 2) *
                             pow(rDOT, 3) -
                         64 * PhiDOT * r * pow(rDOT, 4) -
                         Complex(0, 16) * pow(rDOT, 5)) *
                        (S1z + delta * S1z + S2z - delta * S2z) -
                    3 * Nu *
                        (Complex(0, 20) * pow(PhiDOT, 2) * pow(r, 2) *
                             pow(rDOT, 3) *
                             (Complex(0, 3109) + (9208 + 9384 * delta) * S1z +
                              (9208 - 9384 * delta) * S2z) -
                         32 * PhiDOT * r * pow(rDOT, 4) *
                             (Complex(0, 1004) +
                              5 * (1091 + 869 * delta) * S1z +
                              (5455 - 4345 * delta) * S2z) -
                         Complex(0, 8) * pow(rDOT, 5) *
                             (Complex(0, 753) + 20 * (303 + 215 * delta) * S1z +
                              (6060 - 4300 * delta) * S2z) +
                         20 * pow(PhiDOT, 5) * pow(r, 5) *
                             (Complex(0, 4690) + (652 + 3807 * delta) * S1z +
                              (652 - 3807 * delta) * S2z) +
                         16 * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) *
                             (Complex(0, 1872) +
                              5 * (-591 + 605 * delta) * S1z -
                              5 * (591 + 605 * delta) * S2z) +
                         Complex(0, 1) * pow(PhiDOT, 4) * pow(r, 4) * rDOT *
                             (Complex(0, 78819) +
                              20 * (12535 + 5539 * delta) * S1z -
                              20 * (-12535 + 5539 * delta) * S2z)) +
                    pow(Nu, 2) *
                        (Complex(0, 20) * pow(PhiDOT, 2) * pow(r, 2) *
                             pow(rDOT, 3) *
                             (Complex(0, 30451) +
                              12 * (2831 + 2487 * delta) * S1z +
                              (33972 - 29844 * delta) * S2z) -
                         64 * PhiDOT * r * pow(rDOT, 4) *
                             (Complex(0, 4898) +
                              15 * (1062 + 449 * delta) * S1z +
                              (15930 - 6735 * delta) * S2z) -
                         Complex(0, 24) * pow(rDOT, 5) *
                             (Complex(0, 2449) + 400 * (34 + 11 * delta) * S1z -
                              400 * (-34 + 11 * delta) * S2z) +
                         96 * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) *
                             (Complex(0, 2838) +
                              25 * (-445 + 79 * delta) * S1z -
                              25 * (445 + 79 * delta) * S2z) +
                         Complex(0, 1) * pow(PhiDOT, 4) * pow(r, 4) * rDOT *
                             (Complex(0, 761581) +
                              240 * (9441 + 1247 * delta) * S1z -
                              240 * (-9441 + 1247 * delta) * S2z) +
                         20 * pow(PhiDOT, 5) * pow(r, 5) *
                             (Complex(0, 44900) +
                              3 * (-1858 + 5829 * delta) * S1z -
                              3 * (1858 + 5829 * delta) * S2z))) +
               12 * mass * r *
                   (30 * pow(Nu, 3) *
                        (73 * pow(PhiDOT, 3) * pow(r, 3) -
                         Complex(0, 11024) * pow(PhiDOT, 2) * pow(r, 2) * rDOT +
                         20149 * PhiDOT * r * pow(rDOT, 2) +
                         Complex(0, 8184) * pow(rDOT, 3)) *
                        (S1z + S2z) +
                    42 *
                        (871 * pow(PhiDOT, 3) * pow(r, 3) +
                         Complex(0, 1684) * pow(PhiDOT, 2) * pow(r, 2) * rDOT -
                         1543 * PhiDOT * r * pow(rDOT, 2) -
                         Complex(0, 472) * pow(rDOT, 3)) *
                        (S1z + delta * S1z + S2z - delta * S2z) +
                    Nu * (-2 * pow(PhiDOT, 3) * pow(r, 3) *
                              (Complex(0, 69725) +
                               (89007 + 95170 * delta) * S1z +
                               (89007 - 95170 * delta) * S2z) +
                          Complex(0, 1) * pow(rDOT, 3) *
                              (Complex(0, 24843) +
                               20 * (7293 + 5153 * delta) * S1z -
                               20 * (-7293 + 5153 * delta) * S2z) -
                          Complex(0, 2) * pow(PhiDOT, 2) * pow(r, 2) * rDOT *
                              (Complex(0, 66847) +
                               8 * (24951 + 24245 * delta) * S1z -
                               8 * (-24951 + 24245 * delta) * S2z) +
                          2 * PhiDOT * r * pow(rDOT, 2) *
                              (Complex(0, 50631) +
                               5 * (43338 + 34327 * delta) * S1z -
                               5 * (-43338 + 34327 * delta) * S2z)) +
                    pow(Nu, 2) *
                        (10 * pow(PhiDOT, 3) * pow(r, 3) *
                             (Complex(0, 44862) +
                              (16586 + 19777 * delta) * S1z +
                              (16586 - 19777 * delta) * S2z) +
                         Complex(0, 2) * pow(PhiDOT, 2) * pow(r, 2) * rDOT *
                             (Complex(0, 216321) +
                              100 * (3079 + 2111 * delta) * S1z -
                              100 * (-3079 + 2111 * delta) * S2z) -
                         Complex(0, 1) * pow(rDOT, 3) *
                             (Complex(0, 80079) +
                              40 * (8419 + 2675 * delta) * S1z -
                              40 * (-8419 + 2675 * delta) * S2z) -
                         2 * PhiDOT * r * pow(rDOT, 2) *
                             (Complex(0, 161463) +
                              5 * (89002 + 36251 * delta) * S1z -
                              5 * (-89002 + 36251 * delta) * S2z))))) /
             (166320. * sqrt(78) * pow(r, 4))));
  }

  else {
    return 0;
  }
}

/* static COMPLEX16 hQC_6_m_4(REAL8 Nu, UINT4 vpnorder, REAL8 x){

    if(vpnorder == 4){
        return(((-256*sqrt(0.05128205128205128))/495. +
        (256*sqrt(0.05128205128205128)*Nu)/99. -
        (256*sqrt(0.05128205128205128)*pow(Nu,2))/99.)*pow(x,3));
    }
    else if(vpnorder == 6){
        return(((3968*sqrt(0.05128205128205128))/1155. -
        (9088*sqrt(0.05128205128205128)*Nu)/495. +
        (1024*sqrt(0.05128205128205128)*pow(Nu,2))/45. -
        (2432*sqrt(0.05128205128205128)*pow(Nu,3))/495.)*pow(x,4));
    }
    else{
        return 0;
    }

} */

static COMPLEX16 hl_6_m_4(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_6_m_4: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_6_m_4(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z,
                      S2z) /* +hQC_6_m_4(Nu,vpnorder,x) */) *
           cpolar(1, -4 * Phi);
  }
}

static COMPLEX16 hl_6_m_min4(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_6_m_min4: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) * pow(-1, 6) *
           conj(hGO_6_m_4(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z,
                          S2z) /* +hQC_6_m_4(Nu,vpnorder,x) */) *
           cpolar(1, 4 * Phi);
  }
}

// 63
static COMPLEX16 hGO_6_m_3(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  if (vpnorder == 5) {
    return (Complex(0, -0.00036075036075036075) * delta * mass *
            (1 - 4 * Nu + 3 * pow(Nu, 2)) * PhiDOT *
            (410 * pow(mass, 2) +
             2 * mass * r *
                 (929 * pow(PhiDOT, 2) * pow(r, 2) +
                  Complex(0, 1029) * PhiDOT * r * rDOT - 310 * pow(rDOT, 2)) -
             3 * pow(r, 2) *
                 (513 * pow(PhiDOT, 4) * pow(r, 4) +
                  Complex(0, 28) * pow(PhiDOT, 3) * pow(r, 3) * rDOT +
                  228 * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) +
                  Complex(0, 168) * PhiDOT * r * pow(rDOT, 3) -
                  40 * pow(rDOT, 4)))) /
           (sqrt(65) * pow(r, 2));
  }

  else if (vpnorder == 6) {

    return (/* (Complex(0,0.4602272727272727)*Nu*(S1z + delta*(-1 + Nu)*S1z -
      3*Nu*S1z
      - (1 + delta)*S2z + (3 + delta)*Nu*S2z)*pow(x,4))/sqrt(65)
      + */
            /* Henry et al. ecc spin terms */ (
                (Complex(0, 0.0006313131313131314) * pow(mass, 2) * Nu *
                 (410 * pow(mass, 2) +
                  2 * mass * r *
                      (929 * pow(PhiDOT, 2) * pow(r, 2) +
                       Complex(0, 1029) * PhiDOT * r * rDOT -
                       310 * pow(rDOT, 2)) -
                  3 * pow(r, 2) *
                      (513 * pow(PhiDOT, 4) * pow(r, 4) +
                       Complex(0, 28) * pow(PhiDOT, 3) * pow(r, 3) * rDOT +
                       228 * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) +
                       Complex(0, 168) * PhiDOT * r * pow(rDOT, 3) -
                       40 * pow(rDOT, 4))) *
                 ((1 + delta * (-1 + Nu) - 3 * Nu) * S1z +
                  (-1 + delta * (-1 + Nu) + 3 * Nu) * S2z)) /
                (sqrt(65) * pow(r, 4))));
  }

  else if (vpnorder == 7) {

    return (
        (Complex(0, -0.000015031265031265032) * delta * mass * PhiDOT *
         (16 * pow(mass, 3) *
              (-3893 + 16386 * Nu - 17175 * pow(Nu, 2) + 6922 * pow(Nu, 3)) +
          8 * pow(mass, 2) * r *
              ((-22256 + 111525 * Nu - 170247 * pow(Nu, 2) +
                94453 * pow(Nu, 3)) *
                   pow(PhiDOT, 2) * pow(r, 2) +
               Complex(0, 18) *
                   (1954 - 7090 * Nu + 1431 * pow(Nu, 2) + 5232 * pow(Nu, 3)) *
                   PhiDOT * r * rDOT -
               3 *
                   (6603 - 26446 * Nu + 16325 * pow(Nu, 2) +
                    7138 * pow(Nu, 3)) *
                   pow(rDOT, 2)) -
          12 * mass * pow(r, 2) *
              (2 * (771 + 6004 * Nu - 44896 * pow(Nu, 2) + 48978 * pow(Nu, 3)) *
                   pow(PhiDOT, 4) * pow(r, 4) +
               Complex(0, 1) *
                   (10771 - 40224 * Nu + 14456 * pow(Nu, 2) +
                    21414 * pow(Nu, 3)) *
                   pow(PhiDOT, 3) * pow(r, 3) * rDOT +
               (34437 - 136928 * Nu + 85844 * pow(Nu, 2) + 30834 * pow(Nu, 3)) *
                   pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) +
               Complex(0, 3) *
                   (9713 - 38048 * Nu + 21348 * pow(Nu, 2) +
                    11562 * pow(Nu, 3)) *
                   PhiDOT * r * pow(rDOT, 3) -
               4 * (1811 - 7007 * Nu + 3585 * pow(Nu, 2) + 2511 * pow(Nu, 3)) *
                   pow(rDOT, 4)) +
          9 * pow(r, 3) *
              ((4339 - 4350 * Nu - 50298 * pow(Nu, 2) + 61600 * pow(Nu, 3)) *
                   pow(PhiDOT, 6) * pow(r, 6) +
               Complex(0, 2) *
                   (-12661 + 53118 * Nu - 45494 * pow(Nu, 2) +
                    2652 * pow(Nu, 3)) *
                   pow(PhiDOT, 5) * pow(r, 5) * rDOT +
               12 * (885 - 2346 * Nu - 3253 * pow(Nu, 2) + 5846 * pow(Nu, 3)) *
                   pow(PhiDOT, 4) * pow(r, 4) * pow(rDOT, 2) +
               Complex(0, 8) *
                   (211 - 90 * Nu - 2692 * pow(Nu, 2) + 2880 * pow(Nu, 3)) *
                   pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 3) +
               16 * (181 - 522 * Nu - 493 * pow(Nu, 2) + 1062 * pow(Nu, 3)) *
                   pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 4) +
               Complex(0, 48) *
                   (45 - 102 * Nu - 236 * pow(Nu, 2) + 352 * pow(Nu, 3)) *
                   PhiDOT * r * pow(rDOT, 5) -
               160 * (3 - 6 * Nu - 19 * pow(Nu, 2) + 26 * pow(Nu, 3)) *
                   pow(rDOT, 6)))) /
        (sqrt(65) * pow(r, 3)));
  }

  else {
    return 0;
  }
}

/* static COMPLEX16 hQC_6_m_3(REAL8 Nu, UINT4 vpnorder, REAL8 x){
    REAL8 delta = sqrt(1-4*Nu);

    if(vpnorder == 5){
        return((Complex(0,-0.262987012987013)*delta*(1 - 4*Nu + 3*pow(Nu,2))*
        pow(x,3.5))/sqrt(65));
    }
    else{
        return 0;
    }
} */

static COMPLEX16 hl_6_m_3(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_6_m_3: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_6_m_3(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z,
                      S2z) /* +hQC_6_m_3(Nu,vpnorder,x) */) *
           cpolar(1, -3 * Phi);
  }
}

static COMPLEX16 hl_6_m_min3(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_6_m_min3: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) * pow(-1, 6) *
           conj(hGO_6_m_3(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z,
                          S2z) /* +hQC_6_m_3(Nu,vpnorder,x) */) *
           cpolar(1, 3 * Phi);
  }
}

// 62
static COMPLEX16 hGO_6_m_2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  if (vpnorder == 4) {
    return ((1 - 5 * Nu + 5 * pow(Nu, 2)) *
            (516 * pow(mass, 3) +
             360 * pow(r, 3) * pow(PhiDOT * r - Complex(0, 1) * rDOT, 2) *
                 pow(PhiDOT * r + Complex(0, 1) * rDOT, 4) +
             pow(mass, 2) * r *
                 (-145 * pow(PhiDOT, 2) * pow(r, 2) +
                  Complex(0, 2920) * PhiDOT * r * rDOT - 2418 * pow(rDOT, 2)) -
             3 * mass * pow(r, 2) *
                 (233 * pow(PhiDOT, 4) * pow(r, 4) +
                  Complex(0, 1050) * pow(PhiDOT, 3) * pow(r, 3) * rDOT -
                  252 * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) +
                  Complex(0, 1040) * PhiDOT * r * pow(rDOT, 3) -
                  600 * pow(rDOT, 4)))) /
           (2376. * sqrt(65) * pow(r, 3));
  }

  else if (vpnorder == 6) {
    return (14 * pow(mass, 4) *
                (-5740 + 29361 * Nu - 33348 * pow(Nu, 2) + 7334 * pow(Nu, 3)) -
            7560 * (-1 + 9 * Nu - 26 * pow(Nu, 2) + 23 * pow(Nu, 3)) *
                pow(r, 4) * pow(PhiDOT * r - Complex(0, 1) * rDOT, 3) *
                pow(PhiDOT * r + Complex(0, 1) * rDOT, 5) +
            2 * pow(mass, 3) * r *
                ((9187 - 37401 * Nu + 16548 * pow(Nu, 2) + 2821 * pow(Nu, 3)) *
                     pow(PhiDOT, 2) * pow(r, 2) +
                 Complex(0, 2) *
                     (-47984 + 275121 * Nu - 442540 * pow(Nu, 2) +
                      255850 * pow(Nu, 3)) *
                     PhiDOT * r * rDOT +
                 14 *
                     (3614 - 21621 * Nu + 39684 * pow(Nu, 2) -
                      29332 * pow(Nu, 3)) *
                     pow(rDOT, 2)) +
            pow(mass, 2) * pow(r, 2) *
                ((54699 - 336749 * Nu + 596995 * pow(Nu, 2) -
                  337960 * pow(Nu, 3)) *
                     pow(PhiDOT, 4) * pow(r, 4) -
                 Complex(0, 4) *
                     (-5781 + 89572 * Nu - 379358 * pow(Nu, 2) +
                      444689 * pow(Nu, 3)) *
                     pow(PhiDOT, 3) * pow(r, 3) * rDOT +
                 (-9351 + 101899 * Nu - 419300 * pow(Nu, 2) +
                  566195 * pow(Nu, 3)) *
                     pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) -
                 Complex(0, 32) *
                     (3921 - 11389 * Nu - 27265 * pow(Nu, 2) +
                      58450 * pow(Nu, 3)) *
                     PhiDOT * r * pow(rDOT, 3) +
                 42 *
                     (1867 - 5501 * Nu - 12824 * pow(Nu, 2) +
                      28137 * pow(Nu, 3)) *
                     pow(rDOT, 4)) +
            3 * mass * pow(r, 3) *
                ((-7885 + 64211 * Nu - 170905 * pow(Nu, 2) +
                  146580 * pow(Nu, 3)) *
                     pow(PhiDOT, 6) * pow(r, 6) +
                 Complex(0, 4) *
                     (1438 + 9779 * Nu - 86023 * pow(Nu, 2) +
                      109949 * pow(Nu, 3)) *
                     pow(PhiDOT, 5) * pow(r, 5) * rDOT +
                 (16353 - 68054 * Nu + 10297 * pow(Nu, 2) +
                  77294 * pow(Nu, 3)) *
                     pow(PhiDOT, 4) * pow(r, 4) * pow(rDOT, 2) +
                 Complex(0, 2) *
                     (14341 - 392 * Nu - 316841 * pow(Nu, 2) +
                      452508 * pow(Nu, 3)) *
                     pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 3) -
                 4 *
                     (13 + 12530 * Nu - 68803 * pow(Nu, 2) +
                      80654 * pow(Nu, 3)) *
                     pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 4) +
                 Complex(0, 80) *
                     (228 - 217 * Nu - 3871 * pow(Nu, 2) + 5803 * pow(Nu, 3)) *
                     PhiDOT * r * pow(rDOT, 5) -
                 840 * (9 + 4 * Nu - 221 * pow(Nu, 2) + 308 * pow(Nu, 3)) *
                     pow(rDOT, 6))) /
           (33264. * sqrt(65) * pow(r, 4));
  }

  else if (vpnorder == 7) {

    return (/* (4*(174*pow(Nu,3)*(S1z + S2z) - 14*(S1z + delta*S1z + S2z -
      delta*S2z)
      + Nu*(53*delta*(S1z - S2z) + 115*(S1z + S2z)) - 2*pow(Nu,2)*(18*delta*(S1z
      - S2z)
      + 139*(S1z + S2z)))*pow(x,4.5))/(2079.*sqrt(65))
      + */
            /*Henry et al. ecc+spin terms */
            ((pow(mass, 2) *
              (40 * pow(mass, 2) *
                   (6 * pow(Nu, 3) *
                        (3489 * PhiDOT * r + Complex(0, 24880) * rDOT) *
                        (S1z + S2z) -
                    21 * (683 * PhiDOT * r + Complex(0, 476) * rDOT) *
                        (S1z + delta * S1z + S2z - delta * S2z) -
                    3 * pow(Nu, 2) *
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
                   (30 * pow(Nu, 3) *
                        (9893 * pow(PhiDOT, 3) * pow(r, 3) +
                         Complex(0, 55432) * pow(PhiDOT, 2) * pow(r, 2) * rDOT -
                         5479 * PhiDOT * r * pow(rDOT, 2) +
                         Complex(0, 20460) * pow(rDOT, 3)) *
                        (S1z + S2z) -
                    210 *
                        (67 * pow(PhiDOT, 3) * pow(r, 3) -
                         Complex(0, 122) * pow(PhiDOT, 2) * pow(r, 2) * rDOT +
                         433 * PhiDOT * r * pow(rDOT, 2) +
                         Complex(0, 236) * pow(rDOT, 3)) *
                        (S1z + delta * S1z + S2z - delta * S2z) +
                    Nu * (pow(PhiDOT, 3) * pow(r, 3) *
                              (Complex(0, 285011) +
                               10 * (10608 + 10973 * delta) * S1z +
                               (106080 - 109730 * delta) * S2z) +
                          5 * PhiDOT * r * pow(rDOT, 2) *
                              (Complex(0, 50631) +
                               (80562 + 97252 * delta) * S1z +
                               (80562 - 97252 * delta) * S2z) +
                          Complex(0, 5) * pow(rDOT, 3) *
                              (Complex(0, 24843) +
                               10 * (7293 + 5153 * delta) * S1z +
                               (72930 - 51530 * delta) * S2z) -
                          Complex(0, 2) * pow(PhiDOT, 2) * pow(r, 2) * rDOT *
                              (Complex(0, -12613) +
                               40 * (-2676 + 1801 * delta) * S1z -
                               40 * (2676 + 1801 * delta) * S2z)) -
                    pow(Nu, 2) *
                        (5 * PhiDOT * r * pow(rDOT, 2) *
                             (Complex(0, 161463) +
                              2 * (22706 + 51787 * delta) * S1z +
                              (45412 - 103574 * delta) * S2z) +
                         2 * pow(PhiDOT, 3) * pow(r, 3) *
                             (Complex(0, 460269) +
                              5 * (28838 + 14911 * delta) * S1z +
                              (144190 - 74555 * delta) * S2z) +
                         Complex(0, 5) * pow(rDOT, 3) *
                             (Complex(0, 80079) +
                              20 * (8419 + 2675 * delta) * S1z -
                              20 * (-8419 + 2675 * delta) * S2z) -
                         Complex(0, 2) * pow(PhiDOT, 2) * pow(r, 2) * rDOT *
                             (Complex(0, -36879) +
                              10 * (-78341 + 8003 * delta) * S1z -
                              10 * (78341 + 8003 * delta) * S2z))) +
               3 * pow(r, 2) *
                   (60 * pow(Nu, 3) *
                        (11215 * pow(PhiDOT, 5) * pow(r, 5) +
                         Complex(0, 57242) * pow(PhiDOT, 4) * pow(r, 4) * rDOT +
                         28128 * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) +
                         Complex(0, 57112) * pow(PhiDOT, 2) * pow(r, 2) *
                             pow(rDOT, 3) -
                         6992 * PhiDOT * r * pow(rDOT, 4) +
                         Complex(0, 9600) * pow(rDOT, 5)) *
                        (S1z + S2z) +
                    6300 *
                        (9 * pow(PhiDOT, 5) * pow(r, 5) +
                         Complex(0, 9) * pow(PhiDOT, 4) * pow(r, 4) * rDOT -
                         28 * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) -
                         Complex(0, 6) * pow(PhiDOT, 2) * pow(r, 2) *
                             pow(rDOT, 3) -
                         16 * PhiDOT * r * pow(rDOT, 4) -
                         Complex(0, 8) * pow(rDOT, 5)) *
                        (S1z + delta * S1z + S2z - delta * S2z) +
                    pow(Nu, 2) *
                        (Complex(0, -4) * pow(PhiDOT, 2) * pow(r, 2) *
                             pow(rDOT, 3) *
                             (Complex(0, 37829) +
                              30 * (30617 + 1089 * delta) * S1z +
                              (918510 - 32670 * delta) * S2z) -
                         48 * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) *
                             (Complex(0, 44162) +
                              5 * (10459 + 3923 * delta) * S1z +
                              (52295 - 19615 * delta) * S2z) -
                         Complex(0, 120) * pow(rDOT, 5) *
                             (Complex(0, 2449) + 200 * (34 + 11 * delta) * S1z -
                              200 * (-34 + 11 * delta) * S2z) +
                         4 * pow(PhiDOT, 5) * pow(r, 5) *
                             (Complex(0, -354934) +
                              75 * (-650 + 297 * delta) * S1z -
                              75 * (650 + 297 * delta) * S2z) -
                         320 * PhiDOT * r * pow(rDOT, 4) *
                             (Complex(0, 2449) + 3 * (36 + 577 * delta) * S1z -
                              3 * (-36 + 577 * delta) * S2z) +
                         Complex(0, 1) * pow(PhiDOT, 4) * pow(r, 4) * rDOT *
                             (Complex(0, 880129) +
                              120 * (-28647 + 4751 * delta) * S1z -
                              120 * (28647 + 4751 * delta) * S2z)) -
                    3 * Nu *
                        (Complex(0, -40) * pow(rDOT, 5) *
                             (Complex(0, 753) + 10 * (303 + 215 * delta) * S1z +
                              (3030 - 2150 * delta) * S2z) -
                         320 * PhiDOT * r * pow(rDOT, 4) *
                             (Complex(0, 251) + (422 + 551 * delta) * S1z +
                              (422 - 551 * delta) * S2z) -
                         Complex(0, 4) * pow(PhiDOT, 2) * pow(r, 2) *
                             pow(rDOT, 3) *
                             (Complex(0, 3971) +
                              40 * (1858 + 333 * delta) * S1z -
                              40 * (-1858 + 333 * delta) * S2z) +
                         4 * pow(PhiDOT, 5) * pow(r, 5) *
                             (Complex(0, -37091) +
                              5 * (3454 + 3105 * delta) * S1z -
                              5 * (-3454 + 3105 * delta) * S2z) -
                         16 * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) *
                             (Complex(0, 13464) +
                              5 * (5007 + 3799 * delta) * S1z -
                              5 * (-5007 + 3799 * delta) * S2z) +
                         Complex(0, 1) * pow(PhiDOT, 4) * pow(r, 4) * rDOT *
                             (Complex(0, 94671) +
                              10 * (-16313 + 14227 * delta) * S1z -
                              10 * (16313 + 14227 * delta) * S2z))))) /
             (997920. * sqrt(65) * pow(r, 4))));
  }

  else {
    return 0;
  }
}

/* static COMPLEX16 hQC_6_m_2(REAL8 Nu, UINT4 vpnorder, REAL8 x){

    if(vpnorder == 4){
        return((4/(297.*sqrt(65)) - (4*sqrt(0.38461538461538464)*Nu)/297. +
        (4*sqrt(0.38461538461538464)*pow(Nu,2))/297.)*pow(x,3));
    }
    else if(vpnorder == 6){
        return((-6/(77.*sqrt(65)) + (118*Nu)/(297.*sqrt(65)) -
        (128*pow(Nu,2))/(297.*sqrt(65)) +
        (14*pow(Nu,3))/(297.*sqrt(65)))*pow(x,4));
    }
    else{
        return 0;
    }
} */

static COMPLEX16 hl_6_m_2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_6_m_1: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_6_m_2(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z,
                      S2z) /* +hQC_6_m_2(Nu,vpnorder,x) */) *
           cpolar(1, -2 * Phi);
  }
}

static COMPLEX16 hl_6_m_min2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_6_m_min2: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) * pow(-1, 6) *
           conj(hGO_6_m_2(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z,
                          S2z) /* +hQC_6_m_2(Nu,vpnorder,x) */) *
           cpolar(1, 2 * Phi);
  }
}

// 61
static COMPLEX16 hGO_6_m_1(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  if (vpnorder == 5) {
    return (Complex(0, 0.0002405002405002405) * delta * mass *
            (1 - 4 * Nu + 3 * pow(Nu, 2)) * PhiDOT *
            (410 * pow(mass, 2) -
             2 * mass * r *
                 (359 * pow(PhiDOT, 2) * pow(r, 2) -
                  Complex(0, 343) * PhiDOT * r * rDOT + 310 * pow(rDOT, 2)) +
             3 * pow(r, 2) *
                 (103 * pow(PhiDOT, 4) * pow(r, 4) -
                  Complex(0, 196) * pow(PhiDOT, 3) * pow(r, 3) * rDOT +
                  108 * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) -
                  Complex(0, 56) * PhiDOT * r * pow(rDOT, 3) +
                  40 * pow(rDOT, 4)))) /
           (sqrt(26) * pow(r, 2));
  }

  else if (vpnorder == 6) {

    return (/* (Complex(0,-0.00042087542087542086)*Nu*(S1z + delta*(-1 + Nu)*S1z
      - 3*Nu*S1z
      - (1 + delta)*S2z + (3 + delta)*Nu*S2z)*pow(x,4))/sqrt(26)
      + */
            /* Henry et al. ecc spin terms */ (
                (Complex(0, -0.00042087542087542086) * pow(mass, 2) * Nu *
                 (410 * pow(mass, 2) -
                  2 * mass * r *
                      (359 * pow(PhiDOT, 2) * pow(r, 2) -
                       Complex(0, 343) * PhiDOT * r * rDOT +
                       310 * pow(rDOT, 2)) +
                  3 * pow(r, 2) *
                      (103 * pow(PhiDOT, 4) * pow(r, 4) -
                       Complex(0, 196) * pow(PhiDOT, 3) * pow(r, 3) * rDOT +
                       108 * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) -
                       Complex(0, 56) * PhiDOT * r * pow(rDOT, 3) +
                       40 * pow(rDOT, 4))) *
                 ((1 + delta * (-1 + Nu) - 3 * Nu) * S1z +
                  (-1 + delta * (-1 + Nu) + 3 * Nu) * S2z)) /
                (sqrt(26) * pow(r, 4))));
  }

  else if (vpnorder == 7) {

    return ((delta * mass * PhiDOT *
             (Complex(0, 16) * pow(mass, 3) *
                  (-3893 + Nu * (16386 + Nu * (-17175 + 6922 * Nu))) +
              24 * pow(mass, 2) * r *
                  (Complex(0, -1) *
                       (-1608 + Nu * (7857 + Nu * (-14179 + 11585 * Nu))) *
                       pow(PhiDOT, 2) * pow(r, 2) -
                   2 * (1954 + Nu * (-7090 + 3 * Nu * (477 + 1744 * Nu))) *
                       PhiDOT * r * rDOT -
                   Complex(0, 1) *
                       (6603 + Nu * (-26446 + Nu * (16325 + 7138 * Nu))) *
                       pow(rDOT, 2)) +
              12 * mass * pow(r, 2) *
                  (Complex(0, 2) *
                       (973 + 2 * Nu * (-898 + Nu * (-4168 + 6015 * Nu))) *
                       pow(PhiDOT, 4) * pow(r, 4) +
                   (25393 + 2 * Nu * (-48592 + Nu * (24716 + 15777 * Nu))) *
                       pow(PhiDOT, 3) * pow(r, 3) * rDOT +
                   Complex(0, 3) *
                       (6017 + 2 * Nu * (-11616 + Nu * (5954 + 4053 * Nu))) *
                       pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) +
                   (9713 + 2 * Nu * (-19024 + 3 * Nu * (3558 + 1927 * Nu))) *
                       PhiDOT * r * pow(rDOT, 3) +
                   Complex(0, 4) *
                       (1811 + Nu * (-7007 + 3 * Nu * (1195 + 837 * Nu))) *
                       pow(rDOT, 4)) +
              9 * pow(r, 3) *
                  (Complex(0, -1) *
                       (781 + 2 * Nu * (-81 + Nu * (-5699 + 6336 * Nu))) *
                       pow(PhiDOT, 6) * pow(r, 6) -
                   2 * (8025 + 2 * Nu * (-15003 + Nu * (6791 + 5258 * Nu))) *
                       pow(PhiDOT, 5) * pow(r, 5) * rDOT -
                   Complex(0, 844) * (3 + Nu * (-6 + Nu * (-19 + 26 * Nu))) *
                       pow(PhiDOT, 4) * pow(r, 4) * pow(rDOT, 2) +
                   8 * (-425 + 2 * Nu * (519 + 2 * (481 - 772 * Nu) * Nu)) *
                       pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 3) -
                   Complex(0, 592) * (3 + Nu * (-6 + Nu * (-19 + 26 * Nu))) *
                       pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 4) +
                   16 * (-45 + 2 * Nu * (51 + 2 * (59 - 88 * Nu) * Nu)) *
                       PhiDOT * r * pow(rDOT, 5) -
                   Complex(0, 160) * (3 + Nu * (-6 + Nu * (-19 + 26 * Nu))) *
                       pow(rDOT, 6)))) /
            (99792. * sqrt(26) * pow(r, 3)));
  }

  else {
    return 0;
  }
}

/* static COMPLEX16 hQC_6_m_1(REAL8 Nu, UINT4 vpnorder, REAL8 x){
    REAL8 delta = sqrt(1-4*Nu);

    if(vpnorder == 5){
        return((Complex(0,0.0002405002405002405)*delta*(1 - 4*Nu + 3*pow(Nu,2))*
        pow(x,3.5))/sqrt(26));
    }
    else{
        return 0;
    }
} */

static COMPLEX16 hl_6_m_1(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_6_m_1: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           (hGO_6_m_1(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z,
                      S2z) /* +hQC_6_m_1(Nu,vpnorder,x) */) *
           cpolar(1, -1 * Phi);
  }
}

static COMPLEX16 hl_6_m_min1(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_6_m_min1: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) * pow(-1, 6) *
           conj(hGO_6_m_1(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z,
                          S2z) /* +hQC_6_m_1(Nu,vpnorder,x) */) *
           cpolar(1, 1 * Phi);
  }
}

// 77
static COMPLEX16 hGO_7_m_7(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  if (vpnorder == 5) {
    return (delta * (1 - 4 * Nu + 3 * pow(Nu, 2)) *
            (720 * pow(r, 3) * pow(Complex(0, 1) * PhiDOT * r - rDOT, 7) +
             2 * pow(mass, 3) * (Complex(0, -4559) * PhiDOT * r + 1976 * rDOT) +
             18 * pow(mass, 2) * r *
                 (Complex(0, -3317) * pow(PhiDOT, 3) * pow(r, 3) +
                  4143 * pow(PhiDOT, 2) * pow(r, 2) * rDOT +
                  Complex(0, 2178) * PhiDOT * r * pow(rDOT, 2) -
                  442 * pow(rDOT, 3)) +
             45 * mass * pow(r, 2) *
                 (Complex(0, -1069) * pow(PhiDOT, 5) * pow(r, 5) +
                  2112 * pow(PhiDOT, 4) * pow(r, 4) * rDOT +
                  Complex(0, 2370) * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) -
                  1600 * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 3) -
                  Complex(0, 600) * PhiDOT * r * pow(rDOT, 4) +
                  96 * pow(rDOT, 5)))) /
           (720. * sqrt(6006) * pow(r, 3));
  }

  else if (vpnorder == 7) {

    return (
        (delta *
         (5040 * (-59 + 473 * Nu - 1185 * pow(Nu, 2) + 831 * pow(Nu, 3)) *
              pow(r, 4) * pow(PhiDOT * r + Complex(0, 1) * rDOT, 8) *
              (Complex(0, 1) * PhiDOT * r + rDOT) +
          4 * pow(mass, 4) *
              (Complex(0, -7) *
                   (-845315 + 3567934 * Nu - 3505793 * pow(Nu, 2) +
                    1006326 * pow(Nu, 3)) *
                   PhiDOT * r +
               4 *
                   (-482641 + 2069935 * Nu - 2198091 * pow(Nu, 2) +
                    803481 * pow(Nu, 3)) *
                   rDOT) +
          2 * pow(mass, 3) * r *
              (Complex(0, -7) *
                   (-8046605 + 37481621 * Nu - 46964439 * pow(Nu, 2) +
                    19173243 * pow(Nu, 3)) *
                   pow(PhiDOT, 3) * pow(r, 3) +
               2 *
                   (-17292661 + 90383188 * Nu - 148920723 * pow(Nu, 2) +
                    88022760 * pow(Nu, 3)) *
                   pow(PhiDOT, 2) * pow(r, 2) * rDOT +
               Complex(0, 7) *
                   (-1175695 + 7490560 * Nu - 17116047 * pow(Nu, 2) +
                    13239024 * pow(Nu, 3)) *
                   PhiDOT * r * pow(rDOT, 2) -
               8 *
                   (-74093 + 755216 * Nu - 2568447 * pow(Nu, 2) +
                    2398116 * pow(Nu, 3)) *
                   pow(rDOT, 3)) +
          18 * pow(mass, 2) * pow(r, 2) *
              (Complex(0, -14) *
                   (-270811 + 1455597 * Nu - 2122103 * pow(Nu, 2) +
                    757575 * pow(Nu, 3)) *
                   pow(PhiDOT, 5) * pow(r, 5) +
               (1113940 + 4956889 * Nu - 34515698 * pow(Nu, 2) +
                28651791 * pow(Nu, 3)) *
                   pow(PhiDOT, 4) * pow(r, 4) * rDOT +
               Complex(0, 7) *
                   (715061 - 1359787 * Nu - 4470709 * pow(Nu, 2) +
                    5729499 * pow(Nu, 3)) *
                   pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) +
               (-4243543 + 9349627 * Nu + 22353219 * pow(Nu, 2) -
                32044971 * pow(Nu, 3)) *
                   pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 3) -
               Complex(0, 14) *
                   (119115 - 252749 * Nu - 689509 * pow(Nu, 2) +
                    975153 * pow(Nu, 3)) *
                   PhiDOT * r * pow(rDOT, 4) +
               2 *
                   (131615 - 252847 * Nu - 892507 * pow(Nu, 2) +
                    1206639 * pow(Nu, 3)) *
                   pow(rDOT, 5)) +
          315 * mass * pow(r, 3) *
              (Complex(0, 1) *
                   (5005 + 44124 * Nu - 401470 * pow(Nu, 2) +
                    512250 * pow(Nu, 3)) *
                   pow(PhiDOT, 7) * pow(r, 7) +
               2 *
                   (32339 - 197757 * Nu + 459526 * pow(Nu, 2) -
                    383013 * pow(Nu, 3)) *
                   pow(PhiDOT, 6) * pow(r, 6) * rDOT -
               Complex(0, 1) *
                   (-168255 + 685408 * Nu - 635707 * pow(Nu, 2) +
                    199944 * pow(Nu, 3)) *
                   pow(PhiDOT, 5) * pow(r, 5) * pow(rDOT, 2) +
               4 *
                   (-53106 + 146273 * Nu + 123916 * pow(Nu, 2) -
                    235713 * pow(Nu, 3)) *
                   pow(PhiDOT, 4) * pow(r, 4) * pow(rDOT, 3) -
               Complex(0, 10) *
                   (16287 - 22820 * Nu - 136891 * pow(Nu, 2) +
                    159864 * pow(Nu, 3)) *
                   pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 4) +
               16 *
                   (4761 + 725 * Nu - 72818 * pow(Nu, 2) + 75357 * pow(Nu, 3)) *
                   pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 5) +
               Complex(0, 40) *
                   (499 + 1018 * Nu - 11789 * pow(Nu, 2) + 11502 * pow(Nu, 3)) *
                   PhiDOT * r * pow(rDOT, 6) -
               32 * (70 + 311 * Nu - 2394 * pow(Nu, 2) + 2253 * pow(Nu, 3)) *
                   pow(rDOT, 7)))) /
        (171360. * sqrt(6006) * pow(r, 4)));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hl_7_m_7(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_7_m_7: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           hGO_7_m_7(mass, Nu, r, rDOT, PhiDOT, vpnorder) * cpolar(1, -7 * Phi);
  }
}

static COMPLEX16 hl_7_m_min7(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_7_m_min7: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) * pow(-1, 7) *
           conj(hGO_7_m_7(mass, Nu, r, rDOT, PhiDOT, vpnorder)) *
           cpolar(1, 7 * Phi);
  }
}

// 75
static COMPLEX16 hGO_7_m_5(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  if (vpnorder == 5) {
    return (
        (delta * (1 - 4 * Nu + 3 * pow(Nu, 2)) *
         (2 * pow(mass, 3) * (Complex(0, 22795) * PhiDOT * r - 13832 * rDOT) -
          Complex(0, 5040) * pow(r, 3) * (PhiDOT * r - Complex(0, 1) * rDOT) *
              pow(PhiDOT * r + Complex(0, 1) * rDOT, 6) +
          18 * pow(mass, 2) * r *
              (Complex(0, 5105) * pow(PhiDOT, 3) * pow(r, 3) -
               13041 * pow(PhiDOT, 2) * pow(r, 2) * rDOT -
               Complex(0, 10890) * PhiDOT * r * pow(rDOT, 2) +
               3094 * pow(rDOT, 3)) +
          45 * mass * pow(r, 2) *
              (Complex(0, -1207) * pow(PhiDOT, 5) * pow(r, 5) +
               336 * pow(PhiDOT, 4) * pow(r, 4) * rDOT -
               Complex(0, 3114) * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) +
               4928 * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 3) +
               Complex(0, 3000) * PhiDOT * r * pow(rDOT, 4) -
               672 * pow(rDOT, 5)))) /
        (65520. * sqrt(66) * pow(r, 3)));
  }

  else if (vpnorder == 7) {

    return (
        (delta *
         (Complex(0, 5040) *
              (-59 + 473 * Nu - 1185 * pow(Nu, 2) + 831 * pow(Nu, 3)) *
              pow(r, 4) * pow(PhiDOT * r - Complex(0, 1) * rDOT, 2) *
              pow(PhiDOT * r + Complex(0, 1) * rDOT, 7) +
          4 * pow(mass, 4) *
              (Complex(0, 5) *
                   (-845315 + 3567934 * Nu - 3505793 * pow(Nu, 2) +
                    1006326 * pow(Nu, 3)) *
                   PhiDOT * r +
               4 *
                   (482641 - 2069935 * Nu + 2198091 * pow(Nu, 2) -
                    803481 * pow(Nu, 3)) *
                   rDOT) +
          2 * pow(mass, 3) * r *
              (Complex(0, 5) *
                   (-2480333 + 11838029 * Nu - 15223575 * pow(Nu, 2) +
                    5981667 * pow(Nu, 3)) *
                   pow(PhiDOT, 3) * pow(r, 3) -
               2 *
                   (-7732189 + 40782964 * Nu - 67735947 * pow(Nu, 2) +
                    39807720 * pow(Nu, 3)) *
                   pow(PhiDOT, 2) * pow(r, 2) * rDOT -
               Complex(0, 5) *
                   (-1175695 + 7490560 * Nu - 17116047 * pow(Nu, 2) +
                    13239024 * pow(Nu, 3)) *
                   PhiDOT * r * pow(rDOT, 2) +
               8 *
                   (-74093 + 755216 * Nu - 2568447 * pow(Nu, 2) +
                    2398116 * pow(Nu, 3)) *
                   pow(rDOT, 3)) +
          18 * pow(mass, 2) * pow(r, 2) *
              (Complex(0, -10) *
                   (-58221 + 312011 * Nu - 571313 * pow(Nu, 2) +
                    397665 * pow(Nu, 3)) *
                   pow(PhiDOT, 5) * pow(r, 5) +
               (-117884 + 723063 * Nu - 2639174 * pow(Nu, 2) +
                3313409 * pow(Nu, 3)) *
                   pow(PhiDOT, 4) * pow(r, 4) * rDOT -
               Complex(0, 5) *
                   (236029 - 577139 * Nu - 804413 * pow(Nu, 2) +
                    1190115 * pow(Nu, 3)) *
                   pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) +
               (1940751 - 4452163 * Nu - 9331755 * pow(Nu, 2) +
                13753811 * pow(Nu, 3)) *
                   pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 3) +
               Complex(0, 10) *
                   (119115 - 252749 * Nu - 689509 * pow(Nu, 2) +
                    975153 * pow(Nu, 3)) *
                   PhiDOT * r * pow(rDOT, 4) +
               2 *
                   (-131615 + 252847 * Nu + 892507 * pow(Nu, 2) -
                    1206639 * pow(Nu, 3)) *
                   pow(rDOT, 5)) +
          45 * mass * pow(r, 3) *
              (Complex(0, 1) *
                   (39847 - 121612 * Nu - 163210 * pow(Nu, 2) +
                    376622 * pow(Nu, 3)) *
                   pow(PhiDOT, 7) * pow(r, 7) +
               2 *
                   (79987 - 238709 * Nu - 93106 * pow(Nu, 2) +
                    259939 * pow(Nu, 3)) *
                   pow(PhiDOT, 6) * pow(r, 6) * rDOT +
               Complex(0, 1) *
                   (173997 + 166656 * Nu - 3315983 * pow(Nu, 2) +
                    3362728 * pow(Nu, 3)) *
                   pow(PhiDOT, 5) * pow(r, 5) * pow(rDOT, 2) -
               4 *
                   (-13074 + 266583 * Nu - 998488 * pow(Nu, 2) +
                    847097 * pow(Nu, 3)) *
                   pow(PhiDOT, 4) * pow(r, 4) * pow(rDOT, 3) +
               Complex(0, 2) *
                   (135447 - 440484 * Nu - 25811 * pow(Nu, 2) +
                    357784 * pow(Nu, 3)) *
                   pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 4) -
               16 *
                   (15423 - 4597 * Nu - 205150 * pow(Nu, 2) +
                    217363 * pow(Nu, 3)) *
                   pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 5) -
               Complex(0, 200) *
                   (499 + 1018 * Nu - 11789 * pow(Nu, 2) + 11502 * pow(Nu, 3)) *
                   PhiDOT * r * pow(rDOT, 6) +
               224 * (70 + 311 * Nu - 2394 * pow(Nu, 2) + 2253 * pow(Nu, 3)) *
                   pow(rDOT, 7)))) /
        (2.22768e6 * sqrt(66) * pow(r, 4)));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hl_7_m_5(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_7_m_5: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           hGO_7_m_5(mass, Nu, r, rDOT, PhiDOT, vpnorder) * cpolar(1, -5 * Phi);
  }
}

static COMPLEX16 hl_7_m_min5(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_7_m_min5: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) * pow(-1, 7) *
           conj(hGO_7_m_5(mass, Nu, r, rDOT, PhiDOT, vpnorder)) *
           cpolar(1, 5 * Phi);
  }
}

// 73
static COMPLEX16 hGO_7_m_3(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  if (vpnorder == 5) {
    return (delta * (1 - 4 * Nu + 3 * pow(Nu, 2)) *
            (5040 * pow(r, 3) * pow(PhiDOT * r - Complex(0, 1) * rDOT, 2) *
                 pow(Complex(0, -1) * PhiDOT * r + rDOT, 5) +
             2 * pow(mass, 3) *
                 (Complex(0, -13677) * PhiDOT * r + 13832 * rDOT) +
             18 * pow(mass, 2) * r *
                 (Complex(0, 1529) * pow(PhiDOT, 3) * pow(r, 3) +
                  2401 * pow(PhiDOT, 2) * pow(r, 2) * rDOT +
                  Complex(0, 6534) * PhiDOT * r * pow(rDOT, 2) -
                  3094 * pow(rDOT, 3)) +
             15 * mass * pow(r, 2) *
                 (Complex(0, 179) * pow(PhiDOT, 5) * pow(r, 5) -
                  4368 * pow(PhiDOT, 4) * pow(r, 4) * rDOT -
                  Complex(0, 4878) * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) -
                  2240 * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 3) -
                  Complex(0, 5400) * PhiDOT * r * pow(rDOT, 4) +
                  2016 * pow(rDOT, 5)))) /
           (240240. * sqrt(6) * pow(r, 3));
  }

  else if (vpnorder == 7) {
    return (
        (delta *
         (-5040 * (-59 + 473 * Nu - 1185 * pow(Nu, 2) + 831 * pow(Nu, 3)) *
              pow(r, 4) * pow(PhiDOT * r + Complex(0, 1) * rDOT, 6) *
              pow(Complex(0, 1) * PhiDOT * r + rDOT, 3) +
          4 * pow(mass, 4) *
              (Complex(0, -3) *
                   (-845315 + 3567934 * Nu - 3505793 * pow(Nu, 2) +
                    1006326 * pow(Nu, 3)) *
                   PhiDOT * r +
               4 *
                   (-482641 + 2069935 * Nu - 2198091 * pow(Nu, 2) +
                    803481 * pow(Nu, 3)) *
                   rDOT) +
          2 * pow(mass, 3) * r *
              (Complex(0, 3) *
                   (-1230515 + 5257699 * Nu - 5937001 * pow(Nu, 2) +
                    2812717 * pow(Nu, 3)) *
                   pow(PhiDOT, 3) * pow(r, 3) +
               2 *
                   (-1358541 + 7716148 * Nu - 13612763 * pow(Nu, 2) +
                    7664360 * pow(Nu, 3)) *
                   pow(PhiDOT, 2) * pow(r, 2) * rDOT +
               Complex(0, 3) *
                   (-1175695 + 7490560 * Nu - 17116047 * pow(Nu, 2) +
                    13239024 * pow(Nu, 3)) *
                   PhiDOT * r * pow(rDOT, 2) -
               8 *
                   (-74093 + 755216 * Nu - 2568447 * pow(Nu, 2) +
                    2398116 * pow(Nu, 3)) *
                   pow(rDOT, 3)) +
          6 * pow(mass, 2) * pow(r, 2) *
              (Complex(0, -2) *
                   (76741 - 370227 * Nu + 248553 * pow(Nu, 2) +
                    279655 * pow(Nu, 3)) *
                   pow(PhiDOT, 5) * pow(r, 5) -
               3 *
                   (65300 + 667351 * Nu - 4038422 * pow(Nu, 2) +
                    3825889 * pow(Nu, 3)) *
                   pow(PhiDOT, 4) * pow(r, 4) * rDOT -
               Complex(0, 3) *
                   (249977 + 166121 * Nu - 4919353 * pow(Nu, 2) +
                    5508423 * pow(Nu, 3)) *
                   pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) +
               (-1216669 + 3561561 * Nu + 1952337 * pow(Nu, 2) -
                4679113 * pow(Nu, 3)) *
                   pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 3) -
               Complex(0, 18) *
                   (119115 - 252749 * Nu - 689509 * pow(Nu, 2) +
                    975153 * pow(Nu, 3)) *
                   PhiDOT * r * pow(rDOT, 4) +
               6 *
                   (131615 - 252847 * Nu - 892507 * pow(Nu, 2) +
                    1206639 * pow(Nu, 3)) *
                   pow(rDOT, 5)) +
          15 * mass * pow(r, 3) *
              (Complex(0, -17) *
                   (-3693 + 24516 * Nu - 50930 * pow(Nu, 2) +
                    30982 * pow(Nu, 3)) *
                   pow(PhiDOT, 7) * pow(r, 7) +
               6 *
                   (14141 + 89013 * Nu - 604958 * pow(Nu, 2) +
                    566877 * pow(Nu, 3)) *
                   pow(PhiDOT, 6) * pow(r, 6) * rDOT +
               Complex(0, 1) *
                   (-125697 + 1209408 * Nu - 3554741 * pow(Nu, 2) +
                    2822200 * pow(Nu, 3)) *
                   pow(PhiDOT, 5) * pow(r, 5) * pow(rDOT, 2) +
               4 *
                   (90786 - 8859 * Nu - 1290784 * pow(Nu, 2) +
                    1354859 * pow(Nu, 3)) *
                   pow(PhiDOT, 4) * pow(r, 4) * pow(rDOT, 3) +
               Complex(0, 6) *
                   (27423 + 212284 * Nu - 1343099 * pow(Nu, 2) +
                    1240856 * pow(Nu, 3)) *
                   pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 4) +
               16 *
                   (10461 - 33135 * Nu - 6298 * pow(Nu, 2) +
                    31817 * pow(Nu, 3)) *
                   pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 5) +
               Complex(0, 360) *
                   (499 + 1018 * Nu - 11789 * pow(Nu, 2) + 11502 * pow(Nu, 3)) *
                   PhiDOT * r * pow(rDOT, 6) -
               672 * (70 + 311 * Nu - 2394 * pow(Nu, 2) + 2253 * pow(Nu, 3)) *
                   pow(rDOT, 7)))) /
        (8.16816e6 * sqrt(6) * pow(r, 4)));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hl_7_m_3(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_7_m_3: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           hGO_7_m_3(mass, Nu, r, rDOT, PhiDOT, vpnorder) * cpolar(1, -3 * Phi);
  }
}

static COMPLEX16 hl_7_m_min3(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_7_m_min3: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) * pow(-1, 7) *
           conj(hGO_7_m_3(mass, Nu, r, rDOT, PhiDOT, vpnorder)) *
           cpolar(1, 3 * Phi);
  }
}

// 71
static COMPLEX16 hGO_7_m_1(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  if (vpnorder == 5) {
    return (delta * (1 - 4 * Nu + 3 * pow(Nu, 2)) *
            (2 * pow(mass, 3) * (Complex(0, 4559) * PhiDOT * r - 13832 * rDOT) -
             Complex(0, 5040) * pow(r, 3) *
                 pow(PhiDOT * r - Complex(0, 1) * rDOT, 3) *
                 pow(PhiDOT * r + Complex(0, 1) * rDOT, 4) +
             18 * pow(mass, 2) * r *
                 (Complex(0, -1275) * pow(PhiDOT, 3) * pow(r, 3) +
                  2919 * pow(PhiDOT, 2) * pow(r, 2) * rDOT -
                  Complex(0, 2178) * PhiDOT * r * pow(rDOT, 2) +
                  3094 * pow(rDOT, 3)) +
             Complex(0, 27) * mass * pow(r, 2) *
                 (699 * pow(PhiDOT, 5) * pow(r, 5) +
                  Complex(0, 1120) * pow(PhiDOT, 4) * pow(r, 4) * rDOT +
                  1874 * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) +
                  Complex(0, 2240) * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 3) +
                  1000 * PhiDOT * r * pow(rDOT, 4) +
                  Complex(0, 1120) * pow(rDOT, 5)))) /
           (432432. * sqrt(2) * pow(r, 3));
  }

  else if (vpnorder == 7) {

    return (
        (delta *
         (Complex(0, 5040) *
              (-59 + 473 * Nu - 1185 * pow(Nu, 2) + 831 * pow(Nu, 3)) *
              pow(r, 4) * pow(PhiDOT * r - Complex(0, 1) * rDOT, 4) *
              pow(PhiDOT * r + Complex(0, 1) * rDOT, 5) +
          4 * pow(mass, 4) *
              (Complex(0, 1) *
                   (-845315 + 3567934 * Nu - 3505793 * pow(Nu, 2) +
                    1006326 * pow(Nu, 3)) *
                   PhiDOT * r +
               4 *
                   (482641 - 2069935 * Nu + 2198091 * pow(Nu, 2) -
                    803481 * pow(Nu, 3)) *
                   rDOT) +
          2 * pow(mass, 3) * r *
              (Complex(0, -1) *
                   (-3085939 + 13805563 * Nu - 16517289 * pow(Nu, 2) +
                    7209909 * pow(Nu, 3)) *
                   pow(PhiDOT, 3) * pow(r, 3) +
               2 *
                   (-1828283 + 8817260 * Nu - 13448829 * pow(Nu, 2) +
                    8407320 * pow(Nu, 3)) *
                   pow(PhiDOT, 2) * pow(r, 2) * rDOT -
               Complex(0, 1) *
                   (-1175695 + 7490560 * Nu - 17116047 * pow(Nu, 2) +
                    13239024 * pow(Nu, 3)) *
                   PhiDOT * r * pow(rDOT, 2) +
               8 *
                   (-74093 + 755216 * Nu - 2568447 * pow(Nu, 2) +
                    2398116 * pow(Nu, 3)) *
                   pow(rDOT, 3)) +
          18 * pow(mass, 2) * pow(r, 2) *
              (Complex(0, 2) *
                   (-97035 + 529085 * Nu - 946023 * pow(Nu, 2) +
                    605111 * pow(Nu, 3)) *
                   pow(PhiDOT, 5) * pow(r, 5) +
               (12636 - 513209 * Nu + 2273154 * pow(Nu, 2) -
                2157167 * pow(Nu, 3)) *
                   pow(PhiDOT, 4) * pow(r, 4) * rDOT +
               Complex(0, 3) *
                   (81001 - 68503 * Nu - 953961 * pow(Nu, 2) +
                    1116423 * pow(Nu, 3)) *
                   pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) +
               (-362041 + 445301 * Nu + 3689709 * pow(Nu, 2) -
                4537349 * pow(Nu, 3)) *
                   pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 3) +
               Complex(0, 2) *
                   (119115 - 252749 * Nu - 689509 * pow(Nu, 2) +
                    975153 * pow(Nu, 3)) *
                   PhiDOT * r * pow(rDOT, 4) +
               2 *
                   (-131615 + 252847 * Nu + 892507 * pow(Nu, 2) -
                    1206639 * pow(Nu, 3)) *
                   pow(rDOT, 5)) -
          Complex(0, 9) * mass * pow(r, 3) *
              ((-107407 + 883020 * Nu - 2319286 * pow(Nu, 2) +
                1727170 * pow(Nu, 3)) *
                   pow(PhiDOT, 7) * pow(r, 7) +
               Complex(0, 2) *
                   (-937 + 307287 * Nu - 1356450 * pow(Nu, 2) +
                    1189583 * pow(Nu, 3)) *
                   pow(PhiDOT, 6) * pow(r, 6) * rDOT +
               (216555 + 690656 * Nu - 6235017 * pow(Nu, 2) +
                5984984 * pow(Nu, 3)) *
                   pow(PhiDOT, 5) * pow(r, 5) * pow(rDOT, 2) +
               Complex(0, 4) *
                   (34014 + 335581 * Nu - 1987500 * pow(Nu, 2) +
                    1820899 * pow(Nu, 3)) *
                   pow(PhiDOT, 4) * pow(r, 4) * pow(rDOT, 3) +
               2 *
                   (136281 + 310468 * Nu - 3370653 * pow(Nu, 2) +
                    3281032 * pow(Nu, 3)) *
                   pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 4) +
               Complex(0, 80) *
                   (2481 + 14269 * Nu - 99426 * pow(Nu, 2) +
                    92773 * pow(Nu, 3)) *
                   pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 5) +
               200 *
                   (499 + 1018 * Nu - 11789 * pow(Nu, 2) + 11502 * pow(Nu, 3)) *
                   PhiDOT * r * pow(rDOT, 6) +
               Complex(0, 1120) *
                   (70 + 311 * Nu - 2394 * pow(Nu, 2) + 2253 * pow(Nu, 3)) *
                   pow(rDOT, 7)))) /
        (1.4702688e7 * sqrt(2) * pow(r, 4)));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hl_7_m_1(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_7_m_1: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           hGO_7_m_1(mass, Nu, r, rDOT, PhiDOT, vpnorder) * cpolar(1, -1 * Phi);
  }
}

static COMPLEX16 hl_7_m_min1(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_7_m_min1: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) * pow(-1, 7) *
           conj(hGO_7_m_1(mass, Nu, r, rDOT, PhiDOT, vpnorder)) *
           cpolar(1, 1 * Phi);
  }
}

// 72
static COMPLEX16 hGO_7_m_2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  if (vpnorder == 6) {
    return (-(mass * (-1 + 7 * Nu - 14 * pow(Nu, 2) + 7 * pow(Nu, 3)) * PhiDOT *
              (8 * pow(mass, 2) * (494 * PhiDOT * r + Complex(0, 411) * rDOT) -
               12 * mass * r *
                   (530 * pow(PhiDOT, 3) * pow(r, 3) -
                    Complex(0, 6) * pow(PhiDOT, 2) * pow(r, 2) * rDOT +
                    453 * PhiDOT * r * pow(rDOT, 2) +
                    Complex(0, 197) * pow(rDOT, 3)) +
               3 * pow(r, 2) *
                   (824 * pow(PhiDOT, 5) * pow(r, 5) -
                    Complex(0, 671) * pow(PhiDOT, 4) * pow(r, 4) * rDOT +
                    864 * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) +
                    Complex(0, 44) * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 3) +
                    320 * PhiDOT * r * pow(rDOT, 4) +
                    Complex(0, 120) * pow(rDOT, 5)))) /
            (96096. * sqrt(3) * pow(r, 2)));
  }

  else if (vpnorder == 7) {

    return (/* (4*Nu*(S1z - delta*S1z + 2*(-2 + delta)*Nu*S1z + S2z + delta*S2z
      - 2*(2 + delta)*Nu*S2z + 2*pow(Nu,2)*(S1z +
      S2z))*pow(x,4.5))/(3003.*sqrt(3))
      + */
            /* Henry et al. ecc+spin terms */ (
                (pow(mass, 2) * Nu *
                 (8 * pow(mass, 2) *
                      (494 * PhiDOT * r + Complex(0, 411) * rDOT) -
                  12 * mass * r *
                      (530 * pow(PhiDOT, 3) * pow(r, 3) -
                       Complex(0, 6) * pow(PhiDOT, 2) * pow(r, 2) * rDOT +
                       453 * PhiDOT * r * pow(rDOT, 2) +
                       Complex(0, 197) * pow(rDOT, 3)) +
                  3 * pow(r, 2) *
                      (824 * pow(PhiDOT, 5) * pow(r, 5) -
                       Complex(0, 671) * pow(PhiDOT, 4) * pow(r, 4) * rDOT +
                       864 * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) +
                       Complex(0, 44) * pow(PhiDOT, 2) * pow(r, 2) *
                           pow(rDOT, 3) +
                       320 * PhiDOT * r * pow(rDOT, 4) +
                       Complex(0, 120) * pow(rDOT, 5))) *
                 ((1 - 4 * Nu + 2 * pow(Nu, 2) + delta * (-1 + 2 * Nu)) * S1z +
                  (1 + delta - 4 * Nu - 2 * delta * Nu + 2 * pow(Nu, 2)) *
                      S2z)) /
                (48048. * sqrt(3) * pow(r, 4))));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hl_7_m_2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_7_m_2: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           hGO_7_m_2(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z) *
           cpolar(1, -2 * Phi);
  }
}

static COMPLEX16 hl_7_m_min2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_7_m_min2: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) * pow(-1, 7) *
           conj(hGO_7_m_2(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z)) *
           cpolar(1, 2 * Phi);
  }
}

// 74
static COMPLEX16 hGO_7_m_4(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  if (vpnorder == 6) {
    return ((mass * (-1 + 7 * Nu - 14 * pow(Nu, 2) + 7 * pow(Nu, 3)) * PhiDOT *
             (8 * pow(mass, 2) * (988 * PhiDOT * r + Complex(0, 411) * rDOT) +
              12 * mass * r *
                  (844 * pow(PhiDOT, 3) * pow(r, 3) +
                   Complex(0, 1518) * pow(PhiDOT, 2) * pow(r, 2) * rDOT -
                   906 * PhiDOT * r * pow(rDOT, 2) -
                   Complex(0, 197) * pow(rDOT, 3)) -
              15 * pow(r, 2) *
                  (656 * pow(PhiDOT, 5) * pow(r, 5) +
                   Complex(0, 179) * pow(PhiDOT, 4) * pow(r, 4) * rDOT +
                   192 * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) +
                   Complex(0, 260) * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 3) -
                   128 * PhiDOT * r * pow(rDOT, 4) -
                   Complex(0, 24) * pow(rDOT, 5)))) /
            (21840. * sqrt(66) * pow(r, 2)));
  }

  else if (vpnorder == 7) {

    return (/* (-512*sqrt(0.06060606060606061)*Nu*(S1z - delta*S1z + 2*(-2 +
      delta)*Nu*S1z
      + S2z + delta*S2z - 2*(2 + delta)*Nu*S2z + 2*pow(Nu,2)*(S1z +
      S2z))*pow(x,4.5))/1365.
      + */
            /* Henry et al. ecc+spin terms */ (
                -0.00009157509157509158 *
                (pow(mass, 2) * Nu *
                 (8 * pow(mass, 2) *
                      (988 * PhiDOT * r + Complex(0, 411) * rDOT) +
                  12 * mass * r *
                      (844 * pow(PhiDOT, 3) * pow(r, 3) +
                       Complex(0, 1518) * pow(PhiDOT, 2) * pow(r, 2) * rDOT -
                       906 * PhiDOT * r * pow(rDOT, 2) -
                       Complex(0, 197) * pow(rDOT, 3)) -
                  15 * pow(r, 2) *
                      (656 * pow(PhiDOT, 5) * pow(r, 5) +
                       Complex(0, 179) * pow(PhiDOT, 4) * pow(r, 4) * rDOT +
                       192 * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) +
                       Complex(0, 260) * pow(PhiDOT, 2) * pow(r, 2) *
                           pow(rDOT, 3) -
                       128 * PhiDOT * r * pow(rDOT, 4) -
                       Complex(0, 24) * pow(rDOT, 5))) *
                 ((1 - 4 * Nu + 2 * pow(Nu, 2) + delta * (-1 + 2 * Nu)) * S1z +
                  (1 + delta - 4 * Nu - 2 * delta * Nu + 2 * pow(Nu, 2)) *
                      S2z)) /
                (sqrt(66) * pow(r, 4))));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hl_7_m_4(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_7_m_4: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           hGO_7_m_4(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z) *
           cpolar(1, -4 * Phi);
  }
}

static COMPLEX16 hl_7_m_min4(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_7_m_min4: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) * pow(-1, 7) *
           conj(hGO_7_m_4(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z)) *
           cpolar(1, 4 * Phi);
  }
}

// 76
static COMPLEX16 hGO_7_m_6(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z) {
  REAL8 delta = sqrt(1 - 4 * Nu);

  if (vpnorder == 6) {
    return (
        -(mass * (-1 + 7 * Nu - 14 * pow(Nu, 2) + 7 * pow(Nu, 3)) * PhiDOT *
          (8 * pow(mass, 2) * (494 * PhiDOT * r + Complex(0, 137) * rDOT) +
           4 * mass * r *
               (6026 * pow(PhiDOT, 3) * pow(r, 3) +
                Complex(0, 4038) * pow(PhiDOT, 2) * pow(r, 2) * rDOT -
                1359 * PhiDOT * r * pow(rDOT, 2) -
                Complex(0, 197) * pow(rDOT, 3)) +
           15 * pow(r, 2) *
               (1240 * pow(PhiDOT, 5) * pow(r, 5) +
                Complex(0, 911) * pow(PhiDOT, 4) * pow(r, 4) * rDOT -
                544 * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) -
                Complex(0, 236) * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 3) +
                64 * PhiDOT * r * pow(rDOT, 4) +
                Complex(0, 8) * pow(rDOT, 5)))) /
        (3360. * sqrt(429) * pow(r, 2)));
  } else if (vpnorder == 7) {

    return (/* (324*sqrt(0.02097902097902098)*Nu*(S1z - delta*S1z
      + 2*(-2 + delta)*Nu*S1z + S2z + delta*S2z - 2*(2 + delta)*Nu*S2z
      + 2*pow(Nu,2)*(S1z + S2z))*pow(x,4.5))/35.
      + */
            /* Henry et al. ecc+spin terms */ (
                (pow(mass, 2) * Nu *
                 (8 * pow(mass, 2) *
                      (494 * PhiDOT * r + Complex(0, 137) * rDOT) +
                  4 * mass * r *
                      (6026 * pow(PhiDOT, 3) * pow(r, 3) +
                       Complex(0, 4038) * pow(PhiDOT, 2) * pow(r, 2) * rDOT -
                       1359 * PhiDOT * r * pow(rDOT, 2) -
                       Complex(0, 197) * pow(rDOT, 3)) +
                  15 * pow(r, 2) *
                      (1240 * pow(PhiDOT, 5) * pow(r, 5) +
                       Complex(0, 911) * pow(PhiDOT, 4) * pow(r, 4) * rDOT -
                       544 * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 2) -
                       Complex(0, 236) * pow(PhiDOT, 2) * pow(r, 2) *
                           pow(rDOT, 3) +
                       64 * PhiDOT * r * pow(rDOT, 4) +
                       Complex(0, 8) * pow(rDOT, 5))) *
                 ((1 - 4 * Nu + 2 * pow(Nu, 2) + delta * (-1 + 2 * Nu)) * S1z +
                  (1 + delta - 4 * Nu - 2 * delta * Nu + 2 * pow(Nu, 2)) *
                      S2z)) /
                (1680. * sqrt(429) * pow(r, 4))));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hl_7_m_6(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder, REAL8 S1z,
                          REAL8 S2z) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_7_m_6: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           hGO_7_m_6(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z) *
           cpolar(1, -6 * Phi);
  }
}

static COMPLEX16 hl_7_m_min6(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder,
                             REAL8 S1z, REAL8 S2z) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_7_m_min6: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) * pow(-1, 7) *
           conj(hGO_7_m_6(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z)) *
           cpolar(1, 6 * Phi);
  }
}

// 88
static COMPLEX16 hGO_8_m_8(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder) {

  if (vpnorder == 6) {
    return (
        ((-1 + 7 * Nu - 14 * pow(Nu, 2) + 7 * pow(Nu, 3)) *
         (9118 * pow(mass, 4) +
          5040 * pow(r, 4) * pow(PhiDOT * r + Complex(0, 1) * rDOT, 8) +
          4 * pow(mass, 3) * r *
              (82543 * pow(PhiDOT, 2) * pow(r, 2) +
               Complex(0, 67760) * PhiDOT * r * rDOT - 16717 * pow(rDOT, 2)) +
          9 * pow(mass, 2) * pow(r, 2) *
              (124583 * pow(PhiDOT, 4) * pow(r, 4) +
               Complex(0, 192640) * pow(PhiDOT, 3) * pow(r, 3) * rDOT -
               144024 * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) -
               Complex(0, 56280) * PhiDOT * r * pow(rDOT, 3) +
               9188 * pow(rDOT, 4)) +
          315 * mass * pow(r, 3) *
              (2005 * pow(PhiDOT, 6) * pow(r, 6) +
               Complex(0, 4250) * pow(PhiDOT, 5) * pow(r, 5) * rDOT -
               5538 * pow(PhiDOT, 4) * pow(r, 4) * pow(rDOT, 2) -
               Complex(0, 4760) * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 3) +
               2600 * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 4) +
               Complex(0, 816) * PhiDOT * r * pow(rDOT, 5) -
               112 * pow(rDOT, 6)))) /
        (2016. * sqrt(170170) * pow(r, 4)));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hl_8_m_8(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_8_m_8: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           hGO_8_m_8(mass, Nu, r, rDOT, PhiDOT, vpnorder) * cpolar(1, -8 * Phi);
  }
}

static COMPLEX16 hl_8_m_min8(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_8_m_min8: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) * pow(-1, 8) *
           conj(hGO_8_m_8(mass, Nu, r, rDOT, PhiDOT, vpnorder)) *
           cpolar(1, 8 * Phi);
  }
}

// 86
static COMPLEX16 hGO_8_m_6(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder) {

  if (vpnorder == 6) {
    return (
        -((-1 + 7 * Nu - 14 * pow(Nu, 2) + 7 * pow(Nu, 3)) *
          (18236 * pow(mass, 4) -
           10080 * pow(r, 4) * (PhiDOT * r - Complex(0, 1) * rDOT) *
               pow(PhiDOT * r + Complex(0, 1) * rDOT, 7) +
           8 * pow(mass, 3) * r *
               (42923 * pow(PhiDOT, 2) * pow(r, 2) +
                Complex(0, 50820) * PhiDOT * r * rDOT - 16717 * pow(rDOT, 2)) +
           18 * pow(mass, 2) * pow(r, 2) *
               (15663 * pow(PhiDOT, 4) * pow(r, 4) +
                Complex(0, 58240) * pow(PhiDOT, 3) * pow(r, 3) * rDOT -
                74094 * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) -
                Complex(0, 42210) * PhiDOT * r * pow(rDOT, 3) +
                9188 * pow(rDOT, 4)) -
           315 * mass * pow(r, 3) *
               (678 * pow(PhiDOT, 6) * pow(r, 6) +
                Complex(0, 681) * pow(PhiDOT, 5) * pow(r, 5) * rDOT +
                852 * pow(PhiDOT, 4) * pow(r, 4) * pow(rDOT, 2) +
                Complex(0, 2660) * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 3) -
                2640 * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 4) -
                Complex(0, 1224) * PhiDOT * r * pow(rDOT, 5) +
                224 * pow(rDOT, 6)))) /
        (10080. * sqrt(51051) * pow(r, 4)));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hl_8_m_6(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_8_m_6: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           hGO_8_m_6(mass, Nu, r, rDOT, PhiDOT, vpnorder) * cpolar(1, -6 * Phi);
  }
}

static COMPLEX16 hl_8_m_min6(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_8_m_min6: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) * pow(-1, 8) *
           conj(hGO_8_m_6(mass, Nu, r, rDOT, PhiDOT, vpnorder)) *
           cpolar(1, 6 * Phi);
  }
}

// 84
static COMPLEX16 hGO_8_m_4(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder) {

  if (vpnorder == 6) {
    return (
        ((-1 + 7 * Nu - 14 * pow(Nu, 2) + 7 * pow(Nu, 3)) *
         (9118 * pow(mass, 4) +
          5040 * pow(r, 4) * pow(PhiDOT * r - Complex(0, 1) * rDOT, 2) *
              pow(PhiDOT * r + Complex(0, 1) * rDOT, 6) +
          4 * pow(mass, 3) * r *
              (14623 * pow(PhiDOT, 2) * pow(r, 2) +
               Complex(0, 33880) * PhiDOT * r * rDOT - 16717 * pow(rDOT, 2)) -
          9 * pow(mass, 2) * pow(r, 2) *
              (8377 * pow(PhiDOT, 4) * pow(r, 4) +
               Complex(0, 2240) * pow(PhiDOT, 3) * pow(r, 3) * rDOT +
               24144 * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) +
               Complex(0, 28140) * PhiDOT * r * pow(rDOT, 3) -
               9188 * pow(rDOT, 4)) +
          45 * mass * pow(r, 3) *
              (243 * pow(PhiDOT, 6) * pow(r, 6) -
               Complex(0, 1701) * pow(PhiDOT, 5) * pow(r, 5) * rDOT +
               3762 * pow(PhiDOT, 4) * pow(r, 4) * pow(rDOT, 2) +
               Complex(0, 1260) * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 3) +
               2840 * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 4) +
               Complex(0, 2856) * PhiDOT * r * pow(rDOT, 5) -
               784 * pow(rDOT, 6)))) /
        (65520. * sqrt(374) * pow(r, 4)));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hl_8_m_4(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_8_m_4: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           hGO_8_m_4(mass, Nu, r, rDOT, PhiDOT, vpnorder) * cpolar(1, -4 * Phi);
  }
}

static COMPLEX16 hl_8_m_min4(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_8_m_min4: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) * pow(-1, 8) *
           conj(hGO_8_m_4(mass, Nu, r, rDOT, PhiDOT, vpnorder)) *
           cpolar(1, 4 * Phi);
  }
}

// 82
static COMPLEX16 hGO_8_m_2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                           REAL8 PhiDOT, UINT4 vpnorder) {

  if (vpnorder == 6) {
    return (
        ((-1 + 7 * Nu - 14 * pow(Nu, 2) + 7 * pow(Nu, 3)) *
         (-18236 * pow(mass, 4) +
          10080 * pow(r, 4) * pow(PhiDOT * r - Complex(0, 1) * rDOT, 3) *
              pow(PhiDOT * r + Complex(0, 1) * rDOT, 5) +
          8 * pow(mass, 3) * r *
              (2357 * pow(PhiDOT, 2) * pow(r, 2) -
               Complex(0, 16940) * PhiDOT * r * rDOT + 16717 * pow(rDOT, 2)) +
          18 * pow(mass, 2) * pow(r, 2) *
              (1297 * pow(PhiDOT, 4) * pow(r, 4) +
               Complex(0, 13440) * pow(PhiDOT, 3) * pow(r, 3) * rDOT -
               5826 * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 2) +
               Complex(0, 14070) * PhiDOT * r * pow(rDOT, 3) -
               9188 * pow(rDOT, 4)) -
          45 * mass * pow(r, 3) *
              (758 * pow(PhiDOT, 6) * pow(r, 6) +
               Complex(0, 2891) * pow(PhiDOT, 5) * pow(r, 5) * rDOT +
               564 * pow(PhiDOT, 4) * pow(r, 4) * pow(rDOT, 2) +
               Complex(0, 5740) * pow(PhiDOT, 3) * pow(r, 3) * pow(rDOT, 3) -
               2000 * pow(PhiDOT, 2) * pow(r, 2) * pow(rDOT, 4) +
               Complex(0, 2856) * PhiDOT * r * pow(rDOT, 5) -
               1568 * pow(rDOT, 6)))) /
        (288288. * sqrt(85) * pow(r, 4)));
  }

  else {
    return 0;
  }
}

static COMPLEX16 hl_8_m_2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT, REAL8 Phi,
                          REAL8 PhiDOT, REAL8 R, UINT4 vpnorder) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_8_m_2: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) *
           hGO_8_m_2(mass, Nu, r, rDOT, PhiDOT, vpnorder) * cpolar(1, -2 * Phi);
  }
}

static COMPLEX16 hl_8_m_min2(REAL8 mass, REAL8 Nu, REAL8 r, REAL8 rDOT,
                             REAL8 Phi, REAL8 PhiDOT, REAL8 R, UINT4 vpnorder) {

  if ((vpnorder < 0) || (vpnorder > 8)) {
    XLAL_ERROR(XLAL_EINVAL, "Error in hl_8_m_min2: Input PN order parameter "
                            "should be between [0, 8].");
  }

  else {
    return ((4 * mass * Nu * sqrt(M_PI / 5.)) / R) * pow(-1, 8) *
           conj(hGO_8_m_2(mass, Nu, r, rDOT, PhiDOT, vpnorder)) *
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
        hlm += hl_2_m_2(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x);
      }
      return (hlm);

    case 1:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_2_m_1(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x);
      }
      return (hlm);
    case 0:
      return (0);

    case -1:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_2_m_min1(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x);
      }
      return (hlm);
    case -2:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_2_m_min2(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x);
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
        hlm += hl_3_m_3(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x);
      }
      return (hlm);
    case 2:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_3_m_2(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x);
      }
      return (hlm);
    case 1:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_3_m_1(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x);
      }
      return (hlm);
    case 0:
      return (0);

    case -1:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_3_m_min1(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x);
      }
      return (hlm);
    case -2:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_3_m_min2(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x);
      }
      return (hlm);
    case -3:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_3_m_min3(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x);
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
        hlm += hl_4_m_4(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x);
      }
      return (hlm);
    case 3:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_4_m_3(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x);
      }
      return (hlm);
    case 2:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_4_m_2(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x);
      }
      return (hlm);
    case 1:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_4_m_1(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x);
      }
      return (hlm);
    case 0:
      return (0);

    case -1:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_4_m_min1(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x);
      }
      return (hlm);
    case -2:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_4_m_min2(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x);
      }
      return (hlm);
    case -3:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_4_m_min3(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x);
      }
      return (hlm);
    case -4:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_4_m_min4(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x);
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
        hlm += hl_5_m_5(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x);
      }
      return (hlm);
    case 4:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_5_m_4(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z);
      }
      return (hlm);
    case 3:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_5_m_3(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x);
      }
      return (hlm);
    case 2:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_5_m_2(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z);
      }
      return (hlm);
    case 1:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_5_m_1(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x);
      }
      return (hlm);
    case 0:
      return (0);
    case -1:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_5_m_min1(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x);
      }
      return (hlm);
    case -2:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_5_m_min2(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z);
      }
      return (hlm);
    case -3:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_5_m_min3(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x);
      }
      return (hlm);
    case -4:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_5_m_min4(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z);
      }
      return (hlm);
    case -5:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_5_m_min5(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x);
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
        hlm += hl_6_m_6(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z);
      }
      return (hlm);
    case 5:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_6_m_5(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z);
      }
      return (hlm);
    case 4:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_6_m_4(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z);
      }
      return (hlm);
    case 3:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_6_m_3(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z);
      }
      return (hlm);
    case 2:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_6_m_2(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z);
      }
      return (hlm);
    case 1:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_6_m_1(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z);
      }
      return (hlm);
    case 0:
      return (0);
    case -1:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_6_m_min1(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z);
      }
      return (hlm);
    case -2:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_6_m_min2(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z);
      }
      return (hlm);
    case -3:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_6_m_min3(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z);
      }
      return (hlm);
    case -4:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_6_m_min4(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z);
      }
      return (hlm);
    case -5:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_6_m_min5(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z);
      }
      return (hlm);
    case -6:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_6_m_min6(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z);
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
        hlm += hl_7_m_7(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno);
      }
      return (hlm);
    case 6:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_7_m_6(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z);
      }
      return (hlm);
    case 5:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_7_m_5(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno);
      }
      return (hlm);
    case 4:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_7_m_4(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z);
      }
      return (hlm);
    case 3:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_7_m_3(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno);
      }
      return (hlm);
    case 2:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_7_m_2(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z);
      }
      return (hlm);
    case 1:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_7_m_1(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno);
      }
      return (hlm);
    case 0:
      return (0);
    case -1:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_7_m_min1(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno);
      }
      return (hlm);
    case -2:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_7_m_min2(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z);
      }
      return (hlm);
    case -3:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_7_m_min3(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno);
      }
      return (hlm);
    case -4:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_7_m_min4(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z);
      }
      return (hlm);
    case -5:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_7_m_min5(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno);
      }
      return (hlm);
    case -6:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_7_m_min6(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z);
      }
      return (hlm);
    case -7:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_7_m_min7(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno);
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
        hlm += hl_8_m_8(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno);
      }
      return (hlm);
    case 7:
      return (0);
    case 6:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_8_m_6(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno);
      }
      return (hlm);
    case 5:
      return (0);
    case 4:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_8_m_4(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno);
      }
      return (hlm);
    case 3:
      return (0);
    case 2:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_8_m_2(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno);
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
        hlm += hl_8_m_min2(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno);
      }
      return (hlm);
    case -3:
      return (0);
    case -4:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_8_m_min4(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno);
      }
      return (hlm);
    case -5:
      return (0);
    case -6:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_8_m_min6(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno);
      }
      return (hlm);
    case -7:
      return (0);
    case -8:
      for (INT4 pno = vpnorder; pno >= 0; pno--) {
        hlm += hl_8_m_min8(mass, Nu, r, rDOT, Phi, PhiDOT, R, pno);
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