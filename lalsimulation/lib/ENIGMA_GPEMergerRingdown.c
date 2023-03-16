static int LoadGPEhyperparameterFromFile(LALH5File *file, const char *dsetname,
                                         int *AmpIndex, REAL8 *AmpLength,
                                         int *PhaseIndex, REAL8 *PhaseLength);

// Function to return a new, empty training set with Loaded = 0
static TrainingSet *CreateEmptyTrainingSet() {
  static TrainingSet Dset;
  Dset.Loaded = 0;
  return &Dset;
}

// Load GPE hyperparameter from file
// The current file convention has 6 lines, 1 number per line: AmpIndex,
// AmpLength, AmpSigmaF, PhaseIndex, PhaseLength, PhaseSigmaF
static int LoadGPEhyperparameterFromFile(LALH5File *file, const char *dsetname,
                                         int *AmpIndex, REAL8 *AmpLength,
                                         int *PhaseIndex, REAL8 *PhaseLength) {
  int errorcode = XLAL_SUCCESS;

  LALH5Dataset *cov_dset = NULL;
  REAL8Vector *cov_data = NULL;

  cov_dset = XLALH5DatasetRead(file, dsetname);
  if (cov_dset == NULL) {
    errorcode = XLAL_EFUNC;
    XLAL_ERROR_FAIL(errorcode);
  }
  cov_data = XLALH5DatasetReadREAL8Vector(cov_dset);
  if (cov_dset == NULL) {
    errorcode = XLAL_EFUNC;
    XLAL_ERROR_FAIL(errorcode);
  }
  if (cov_data->length != 6) {
    errorcode = XLAL_EIO;
    XLAL_ERROR_FAIL(errorcode,
                    "Incorrect File Format. Expect 6 elements but got %d",
                    (int)cov_data->length);
  }

  // AmpIndex and PhaseIndex are actually int but are stored as double in the
  // HDF5 file, which is fine as long as the int is not too big and can be
  // exactly represented.
  // Could have also stored these as attributes of the root group "/" in the
  // file.
  *AmpIndex = (int)cov_data->data[0];
  *AmpLength = cov_data->data[1];
  *PhaseIndex = (int)cov_data->data[3];
  *PhaseLength = cov_data->data[4];

XLAL_FAIL:
  XLALDestroyREAL8Vector(cov_data);
  XLALH5DatasetFree(cov_dset);

  return errorcode;
}

// Load Surrogate waveform from file
// Include period -1000M < t < 100M around peak amplitude
// Resample to required time values
// Renormalise s.t. peak amplitude is 1.0
static int LoadSurrogateWaveformFile(REAL8 **hp, REAL8 **hc, LALH5File *file,
                                     const char *dsetname) {

  int errorcode = XLAL_SUCCESS;

  assert(hp && hc);
  assert(!*hp && !*hc);

  LALH5Dataset *hp_dset = NULL;
  REAL8Vector *hp_data = NULL;
  LALH5Dataset *hc_dset = NULL;
  REAL8Vector *hc_data = NULL;

  (*hp) = (REAL8 *)LALMalloc(4300 * sizeof(REAL8));
  (*hc) = (REAL8 *)LALMalloc(4300 * sizeof(REAL8));
  if (*hp == NULL || *hc == NULL) {
    errorcode = XLAL_EFUNC;
    XLAL_ERROR_FAIL(errorcode);
  }

  char buf[1024];
  size_t written;

  written = snprintf(buf, sizeof(buf), "%s.hp", dsetname);
  if (written >= sizeof(buf)) {
    errorcode = XLAL_EFAILED;
    XLAL_ERROR_FAIL(errorcode);
  }
  hp_dset = XLALH5DatasetRead(file, buf);
  hp_data = XLALH5DatasetReadREAL8Vector(hp_dset);
  if (hp_data->length != 4300) {
    errorcode = XLAL_EIO;
    XLAL_ERROR_FAIL(errorcode,
                    "Incorrect File Format. Expect 4300 elements but got %d",
                    (int)hp_data->length);
  }
  memcpy(*hp, hp_data->data, sizeof(hp_data->data[0]) * hp_data->length);

  written = snprintf(buf, sizeof(buf), "%s.hc", dsetname);
  if (written >= sizeof(buf)) {
    errorcode = XLAL_EFAILED;
    XLAL_ERROR_FAIL(errorcode);
  }
  hc_dset = XLALH5DatasetRead(file, buf);
  hc_data = XLALH5DatasetReadREAL8Vector(hc_dset);
  if (hc_data->length != 4300) {
    errorcode = XLAL_EIO;
    XLAL_ERROR_FAIL(errorcode,
                    "Incorrect File Format. Expect 4300 elements but got %d",
                    (int)hc_data->length);
  }
  memcpy(*hc, hc_data->data, sizeof(hc_data->data[0]) * hc_data->length);

XLAL_FAIL:
  XLALDestroyREAL8Vector(hp_data);
  XLALH5DatasetFree(hp_dset);

  XLALDestroyREAL8Vector(hc_data);
  XLALH5DatasetFree(hc_dset);

  return errorcode;
}

// Convert waveform to amplitude-phase representation
// On input the arrays T is time, and A and B are assumed to be the plus and
// cross polarisation states of the waveform On output the arrays A and B are
// the amplitude and phase of the waveforms After merger, when the amplitude
// drops below 0.01 we set the phase to a constant
static void WaveToAmpPhase(REAL8 *T, REAL8 *A, REAL8 *B, int length,
                           REAL8 *const_phase_array) {
  REAL8 t, amp, phi = 0.0, phi_previous = -10.0;
  int i, branch = 0;
  for (i = 0; i < length; i++) {
    t = T[i];
    amp = sqrt(A[i] * A[i] + B[i] * B[i]);
    phi = atan2(B[i], A[i]);
    if (phi < phi_previous) {
      branch = branch + 1;
    }
    A[i] = amp;
    if (t > 0 && amp < 0.0001) {
      B[i] = B[i - 1];
    } else {
      B[i] = phi + 2 * branch * M_PI;
    }
    phi_previous = phi;
  }

  // Subtract off the phase at merger, so that phi(t=0)=0
  REAL8 maxamp = 0.0;
  for (i = 0; i < length; i++) {
    if (A[i] > maxamp) {
      phi = B[i];
      maxamp = A[i];
    }
  } // phi is now the (unwrapped) phase at merger

  for (i = 0; i < length; i++) {
    B[i] = B[i] - phi;
  }

  // Now subtract off the reference phase
  if (length != 4300) {
    XLAL_ERROR_VOID(XLAL_EINVAL);
  }
  for (i = 0; i < length; i++) {
    B[i] = (B[i] - const_phase_array[i]);
  }
}

// Convert amplitude-phase representation to waveform
// On input the arrays A and B are assumed to be the amplitude and phase of the
// waveforms On output the arrays A and B are the plus and cross polarisation
// states of the waveform
static void AmpPhaseToWave(REAL8 *A, REAL8 *B, int length,
                           REAL8 *const_phase_array) {
  int i;

  // Now add on the reference phase
  if (length != 4300) {
    XLAL_ERROR_VOID(XLAL_EINVAL);
  }
  for (i = 0; i < length; i++) {
    B[i] = (B[i] + const_phase_array[i]);
  }

  REAL8 hp, hc;
  for (i = 0; i < length; i++) {
    hp = A[i] * cos(B[i]);
    hc = A[i] * sin(B[i]);
    A[i] = hp;
    B[i] = hc;
  }
}

// The first time, and only the first time, the function void interp(REAL8
// *tsamples, TrainingSet *Dset) is called it should run this function which
// loads the training set into memeory and prepares other quantities such as the
// inverse covariance matrix which will be needed for all subsequent
// interpolation calls
static int setup(TrainingSet *Dset) {

  int errorcode = XLAL_SUCCESS;

  assert(Dset != NULL);
  assert(Dset->time == NULL);
  assert(Dset->xvals == NULL);
  assert(Dset->Waves == NULL);
  assert(Dset->const_phase_array == NULL);
  LALH5File *file = NULL;

  // Number of training set points: 19 NR + 36 Surrogate + 14 Surrogate Cluster
  Dset->N_Dset_points = 70;

  // The root file name for the GW_data folder
  char *enigma = XLALFileResolvePathLong("ENIGMA.h5", PKG_DATA_DIR);
  if (enigma == NULL)
    XLAL_ERROR_FAIL(XLAL_EIO,
                    "Unable to resolve data file %s in $LAL_DATA_PATH",
                    "ENIGMA.h5");
  file = XLALH5FileOpen(enigma, "r");
  if (file == NULL)
    XLAL_ERROR_FAIL(XLAL_EFUNC);

  // Set the time sampling values
  Dset->time = (REAL8 *)LALMalloc(4300 * sizeof(REAL8));
  Dset->length = 4300;
  for (int i = 0; i < 1500; i++) {
    Dset->time[i] = -2500.0 + 1.0 * i;
  } // for the early inspiral use steps of deltaT = 0.5M
  for (int i = 1500; i < 3300; i++) {
    Dset->time[i] = -1000.0 + 0.5 * (i - 1500);
  } // for the late inspiral use steps of deltaT = 0.5M
  for (int i = 3300; i < 4300; i++) {
    Dset->time[i] = -100.0 + 0.2 * (i - 3300);
  } // for the merger use steps of deltaT = 0.2M

  // Set the training set parameter values
  Dset->xvals = (REAL8 *)LALMalloc(Dset->N_Dset_points * sizeof(REAL8));
  if (Dset->xvals == NULL) {
    errorcode = XLAL_ENOMEM;
    XLAL_ERROR_FAIL(errorcode);
  }

  // Malloc space for the training set waveforms
  Dset->Waves = (wave *)LALCalloc(Dset->N_Dset_points, sizeof(wave));
  if (Dset->Waves == NULL) {
    errorcode = XLAL_ENOMEM;
    XLAL_ERROR_FAIL(errorcode);
  }

  // Reference waveform phase
  {
    Dset->const_phase_array = (REAL8 *)LALCalloc(4300, sizeof(REAL8));
    if (Dset->const_phase_array == NULL) {
      errorcode = XLAL_ENOMEM;
      XLAL_ERROR_FAIL(errorcode);
    }
    REAL8 *const_hp = NULL;
    REAL8 *const_hc = NULL;
    errorcode = LoadSurrogateWaveformFile(&const_hp, &const_hc, file,
                                          "m_2_m_2_q1.0_f15");
    if (errorcode != XLAL_SUCCESS)
      XLAL_ERROR_FAIL(XLAL_EFUNC);
    WaveToAmpPhase(Dset->time, const_hp, const_hc, Dset->length,
                   Dset->const_phase_array);
    for (int i = 0; i < 4300; i++) {
      Dset->const_phase_array[i] = const_hc[i];
    }
    LALFree(const_hp);
    LALFree(const_hc);
  }

  // Surrogate waveforms - This is a denser set of simulations evenly placed in
  // log(q)
  REAL8 qValuesSurrogateLogq[70] = {
      10.00, 9.67, 9.35, 9.05, 8.75, 8.46, 8.19, 7.92, 7.66, 7.41, 7.16, 6.93,
      6.70,  6.48, 6.27, 6.06, 5.86, 5.67, 5.48, 5.30, 5.13, 4.96, 4.80, 4.64,
      4.49,  4.34, 4.20, 4.06, 3.93, 3.80, 3.67, 3.55, 3.44, 3.32, 3.22, 3.11,
      3.01,  2.91, 2.81, 2.72, 2.63, 2.55, 2.46, 2.38, 2.30, 2.23, 2.15, 2.08,
      2.02,  1.95, 1.89, 1.82, 1.76, 1.71, 1.65, 1.60, 1.54, 1.49, 1.44, 1.40,
      1.35,  1.31, 1.26, 1.22, 1.18, 1.14, 1.11, 1.07, 1.03, 1.00};
  for (int i = 0; i < 70; i++) {
    Dset->xvals[i] = qValuesSurrogateLogq[i];
  }
  const char *SurrogateLogqfiles[70];
  SurrogateLogqfiles[0] = "m_2_m_2_q10.0_f15";
  SurrogateLogqfiles[1] = "m_2_m_2_q9.67_f15";
  SurrogateLogqfiles[2] = "m_2_m_2_q9.35_f15";
  SurrogateLogqfiles[3] = "m_2_m_2_q9.05_f15";
  SurrogateLogqfiles[4] = "m_2_m_2_q8.75_f15";
  SurrogateLogqfiles[5] = "m_2_m_2_q8.46_f15";
  SurrogateLogqfiles[6] = "m_2_m_2_q8.19_f15";
  SurrogateLogqfiles[7] = "m_2_m_2_q7.92_f15";
  SurrogateLogqfiles[8] = "m_2_m_2_q7.66_f15";
  SurrogateLogqfiles[9] = "m_2_m_2_q7.41_f15";

  SurrogateLogqfiles[10] = "m_2_m_2_q7.16_f15";
  SurrogateLogqfiles[11] = "m_2_m_2_q6.93_f15";
  SurrogateLogqfiles[12] = "m_2_m_2_q6.7_f15";
  SurrogateLogqfiles[13] = "m_2_m_2_q6.48_f15";
  SurrogateLogqfiles[14] = "m_2_m_2_q6.27_f15";
  SurrogateLogqfiles[15] = "m_2_m_2_q6.06_f15";
  SurrogateLogqfiles[16] = "m_2_m_2_q5.86_f15";
  SurrogateLogqfiles[17] = "m_2_m_2_q5.67_f15";
  SurrogateLogqfiles[18] = "m_2_m_2_q5.48_f15";
  SurrogateLogqfiles[19] = "m_2_m_2_q5.3_f15";

  SurrogateLogqfiles[20] = "m_2_m_2_q5.13_f15";
  SurrogateLogqfiles[21] = "m_2_m_2_q4.96_f15";
  SurrogateLogqfiles[22] = "m_2_m_2_q4.8_f15";
  SurrogateLogqfiles[23] = "m_2_m_2_q4.64_f15";
  SurrogateLogqfiles[24] = "m_2_m_2_q4.49_f15";
  SurrogateLogqfiles[25] = "m_2_m_2_q4.34_f15";
  SurrogateLogqfiles[26] = "m_2_m_2_q4.2_f15";
  SurrogateLogqfiles[27] = "m_2_m_2_q4.06_f15";
  SurrogateLogqfiles[28] = "m_2_m_2_q3.93_f15";
  SurrogateLogqfiles[29] = "m_2_m_2_q3.8_f15";

  SurrogateLogqfiles[30] = "m_2_m_2_q3.67_f15";
  SurrogateLogqfiles[31] = "m_2_m_2_q3.55_f15";
  SurrogateLogqfiles[32] = "m_2_m_2_q3.44_f15";
  SurrogateLogqfiles[33] = "m_2_m_2_q3.32_f15";
  SurrogateLogqfiles[34] = "m_2_m_2_q3.22_f15";
  SurrogateLogqfiles[35] = "m_2_m_2_q3.11_f15";
  SurrogateLogqfiles[36] = "m_2_m_2_q3.01_f15";
  SurrogateLogqfiles[37] = "m_2_m_2_q2.91_f15";
  SurrogateLogqfiles[38] = "m_2_m_2_q2.81_f15";
  SurrogateLogqfiles[39] = "m_2_m_2_q2.72_f15";

  SurrogateLogqfiles[40] = "m_2_m_2_q2.63_f15";
  SurrogateLogqfiles[41] = "m_2_m_2_q2.55_f15";
  SurrogateLogqfiles[42] = "m_2_m_2_q2.46_f15";
  SurrogateLogqfiles[43] = "m_2_m_2_q2.38_f15";
  SurrogateLogqfiles[44] = "m_2_m_2_q2.3_f15";
  SurrogateLogqfiles[45] = "m_2_m_2_q2.23_f15";
  SurrogateLogqfiles[46] = "m_2_m_2_q2.15_f15";
  SurrogateLogqfiles[47] = "m_2_m_2_q2.08_f15";
  SurrogateLogqfiles[48] = "m_2_m_2_q2.02_f15";
  SurrogateLogqfiles[49] = "m_2_m_2_q1.95_f15";

  SurrogateLogqfiles[50] = "m_2_m_2_q1.89_f15";
  SurrogateLogqfiles[51] = "m_2_m_2_q1.82_f15";
  SurrogateLogqfiles[52] = "m_2_m_2_q1.76_f15";
  SurrogateLogqfiles[53] = "m_2_m_2_q1.71_f15";
  SurrogateLogqfiles[54] = "m_2_m_2_q1.65_f15";
  SurrogateLogqfiles[55] = "m_2_m_2_q1.6_f15";
  SurrogateLogqfiles[56] = "m_2_m_2_q1.54_f15";
  SurrogateLogqfiles[57] = "m_2_m_2_q1.49_f15";
  SurrogateLogqfiles[58] = "m_2_m_2_q1.44_f15";
  SurrogateLogqfiles[59] = "m_2_m_2_q1.4_f15";

  SurrogateLogqfiles[60] = "m_2_m_2_q1.35_f15";
  SurrogateLogqfiles[61] = "m_2_m_2_q1.31_f15";
  SurrogateLogqfiles[62] = "m_2_m_2_q1.26_f15";
  SurrogateLogqfiles[63] = "m_2_m_2_q1.22_f15";
  SurrogateLogqfiles[64] = "m_2_m_2_q1.18_f15";
  SurrogateLogqfiles[65] = "m_2_m_2_q1.14_f15";
  SurrogateLogqfiles[66] = "m_2_m_2_q1.11_f15";
  SurrogateLogqfiles[67] = "m_2_m_2_q1.07_f15";
  SurrogateLogqfiles[68] = "m_2_m_2_q1.03_f15";
  SurrogateLogqfiles[69] = "m_2_m_2_q1.0_f15";

  for (int i = 0; i < 70; i++) {
    errorcode = LoadSurrogateWaveformFile(&(((Dset->Waves)[i]).hp),
                                          &(((Dset->Waves)[i]).hc), file,
                                          SurrogateLogqfiles[i]);
    if (errorcode != XLAL_SUCCESS)
      XLAL_ERROR_FAIL(XLAL_EFUNC);
    WaveToAmpPhase(Dset->time, Dset->Waves[i].hp, Dset->Waves[i].hc,
                   Dset->length, Dset->const_phase_array);
  }

  // Convert mass ratios into log q
  for (int i = 0; i < Dset->N_Dset_points; i++) {
    Dset->xvals[i] = log(Dset->xvals[i]);
  }

  // Covariance functions and hyperparameters
  errorcode = LoadGPEhyperparameterFromFile(
      file, "CovHyp.dat", &(Dset->cov_index_amp), &(Dset->Log10DeltaAmp),
      &(Dset->cov_index_phase), &(Dset->Log10DeltaPhase));
  if (errorcode != XLAL_SUCCESS)
    XLAL_ERROR_FAIL(XLAL_EFUNC);

  // Setup inverse covariance matrix and other args for interp function
  GPR_InverseCovarianceMatrix(&(Dset->inv_cov_amp), Dset->xvals,
                              Dset->N_Dset_points, Dset->Log10DeltaAmp,
                              Dset->cov_index_amp);
  GPR_InverseCovarianceMatrix(&(Dset->inv_cov_phase), Dset->xvals,
                              Dset->N_Dset_points, Dset->Log10DeltaPhase,
                              Dset->cov_index_phase);

  // Set the loaded int to one
  (*Dset).Loaded = 1;

XLAL_FAIL:
  if (errorcode != XLAL_SUCCESS) {
    for (int i = 0; i < 70; i++) {
      XLALFree(Dset->Waves[i].hp);
      XLALFree(Dset->Waves[i].hc);
    }
    XLALFree(Dset->Waves);
    Dset->Waves = NULL;
    XLALFree(Dset->time);
    Dset->time = NULL;
    XLALFree(Dset->xvals);
    Dset->xvals = NULL;
    XLALFree(Dset->const_phase_array);
    Dset->const_phase_array = NULL;
  }
  XLALH5FileClose(file);
  XLALFree(enigma);

  return errorcode;
}

// Interpolation
static int interp(REAL8 eta, REAL8 **hp, REAL8 **hc, REAL8 **tsamples,
                  int *len_tsamples, TrainingSet *Dset) {
  int errorcode = XLAL_SUCCESS;
  REAL8 *kstar_amp = NULL;
  REAL8 *kstar_phase = NULL;
  REAL8 *coeff_amp = NULL;
  REAL8 *coeff_phase = NULL;
  REAL8 *Hp = NULL;
  REAL8 *Hc = NULL;
  REAL8 *t = NULL;

  assert(hp && hc);
  assert(*hp == NULL && *hc == NULL);
  assert(tsamples);
  assert(len_tsamples);
  assert(Dset);

  // The function setup() is a time consuming step. The function setup() changes
  // the value of Dset.Loaded to 1, so setup() is only called the first time
  // this function is run.
  if ((*Dset).Loaded == 0) {
    errorcode = setup(Dset);
    if (errorcode != XLAL_SUCCESS)
      XLAL_ERROR_FAIL(XLAL_EFUNC);
  }

  // Compute the covariance between q and and all the Training set points
  kstar_amp = (REAL8 *)LALMalloc(Dset->N_Dset_points * sizeof(REAL8));
  kstar_phase = (REAL8 *)LALMalloc(Dset->N_Dset_points * sizeof(REAL8));
  if (kstar_amp == NULL || kstar_phase == NULL) {
    errorcode = XLAL_ENOMEM;
    XLAL_ERROR_FAIL(errorcode);
  }
  for (int i = 0; i < Dset->N_Dset_points; i++) {
    kstar_amp[i] = Covariance(-log(SmallMassRatio(eta)), (Dset->xvals)[i],
                              Dset->Log10DeltaAmp, 0.0, Dset->cov_index_amp);
    kstar_phase[i] =
        Covariance(-log(SmallMassRatio(eta)), (Dset->xvals)[i],
                   Dset->Log10DeltaPhase, 0.0, Dset->cov_index_phase);
    if (kstar_amp[i] == XLAL_REAL8_FAIL_NAN ||
        kstar_phase[i] == XLAL_REAL8_FAIL_NAN) {
      errorcode = XLAL_EFUNC;
      XLAL_ERROR_FAIL(errorcode);
    }
  }

  // Compute the coefficients kstar.inv_cov
  coeff_amp = (REAL8 *)LALCalloc(Dset->N_Dset_points, sizeof(REAL8));
  coeff_phase = (REAL8 *)LALCalloc(Dset->N_Dset_points, sizeof(REAL8));
  if (coeff_amp == NULL || coeff_phase == NULL) {
    errorcode = XLAL_ENOMEM;
    XLAL_ERROR_FAIL(errorcode);
  }
  for (int i = 0; i < Dset->N_Dset_points; i++) {
    for (int j = 0; j < Dset->N_Dset_points; j++) {
      coeff_amp[i] +=
          (Dset->inv_cov_amp)[i * (Dset->N_Dset_points) + j] * kstar_amp[j];
      coeff_phase[i] +=
          (Dset->inv_cov_phase)[i * (Dset->N_Dset_points) + j] * kstar_phase[j];
    }
  }

  // Take linear combination of training set waveforms
  Hp = LALCalloc(Dset->length, sizeof(REAL8));
  Hc = LALCalloc(Dset->length, sizeof(REAL8));
  if (Hp == NULL || Hc == NULL) {
    errorcode = XLAL_ENOMEM;
    XLAL_ERROR_FAIL(errorcode);
  }
  for (int i = 0; i < Dset->length; i++) {
    for (int j = 0; j < Dset->N_Dset_points; j++) {
      Hp[i] += coeff_amp[j] * (Dset->Waves[j].hp)[i];
      Hc[i] += coeff_phase[j] * (Dset->Waves[j].hc)[i];
    }
  }

  // Training set waveforms are stored as amplitude and phase. Convert to h_plus
  // and h_cross
  AmpPhaseToWave(Hp, Hc, Dset->length, Dset->const_phase_array);

  t = LALMalloc(sizeof(*t) * Dset->length);
  if (t == NULL) {
    errorcode = XLAL_ENOMEM;
    XLAL_ERROR_FAIL(errorcode);
  }
  memcpy(t, Dset->time, sizeof(*t) * Dset->length);

  *hp = Hp;
  *hc = Hc;
  *tsamples = t;
  *len_tsamples = Dset->length;

XLAL_FAIL:
  // Free memory and exit
  LALFree(kstar_amp);
  LALFree(kstar_phase);
  LALFree(coeff_amp);
  LALFree(coeff_phase);
  if (errorcode != XLAL_SUCCESS) {
    LALFree(Hp);
    LALFree(Hc);
    LALFree(t);
  }

  return errorcode;
}

static REAL8 interp_for_root_finding(REAL8 x, void *params) {
  interp_params *P = (interp_params *)params;
  gsl_interp_accel *accHp = P->accHp, *accHc = P->accHc;
  gsl_interp *interpHp = P->interpHp, *interpHc = P->interpHc;
  REAL8 omega_attach = P->omega_attach;
  REAL8 *tsamples = P->tsamples;
  REAL8 *hp_merger = P->hp_merger;
  REAL8 *hc_merger = P->hc_merger;

  if (x > 0)
    x = 0.; // During the ringdown set frequency to a constant

  REAL8 hp_merger_dot_val =
      gsl_interp_eval_deriv(interpHp, tsamples, hp_merger, x, accHp);
  REAL8 hp_merger_val =
      gsl_interp_eval(interpHp, tsamples, hp_merger, x, accHp);
  REAL8 hc_merger_dot_val =
      gsl_interp_eval_deriv(interpHc, tsamples, hc_merger, x, accHc);
  REAL8 hc_merger_val =
      gsl_interp_eval(interpHc, tsamples, hc_merger, x, accHc);

  REAL8 freq =
      -omega_attach -
      (hp_merger_dot_val * hc_merger_val - hc_merger_dot_val * hp_merger_val) /
          (2.0 * (pow2(hp_merger_val) + pow2(hc_merger_val)));

  return freq;
}

int Attach_GPE_Merger_Ringdown(REAL8 dt, REAL8 **h_plus, REAL8 **h_cross,
                               REAL8 inspiral_matching_time,
                               REAL8 inspiral_matching_Hp,
                               REAL8 inspiral_matching_Hc, int *Length,
                               TrainingSet *Dset, REAL8 m1, REAL8 m2, REAL8 UNUSED S1z, REAL8 UNUSED S2z,
                               REAL8 inc, REAL8 *time_of_merger) {

  int errorcode = XLAL_SUCCESS;

  assert(h_plus && h_cross);
  assert(*h_plus && *h_cross);
  assert(Dset);

  gsl_interp_accel *accHp = NULL, *accHc = NULL;
  gsl_interp *interpHp = NULL, *interpHc = NULL;
  gsl_root_fsolver *s = NULL;

  // Generate GPE interpolation for merger ringdown
  REAL8 *tsamples = NULL;
  REAL8 *hp_merger = NULL;
  REAL8 *hc_merger = NULL;
  int GPELength = 0;
  errorcode = interp(SymMassRatio(m1 / m2), &hp_merger, &hc_merger, &tsamples,
                     &GPELength, Dset);
  if (errorcode != XLAL_SUCCESS)
    XLAL_ERROR_FAIL(errorcode);

  // Get the matching frequency
  REAL8 mr = m1 / m2;
  REAL8 omega_attach;
  errorcode = PN_Omega(mr, m1 + m2, NULL, &omega_attach);
  if (errorcode != XLAL_SUCCESS)
    XLAL_ERROR_FAIL(errorcode);

  int LengthOfInspiral = *Length;

  // Create interpolants for h_plus and h_cross
  accHp = gsl_interp_accel_alloc();
  accHc = gsl_interp_accel_alloc();
  interpHp = gsl_interp_alloc(gsl_interp_steffen, GPELength);
  interpHc = gsl_interp_alloc(gsl_interp_steffen, GPELength);
  gsl_interp_init(interpHp, tsamples, hp_merger, GPELength);
  gsl_interp_init(interpHc, tsamples, hc_merger, GPELength);

  // Create interpolation for frequency and store in interp_params structure
  interp_params InterpParams;
  InterpParams.tsamples = tsamples;
  InterpParams.hp_merger = hp_merger;
  InterpParams.hc_merger = hc_merger;
  InterpParams.interpHp = interpHp;
  InterpParams.accHp = accHp;
  InterpParams.interpHc = interpHc;
  InterpParams.accHc = accHc;
  InterpParams.omega_attach = omega_attach;

  // Find time in GPE waveform where f = fmatch
  gsl_function F;
  F.function = &interp_for_root_finding;
  F.params = &InterpParams;
  const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc(T);
  REAL8 r = 0, x_lo = tsamples[0], x_hi = tsamples[GPELength - 1];
  gsl_root_fsolver_set(s, &F, x_lo, x_hi);
  int status, iter = 0, max_iter = 200;
  do {
    iter++;
    status = gsl_root_fsolver_iterate(s);
    r = gsl_root_fsolver_root(s);
    x_lo = gsl_root_fsolver_x_lower(s);
    x_hi = gsl_root_fsolver_x_upper(s);
    status = gsl_root_test_interval(x_lo, x_hi, 0, 0.0001);
  } while (status == GSL_CONTINUE && iter < max_iter);
  REAL8 matching_time = r;

  REAL8 Hp_MERGER, Hc_MERGER;

  // The inclination factors
  REAL8 cos_inc = cos(inc);
  REAL8 one_plus_cos_squared_inc = 0.5 * (1. + cos_inc * cos_inc);

  // Phase difference and amplitude ratio
  int end = LengthOfInspiral - 1;
  Hp_MERGER =
      gsl_interp_eval(interpHp, tsamples, hp_merger, matching_time, accHp);
  Hc_MERGER =
      gsl_interp_eval(interpHc, tsamples, hc_merger, matching_time, accHc);
  REAL8 delta_phi =
      atan2(inspiral_matching_Hp, inspiral_matching_Hc) -
      atan2(one_plus_cos_squared_inc * Hp_MERGER, cos_inc * Hc_MERGER);
  REAL8 ratioAmplitude =
      sqrt(pow2(inspiral_matching_Hp) + pow2(inspiral_matching_Hc)) /
      sqrt(pow2(one_plus_cos_squared_inc * Hp_MERGER) +
           pow2(cos_inc * Hc_MERGER));

  // Join the Inspiral and Merger together
  int LengthOfMerger = floor((tsamples[GPELength - 1] - matching_time) / dt);
  if (LengthOfMerger <= 0) {
    errorcode = XLAL_EINVAL;
    XLAL_ERROR_FAIL(errorcode, "Invalid value for LengthOfMerger");
  }
  *Length += LengthOfMerger - 1; // -1 b/c I overwrite the last inspiral entry
  *h_plus = my_realloc(*h_plus, *Length * sizeof(**h_plus));
  *h_cross = my_realloc(*h_cross, *Length * sizeof(**h_cross));
  if (h_plus == NULL || h_cross == NULL) {
    errorcode = XLAL_ENOMEM;
    XLAL_ERROR_FAIL(errorcode);
  }
  if (inspiral_matching_time <= (end - 1) * dt ||
      inspiral_matching_time > end * dt) {
    errorcode = XLAL_EINVAL;
    XLAL_ERROR_FAIL(errorcode, "inspiral_matching_time <= (end - 1)*dt || "
                               "inspiral_matching_time > end*dt, failed");
  }

  for (int i = 0; i < LengthOfMerger; i++) {
    /* the waveforms are glued together at the times where omega =
     * omega_match corresponding to inspiral_matching_time and matching_time
     * respectively */
    REAL8 t_inspiral =
        (end + i) *
        dt; /* this is the time of the uniform sampling of our output */
    REAL8 t_merger = t_inspiral - inspiral_matching_time +
                     matching_time; /* correct for different times associated
                                       with omega_match */
    if (t_merger < matching_time) {
      errorcode = XLAL_EINVAL;
      XLAL_ERROR_FAIL(errorcode, "t_merger >= matching_time, failed");
    }

    Hp_MERGER = one_plus_cos_squared_inc *
                gsl_interp_eval(interpHp, tsamples, hp_merger, t_merger, accHp);
    Hc_MERGER = cos_inc *
                gsl_interp_eval(interpHc, tsamples, hc_merger, t_merger, accHc);

    (*h_plus)[end + i] = ratioAmplitude * (cos(delta_phi) * Hp_MERGER +
                                           sin(delta_phi) * Hc_MERGER);
    (*h_cross)[end + i] = ratioAmplitude * (cos(delta_phi) * Hc_MERGER -
                                            sin(delta_phi) * Hp_MERGER);
    if (i > 2) {
      REAL8 amp_i = ((*h_plus)[end + i] * (*h_plus)[end + i] +
                     (*h_cross)[end + i] * (*h_cross)[end + i]);
      REAL8 amp_im1 = ((*h_plus)[end + i - 1] * (*h_plus)[end + i - 1] +
                       (*h_cross)[end + i - 1] * (*h_cross)[end + i - 1]);
      REAL8 amp_im2 = ((*h_plus)[end + i - 2] * (*h_plus)[end + i - 2] +
                       (*h_cross)[end + i - 2] * (*h_cross)[end + i - 2]);

      if (amp_i < amp_im1 && amp_im1 > amp_im2) {
        *time_of_merger = dt * (REAL8) (end + i - 1);
      }
    }
  }

XLAL_FAIL:
  // Free up everything before exit
  LALFree(hp_merger);
  LALFree(hc_merger);
  LALFree(tsamples);
  if (interpHp)
    gsl_interp_free(interpHp);
  if (accHp)
    gsl_interp_accel_free(accHp);
  if (interpHc)
    gsl_interp_free(interpHc);
  if (accHc)
    gsl_interp_accel_free(accHc);
  if (s)
    gsl_root_fsolver_free(s);

  return errorcode;
}

static REAL8 GPR_CovarianceFunction_SE(REAL8 x1, REAL8 x2, REAL8 *theta,
                                       REAL8 Jitter) {
  REAL8 SigmaF = theta[0];
  REAL8 Delta = theta[1];
  REAL8 ans = pow2(SigmaF) * exp(-0.5 * pow2((x1 - x2) / Delta));
  if (fabs(x1 - x2) < 0.0001) {
    ans = ans + Jitter * SigmaF;
  }
  return ans;
}

static REAL8 GPR_CovarianceFunction_SE_BC(REAL8 x1, REAL8 x2, REAL8 *theta,
                                          REAL8 Jitter) {
  return (GPR_CovarianceFunction_SE(x1, x2, theta, Jitter) +
          GPR_CovarianceFunction_SE(x1, -x2, theta, Jitter));
}

static REAL8 GPR_CovarianceFunction_Wendland0(REAL8 x1, REAL8 x2, REAL8 *theta,
                                              REAL8 Jitter) {
  REAL8 SigmaF = theta[0];
  REAL8 Delta = theta[1];
  REAL8 x = fabs(x1 - x2) / Delta;
  REAL8 ans;
  int j = 0 + 1;
  if (x < 1.0) {
    ans = pow(1.0 - x, 1.0 * j);
  } else {
    ans = 0.0;
  }
  if (fabs(x1 - x2) < 0.0001) {
    ans = ans + Jitter * SigmaF;
  }
  return ans;
}

static REAL8 GPR_CovarianceFunction_Wendland0_BC(REAL8 x1, REAL8 x2,
                                                 REAL8 *theta, REAL8 Jitter) {
  return GPR_CovarianceFunction_Wendland0(x1, x2, theta, Jitter) +
         GPR_CovarianceFunction_Wendland0(x1, -x2, theta, Jitter);
}

static REAL8 GPR_CovarianceFunction_Wendland1(REAL8 x1, REAL8 x2, REAL8 *theta,
                                              REAL8 Jitter) {
  REAL8 SigmaF = theta[0];
  REAL8 Delta = theta[1];
  REAL8 x = fabs(x1 - x2) / Delta;
  REAL8 ans;
  int j = 1 + 1;
  if (x < 1.0) {
    ans = pow(1.0 - x, j + 1.0) * ((j + 1.0) * x + 1.0);
  } else {
    ans = 0.0;
  }
  if (fabs(x1 - x2) < 0.0001) {
    ans = ans + Jitter * SigmaF;
  }
  return ans;
}

static REAL8 GPR_CovarianceFunction_Wendland1_BC(REAL8 x1, REAL8 x2,
                                                 REAL8 *theta, REAL8 Jitter) {
  return GPR_CovarianceFunction_Wendland1(x1, x2, theta, Jitter) +
         GPR_CovarianceFunction_Wendland1(x1, -x2, theta, Jitter);
}

static REAL8 GPR_CovarianceFunction_Wendland2(REAL8 x1, REAL8 x2, REAL8 *theta,
                                              REAL8 Jitter) {
  REAL8 SigmaF = theta[0];
  REAL8 Delta = theta[1];
  REAL8 x = fabs(x1 - x2) / Delta;
  REAL8 ans;
  int j = 2 + 1;
  if (x < 1.0) {
    ans = pow(1.0 - x, j + 2.0) *
          ((j * j + 4.0 * j + 3.0) * x * x + (3.0 * j + 6.0) * x + 3.0) / 3.0;
  } else {
    ans = 0.0;
  }
  if (fabs(x1 - x2) < 0.0001) {
    ans = ans + Jitter * SigmaF;
  }
  return ans;
}

static REAL8 GPR_CovarianceFunction_Wendland2_BC(REAL8 x1, REAL8 x2,
                                                 REAL8 *theta, REAL8 Jitter) {
  return GPR_CovarianceFunction_Wendland2(x1, x2, theta, Jitter) +
         GPR_CovarianceFunction_Wendland2(x1, -x2, theta, Jitter);
}

static REAL8 GPR_CovarianceFunction_Wendland3(REAL8 x1, REAL8 x2, REAL8 *theta,
                                              REAL8 Jitter) {
  REAL8 SigmaF = theta[0];
  REAL8 Delta = theta[1];
  REAL8 x = fabs(x1 - x2) / Delta;
  REAL8 ans;
  int j = 3 + 1;
  if (x < 1.0) {
    ans = pow(1.0 - x, j + 3.0) *
          ((j * j * j + 9.0 * j * j + 23.0 * j + 15.0) * x * x * x +
           (6.0 * j * j + 36.0 * j + 45.0) * x * x + (15.0 * j + 45.0) * x +
           15.0) /
          15.0;
  } else {
    ans = 0.0;
  }
  if (fabs(x1 - x2) < 0.0001) {
    ans = ans + Jitter * SigmaF;
  }
  return ans;
}

static REAL8 GPR_CovarianceFunction_Wendland3_BC(REAL8 x1, REAL8 x2,
                                                 REAL8 *theta, REAL8 Jitter) {
  return GPR_CovarianceFunction_Wendland3(x1, x2, theta, Jitter) +
         GPR_CovarianceFunction_Wendland3(x1, -x2, theta, Jitter);
}

// Covariance function
static REAL8 Covariance(REAL8 x1, REAL8 x2, REAL8 log10Delta, REAL8 Jitter,
                        int covindex) {
  REAL8 cov = XLAL_REAL8_FAIL_NAN;

  REAL8 Theta[2];
  Theta[0] = 1.0;
  Theta[1] = pow(10.0, log10Delta);

  if (covindex == -1) {
    cov = GPR_CovarianceFunction_SE_BC(x1, x2, Theta, Jitter);
  } else if (covindex == 0) {
    cov = GPR_CovarianceFunction_Wendland0_BC(x1, x2, Theta, Jitter);
  } else if (covindex == 1) {
    cov = GPR_CovarianceFunction_Wendland1_BC(x1, x2, Theta, Jitter);
  } else if (covindex == 2) {
    cov = GPR_CovarianceFunction_Wendland2_BC(x1, x2, Theta, Jitter);
  } else if (covindex == 3) {
    cov = GPR_CovarianceFunction_Wendland3_BC(x1, x2, Theta, Jitter);
  }
  // Invalid covariance function
  else {
    XLAL_ERROR_REAL8(XLAL_EINVAL, "Invalid covariance function.");
  }
  return cov;
}

// Covariance matrix
static void GPR_CovarianceMatrix(REAL8 **covMatrix, REAL8 *xvals, int Npoints,
                                 REAL8 log10Delta, int covindex) {

  REAL8 Jitter = 0.0001;

  *covMatrix = (REAL8 *)LALCalloc(sizeof(REAL8), Npoints * Npoints);

  int i, j;
  for (i = 0; i < Npoints; i++) {
    for (j = 0; j < Npoints; j++) {
      if (i == j) {
        (*covMatrix)[i * Npoints + i] =
            Covariance(xvals[i], xvals[i], log10Delta, Jitter, covindex);
      } else {
        (*covMatrix)[i * Npoints + j] =
            Covariance(xvals[i], xvals[j], log10Delta, Jitter, covindex);
        (*covMatrix)[j * Npoints + i] = (*covMatrix)[i * Npoints + j];
      }
    }
  }

  return;
}

// Covariance matrix
static void GPR_InverseCovarianceMatrix(REAL8 **covMatrix, REAL8 *xvals,
                                        int Npoints, REAL8 log10Delta,
                                        int covindex) {

  GPR_CovarianceMatrix(covMatrix, xvals, Npoints, log10Delta, covindex);

  int i, j;
  gsl_matrix *m = gsl_matrix_alloc(Npoints, Npoints);
  for (i = 0; i < Npoints; i++) {
    for (j = 0; j < Npoints; j++) {
      gsl_matrix_set(m, i, j, (*covMatrix)[i * Npoints + j]);
    }
  }

  gsl_linalg_cholesky_decomp(m);
  gsl_linalg_cholesky_invert(m);

  for (i = 0; i < Npoints; i++) {
    for (j = 0; j < Npoints; j++) {
      (*covMatrix)[i * Npoints + j] = gsl_matrix_get(m, i, j);
    }
  }

  gsl_matrix_free(m);
}
