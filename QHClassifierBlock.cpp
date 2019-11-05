/*
   Copyright 2019 Anton Runov <anton.runov.oca@gmail.com>

   This file is part of qhear.

   Qhear is free software: you can redistribute it and/or modify it under the
   terms of the GNU Affero General Public License as published by the Free
   Software Foundation, either version 3 of the License, or (at your option)
   any later version.

   Qhear is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
   FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for
   more details.

   You should have received a copy of the GNU Affero General Public License
   along with qhear.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "QHClassifierBlock.h"
#include <string.h>
#include <cmath>
#include "qhlog.h"

// ----------------------------------------------------------------------------

QHClassifierBlock::QHClassifierBlock()
{
  m_prmVecLen = 0;
  m_order = 0;
  m_minDt = 0;
  m_maxDt = 0;
  m_maxImgIdx = 0;
  m_A = NAN;
}

// ----------------------------------------------------------------------------

bool QHClassifierBlock::init(const double* H, int prmVecLen, int len, int bandwidth, double min_f, double max_f)
{
  if (14 != prmVecLen) {
    qhlog_error("error, prmVecLen should be 14\n");
    return false;
  }
  double A = 0;
  int max_hrm = 0;
  double dt_min = 0;
  double dt_max = 0;
  double df_min = 0;
  double df_max = 0;
  int img_idx_max = 0;

  for (int jh=0; jh<len; ++jh) {
    const double* prm = H + prmVecLen*jh;
    A += prm[11];
    max_hrm = std::max(max_hrm, (int)prm[12]);
    img_idx_max = std::max(img_idx_max, (int)prm[8]);
    
    dt_min = std::min(dt_min, prm[0]);
    dt_max = std::max(dt_max, prm[1]);
    df_min = std::min(df_min, prm[2]);
    df_max = std::max(df_max, prm[3]);
    dt_min = std::min(dt_min, prm[4]);
    dt_max = std::max(dt_max, prm[5]);
    df_min = std::min(df_min, prm[6]);
    df_max = std::max(df_max, prm[7]);
  }
  if (0 == bandwidth) {
    // Step 2
  }
  else if (df_min < -bandwidth || df_max > bandwidth) {
    qhlog_error("error, df out of bandwidth\n");
    return false;
  }
  m_H = std::vector<double>(len*prmVecLen);
  memcpy(m_H.data(), H, sizeof(double)*len*prmVecLen);
  
  for (int jh=0; jh<len; ++jh) {
    double* prm = m_H.data() + prmVecLen*jh;
    prm[11] /= A;
  }

  m_qhFeature.init(bandwidth);
  m_prmVecLen = prmVecLen;
  m_order = max_hrm+1;
  m_minDt = std::round(dt_min);
  m_maxDt = std::round(dt_max);
  m_minDf = std::round(df_min);
  m_maxDf = std::round(df_max);
  m_maxImgIdx = img_idx_max;
  m_minF = min_f;
  m_maxF = max_f;
  m_A = A;
  //qhlog_debug("A = %f\n", m_A);
  return true;
}

// ----------------------------------------------------------------------------

bool QHClassifierBlock::checkFreq(double f, const QHIimg& iimg) const
{
  double minF;
  double maxF;

  if (isStep2()) {
    minF = std::max(m_minF, (double)-m_minDf);
    maxF = std::min(m_maxF, iimg.lenF() - m_maxDf-1.0);
  }
  else {
    minF = std::max(m_minF, 5.0 * (m_qhFeature.bandwidth()+1) );
    maxF = std::min(m_maxF, 5.0 * ( (iimg.lenF()-2.0)/(m_order+1) ) );
  }
  
  if (f < minF || f > maxF) {
    return false;
  }
  return true;
}

// ----------------------------------------------------------------------------

double QHClassifierBlock::classify(const QHIimg& iimg, double t, double f, int out_scale_mode) const
{
  int minT = -m_minDt + 1;
  int maxT = iimg.lenT() - m_maxDt -1;
  if (t < minT || t > maxT) {
    //qhlog_debug("skipping: %d %f\n", t, f);
    if (1 == out_scale_mode) {
      return -6.0;
    }
    else if (2 == out_scale_mode) {
      return -3.0;
    }
    return 0.0;
  }
  int NH = m_H.size() / m_prmVecLen;
  const double* prm_base = m_H.data();

  double x = 0;
  for (int jh=0; jh < NH; ++jh) {
    const double* prm = prm_base + m_prmVecLen*jh;
    double alpha = prm[11];
    double p = prm[9];
    double theta = prm[10];

    double v = m_qhFeature.calc(iimg, t, f, prm);
    v = (v*p < theta*p) ? 1 : 0;
    x += v*alpha;
  }
  if (1 == out_scale_mode) {
    x = (x-0.5)*m_A;
    x = std::max(-6.0, std::min(x, 16.0));
  }
  else if (2 == out_scale_mode) {
    x = (x-0.5)*sqrt(m_A);
    x = std::max(-3.0, std::min(x, 8.0));
  }
  return x;
}

// ----------------------------------------------------------------------------

