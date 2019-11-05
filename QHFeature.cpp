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

#include "QHFeature.h"
#include "QHIimg.h"
#include <cmath>

// ----------------------------------------------------------------------------

QHFeature::QHFeature()
{
  m_bandwidth = 0;
  m_k0s = 0;
  m_k1s = 0;
}

// ----------------------------------------------------------------------------

void QHFeature::init(int bandwidth)
{
  m_bandwidth = bandwidth;
  if (0 < bandwidth) {
    m_k0s = std::pow(1.0/m_bandwidth, 2);
    m_k1s = 1.0/5/m_bandwidth;
  }
  else {
    m_k0s = 0;
    m_k1s = 0;
  }
}

// ---------------------------------------------------------------------------

inline static long stretch_bounds0(double f, double df, double k)
{
  return round(f + df*k);
}

// ---------------------------------------------------------------------------

inline static long stretch_bounds1(double f, double df, double k)
{
  double tmp = df*df;
  return round( f + df + k*df*tmp );
}

// ----------------------------------------------------------------------------

double  QHFeature::calc(const QHIimg& iimg, double t, double f, const double* prm) const
{
  long frq[4];
  if (0 < m_bandwidth) {
    if (0.0 > prm[13]) {
      double k = m_k0s * (f*m_k1s - 1);
      f = f*(1+prm[12])/5;
      frq[0] = stretch_bounds1(f,prm[2],k);
      frq[1] = stretch_bounds1(f,prm[3],k);
      frq[2] = stretch_bounds1(f,prm[6],k);
      frq[3] = stretch_bounds1(f,prm[7],k);
    }
    else {
      double k = 1.0 + (f*m_k1s-1.0)*prm[13];
      f = f*(1+prm[12])/5;
      frq[0] = stretch_bounds0(f,prm[2],k);
      frq[1] = stretch_bounds0(f,prm[3],k);
      frq[2] = stretch_bounds0(f,prm[6],k);
      frq[3] = stretch_bounds0(f,prm[7],k);
    }
  }
  else {
      frq[0] = f + prm[2];
      frq[1] = f + prm[3];
      frq[2] = f + prm[6];
      frq[3] = f + prm[7];
  }
  
  t -= 1;
  int img_idx = prm[8];
  if (iimg.numLayers() > img_idx) {
    double x = iimg.rectSum(img_idx, t+prm[0], t+prm[1], frq[0], frq[1] );
    x -= iimg.rectSum(img_idx, t+prm[4], t+prm[5], frq[2], frq[3] );
    return iimg.unscaleDiff(img_idx, x);
  }
  else if (2*iimg.numLayers() > img_idx){
    img_idx -= iimg.numLayers();
    double x = iimg.rectSum(img_idx, t+prm[0], t+prm[1], frq[0], frq[1] );
    return iimg.unscale(img_idx, x);
  }
  return NAN;
}

// ----------------------------------------------------------------------------

double QHFeature::calcPlain(const QHIimg& iimg, double t, double f, const double* prm)
{
  t -= 1;
  int img_idx = prm[8];
  if (iimg.numLayers() > img_idx) {
    double x = iimg.rectSum(img_idx, t+prm[0], t+prm[1], f+prm[2], f+prm[3] );
    x -= iimg.rectSum(img_idx, t+prm[4], t+prm[5], f+prm[6], f+prm[7] );
    return iimg.unscaleDiff(img_idx, x);
  }
  else if (2*iimg.numLayers() > img_idx){
    img_idx -= iimg.numLayers();
    double x = iimg.rectSum(img_idx, t+prm[0], t+prm[1], f+prm[2], f+prm[3] );
    return iimg.unscale(img_idx, x);
  }
  return NAN;
}

// ----------------------------------------------------------------------------

