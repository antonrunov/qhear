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

#include "QHIimg.h"
#include "DataChunk.h"

#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <stdint.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <string.h>
#include <fftw3.h>
#include "qhlog.h"

// ----------------------------------------------------------------------------

QHIimg::QHIimg()
{
  m_L1 = 800;
  m_L2 = 0;
  m_extraMarginLeft = 0.0;
  m_extraMarginRight = 1.0;
}

// ----------------------------------------------------------------------------

QHIimg::QHIimg(long len_F, double margin_left, double margin_right)
{
  m_L1 = len_F;
  m_L2 = 0;
  m_extraMarginLeft = margin_left;
  m_extraMarginRight = margin_right;
}

// ----------------------------------------------------------------------------

QHIimg::~QHIimg()
{
  clear();
}

// ----------------------------------------------------------------------------

static void calc_iimg(uint32_t* dst, double* src, long M, long N, double sh, double tol)
{
  //qhlog_debug("calc_iimg (%ldx%ld) sh=%f, tol=%f\n", M, N, sh, tol);
  int cnt = 0;
  double s = 0;
  double s_min = 0;
  double s_max = 0;
  double buf[M];
  for (long j2=0; j2<M; ++j2) {
    s += ( (src[j2] - sh) / tol );
    buf[j2] = s;
    dst[j2] = std::round(s);
    s_min = std::min(s_min, s);
    s_max = std::max(s_max, s);
    ++cnt;
  }
  for (long j1=1; j1<N; ++j1) {
    s = 0;
    double*   p_src = src + j1*M;
    uint32_t*  p_dst = dst + (j1-1)*M;
    for (long j2=0; j2<M; ++j2) {
      s += ( (p_src[j2]-sh) / tol );
      buf[j2] += s;
      p_dst[j2+M] = std::round(buf[j2]);
      s_min = std::min(s_min, buf[j2]);
      s_max = std::max(s_max, buf[j2]);
      ++cnt;
    }
  }

  qhlog_info("calc_iimg: %.3f - %.3f\n", s_min/0x7fffffff, s_max/0x7fffffff);
}

// ----------------------------------------------------------------------------

bool QHIimg::create(std::vector<double> data, int sr, const char* wnd_file, const char* avg_file)
{
  clear();
  DataChunk wnd;
  if (8000 > sr) {
    qhlog_error("ERROR: invalid sample rate\n");
    return false;
  }
  if (! wnd.read(wnd_file)) {
    qhlog_error("ERROR: failed to read wnd file\n");
    return false;
  }
  int wnd_ver = wnd.getLong("qH_STFT_wnd_version", -1);
  if (-1 == wnd_ver) {
    qhlog_error("ERROR: no version\n");
    return false;
  }

  double d = wnd.getDouble("d", 0.0);
  double ta = wnd.getDouble("ta", 0.0);
  if (0.0 >= ta || ta >= d) {
    qhlog_error("ERROR: wrong ta/d %f %f\n", ta, d);
    return false;
  }
  
  long win_len = std::round(d*sr);
  if (win_len/2+1 < m_L1) {
    qhlog_error("ERROR: insufficient sample rate\n");
    return false;
  }

  DataChunk spectr_avg;
  if (! spectr_avg.read(avg_file)) {
    qhlog_error("ERROR: failed to read spectr_avg file\n");
    return false;
  }
  std::vector<double> avg = spectr_avg.getDoubleData();
  if (5 > avg.size()) {
    qhlog_error("ERROR: invalid spectr_avg file\n");
    return false;
  }

  long M = m_L1;
  long margin_left = std::round(sr*(m_extraMarginLeft+ta)) ; //std::round(ta*sr);
  long margin_right = std::round(sr*m_extraMarginRight);
  long win_step = std::round(0.01*sr);
  //qhlog_debug("win_len = %d, step = %d\n", win_len, win_step);

  std::vector<double> w0 = wnd.getDoubleData();
  if (w0.empty()) {
    qhlog_error("ERROR: empty w0\n");
    return false;
  }
  std::vector<double> w(win_len);
  double interp_c = (w0.size()-1.0)/win_len;
  for (long j=0; j<win_len; ++j) {
    double dx = j*interp_c;
    long j1 = floor(dx);
    dx -= j1;
    w[j] = w0[j1] + dx*(w0[j1+1] - w0[j1]);
  }

  long len = data.size();
  const double* inp = data.data();
  double* in;
  fftw_complex* out;
  fftw_plan p;
  in = (double*) fftw_malloc(sizeof(double) * win_len);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (win_len/2+1));
  p = fftw_plan_dft_r2c_1d(win_len, in, out, FFTW_ESTIMATE);

  long N = (1+(len-win_len+margin_left+margin_right)/win_step); 
  std::vector<double> z(M*N);
  double* zp = z.data();
  
  std::vector<double> phi(M*N);
  double* phip = phi.data();
  double p0_sum = 0.0;
  std::vector<double> spectrum(M);
  double* spp = spectrum.data();

  for (long j=-margin_left; j+win_len < len+margin_right; j+=win_step) {
    if (0 > j) {
      long dj = std::min(-j, win_len);
      memset(in, 0, sizeof(double)*dj);
      if (0 < dj) {
        memcpy(in+dj, inp+j+dj, sizeof(double)*(win_len-dj));
      }
    }
    else if(j+win_len > len) {
      long dj = std::min(win_len, j+win_len - len);
      if (0 < dj) {
        memcpy(in, inp+j, sizeof(double)*(win_len-dj));
      }
      memset(in+win_len-dj, 0, sizeof(double)*dj);
    }
    else {
      memcpy(in, inp+j, sizeof(double)*win_len);
    }
    for (long j1=0; j1<win_len; ++j1) {
      in[j1] *= w[j1];
    }
    fftw_execute(p);
    double* c = (double*)out;
    for (long k=0; k<M; ++k) {
      double tmp = c[0]*c[0] + c[1]*c[1];
      tmp = std::max(1e-5, tmp);
      *zp = 10*std::log10(tmp);
      *phip++ = std::atan2(c[1], c[0]) / M_PI;
      //p0_sum += *zp++;
      if (0 <= j && j+win_len < len) {
        spp[k] += *zp++;
      }
      c += 2;
    }
  }
  fftw_destroy_plan(p);
  fftw_free(in);
  fftw_free(out);

  double p1_sum = 0.0;
  double p2_sum = 0.0;
  double delta_phi2 = -0.1;
  double Kf = d/0.01/2;
  std::vector<double> phi2(M*N, 0);
  for (long j1=0; j1<N-1; ++j1) {
    double *p1 = phi.data() + j1*M;
    double *p2 = phi2.data() + j1*M;
    for (long j2=0; j2<M-1; ++j2) {
      p2[j2] = p1[j2+M] - p1[j2] - std::fmod(j2,Kf)+0.1; // 0.1
      p2[j2] -= floor(0.5*p2[j2]+0.5)*2;
      p2_sum += p2[j2];

      p1[j2] = p1[j2+1] - p1[j2] + delta_phi2;
      p1[j2] -= floor(0.5*p1[j2]+0.5)*2;
      p1_sum += p1[j2];
    }
    p1[M-1] = 0.0;
    p2[M-1] = 0.0;
  }
  {
    double *p1 = phi.data() + (N-1)*M;
    double *p2 = phi2.data() + (N-1)*M;
    for (long j2=0; j2<M; ++j2) {
      p2[j2] = 0.0;
      p1[j2] = 0.0;
    }
  }

  int sp_dt = avg.size()/2;
  int sp_n = M / sp_dt - 1;
  qhlog_debug("sp: dt = %d  n = %d\n", sp_dt, sp_n);
  std::vector<double> sp0(sp_n+2);
  double k_sp = 0;
  for (auto x : avg) {
    k_sp += x;
  }
  k_sp = 1.0/(k_sp*N);
  for (int j=0; j<sp_n; ++j) {
    double* sp = spp + j*sp_dt;
    double tmp = 0;
    for (int j1=0; j1<avg.size(); ++j1) {
      tmp += sp[j1]*avg[j1];
    }
    sp0[j+1] = tmp * k_sp;
    //qhlog_debug( " %.0f ", sp0[j+1] );
  }
  //qhlog_debug("\n");
  sp0[0] = sp0[1];
  sp0[sp_n+1] = sp0[sp_n];

  // Interpolate and normalize
  interp_c = 1.0/sp_dt;
  for (long k=0; k<M; ++k) {
    double dx = k*interp_c;
    long j1 = floor(dx);
    dx -= j1;
    double tmp = sp0[j1] + dx*(sp0[j1+1] - sp0[j1]);

    for (long j=0; j<N; ++j) {
      z[k + j*M] -= tmp;
    }
  }

  double z_min = 1e10;
  double z_max = -1e10;
  double z_avg = 0.0;
  for (auto x : z) {
    z_min = std::min(z_min, x);
    z_max = std::max(z_max, x);
    z_avg += x;
  }
  z_avg /= M*(len-win_len)/win_step;
  qhlog_debug("z: %g %g %g\n", z_min, z_avg, z_max);

  p0_sum /= M*N;
  p1_sum /= (M-1)*(N-1);
  p2_sum /= (M-1)*(N-1);

  qhlog_debug("avgs: %g %g %g\n", p0_sum, p1_sum, p2_sum);

# if 0
  {
    double *p1 = phi.data() + (N-1)*M;
    double *p2 = phi2.data() + (N-1)*M;
    for (int j2=0; j2<M-1; ++j2) {
      p2[j2] = 0.0;
      p1[j2] -= p1[j2+1];
    }
  }
  
  for (int j1=0; j1<N-1; ++j1) {
    double *p1 = phi.data() + j1*M;
    double *p2 = phi2.data() + j1*M;
    int j2 = M-1;
    p2[j2] = p1[j2+M] - p1[j2];
    p1[j2] = 0.0;
  }
#endif
  m_iimgs.resize(3);
  m_iimgs[0] = new uint32_t[M*N];
  m_iimgs[1] = new uint32_t[M*N];
  m_iimgs[2] = new uint32_t[M*N];

  m_shifts.resize(3);
  m_shifts[0] = p0_sum;
  m_shifts[1] = p2_sum;
  m_shifts[2] = p1_sum;

  m_tols.resize(3);
  m_tols[0] = 0.02;
  m_tols[1] = 0.0005;
  m_tols[2] = 0.0005;

  calc_iimg(m_iimgs[0], z.data(), M, N, m_shifts[0], m_tols[0]);
  calc_iimg(m_iimgs[1], phi2.data(), M, N, m_shifts[1], m_tols[1]/M_PI);
  calc_iimg(m_iimgs[2], phi.data(), M, N, m_shifts[2], m_tols[2]/M_PI);

  m_L2 = N;
  
  return true;
}

// ----------------------------------------------------------------------------

std::vector<double> QHIimg::transpose(const std::vector<double> src, int M, int N)
{
  if (src.size() != N*M) {
    qhlog_error("transpose error, wrong dimensions %d != %d * %d\n", src.size(), M, N);
    return src;
  }
  std::vector<double> dst(M*N);
  for (int j2=0; j2<M; ++j2) {
    for (int j1=0; j1<N; ++j1) {
      dst[j2*N + j1] = src[j1*M + j2];
    }
  }
  return dst;
}

// ----------------------------------------------------------------------------

bool QHIimg::createStep2(std::vector<double> dataS, std::vector<double> dataL)
{
  clear();
  if (0 > m_L1 || dataS.size() != dataL.size()
      || 0 != (dataS.size() % m_L1) || dataS.size() <= m_L1) {
    qhlog_error("L1=%d sizes=%d,%d  (%d)\n", m_L1, dataS.size(), dataL.size(), dataS.size() % m_L1);
    return false;
  }

  m_L2 = dataS.size() / m_L1;
  long M = m_L1;
  long N = m_L2;
  dataS = transpose(dataS, N, M);
  dataL = transpose(dataL, N, M);
  qhlog_info("size: %d x %d\n", M, N);
  double v_min = dataS[0];
  double v_max = dataS[0];
  double sum_s = 0;
  for (double x : dataS) {
    v_min = std::min(v_min, x);
    v_max = std::max(v_max, x);
    sum_s += x;
  }
  qhlog_debug("values S: %f - %f  => %f %f\n", v_min, v_max, sum_s, sum_s/M/N);

  v_min = dataL[0];
  v_max = dataL[0];
  double sum_l = 0;
  for (double x : dataL) {
    v_min = std::min(v_min, x);
    v_max = std::max(v_max, x);
    sum_l += x;
  }
  qhlog_debug("values L: %f - %f  => %f %f\n", v_min, v_max, sum_l, sum_l/M/N);

  m_iimgs.resize(2);
  m_iimgs[0] = new uint32_t[M*N];
  m_iimgs[1] = new uint32_t[M*N];

  m_shifts.resize(2);
  m_shifts[0] = sum_s/M/N;
  m_shifts[1] = sum_l/M/N;

  m_tols.resize(2);
  m_tols[0] = 0.0001;
  m_tols[1] = 0.0001;

  calc_iimg(m_iimgs[0], dataS.data(), M, N, m_shifts[0], m_tols[0]);
  calc_iimg(m_iimgs[1], dataL.data(), M, N, m_shifts[1], m_tols[1]);

  return true;
}

// ----------------------------------------------------------------------------

double QHIimg::rectSum(int idx, long j1, long j2, long i1, long i2) const
{
#if 1
  if (0 >j1 || 0 >j2 || m_L2<=j1 || m_L2<=j2) {
    qhlog_error("j out of range 0 <= %d %d < %d\n", j1, j2, m_L2);
    return NAN;
  }
  if (0 >i1 || 0 >i2 || m_L1<=i1 || m_L1<=i2) {
    qhlog_error("i out of range 0 <= %d %d < %d\n", i1, i2, m_L1);
    return NAN;
  }
#endif
  uint32_t* p = m_iimgs[idx];
  int32_t sum = p[i2+j2*m_L1]+p[i1+j1*m_L1]-p[i1+j2*m_L1]-p[i2+j1*m_L1];
  double s = sum * 1.0 / ( (i2-i1)*(j2-j1) );
  return s;
}

// ----------------------------------------------------------------------------

void QHIimg::clear()
{
  for (uint32_t* p : m_iimgs) {
    delete [] p;
  }
  m_iimgs.clear();
  m_shifts.clear();
  m_tols.clear();
  m_L2 = 0;
}

// ----------------------------------------------------------------------------

