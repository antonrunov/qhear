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

#include "QHEar.h"
#include "QHIimg.h"
#include <cmath>
#include <glob.h>
#include "DataChunk.h"
#include "qhlog.h"

// ---------------------------------------------------------------------------

QHEar::QHEar()
{
  m_KF = 12*4;
  m_minF = 12.792;
  m_maxF = 5124.0;
  m_refF = 440.0;
  m_jf0 = 0;
}

// ---------------------------------------------------------------------------

QHEar::QHEar(const char* wnd_file, std::vector<QHClassifier> classifiers)
{
  m_wndFile = wnd_file;
  m_avgFile = "qh_sp_avg.txt";
  m_classifiers = classifiers;
  m_KF = 12*4;
  m_minF = 12.792;
  m_maxF = 5124.0;
  m_refF = 440.0;
  m_jf0 = 0;
}

// ---------------------------------------------------------------------------

static bool load_classifier_blocks(std::vector<QHClassifierBlock>& cb, const glob_t* gbuf)
{
  cb.resize(gbuf->gl_pathc);
  for (int j=0; j<cb.size(); ++j) {
    DataChunk cls;
    if (!cls.read(gbuf->gl_pathv[j])) {
      qhlog_error("failed to read data block from %s\n", gbuf->gl_pathv[j]);
      return false;
    }
    int bw = cls.getLong("bw", 0);
    double min_f = cls.getDouble("min_f", 0.0);
    double max_f = cls.getDouble("max_f", 0.0);
    int L = cls.getLong("L", 0);
    int PN = cls.getLong("PN", 0);
    std::vector<double> H = cls.getDoubleData();
    if (H.empty() || H.size() != PN*L) {
      qhlog_error("incorrect data in %s\n", gbuf->gl_pathv[j]);
      return false;
    }

    if (! cb[j].init(H.data(), PN, L, bw, min_f, max_f)) {
      qhlog_error("failed to init qh block from %s\n", gbuf->gl_pathv[j]);
      return false;
    }
  }

  return true;
}

// ---------------------------------------------------------------------------

bool QHEar::init(const char* data_dir)
{
  std::string base(data_dir);
  glob_t globbuf;

  DataChunk config;
  if (!config.read((base+"/qh.ini").c_str())) {
    qhlog_error("QHEar::init failed, no ini file\n");
    return false;
  }
  m_KF = 12 * config.getDouble("KF12", 4);
  m_minF = config.getDouble("min_f", m_minF);
  m_maxF = config.getDouble("max_f", m_maxF);
  m_refF = config.getDouble("ref_f", m_refF);

  // line tracking params
  m_step2OutScaleMode = config.getLong("step2_OutScaleMode", m_step2OutScaleMode);

  m_avgV = config.getDoubleVector("avg_v", m_avgV);
  m_avgF = config.getDoubleVector("avg_f", m_avgF);

  m_jfDelta = config.getLong("f_delta", m_jfDelta);
  m_jlDelta = config.getLong("l_delta", m_jlDelta);
  m_scaleF = config.getDouble("scale_f", m_scaleF);

  m_maxLineIdle = config.getLong("max_line_idle", m_maxLineIdle);
  m_eventLenThreshold = config.getLong("event_len_threshold", m_eventLenThreshold);

  m_freqThreshold = config.getDouble("freq_threshold", m_freqThreshold);
  m_freqThresholdK = config.getDouble("freq_threshold_K", m_freqThresholdK);

  m_thresholdStart = config.getDouble("threshold_start", m_thresholdStart);
  m_thresholdContinue = config.getDouble("threshold_continue", m_thresholdContinue);
  m_maxLineWeight = config.getLong("max_line_weight", m_maxLineWeight);
  double threshold_hist = config.getDouble("threshold_hist", NAN);
  if (!std::isnan(threshold_hist)) {
    m_thresholdContinue = m_thresholdStart - threshold_hist;
  }

  m_skipFirst = config.getLong("skip_first", m_skipFirst);
  m_skipLast = config.getLong("skip_last", m_skipLast);
  m_durMin = config.getLong("dur_min", m_durMin);

  m_avgV.resize(m_jfDelta+1, 0.0);
  m_avgF.resize(m_jfDelta+1, 0.0);

  /*
  qhlog_debug("avg_v:");
  for (auto x : m_avgV) {
    qhlog_debug(" %f", x);
  }
  qhlog_debug("\n");
  */

  if (0 != glob( (base+"/wnd1.chnk").c_str(), 0, NULL, &globbuf)) {
    qhlog_error("QHEar::init failed, no wnd file\n");
    globfree(&globbuf);
    return false;
  }
  m_wndFile = globbuf.gl_pathv[0];
  globfree(&globbuf);

  if (0 != glob( (base+"/avg1.chnk").c_str(), 0, NULL, &globbuf)) {
    qhlog_error("QHEar::init failed, no avg file\n");
    globfree(&globbuf);
    return false;
  }
  m_avgFile = globbuf.gl_pathv[0];
  globfree(&globbuf);

  m_classifiers.reserve(3);
  if (0 != glob( (base+"/cls_f_*.chnk").c_str(), 0, NULL, &globbuf)) {
    qhlog_error("QHEar::init failed, no cls_f\n");
    globfree(&globbuf);
    return false;
  }
  else {
    std::vector<QHClassifierBlock> cb;
    if (!load_classifier_blocks(cb, &globbuf)) {
      qhlog_error("QHEar::init failed, invalid cls_f");
      globfree(&globbuf);
      return false;
    }
    m_classifiers.push_back( QHClassifier(cb) );
  }
  globfree(&globbuf);

  if (0 != glob( (base+"/cls_l_*.chnk").c_str(), 0, NULL, &globbuf)) {
    qhlog_error("QHEar::init failed, no cls_l\n");
    globfree(&globbuf);
    return false;
  }
  else {
    std::vector<QHClassifierBlock> cb;
    if (!load_classifier_blocks(cb, &globbuf)) {
      qhlog_error("QHEar::init failed, invalid cls_f");
      globfree(&globbuf);
      return false;
    }
    m_classifiers.push_back( QHClassifier(cb) );
  }
  globfree(&globbuf);

  if (0 != glob( (base+"/cls_step2_*.chnk").c_str(), 0, NULL, &globbuf)) {
    qhlog_error("QHEar::init failed, no cls_step2\n");
    globfree(&globbuf);
    return false;
  }
  else {
    std::vector<QHClassifierBlock> cb;
    if (!load_classifier_blocks(cb, &globbuf)) {
      qhlog_error("QHEar::init failed, invalid cls_f");
      globfree(&globbuf);
      return false;
    }
    m_classifiers.push_back( QHClassifier(cb) );
  }
  globfree(&globbuf);

  return true;
}

// ---------------------------------------------------------------------------

struct QHEar::FreqState
{
  bool    active() const {return 0 < idleCount;}

  double  f = NAN;
  int     idleCount = 0;
  int     start = 0;
  int     last = 0;
  int     len = 0;
  int     touch = 0;
  double  v_max = 0.0;
};

// ---------------------------------------------------------------------------

void QHEar::releaseFreq(FreqState* st, bool force) {
  if (st->active()) {
    if (force) {
      st->idleCount = 0;
    }
    else {
      --st->idleCount;
    }
    if (0 == st->idleCount && m_eventLenThreshold <= st->len) {
      int dur = st->last - st->start;
      if (dur < m_durMin) {
        FreqInfo fi = {freqHz(st->f), st->v_max, 0};
        for(int k = st->last+1; k <= st->start+m_durMin; ++k) {
          m_freqInfs[k].push_back(fi);
        }
        dur = m_durMin;
      }
      EvInfo ei = {freqHz(st->f), dur};
      m_events.insert(std::make_pair(st->start, ei));
    }
  }
}

// ---------------------------------------------------------------------------

void QHEar::logFreq(FreqState* st, int jt, double v) {
  st->touch = jt;
  if (m_eventLenThreshold > st->len) {
    st->v_max = std::max(v, st->v_max);
    return;
  }

  FreqInfo fi = {freqHz(st->f), st->v_max, 0};
  int next = std::max(st->start-1, jt - m_skipLast);
  for(int k = st->last+1; k <= next; ++k) {
    m_freqInfs[k].push_back(fi);
  }
  st->v_max = v;
  st->last = next;
}

// ---------------------------------------------------------------------------

int QHEar::scan(std::vector<double> data, int sr, std::vector<double>& out)
{
  m_freqInfs.clear();
  m_events.clear();

  QHIimg iimg(1600, 0.0, 2.0);
  if (!iimg.create(data, sr, m_wndFile.c_str(), m_avgFile.c_str())) {
    qhlog_error("QHEar::scan: failed to create first iimg\n");
    return -1;
  }
  int NT = iimg.lenT();
  std::vector<double> freqs;
  std::vector<double> freqs_s2;

  m_jf0 = std::ceil(std::log2(m_minF/m_refF)*m_KF);
  int jf1 = std::floor(std::log2(m_maxF/m_refF)*m_KF);
  freqs.resize(jf1-m_jf0+1);
  freqs_s2.resize(jf1-m_jf0+1);
  for (int j=m_jf0; j<=jf1; ++j) {
    freqs[j-m_jf0] = m_refF * std::pow(2, 1.0*j/m_KF);
    freqs_s2[j-m_jf0] = j-m_jf0;
  }
  qhlog_debug("freqs: %f %f %f - %f (%d)\n", freqs[0], freqs[1], freqs[2], freqs[jf1-m_jf0], jf1-m_jf0+1);

  qhlog_message("scanning: fronts\n");
  std::vector<double> out_s;
  if (!m_classifiers[0].scan(iimg, 0, NT-1, freqs, out_s)) {
    qhlog_error("QHEar::scan: failed to scan (0)\n");
    return -1;
  }

  qhlog_message("scanning: lines\n");
  std::vector<double> out_l;
  if (!m_classifiers[1].scan(iimg, 0, NT-1, freqs, out_l)) {
    qhlog_error("QHEar::scan: failed to scan (1)\n");
    return -1;
  }

  for (double& x : out_s) {
    x = (x-0.5)*2;
  }

  for (double& x : out_l) {
    x = (x-0.5)*2;
  }

  QHIimg iimg_step2(freqs.size(), 0.0, 0.0);
  if (!iimg_step2.createStep2(out_s, out_l)) {
    qhlog_error("QHEar::scan: failed to create step 2 iimg\n");
    return -1;
  }

  qhlog_message("scanning: step2\n");
  std::vector<double> out_step2;
  m_classifiers[2].setOutScaleMode(m_step2OutScaleMode);
  if (!m_classifiers[2].scan(iimg_step2, 0, NT-1, freqs_s2, out_step2)) {
    qhlog_error("QHEar::scan: failed to scan (1)\n");
    return -1;
  }
  out = out_step2;
  return NT;
}

// ---------------------------------------------------------------------------

int QHEar::postproc(std::vector<double> out_step2, int NT)
{
  m_freqInfs.clear();
  m_events.clear();
  
  std::vector<double> freqs;

  m_jf0 = std::ceil(std::log2(m_minF/m_refF)*m_KF);
  int jf1 = std::floor(std::log2(m_maxF/m_refF)*m_KF);
  freqs.resize(jf1-m_jf0+1);
  for (int j=m_jf0; j<=jf1; ++j) {
    freqs[j-m_jf0] = m_refF * std::pow(2, 1.0*j/m_KF);
  }

  std::vector<FreqState> state(freqs.size());
  m_freqInfs.resize(NT);

  qhlog_message("postprocessing\n");
  for (int jt=0; jt<NT; ++jt) {
    for (int jf=m_jfDelta; jf<freqs.size()-m_jfDelta; ++jf) {
      int n = jt + jf*NT;
      releaseFreq(&state[jf-m_jfDelta]);
      bool extr = true;
      for (int dn = NT; dn <= m_jfDelta*NT; dn += NT) {
        if (out_step2[n] < out_step2[n-dn] || out_step2[n] <= out_step2[n+dn]) {
          extr = false;
          break;
        }
      }
      if (extr) {
        double wv = m_avgV[0];
        double v = out_step2[n] * wv;
        double wf = exp(m_scaleF * out_step2[n]) * m_avgF[0];
        double f = jf*wf;
        for (int dn = 1; dn <= m_jfDelta; ++dn) {
          v += m_avgV[dn] * (out_step2[n-NT*dn] + out_step2[n+NT*dn]);
          wv += 2*m_avgV[dn];
          double w_tmp = exp(m_scaleF * out_step2[n-NT*dn]) * m_avgF[dn];
          f += (jf-dn) * w_tmp;
          wf += w_tmp;
          w_tmp = exp(m_scaleF * out_step2[n+NT*dn]) * m_avgF[dn];
          f += (jf+dn) * w_tmp;
          wf += w_tmp;
        }
        v /= wv;
        f /= wf;

        int last_idx = -1;
        double df = NAN;
        for (int k=-m_jlDelta; k<=m_jlDelta; ++k) {
          if (state[jf+k].active() && state[jf+k].touch < jt) {
            double df_tmp = f - state[jf+k].f;
            if ( -1 == last_idx || fabs(df_tmp) < fabs(df) ) {
              if (fabs(df_tmp) < m_freqThreshold*(1 + m_freqThresholdK/state[jf+k].len)) {
                last_idx = jf+k;
                df = df_tmp;
              }
            }
          }
        }

        FreqState* st = NULL;
        if (-1 == last_idx && m_thresholdStart < v) {
          // new freq
          int new_idx = std::round(f);
          st = &state[new_idx];
          releaseFreq(st, true);
          st->idleCount = m_maxLineIdle;
          st->start = jt + m_skipFirst;
          st->last = st->start - 1;
          st->len = 1;
          st->f = f;
          st->v_max = v;
        }
        else if (-1 < last_idx && m_thresholdContinue < v) {
          // continuation
          st = &state[last_idx];
          if (m_maxLineWeight > st->len) {
            ++st->len;
          }
          st->idleCount = m_maxLineIdle;
          st->f += df / st->len;
          //*
          if (std::round(st->f) != last_idx) {
            // move
            last_idx = std::round(st->f);
            releaseFreq(&state[last_idx], true);
            state[last_idx] = *st;
            st->idleCount = 0;
            st = &state[last_idx];
          }
          // */
        }

        if (NULL != st) {
          logFreq(st, jt, v);
        }

      }
    }
  }

  return NT;
}

// ---------------------------------------------------------------------------

