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

#include "QHClassifier.h"

#include <vector>
#include <list>
#include <string>
#include <cmath>

class QHEar
{
  public:
    QHEar();
    QHEar(const char* wnd_file, std::vector<QHClassifier> classifiers);

  public:
    bool init(const char* data_dir);
    int scan(std::vector<double> data, int sr, std::vector<double>& out);
    int postproc(std::vector<double> out_step2, int NT);

    void setKF(int kf) {m_KF = kf;}
    void setMinF(double min_f) {m_minF = min_f;}
    void setMaxF(double max_f) {m_maxF = max_f;}
    void setRefF(double ref_f) {m_refF = ref_f;}

  public:
    struct FreqInfo
    {
      double    f;
      double    v;
      uint32_t  flags;
    };
    std::vector<std::list<FreqInfo>> getFreqs() const {return m_freqInfs;}

    struct EvInfo
    {
      double    f;
      int       dur;
    };
    std::multimap<int,EvInfo> getEvents() const {return m_events;}

  protected:
    struct FreqState;

    double  freqHz(double f) const {return m_refF * std::pow(2, (f+m_jf0-1.0)/m_KF);}
    void    releaseFreq(FreqState* st, bool force = false);
    void    logFreq(FreqState* st, int jt, double v);

  protected:
    std::vector<QHClassifier> m_classifiers;
    std::string   m_wndFile;
    std::string   m_avgFile;

    int     m_KF;
    double  m_minF;
    double  m_maxF;
    double  m_refF;
    int     m_jf0;

    std::vector<std::list<FreqInfo>>  m_freqInfs;
    std::multimap<int,EvInfo>              m_events;
    
  // line tracking params
  protected:
    int     m_step2OutScaleMode = 2;

    std::vector<double> m_avgV = {1.0, 0.7, 0.3};
    std::vector<double> m_avgF = {1.0, 0.8, 0.4};

    int     m_jfDelta = 2;
    int     m_jlDelta = 2;
    double  m_scaleF = 1.0;

    int     m_maxLineIdle = 6;
    int     m_eventLenThreshold = 10;

    double  m_freqThreshold = 1.4;
    double  m_freqThresholdK = 3.0;

    double  m_thresholdStart = -0.0;
    double  m_thresholdContinue = -0.0;
    int     m_maxLineWeight = 20;

    int     m_skipFirst = 0;
    int     m_skipLast = 0;
    int     m_durMin = 0;
};
