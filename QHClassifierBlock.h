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

#pragma once

#include "QHIimg.h"
#include "DataChunk.h"
#include "QHFeature.h"
#include <vector>
#include <atomic>

class QHClassifierBlock
{
  public:
    QHClassifierBlock();

  public:
    bool init(const double* H, int prmVecLen, int len, int bandwidth, double min_f, double max_f);
    bool checkFreq(double f, const QHIimg& iimg) const;
    double classify(const QHIimg& iimg, double t, double f, int out_scale_mode) const;

  protected:
    bool isStep2() const {return 0 == m_qhFeature.bandwidth();}

  protected:
    std::vector<double> m_H;
    QHFeature           m_qhFeature;

    int     m_prmVecLen;
    int     m_order;
    int     m_minDt;
    int     m_maxDt;
    int     m_minDf;
    int     m_maxDf;
    int     m_maxImgIdx;

    double  m_minF;
    double  m_maxF;

    double  m_A;
};
