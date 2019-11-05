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

#include "QHClassifierBlock.h"
#include "QHIimg.h"
#include "DataChunk.h"
#include "QHFeature.h"
#include <vector>
#include <atomic>

class QHClassifier
{
  public:
    QHClassifier(std::vector<QHClassifierBlock>& m_blocks);

  public:
    bool scan(const QHIimg& iimg, int t_min, int t_max,
              const std::vector<double> freqs, std::vector<double>& out);
    void setOutScaleMode(int mode) {m_outScaleMode = mode;}

  public:
    static int NTHREADS;

  protected:
    void classifyPoint( std::atomic<int>& p_idx, int N, int NT, const double* p_freqs,
                        const QHIimg& iimg, double* data, int t0, QHClassifierBlock** blocks );
  protected:
    std::vector<QHClassifierBlock> m_blocks;
    int   m_outScaleMode = 0;
};
