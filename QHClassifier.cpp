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
#include <string.h>
#include <cmath>
#include <thread>
#include <mutex>
#include "qhlog.h"

int QHClassifier::NTHREADS = 1;

// ----------------------------------------------------------------------------

QHClassifier::QHClassifier(std::vector<QHClassifierBlock>& blocks)
{
  m_blocks = blocks;
}

// ----------------------------------------------------------------------------

bool QHClassifier::scan(const QHIimg& iimg, int t_min, int t_max,
          const std::vector<double> freqs, std::vector<double>& out)
{
  std::vector<QHClassifierBlock*> block_map(freqs.size(), NULL);
  for (int j=0; j<freqs.size(); ++j) {
    for(int k=0; k<m_blocks.size(); ++k) {
      if (m_blocks[k].checkFreq(freqs[j], iimg)) {
        block_map[j] = &m_blocks[k];
        break;
      }
    }
    /*
    if (NULL == block_map[j]) {
      qhlog_debug("no classifiers for %f\n", freqs[j]);
    }
    */
  }

  int NT = t_max-t_min+1;
  int N = NT*freqs.size();
  double min_value = 0.0;
  if (1 == m_outScaleMode) {
    min_value = -6.0;
  }
  else if (2 == m_outScaleMode) {
    min_value = -3.0;
  }

  if (out.size() != N) {
    out = std::vector<double>(N, min_value);
  }
  double* data = out.data();

  std::atomic<int> p_idx;
  p_idx = 0;
  std::thread th[NTHREADS];
  const double* p_freqs = freqs.data();
  QHClassifierBlock** blocks = block_map.data();
  for (int j=0; j<NTHREADS; ++j) {
    th[j] = std::thread(  &QHClassifier::classifyPoint, this,
                          std::ref(p_idx), N, NT, p_freqs, std::ref(iimg), data, t_min, blocks);
  }

  for (int j=0; j<NTHREADS; ++j) {
    th[j].join();
  }

  return true;
}

// ----------------------------------------------------------------------------

void QHClassifier::classifyPoint( std::atomic<int>& p_idx, int N, int NT, const double* p_freqs,
                                  const QHIimg& iimg, double* data, int t0, QHClassifierBlock** blocks)
{
  while (true) {
    int j = p_idx++;
    if (j>=N) {
      break;
    }
    int jf = j / NT;
    const QHClassifierBlock* block = blocks[jf];
    if (NULL == block) {
      continue;
    }
    int t = t0 + j-jf*NT;
    double f = p_freqs[jf];
    data[j] = block->classify(iimg, t, f, m_outScaleMode);
  }
}

// ----------------------------------------------------------------------------

