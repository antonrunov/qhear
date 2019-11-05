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

# pragma once

#include <vector>
#include <stdint.h>

class QHIimg
{
  public:
    QHIimg();
    QHIimg(long len_F, double margin_left, double margin_right);
    ~QHIimg();

  public:
    bool create(std::vector<double> data, int sr, const char* wnd_file, const char* avg_file);
    double rectSum(int idx, long j1, long j2, long i1, long i2) const;
    double unscale(int idx, double v) const {return v*m_tols[idx] + m_shifts[idx];}
    double unscaleDiff(int idx, double v) const {return v*m_tols[idx];}
    bool createStep2(std::vector<double> dataS, std::vector<double> dataL);

    long  numLayers() const {return m_iimgs.size();}
    long  lenT() const {return m_L2;}
    long  lenF() const {return m_L1;}

    static std::vector<double> transpose(const std::vector<double> src, int M, int N);

  protected:
    void clear();

  protected:
    long   m_L1;
    long   m_L2;
    double  m_extraMarginLeft;
    double  m_extraMarginRight;
    std::vector<uint32_t*>   m_iimgs;
    std::vector<double>     m_shifts;
    std::vector<double>     m_tols;
};
