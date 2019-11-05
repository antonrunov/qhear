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

class QHIimg;

class QHFeature
{
  public:
    QHFeature();

  public:
    void    init(int bandwidth);
    double  calc(const QHIimg& iimg, double t, double f, const double* prm) const;
    int     bandwidth() const { return m_bandwidth;}

    static double calcPlain(const QHIimg& iimg, double t, double f, const double* prm);

  protected:
    int m_bandwidth;
    double m_k0s;
    double m_k1s;
};
