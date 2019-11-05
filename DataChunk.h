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

#include <string>
#include <vector>
#include <map>

class DataChunk
{
  public:
    DataChunk();
    ~DataChunk();

  public:
    bool read(const char* fname);

    bool        isField(const char* name);
    long        getLong(const char* name, long def);
    double      getDouble(const char* name, double def);
    std::string getString(const char* name, const char* def);

    std::vector<double> getDoubleVector(const char* name, std::vector<double> def);

    std::vector<double> getDoubleData();
    std::vector<char>   getCharData();

  protected:
    std::map<std::string,std::string> m_fields;

    std::string   m_file;
    size_t        m_dataOffset;
    size_t        m_dataSize;
};
