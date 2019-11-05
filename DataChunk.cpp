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

#include "DataChunk.h"
#include <string.h>
#include <limits.h>
#include "qhlog.h"

// ----------------------------------------------------------------------------

DataChunk::DataChunk()
{
  m_dataOffset = 0;
  m_dataSize = 0;
}

// ----------------------------------------------------------------------------

DataChunk::~DataChunk()
{
}

// ----------------------------------------------------------------------------

bool DataChunk::read(const char* fname)
{
  m_fields.clear();
  m_file.clear();
  m_dataOffset = 0;
  m_dataSize = 0;

  FILE* f = fopen(fname, "rb");
  if (NULL == f) {
    return false;
  }

  char buf[4096];
  bool ok = true;
  while (! feof(f)) {
    if (NULL == fgets(buf, sizeof(buf), f)) {
      ok = false;
      break;
    }
    char* s = buf + strspn(buf, " \t");
    if ('\r' == s[0] || '\n' == s[0] ) {
      break;
    }
    if ('#' == s[0] || ';' == s[0]) {
      continue;
    }
    char* s1 = strchr(s, '=');
    if (NULL == s1) {
      continue;
    }
    *s1++ = 0;
    s1 += strspn(s1, " \t");
    // trimming
    for (int j=strlen(s)-1; 0<=j; --j) {
      if (' ' != s[j] && '\t' != s[j]) {
        break;
      }
      s[j] = 0;
    }
    for (int j=strlen(s1)-1; 0<=j; --j) {
      if (' ' != s1[j] && '\t' != s1[j] && '\r' != s1[j] && '\n' != s1[j]) {
        break;
      }
      s1[j] = 0;
    }
    if (0 == s[0]) {
      continue;
    }
    m_fields[s] = s1;
  }

  /*
  qhlog_debug("fields parsing: %s\n", ok ? "ok" : "failed");
  for (auto p : m_fields) {
    qhlog_debug("'%s' => '%s'\n", p.first.c_str(), p.second.c_str());
  }
  */
  if (!ok) {
    m_fields.clear();
    return false;
  }
  m_dataOffset = ftell(f);
  fseek(f, 0, SEEK_END);
  m_dataSize = ftell(f) - m_dataOffset;
  char full_path[PATH_MAX];
  realpath(fname, full_path);
  m_file = full_path;
  fclose(f);
  //qhlog_debug("data: %d bytes, ofs=%d file=%s.\n", m_dataSize, m_dataOffset, m_file.c_str());

  return true;
}

// ----------------------------------------------------------------------------


bool DataChunk::isField(const char* name)
{
  return 0 != m_fields.count(name);
}

// ----------------------------------------------------------------------------

long DataChunk::getLong(const char* name, long def)
{
  if (! isField(name)) {
    return def;
  }
  return std::stol( m_fields.at(name) );
}

// ----------------------------------------------------------------------------

double DataChunk::getDouble(const char* name, double def)
{
  if (! isField(name)) {
    return def;
  }
  return std::stod( m_fields.at(name) );
}

// ----------------------------------------------------------------------------

std::string DataChunk::getString(const char* name, const char* def)
{
  if (! isField(name)) {
    return def;
  }
  return m_fields.at(name);
}

// ----------------------------------------------------------------------------

std::vector<double> DataChunk::getDoubleVector(const char* name, std::vector<double> def)
{
  if (! isField(name)) {
    return def;
  }
  std::string s = m_fields.at(name);
  std::vector<double> ret;
  size_t pos;
  while (! s.empty()) {
    double v = std::stod(s, &pos);
    if (0 == pos) {
      break;
    }
    ret.push_back(v);
    s = s.substr(pos);
  }

  return ret;
}

// ----------------------------------------------------------------------------

std::vector<double> DataChunk::getDoubleData()
{
  std::vector<double> res;
  if (0 < m_dataSize) {
    size_t len = m_dataSize/sizeof(double);
    FILE* f = fopen(m_file.c_str(), "rb");
    if (NULL != f) {
      res.resize(len);
      fseek(f, m_dataOffset, SEEK_SET);
      fread(res.data(), sizeof(double), len, f);
      fclose(f);
    }
  }

  return res;
}

// ----------------------------------------------------------------------------

std::vector<char> DataChunk::getCharData()
{
  std::vector<char> res;
  if (0 < m_dataSize) {
    FILE* f = fopen(m_file.c_str(), "rb");
    if (NULL != f) {
      res.resize(m_dataSize);
      fseek(f, m_dataOffset, SEEK_SET);
      fread(res.data(), 1, m_dataSize, f);
      fclose(f);
    }
  }

  return res;
}

// ----------------------------------------------------------------------------

