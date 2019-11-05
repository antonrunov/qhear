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

#include "qhlog.h"

int qhlog_log_level = 10;
FILE* qhlog_err_file  = stdout;
FILE* qhlog_out_file  = stdout;

// ----------------------------------------------------------------------------

int qhlog_error(const char* fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  int r = vfprintf(qhlog_err_file, fmt, ap);
  va_end(ap);
  return r;
}

// ----------------------------------------------------------------------------

int qhlog_message(const char* fmt, ...)
{
  if (1 > qhlog_log_level) {
    return 0;
  }
  va_list ap;
  va_start(ap, fmt);
  int r = vfprintf(qhlog_out_file, fmt, ap);
  va_end(ap);
  return r;
}

// ----------------------------------------------------------------------------

int qhlog_info(const char* fmt, ...)
{
  if (2 > qhlog_log_level) {
    return 0;
  }
  va_list ap;
  va_start(ap, fmt);
  int r = vfprintf(qhlog_out_file, fmt, ap);
  va_end(ap);
  return r;
}

// ----------------------------------------------------------------------------

int qhlog_debug(const char* fmt, ...)
{
  if (3 > qhlog_log_level) {
    return 0;
  }
  va_list ap;
  va_start(ap, fmt);
  int r = vfprintf(qhlog_out_file, fmt, ap);
  va_end(ap);
  return r;
}

// ----------------------------------------------------------------------------

