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

#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <stdint.h>
#include <sndfile.h>
#include <algorithm>
#include <vector>
#include <cmath>
#include <fftw3.h>
#include <string.h>
#include <getopt.h>
#include <limits.h>

#include "QHIimg.h"
#include "DataChunk.h"
#include "QHEar.h"
#include "qhlog.h"

// ----------------------------------------------------------------------------

static int print_usage()
{
  fprintf(stdout,
      "Usage: qhear [OPTIONS] path/to/file.wav  path/to/output/file.F0\n"
      "\n"
      "OPTIONS\n"
      "  -m MODE          - output mode, 1 for Task 1 (frame level evaluation, default), 2 for Task 2 (note tracking)\n"
      "  -t NUMTHREADS    - number of threads to be used for parallel calculations (use all available cores by default)\n"
      "  -l MAX_LEN       - maximum processed sample length in seconds (300 by default)\n"
      "  -c PATH          - specify an alternative location for qh_config directory\n"
      "  -v               - increase verbosity; use -vv for enabling debug output\n"
      "  -q               - suppress output\n"
      "  -h               - print this message\n"
      );
  return -1;
}

// ----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  int nthreads = 0;
  int out_mode = 1;
  int max_len = 60*5;
  qhlog_err_file = stderr;
  qhlog_log_level = 1;
  std::string qh_config_path;

  while(true) {
    int opt =  getopt(argc, argv, "m:t:l:c:vqh");
    if (-1 == opt) {
      break;
    }
    switch (opt) {
      case 'm':
        out_mode = atoi(optarg);
        break;
      case 't':
        nthreads = atoi(optarg);
        break;
      case 'l':
        max_len = atoi(optarg);
        break;
      case 'c':
        qh_config_path = optarg;
        break;
      case 'v':
        if (0 != qhlog_log_level) ++qhlog_log_level;
        break;
      case 'q':
        qhlog_log_level = 0;
        break;
      default:
        return print_usage();
    }
  }

  if (optind + (0 == out_mode? 0 : 1) >= argc) {
    return print_usage();
  }
  const char* in_file = argv[optind];

  const char* out_bin_file = NULL;
  const char* out_freqs_file = NULL;
  const char* out_events_file = NULL;

  switch (out_mode) {
    case 0:
      out_bin_file = "out.bin";
      out_freqs_file = "out_freqs.txt";
      out_events_file = "out_events.txt";
      break;
    case 1:
      out_freqs_file = argv[optind+1];
      break;
    case 2:
      out_events_file = argv[optind+1];
      break;
    default:
      return print_usage();
  }

  if (1 > nthreads) {
    nthreads = sysconf(_SC_NPROCESSORS_ONLN);
  }
  qhlog_info("using %d threads\n", nthreads);
  QHClassifier::NTHREADS = nthreads;

  SF_INFO info = {0};
  SNDFILE* sf = sf_open(in_file, SFM_READ, &info);
  if (NULL == sf) {
    qhlog_error("failed to open %s (%s)\n", in_file, sf_strerror(sf));
    return -1;
  }
  qhlog_message("%s: %d channels, %d Hz, %d frames\n", in_file, info.channels, info.samplerate, info.frames);
  long len = std::min(info.frames, (sf_count_t)info.samplerate*max_len);

  std::vector<double> inp(len*info.channels);
  long r = sf_readf_double(sf, inp.data(), len);
  if (r != len) {
    qhlog_error("failed to read from file\n");
    return -1;
  }
  if (1 < info.channels) {
    std::vector<double> buf = inp;
    double* p = buf.data();
    inp = std::vector<double>(len);
    for (long j=0; j<len; ++j) {
      double tmp = 0;
      for (int k=0; k<info.channels; ++k) {
        tmp += *p++;
      }
      inp[j] = tmp / info.channels;
    }
  }
  sf_close(sf);
  sf = NULL;

  if (qh_config_path.empty()) {
    char buf[PATH_MAX];
    if (NULL == realpath(argv[0], buf)) {
      qhlog_error("failed to get executable path, using the current dir\n");
      qh_config_path = "./qh_config";
    }
    else {
      char *c = strrchr(buf, '/');
      if (NULL == c) {
        qh_config_path = "./qh_config";
      }
      else {
        *c = 0;
        qh_config_path = buf;
        qh_config_path += "/qh_config";
      }
    }
  }
  qhlog_info("using qh config from %s\n", qh_config_path.c_str());
  QHEar ear;
  if (!ear.init(qh_config_path.c_str())) {
    return -1;
  }
  std::vector<double> out;
  int NT = ear.scan(inp, info.samplerate, out);
  if (0 >= NT) {
    qhlog_error("qh ear scan failed\n");
    return -1;
  }
  ear.postproc(out, NT);

  if (NULL != out_bin_file) {
    FILE* f = fopen(out_bin_file, "w");
    fwrite(out.data(), sizeof(double), out.size(), f);
    fclose(f);
  }

  if (NULL != out_freqs_file) {
    FILE* f = fopen(out_freqs_file, "w");
    const auto out_freqs =  ear.getFreqs();
    for (int j = 0; j < out_freqs.size(); ++j) {
      fprintf(f, "%.2f", j*0.01);
      auto el = out_freqs[j];
      std::vector<double> fr;
      fr.reserve(el.size());
      for (auto fi : el) {
        fr.push_back(fi.f);
      }

      std::sort(fr.begin(), fr.end());
      for (auto fi : fr) {
        fprintf(f, "\t%.2f", fi);
      }
      fprintf(f, "\n");
    }
    fclose(f);
  }

  if (NULL != out_events_file) {
    FILE* f = fopen(out_events_file, "w");
    const auto out_evs =  ear.getEvents();
    for (auto it = out_evs.begin(); it != out_evs.end(); ++it) {
      fprintf(f, "%.2f\t%.2f\t%.2f\n", it->first*0.01, (it->first+it->second.dur)*0.01, it->second.f);
    }
    fclose(f);
  }

  return 0;
}

// ----------------------------------------------------------------------------

