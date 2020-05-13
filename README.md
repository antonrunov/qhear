qHEar is a multiple fundamental frequency tracking tool
presented at [MIREX 2019](https://www.music-ir.org/mirex/wiki/2019:Main_Page). You
can find the evaluation results [here](https://www.music-ir.org/mirex/wiki/2019:Multiple_Fundamental_Frequency_Estimation_%26_Tracking_Results_-_MIREX_Dataset).
The utility depends on the configuration data and classifiers
located in the `qh_config` subfolder. By default it expects this folder to be located
near the binary. An alternative config folder path can be specified with `-c` option.

```
Usage: ./qhear [OPTIONS] path/to/file.wav  path/to/output/file.F0

OPTIONS
  -m MODE          - output mode, 1 for Task 1 (frame level evaluation, default), 2 for Task 2 (note tracking)
  -t NUMTHREADS    - number of threads to be used for parallel calculations (use all available cores by default)
  -l MAX_LEN       - maximum processed sample length in seconds (300 by default)
  -c PATH          - specify an alternative location for qh_config directory
  -v               - increase verbosity; use -vv for enabling debug output
  -q               - suppress output
  -h               - print this message
```

##### Command line calling format for Multiple F0 tasks

- Task 1, frame level evaluation
```
    ./qhear -m1 /path/to/input/file/sample.wav /path/to/output/file/f0_task1_result.txt
```

- Task 2, note tracking
```
    ./qhear -m2 /path/to/input/file/sample.wav /path/to/output/file/f0_task2_result.txt
```

##### Additional information.

The program uses multiple threads for parallel computing. By default it creates `sysconf(_SC_NPROCESSORS_ONLN`)
threads. Use `-t` command line option to override this.

The program performs the entire analysis in memory. This is not an algorithm constraint
but just an implementation peculiarity. The length of the audio data to be analyzed is limited
to 300 seconds (5 minutes) by default. The limit can be changed with `-l` option. Typical memory
footprint for 1 minute sample is ~300 Mb.

On a medium laptop i5 CPU with 2 cores it runs about 2 times slower than real-time. More
advanced 8 Core i7 CPUs can provide better than real-time performance.
