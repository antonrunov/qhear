cmake_minimum_required(VERSION 3.0)
project(qhear)

add_executable(qhear qh_proc_main.cpp QHIimg.cpp DataChunk.cpp QHClassifier.cpp QHClassifierBlock.cpp QHFeature.cpp QHEar.cpp qhlog.cpp)
target_link_libraries(qhear -lsndfile -lFLAC -lvorbis -lvorbisenc -logg -lfftw3 -lpthread)
