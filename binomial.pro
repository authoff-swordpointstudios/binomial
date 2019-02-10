TEMPLATE = app
CONFIG += console c++14
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        main.cpp

QMAKE_CXXFLAGS += \
            -std=c++17 \
#            -std=gnu++17 \
#            -std=c++2a \
#            -std=gnu++2a \
            -fconcepts \
            -fopenmp \
#            -pedantic-errors \
#            -pedantic \
#            -Weffc++ \
            -Wall \
            -Wno-unused-parameter

QMAKE_CXXFLAGS += -DBOOST_LOG_DYN_LINK  # on compile line not link line.

# -DBOOST_SYSTEM_NO_DEPRECATED

LIBS += -pthread  # to avoid some link error on Fedora26

HEADERS += \
    main.h \
    ../grostig_tools/cpp_headers.h \
    ../grostig_tools/boost_headers.h \
    ../grostig_tools/gnuplot-iostream.h \

LIBS += \
        -ltbb \


# boost on fedora26 parial list only derived from Unbuntu list
LIBS += \
#        -lboost_atomic \
#        -lboost_chrono \
#        -lboost_context \
#        -lboost_coroutine \
#        -lboost_date_time \
#        -lboost_fiber \
#        -lboost_filesystem \
#        -lboost_graph_parallel \
#        -lboost_graph \
        -lboost_iostreams \
#        -lboost_locale \
#        -lboost_log_setup \
        -lboost_log \
#        -lboost_math_c99f \
#        -lboost_math_c99l \
#        -lboost_math_c99 \
#        -lboost_math_tr1f \
#        -lboost_math_tr1l \
#        -lboost_math_tr1 \
#        -lboost_mpi_python-py27 \
#        -lboost_mpi_python-py35 \
#        -lboost_mpi_python \
#        -lboost_mpi \
#        -lboost_prg_exec_monitor \
        -lboost_program_options \
#        -lboost_python-py27 \
#        -lboost_python-py35 \
#        -lboost_python \
#        -lboost_random \
#        -lboost_regex \
#        -lboost_serialization \
#        -lboost_signals \
        -lboost_system \
#        -lboost_thread \
#        -lboost_timer \
#        -lboost_type_erasure \
#        -lboost_unit_test_framework \
#        -lboost_wave \
#        -lboost_wserialization

# boost on ubuntu 17.04
#LIBS += \
#        -lboost_atomic \
#        -lboost_chrono \
#        -lboost_context \
#        -lboost_coroutine \
#        -lboost_date_time \
#        -lboost_fiber \
#        -lboost_filesystem \
#        -lboost_graph_parallel \
#        -lboost_graph \
#        -lboost_iostreams \
#        -lboost_locale \
#        -lboost_log_setup \
#        -lboost_log \
#        -lboost_math_c99f \
#        -lboost_math_c99l \
#        -lboost_math_c99 \
#        -lboost_math_tr1f \
#        -lboost_math_tr1l \
#        -lboost_math_tr1 \
#        -lboost_mpi_python-py27 \
#        -lboost_mpi_python-py35 \
#        -lboost_mpi_python \
#        -lboost_mpi \
#        -lboost_prg_exec_monitor \
#        -lboost_program_options \
#        -lboost_python-py27 \
#        -lboost_python-py35 \
#        -lboost_python \
#        -lboost_random \
#        -lboost_regex \
#        -lboost_serialization \
#        -lboost_signals \
#        -lboost_system \
#        -lboost_thread \
#        -lboost_timer \
#        -lboost_type_erasure \
#        -lboost_unit_test_framework \
#        -lboost_wave \
#        -lboost_wserialization
