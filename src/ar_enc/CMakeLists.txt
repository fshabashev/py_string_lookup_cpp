set (CMAKE_CXX_STANDARD 11)

set(source_files
        AdaptiveArithmeticCompress.cpp
        AdaptiveArithmeticDecompress.cpp
        ArithmeticCompress.cpp
        ArithmeticCoder.cpp
        ArithmeticCoder.hpp
        ArithmeticDecompress.cpp
        BitIoStream.cpp
        BitIoStream.hpp
        FrequencyTable.cpp
        FrequencyTable.hpp
        PpmCompress.cpp
        PpmDecompress.cpp
        PpmModel.cpp
        PpmModel.hpp
        )

file(GLOB SOURCES "./*.cpp")
add_library(ArCompress ${SOURCES})
