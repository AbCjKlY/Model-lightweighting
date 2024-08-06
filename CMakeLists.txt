cmake_minimum_required(VERSION 3.21)
project(mesh_simplification)

include_directories(thirdparties/rapidjson/include thirdparties/cpp-httplib)

set(CMAKE_CXX_STANDARD 17)

add_executable(mesh_simplification main.cpp
        src/Vec3f.h
        src/SymmetricMatrix.h
        src/Simplify.cpp
        src/Simplify.h
        src/Obj.cpp
        src/Obj.h
        src/Voxel.cpp
        src/Voxel.h
)

target_link_libraries(mesh_simplification ws2_32)
