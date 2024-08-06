//
// Created by 123 on 24-7-30.
//

#ifndef MESH_SIMPLIFICATION_VOXEL_H
#define MESH_SIMPLIFICATION_VOXEL_H


#include "Obj.h"

class Voxel {
public:
    Obj& obj;
    vector<int> triangle_index;

    Voxel(Obj& obj) : obj(obj) {}
};


#endif //MESH_SIMPLIFICATION_VOXEL_H
