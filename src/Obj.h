//
// Created by 123 on 24-7-23.
//

#ifndef MESH_SIMPLIFICATION_OBJ_H
#define MESH_SIMPLIFICATION_OBJ_H

#include <vector>
#include <string>
#include "Vec3f.h"
using namespace std;

struct Triangle {
    int v[3];
    int vt[3];
    int vn[3];
    bool has_vt, has_vn;
    int mtl_index;
};

struct Tex {
    Vec3f uvw;
    bool has_w;
};

class Obj {
public:
    bool has_vt, has_vn;

    vector<Vec3f> v;
    vector<Tex> vt;
    vector<Vec3f> vn;
    vector<Triangle> f;
    vector<string> materials;
    string mtl;

    void clear();
    void load(const string& content);
    void write(const char *filename);
    string write_to_string();
};


#endif //MESH_SIMPLIFICATION_OBJ_H
