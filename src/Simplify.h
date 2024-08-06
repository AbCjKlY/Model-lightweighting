//
// Created by 123 on 24-7-14.
//

#ifndef MESH_SIMPLIFICATION_SIMPLIFY_H
#define MESH_SIMPLIFICATION_SIMPLIFY_H

#include <vector>
#include <string>
#include "Vec3f.h"
#include "SymmetricMatrix.h"
#include "Obj.h"

enum Attributes {
    NONE,
    NORMAL = 2,
    TEXCOORD = 4,
    COLOR = 8
};

struct ExTriangle {
    double err[4];
    int deleted, dirty;
    Vec3f normal;
    Vec3f uvs[3];
};
struct ExVertex {
    vector<int> ref_triangle;
    SymmetricMatrix q;
    int border;
};

Vec3f barycentric(const Vec3f &p, const Vec3f &a, const Vec3f &b, const Vec3f &c);

Vec3f interpolate(const Vec3f &p, const Vec3f &a, const Vec3f &b, const Vec3f &c, const Vec3f attrs[3]);

class Simplify {
public:
    Obj& obj;
    std::vector<ExTriangle> exTriangles;
    std::vector<ExVertex> exVertices;

    Simplify(Obj &obj) : obj(obj) {}

    double vertex_error(SymmetricMatrix q, double x, double y, double z);
    double calculate_error(int id_v1, int id_v2, Vec3f &p_result);
    bool flipped(Vec3f p, int i0, int i1, ExVertex &v0, ExVertex &v1, std::vector<int> &deleted);
    void update_uvs(int i0, const ExVertex &v, const Vec3f &p, std::vector<int> &deleted);
    void update_triangles(int i0, int i1, std::vector<int> &deleted, int &deleted_triangles);
    void update_mesh();
    void compact_mesh();

    void simplify_mesh(double rate, double agressiveness, bool verbose);

    void load_obj();

    void computeNormals();

    void check_size(string func) const;
};


#endif //MESH_SIMPLIFICATION_SIMPLIFY_H
