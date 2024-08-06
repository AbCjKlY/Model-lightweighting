//
// Created by 123 on 24-7-14.
//

#include <cstring>
#include <cfloat>
#include <map>
#include <sstream>
#include <iostream>
#include "Simplify.h"

#define loopi(start_l,end_l) for ( int i=start_l;i<end_l;++i )
#define loopi(start_l,end_l) for ( int i=start_l;i<end_l;++i )
#define loopj(start_l,end_l) for ( int j=start_l;j<end_l;++j )
#define loopk(start_l,end_l) for ( int k=start_l;k<end_l;++k )

Vec3f barycentric(const Vec3f &p, const Vec3f &a, const Vec3f &b, const Vec3f &c){
    Vec3f v0 = b - a;
    Vec3f v1 = c - a;
    Vec3f v2 = p - a;
    double d00 = v0.dot(v0);
    double d01 = v0.dot(v1);
    double d11 = v1.dot(v1);
    double d20 = v2.dot(v0);
    double d21 = v2.dot(v1);
    double denominator = d00*d11-d01*d01;
    double v = (d11 * d20 - d01 * d21) / denominator;
    double w = (d00 * d21 - d01 * d20) / denominator;
    double u = 1.0 - v - w;
    return {u,v,w};
}//计算点p的重心坐标

Vec3f interpolate(const Vec3f &p, const Vec3f &a, const Vec3f &b, const Vec3f &c, const Vec3f attrs[3])
{
    Vec3f bary = barycentric(p, a, b, c);
    Vec3f out = Vec3f(0, 0, 0);
    out = out + attrs[0] * bary.x;
    out = out + attrs[1] * bary.y;
    out = out + attrs[2] * bary.z;
    return out;
}

// Check if a triangle flips when this edge is removed
bool Simplify::flipped(Vec3f p, int i0, int i1, ExVertex &v0, ExVertex &v1, std::vector<int> &deleted) {
    loopk(0, v0.ref_triangle.size())
    {
        int t_index = v0.ref_triangle[k];
        ExTriangle &et = exTriangles[t_index];
        Triangle &tri = obj.f[t_index];
        if(et.deleted) continue;

        int v_index_in_triangle = 0;
        for (; v_index_in_triangle < 3; ++v_index_in_triangle)
            if (tri.v[v_index_in_triangle] == i0) break;
        int id1 = tri.v[(v_index_in_triangle + 1) % 3];
        int id2 = tri.v[(v_index_in_triangle + 2) % 3];

        if(id1 == i1 || id2 == i1) // delete ?
        {
            deleted[k] = 1;
            continue;
        }
        Vec3f d1 = obj.v[id1] - p;  d1.normalize();
        Vec3f d2 = obj.v[id2] - p;  d2.normalize();
        if (fabs(d1.dot(d2)) > 0.999) return true;

        Vec3f n(0, 0, 0);
        n.cross(d1, d2);
        n.normalize();
        deleted[k] = 0;
        if(n.dot(et.normal) < 0.2) return true;
    }
    return false;
}

// update_uvs
void Simplify::update_uvs(int i0, const ExVertex &ev, const Vec3f &p, vector<int> &deleted)
{
    loopk (0, ev.ref_triangle.size())
    {
        int t_index = ev.ref_triangle[k];
        ExTriangle &et = exTriangles[t_index];
        Triangle &tri = obj.f[t_index];
        if (et.deleted) continue;
        if (deleted[k]) continue;
        int v_index_in_triangle = 0;
        for (; v_index_in_triangle < 3; ++v_index_in_triangle)
        {
            if (tri.v[v_index_in_triangle] == i0)
                break;
        }
        Vec3f p1 = obj.v[tri.v[0]];
        Vec3f p2 = obj.v[tri.v[1]];
        Vec3f p3 = obj.v[tri.v[2]];
        et.uvs[v_index_in_triangle] = interpolate(p, p1, p2, p3, et.uvs);
    }
}

// Update triangle connections and edge error after a edge is collapsed
void Simplify::update_triangles(int i0, int i1, vector<int> &deleted, int &deleted_triangles)
{
    ExVertex &ev0 = exVertices[i0];
    ExVertex &ev = exVertices[i1];
    Vec3f p(0, 0, 0);
    loopk(0, ev.ref_triangle.size())
    {
        int t_index = ev.ref_triangle[k];
        Triangle &tri = obj.f[t_index];
        ExTriangle &et = exTriangles[t_index];
        if (et.deleted) continue;
        if (deleted[k])
        {
            et.deleted = 1;
            deleted_triangles ++;
            continue;
        }
        int v_index_in_triangle = 0;
        for (; v_index_in_triangle < 3; ++v_index_in_triangle)
        {
            if (tri.v[v_index_in_triangle] == i1)
                break;
        }
        tri.v[v_index_in_triangle] = i0;
        if (i0 != i1)
        {
            ev0.ref_triangle.push_back(t_index);
        }
        et.dirty = 1;
        et.err[0] = calculate_error(tri.v[0], tri.v[1], p);
        et.err[1] = calculate_error(tri.v[1], tri.v[2], p);
        et.err[2] = calculate_error(tri.v[2], tri.v[0], p);
        et.err[3] = fmin(et.err[0], fmin(et.err[1], et.err[2]));
    }
}

// compact triangles, compute edge error and build reference list
void Simplify::update_mesh()
{
    int dst = 0;
    loopi (0, exTriangles.size())
    {
        ExTriangle &et = exTriangles[i];
        if (et.deleted == false)
        {
            Triangle &tri = obj.f[i];
            exTriangles[dst] = et;
            obj.f[dst] = tri;
            dst ++;
        }
    }
    exTriangles.resize(dst);
    obj.f.resize(dst);

    // Init Reference TriangleID list
    for (auto & ev : exVertices)
        ev.ref_triangle.clear();
    loopi (0, exTriangles.size())
    {
        loopj (0, 3)
        {
            int v_index = obj.f[i].v[j];
            exVertices[v_index].ref_triangle.push_back(i);
        }
    }
    check_size("Update_mesh");
}

// Finally compact mesh before exiting

void Simplify::compact_mesh()
{
    int dst = 0;

    update_mesh();

    loopi (0, exVertices.size())
    {
        ExVertex &ev = exVertices[i];
        if (!ev.ref_triangle.empty())
        {
            Vec3f &ver = obj.v[i];
            exVertices[dst] = ev;
            obj.v[dst] = ver;
            for (auto & ref_t_index : ev.ref_triangle)
            {
                for (auto & v_index : obj.f[ref_t_index].v)
                    if (v_index == i) v_index = dst;
            }
            dst ++;
        }
    }
    exVertices.resize(dst);
    obj.v.resize(dst);

    check_size("Compact_mesh");
}

// Error between vertex and Quadric
double Simplify::vertex_error(SymmetricMatrix q, double x, double y, double z){
    return q[0]*x*x + 2*q[1]*x*y + 2*q[2]*x*z + 2*q[3]*x + q[4]*y*y
            + 2*q[5]*y*z + 2*q[6]*y + q[7]*z*z + 2*q[8]*z + q[9];
}

// Error for one edge
double Simplify::calculate_error(int id_v1, int id_v2, Vec3f &p_result)
{
    // compute interpolated vertex
    SymmetricMatrix q = exVertices[id_v1].q + exVertices[id_v2].q;
    bool border = exVertices[id_v1].border & exVertices[id_v2].border;
    double error = 0;
    double det = q.det(0, 1, 2, 1, 4, 5, 2, 5, 7);
    if ( det != 0 && !border )
    {
        // q_delta is invertible
        p_result.x = -1/det*(q.det(1, 2, 3, 4, 5, 6, 5, 7 , 8));
        // vx = A41/det(q_delta)
        p_result.y =  1/det*(q.det(0, 2, 3, 1, 5, 6, 2, 7 , 8));
        // vy = A42/det(q_delta)
        p_result.z = -1/det*(q.det(0, 1, 3, 1, 4, 6, 2, 5,  8));
        // vz = A43/det(q_delta)
        error = vertex_error(q, p_result.x, p_result.y, p_result.z);
    }
    else
    {
        // det = 0 -> try to find best result
        Vec3f p1 = obj.v[id_v1];
        Vec3f p2 = obj.v[id_v2];
        Vec3f p3 = (p1 + p2) / 2;
        double error1 = vertex_error(q, p1.x, p1.y, p1.z);
        double error2 = vertex_error(q, p2.x, p2.y, p2.z);
        double error3 = vertex_error(q, p3.x, p3.y, p3.z);
        error = fmin(error1, fmin(error2, error3));
        if (error1 == error) p_result = p1;
        else if (error2 == error) p_result = p2;
        else if (error3 == error) p_result = p3;
    }
    return error;
}

void Simplify::simplify_mesh(double rate, double agressiveness=7, bool verbose=false)
{
    int triangle_count = (int) obj.f.size();
    int target_count = (int) (triangle_count * rate);
    if (target_count < 500) target_count = min(500, triangle_count);
    int deleted_triangles = 0;

    vector<int> deleted0, deleted1;

    for (int iteration = 0; iteration < 120; iteration ++)
    {
        if (triangle_count - deleted_triangles <= target_count) break;

        //Update_mesh
        if (iteration % 5 == 0 && iteration > 0)
        {
            update_mesh();
        }

        for (auto & et : exTriangles) et.dirty = 0;

        // All triangles with edges below the threshold will be removed
        //
        // The following numbers works well for most models.
        // If it does not, try to adjust the 3 parameters
        double threshold = 0.000000001 * pow((iteration + 3), agressiveness);

        if (verbose && iteration % 5 == 0)
        {
            printf("iteration %d - triangles %d threshold %g\n", iteration, triangle_count-deleted_triangles, threshold);
        }

        // remove vertices & mark deleted triangles
        loopi(0, exTriangles.size())
        {
            ExTriangle &et = exTriangles[i];
            Triangle &tri = obj.f[i];
            if(et.err[3] > threshold) continue;
            if(et.deleted || et.dirty) continue;

            loopj(0, 3)
            {
                if (et.err[j] < threshold)
                {
                    int v_index_0 = obj.f[i].v[j];
                    int v_index_1 = obj.f[i].v[(j + 1) % 3];
                    ExVertex &ev0 = exVertices[v_index_0];
                    ExVertex &ev1 = exVertices[v_index_1];
                    // Border check
                    if (ev0.border != ev1.border) continue;

                    // Compute vertex to collapse to
                    Vec3f p;
                    calculate_error(v_index_0, v_index_1, p);

                    deleted0.resize(ev0.ref_triangle.size()); // normals temporarily
                    deleted1.resize(ev1.ref_triangle.size()); // normals temporarily

                    // don't remove if flipped
                    if (flipped(p, v_index_0, v_index_1, ev0, ev1, deleted0)) continue;
                    if (flipped(p, v_index_1, v_index_0, ev1, ev0, deleted1)) continue;

                    if (tri.has_vt)
                    {
                        update_uvs(v_index_0, ev0, p, deleted0);
                        update_uvs(v_index_0, ev1, p, deleted1);
                    }

                    // not flipped, so remove edge
                    obj.v[v_index_0] = p;
                    ev0.q += ev1.q;

                    update_triangles(v_index_0, v_index_0, deleted0, deleted_triangles);
                    update_triangles(v_index_0, v_index_1, deleted1, deleted_triangles);

                    break;
                }
            }
            // done?
            if (triangle_count-deleted_triangles <= target_count) break;
        }
    }
    compact_mesh();
    if (obj.has_vn) computeNormals();
}

void Simplify::check_size(string func) const {
//    cout << func << "   Triangle size: " << obj.f.size() << "   ExTriangle size: " << exTriangles.size() << endl;
//    cout << func << "   Vertices size: " << obj.v.size() << "   ExVertices size: " << exVertices.size() << endl;
}

void Simplify::computeNormals() {
    obj.vn.resize(exVertices.size(), {0, 0, 0});

    for (size_t i = 0; i < exTriangles.size(); ++i) {
        const Triangle &tri = obj.f[i];
        const ExTriangle &et = exTriangles[i];

        obj.vn[tri.v[0]] += et.normal;
        obj.vn[tri.v[1]] += et.normal;
        obj.vn[tri.v[2]] += et.normal;
    }

    for (auto& normal : obj.vn) {
        normal.normalize();
    }
}

void Simplify::load_obj() {
    exVertices.clear();
    exTriangles.clear();

    //initialize
    std::vector<int> v_count, v_ids;

    loopi (0, obj.v.size())
    {
        ExVertex ev;
        ev.border = 0;
        ev.q = SymmetricMatrix(0.0);
        exVertices.push_back(ev);
    }
    loopi (0, obj.f.size())
    {
        ExTriangle et;
        et.deleted = false;
        et.dirty = false;
        et.normal = {0.0, 0.0, 0.0};
        exTriangles.push_back(et);
        for (auto & v_index : obj.f[i].v)
            exVertices[v_index].ref_triangle.push_back(i);
    }

    loopi (0, exVertices.size())
    {
        ExVertex &v = exVertices[i];
        v_count.clear();
        v_ids.clear();
        for (auto & tri_index : v.ref_triangle)
        {
            ExTriangle &t = exTriangles[tri_index];
            loopk (0, 3)
            {
                int ofs = 0, v_index = obj.f[tri_index].v[k];

                for ( ; ofs < v_count.size(); ++ofs)
                    if (v_ids[ofs] == v_index) break;

                if (ofs == v_count.size())
                {
                    v_count.push_back(1);
                    v_ids.push_back(v_index);
                }
                else
                {
                    v_count[ofs]++;
                }
            }
        }
        loopj (0, v_count.size())
            if(v_count[j] == 1)
                exVertices[v_ids[j]].border = 1;
    }

    loopi (0, exTriangles.size())
    {
        ExTriangle &et = exTriangles[i];
        Triangle &tri = obj.f[i];
        Vec3f n, p[3];
        loopj (0, 3) p[j] = obj.v[tri.v[j]];
        n.cross(p[1] - p[0],p[2] - p[0]);
        n.normalize();
        et.normal = n;
        loopj (0, 3) exVertices[tri.v[j]].q += SymmetricMatrix(n.x, n.y, n.z, -n.dot(p[0]));
    }
    loopi (0, exTriangles.size())
    {
        // Calc Edge Error
        ExTriangle &et = exTriangles[i];
        Triangle &tri = obj.f[i];
        Vec3f p;
        loopj (0, 3)
        {
            et.err[j] = calculate_error(tri.v[j], tri.v[(j+1)%3], p);
        }
        et.err[3] = fmin(et.err[0], fmin(et.err[1], et.err[2]));
    }
    check_size("Load");
}
