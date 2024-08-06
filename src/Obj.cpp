//
// Created by 123 on 24-7-23.
//
#include <sstream>
#include <map>
#include <iostream>
#include "Obj.h"

void Obj::clear() {
    v.clear();
    vt.clear();
    vn.clear();
    f.clear();
    has_vn = false;
    has_vt = false;
    mtl = "";
}

void Obj::load(const string &content) {
    clear();
    istringstream file_stream(content);

    string line;
    string cur_mtl;
    int cur_mtl_index = -1;
    map<string, int> material_map;

    while(getline(file_stream, line))
    {
        if (line.empty()) continue;

        std::istringstream line_stream(line);
        string keyword;
        line_stream >> keyword;

        if (keyword == "mtllib")
        {
            line_stream >> mtl;
            continue;
        }
        if (keyword == "usemtl")
        {
            line_stream >> cur_mtl;
            if (material_map.find(cur_mtl) == material_map.end())
            {
                material_map[cur_mtl] = (int) materials.size();
                materials.push_back(cur_mtl);
            }
            cur_mtl_index = material_map[cur_mtl];
            continue;
        }
        if (keyword == "vt")
        {
            has_vt = true;
            Tex now;
            double a, b, c = 0.0;
            int count = sscanf(line.c_str(), "vt %lf %lf %lf", &a, &b, &c);
            if (count == 2)
                now.has_w = false;
            else
                now.has_w = true;
            now.uvw = {a, b, c};
            vt.push_back(now);
            continue;
        }
        if (keyword == "v")
        {
            Vec3f now;
            line_stream >> now.x >> now.y >> now.z;
            v.push_back(now);
            continue;
        }
        if (keyword == "vn")
        {
            has_vn = true;
            continue;
        }
        if (keyword == "f")
        {
            Triangle now{};

            now.mtl_index = cur_mtl_index;

            if (sscanf(line.c_str(), "f %d %d %d", &now.v[0], &now.v[1], &now.v[2]) == 3)
            {
                now.has_vn = false;
                now.has_vt = false;
            }
            else if(sscanf(line.c_str(), "f %d//%d %d//%d %d//%d",
                      &now.v[0], &now.vn[0],
                      &now.v[1], &now.vn[1],
                      &now.v[2], &now.vn[2]) == 6)
            {
                now.has_vn = true;
                now.has_vt = false;
            }
            else if(sscanf(line.c_str(), "f %d/%d/%d %d/%d/%d %d/%d/%d",
                      &now.v[0], &now.vt[0], &now.vn[0],
                      &now.v[1], &now.vt[1], &now.vn[1],
                      &now.v[2], &now.vt[2], &now.vn[2]) == 9)
            {
                now.has_vn = true;
                now.has_vt = true;
            }
            else if (sscanf(line.c_str(), "f %d/%d %d/%d %d/%d",
                            &now.v[0], &now.vt[0],
                            &now.v[1], &now.vt[1],
                            &now.v[2], &now.vt[2]) == 6)
            {
                now.has_vn = false;
                now.has_vt = true;
            }
            for (int i = 0; i < 3; ++i)
            {
                now.v[i] --;  now.vt[i] --; now.vn[i] --;
            }
            f.push_back(now);
        }
    }
}

void Obj::write(const char *filename) {
    FILE *file=fopen(filename, "w");
    int cur_mtl_index = -1;

    if (!file)
    {
        printf("write_obj: can't write data file \"%s\".\n", filename);
        exit(0);
    }

    if (!mtl.empty())
    {
        fprintf(file, "mtllib %s\n", mtl.c_str());
    }

    for (auto & i : v)
        fprintf(file, "v %g %g %g\n", i.x, i.y, i.z);
    if (has_vt)
        for (auto & i : vt)
        {
            if (i.has_w)
                fprintf(file, "vt %g %g %g\n", i.uvw.x, i.uvw.y, i.uvw.z);
            else
                fprintf(file, "vt %g %g\n", i.uvw.x, i.uvw.y);
        }
    if (has_vn)
        for (auto & i : vn)
            fprintf(file, "vn %g %g %g\n", i.x, i.y, i.z);

    for (auto & i : f)
    {
        if (i.mtl_index != cur_mtl_index)
        {
            cur_mtl_index = i.mtl_index;
            fprintf(file, "usemtl %s\n", materials[i.mtl_index].c_str());
        }
        if (i.has_vn && i.has_vt)
        {
            fprintf(file, "f %d/%d/%d %d/%d/%d %d/%d/%d\n",
                    i.v[0]+1, i.vt[0]+1, i.v[0]+1,
                    i.v[1]+1, i.vt[1]+1, i.v[1]+1,
                    i.v[2]+1, i.vt[2]+1, i.v[2]+1);
        }
        else if (i.has_vn)
        {
            fprintf(file, "f %d//%d %d//%d %d//%d\n",
                    i.v[0]+1, i.v[0]+1,
                    i.v[1]+1, i.v[1]+1,
                    i.v[2]+1, i.v[2]+1);
        }
        else if (i.has_vt)
        {
            fprintf(file, "f %d/%d %d/%d %d/%d\n",
                    i.v[0]+1, i.vt[0]+1,
                    i.v[1]+1, i.vt[1]+1,
                    i.v[2]+1, i.vt[2]+1);
        }
        else
        {
            fprintf(file, "f %d %d %d\n",
                    i.v[0]+1, i.v[1]+1, i.v[2]+1);
        }
    }
}

string Obj::write_to_string() {
    ostringstream stream;
    int cur_mtl_index = -1;

    if (!mtl.empty())
    {
        stream << "mtllib " << mtl << "\n";
    }

    for (auto & i : v)
        stream << "v " << i.x << " " << i.y << " " << i.z << "\n";
    if (has_vt)
        for (auto & i : vt)
        {
            if (i.has_w)
                stream << "vt " << i.uvw.x << " " << i.uvw.y << " " << i.uvw.z << "\n";
            else
                stream << "vt " << i.uvw.x << " " << i.uvw.y << "\n";
        }
    if (has_vn)
        for (auto & i : vn)
            stream << "vn " << i.x << " " << i.y << " " << i.z << "\n";

    for (auto & i : f)
    {
        if (i.mtl_index != cur_mtl_index)
        {
            cur_mtl_index = i.mtl_index;
            stream << "usemtl " << materials[i.mtl_index] << "\n";
        }
        if (i.has_vn && i.has_vt)
        {
            stream << "f " << i.v[0] + 1 << "/" << i.vt[0] + 1 << "/" << i.v[0] + 1 << " "
                   << i.v[1] + 1 << "/" << i.vt[1] + 1 << "/" << i.v[1] + 1 << " "
                   << i.v[2] + 1 << "/" << i.vt[2] + 1 << "/" << i.v[2] + 1 << "\n";
        }
        else if (i.has_vn)
        {
            stream << "f " << i.v[0] + 1 << "//" << i.v[0] + 1 << " "
                   << i.v[1] + 1 << "//" << i.v[1] + 1 << " "
                   << i.v[2] + 1 << "//" << i.v[2] + 1 << "\n";
        }
        else if (i.has_vt)
        {
            stream << "f " << i.v[0] + 1 << "/" << i.vt[0] + 1 << " "
                   << i.v[1] + 1 << "/" << i.vt[1] + 1 << " "
                   << i.v[2] + 1 << "/" << i.vt[2] + 1 << "\n";
        }
        else
        {
            stream << "f " << i.v[0] + 1 << " " << i.v[1] + 1 << " " << i.v[2] + 1 << "\n";
        }
    }
    return stream.str();
}