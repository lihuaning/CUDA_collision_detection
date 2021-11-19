#ifndef OBJREAD_H
#define OBJREAD_H

/*
----read the obj file 
----sort the montor


*/
#include <iostream>
#include <fstream>
#include <string>
#include "tri_contact.cuh"
#include "object.cuh"
#include "definition.cuh"

#define MIN 1
#define MAX 0

typedef struct CentroidValue
{
    double min_x, max_x;
    double min_y, max_y;
    double min_z, max_z;
}CentroidValue;



void objread(const std::string objfile_path, 
            std::vector<vec3f> &vertices, 
            std::vector<Object> &triangles, 
            CentroidValue& alltri)
{
    std::ifstream in(objfile_path.c_str());

    if (!in.good())
    {
        std::cout << "open objfile failed!!" << std::endl;
    }

    char buffer[256];
    double x, y, z;
    Object temp_tri;
    unsigned int t1, t2, t3;
    
    //initial 
    alltri.min_x = MIN;
    alltri.min_y = MIN;
    alltri.min_z = MIN;
    alltri.max_x = MAX;
    alltri.max_y = MAX;
    alltri.max_z = MAX;

    while (!in.getline(buffer, 255).eof())
    {
        //read the vertices
        if (buffer[0] == 'v' && buffer[1] == ' ')
        {
            //sscanf return 读取成功的number
            if (sscanf(buffer, "v %lf %lf %lf", &x, &y, &z) == 3)
            {
                vertices.push_back(vec3f(x, y, z));

#ifdef CHAI_HAN
                alltri.min_x = x < alltri.min_x? x : alltri.min_x;
                alltri.min_y = y < alltri.min_y? y : alltri.min_y;
                alltri.min_z = z < alltri.min_z? z : alltri.min_z;

                alltri.max_x = x > alltri.max_x? x : alltri.max_x;
                alltri.max_y = y > alltri.max_y? y : alltri.max_y;
                alltri.max_z = z > alltri.max_z? z : alltri.max_z;
#endif
            }
            else
                std::cout << "ERROR: v number not wanted!!!!" << std::endl;
        }

        else if (buffer[0] == 'f' && buffer[1] == ' ')
        {
            if (sscanf(buffer, "f %d/%d %d/%d %d/%d", &temp_tri.idx[0], &t1, &temp_tri.idx[1], &t2, &temp_tri.idx[2], &t3) != 6)
            {
                std::cout << "ERROR: f number not wanted!!!!" << std::endl;
                exit(-1);
            }
                
            //avoid out of the border
            unsigned int versize = (unsigned int)vertices.size() + 1;
            if (temp_tri.idx[0] > versize || temp_tri.idx[1] > versize || temp_tri.idx[2] > versize)
                std::cout << "ERROR: f vertice index out of the border!!!!" << std::endl;

            //obj文件都是从1开始索引,因此减一之后才和vertice匹配
            temp_tri.idx[0]--;
            temp_tri.idx[1]--;
            temp_tri.idx[2]--;

            //计算centroid
            temp_tri.cx = (vertices[temp_tri.idx[0]].x + vertices[temp_tri.idx[1]].x + vertices[temp_tri.idx[2]].x) / 3;
            temp_tri.cy = (vertices[temp_tri.idx[0]].y + vertices[temp_tri.idx[1]].y + vertices[temp_tri.idx[2]].y) / 3;
            temp_tri.cz = (vertices[temp_tri.idx[0]].z + vertices[temp_tri.idx[1]].z + vertices[temp_tri.idx[2]].z) / 3;
            
#ifndef CHAI_HAN
            //update the minmax
            alltri.min_x = temp_tri.cx < alltri.min_x? temp_tri.cx : alltri.min_x;
            alltri.min_y = temp_tri.cy < alltri.min_y? temp_tri.cy : alltri.min_y;
            alltri.min_z = temp_tri.cz < alltri.min_z? temp_tri.cz : alltri.min_z;

            alltri.max_x = temp_tri.cx > alltri.max_x? temp_tri.cx : alltri.max_x;
            alltri.max_y = temp_tri.cy > alltri.max_y? temp_tri.cy : alltri.max_y;
            alltri.max_z = temp_tri.cz > alltri.max_z? temp_tri.cz : alltri.max_z;
#endif
            //Index of the object
            temp_tri.Index = (unsigned int)triangles.size() + 1;

            triangles.push_back(temp_tri);
        }
    }
}

#endif  //OBJREAD_H