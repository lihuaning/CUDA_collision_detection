#pragma once
#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "definition.cuh"

/*形状为三角形的结构体  */
typedef struct Object
{
    double cx,cy,cz;  //centroid of object
    unsigned int idx[3]; //三个顶点各自的index
    unsigned int Index; //自己的index  start from 1
    mortonCode morton; // morton code of the triangle

    __host__ __device__ bool operator<(const Object& t2)const{
        // printf("%d\n",initialID);
        return morton<t2.morton;//从小到大排则<
    }
}Object;
#endif