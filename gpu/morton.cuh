#pragma once
#ifndef MORTON_H
#define MORTON_H

#include <thrust/sort.h>
#include <thrust/device_ptr.h>
#include  <algorithm>
#include <vector>
#include "object.cuh"
#include "definition.cuh"   
#include "objread.cuh"



__device__ mortonCode expandBits(mortonCode v)
{
#ifdef MORTON_CODE_64
    v = (v | (v << 32)) & 0xFFFF00000000FFFFull;
    v = (v | (v << 16)) & 0x00FF0000FF0000FFull;
    v = (v | (v << 8)) & 0xF00F00F00F00F00Full;
    v = (v | (v << 4)) & 0x30C30C30C30C30C3ull;
    v = (v | (v << 2)) & 0x9249249249249249ull;

#else
    v = (v * 0x00010001u) & 0xFF0000FFu;
    v = (v * 0x00000101u) & 0x0F00F00Fu;
    v = (v * 0x00000011u) & 0xC30C30C3u;
    v = (v * 0x00000005u) & 0x49249249u;

#endif
    return v;
}

//将质心坐标norm到[0,1]
__device__ double Norm(double element, double minval, double maxval)
{
    element = (element - minval)/(maxval - minval);

    return element;
}

__device__ mortonCode morton3D(double x, double y, double z, int tid)
{     

#ifdef MORTON_CODE_64
    x = min(max(x * 1048576.0f, 0.0f), 1048575.0f);
    y = min(max(y * 1048576.0f, 0.0f), 1048575.0f);
    z = min(max(z * 1048576.0f, 0.0f), 1048575.0f);
#else
    x = min(max(x * 1024.0f, 0.0f), 1023.0f);
    y = min(max(y * 1024.0f, 0.0f), 1023.0f);
    z = min(max(z * 1024.0f, 0.0f), 1023.0f);
#endif

    mortonCode xx = expandBits((mortonCode)x);
    mortonCode yy = expandBits((mortonCode)y);
    mortonCode zz = expandBits((mortonCode)z);

    // if(tid == 0){
    //     printf("morton3D code: \n");
    //     printf("%d \n", xx * 4 + yy * 2 + zz);
    // }

    return xx * 4 + yy * 2 + zz;

}

/*
* kernal of calculate Morton codes
*/
__global__ void MortonCode(Object* triangles, 
                            const CentroidValue &alltri, 
                            unsigned int numObjects)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;// 当前线程的坐标
    int skip = blockDim.x * gridDim.x; // 一个grid内的所有线程数量，作为递增量
    /*
    norm the centriod into [0,1]
    calculate the morton code 
    */
    double norm_cx,norm_cy,norm_cz; 
    
    for(unsigned int i = tid; i < numObjects; i += skip){
        norm_cx = Norm(triangles[i].cx, alltri.min_x, alltri.max_x);
        norm_cy = Norm(triangles[i].cy, alltri.min_y, alltri.max_y);
        norm_cz = Norm(triangles[i].cz, alltri.min_z, alltri.max_z);

        triangles[i].morton = morton3D(norm_cx, norm_cy, norm_cz, i);
    }
}

__global__ void AssignMortonWithObjects(Object* triangles, 
                                        mortonCode* MortonCodes, 
                                        int numObjects)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;// 当前线程的坐标
    int skip = blockDim.x * gridDim.x; // 一个grid内的所有线程数量，作为递增量

    for(unsigned int i = tid; i < numObjects; i += skip){
        MortonCodes[i] = triangles[i].morton;
    }
}


#endif 


