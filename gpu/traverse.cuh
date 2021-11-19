#ifndef TRAVERSE_H
#define TRAVERSE_H

#include <thrust/device_vector.h>
#include "definition.cuh"
#include "box.cuh"
#include "Node.cuh"
#include "tri_contact.cuh"
#include <algorithm>

__device__ bool share_vertice(Object* triangle_1,
                              Object* triangle_2)
{
    if((triangle_1->idx[0] == triangle_2->idx[0]) ||
            triangle_1->idx[0] == triangle_2->idx[1] ||
            triangle_1->idx[0] == triangle_2->idx[2] ||
            triangle_1->idx[1] == triangle_2->idx[0] ||
            triangle_1->idx[1] == triangle_2->idx[1] ||
            triangle_1->idx[1] == triangle_2->idx[2] ||
            triangle_1->idx[2] == triangle_2->idx[0] ||
            triangle_1->idx[2] == triangle_2->idx[1] ||
            triangle_1->idx[2] == triangle_2->idx[2])
    {
        return true;
    }
    else{
        return false;
    }
}

/*
* 确认box中的object是否overlap
*/
__device__ bool VerifyCollisions(int2 *pairs,
                                 vec3f *vertices,
                                 Object *triangles)
{
    int total_num = 0;

    // object index start from 1, but vector subscript start from 0
    unsigned int triangle_1 = (*pairs).x;
    unsigned int triangle_2 = (*pairs).y;

    vec3f P1 = vertices[triangles[triangle_1].idx[0]];
    vec3f P2 = vertices[triangles[triangle_1].idx[1]];
    vec3f P3 = vertices[triangles[triangle_1].idx[2]];

    vec3f Q1 = vertices[triangles[triangle_2].idx[0]];
    vec3f Q2 = vertices[triangles[triangle_2].idx[1]];
    vec3f Q3 = vertices[triangles[triangle_2].idx[2]];

    if (!share_vertice(&triangles[triangle_1], &triangles[triangle_2]) && tri_contact(P1, P2, P3, Q1, Q2, Q3))
    {
        return true;
    }
    else
    {
        return false;
    }
    return true;
}


/*
* 检查两个box是否存在覆盖情况
* 将长方体映射到三个维度，只有三个维度都相交才说明长方体是相交的
*/
__device__ bool checkOverlap(Box query,
                             Box checked)
{
    if (query.LeftDown.x > checked.RightUp.x || checked.LeftDown.x > query.RightUp.x)
    {
        return false;
    }
    if (query.LeftDown.y > checked.RightUp.y || checked.LeftDown.y > query.RightUp.y)
    {
        return false;
    }
    if (query.LeftDown.z > checked.RightUp.z || checked.LeftDown.z > query.RightUp.z)
    {
        return false;
    }
    return true;
}

/*
count: list 的原子操作下标
*/
__device__ void traverseIteractive(int2 *list,
                                   Node *BVHTree,
                                   Node *LeafNodes,
                                   Box query,
                                   vec3f *vertices,
                                   Object *triangles,
                                   unsigned int queryObjectIdx,
                                   int *count)
{
    // Allocate traversal stack from thread-local memory,
    // and push NULL to indicate that there are no postponed nodes. 
    Node *stack[64];
    Node **stackPtr = stack;
    *stackPtr++ = NULL; // push

    // Traverse nodes starting from the root.
    Node *node = BVHTree;
    //int over = 0;
    do
    {
        // Check each child node for overlap.
        Node *childL = node->LeftChild;
        Node *childR = node->RightChild;

        bool overlapL = checkOverlap(query,
                                     childL->bounding);
        bool overlapR = checkOverlap(query,
                                     childR->bounding);

        //over += overlapL+overlapR;

        // Query overlaps a leaf node => report collision.
        // avoid pair repeat
        // tips: set the internal ObjectIndex = 0
        if (overlapL && LeafNodes[queryObjectIdx].ObjectIndex < childL->ObjectIndex)
        {
            int2 potential;
            potential.x = queryObjectIdx;
            potential.y = childL->SortedObjectIndex;
            if(VerifyCollisions(&potential, vertices, triangles)){
                int2 temp;
                temp.x = triangles[potential.x].Index;
                temp.y = triangles[potential.y].Index;

                list[atomicAdd(count, 1)] = temp;
            }
        }

        if (overlapR && LeafNodes[queryObjectIdx].ObjectIndex < childR->ObjectIndex)
        {
            int2 potential;
            potential.x = queryObjectIdx;
            potential.y = childR->SortedObjectIndex;
            if(VerifyCollisions(&potential, vertices, triangles)){
                int2 temp;
                temp.x = triangles[potential.x].Index;
                temp.y = triangles[potential.y].Index;

                list[atomicAdd(count, 1)] = temp;
            }
        }

        // Query overlaps an internal node => traverse.
        bool traverseL = (overlapL && !childL->isLeaf);
        bool traverseR = (overlapR && !childR->isLeaf);

        if (!traverseL && !traverseR)
            node = *--stackPtr; // pop
        else
        {
            node = (traverseL) ? childL : childR;
            if (traverseL && traverseR)
                *stackPtr++ = childR; // push
        }
    } while (node != NULL);

    //cout<<"overlap num: "<<over<<endl;
}

/*
* 遍历所有的objects，find collision
* list是三角形真正碰撞的数组
*/
__global__ void FindPotentialCollisions(Node *BVHTree,
                                        Node *LeafNodes,
                                        int2 *list,
                                        vec3f *vertices,
                                        Object *triangles,
                                        unsigned int numObjects,
                                        int *count)
{   
    int tid = threadIdx.x + blockIdx.x * blockDim.x; // 当前线程的坐标
    int skip = blockDim.x * gridDim.x;               // 一个grid内的所有线程数量，作为递增量

    for (unsigned int i = tid; i < numObjects; i += skip)
    {	
        traverseIteractive(list, BVHTree, LeafNodes, LeafNodes[i].bounding, vertices, triangles, i, count);	
    }
}

#endif