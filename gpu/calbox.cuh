#ifndef CALBOX
#define CALBOX

#include "Node.cuh"
#include "tri_contact.cuh"

// #define Max(a, b, c) (a) > (b)? ((a) > (c)? (a) : (c)) : ((b) > (c)? (b):(c))
//#define Min(a, b, c) (a) < (b)? ((a) < (c)? (a) : (c)) : ((b) < (c)? (b):(c))

__device__ double Max(double a, double b, double c)
{
    return a > b ? (a > c ? (a) : (c)) : ((b) > (c) ? (b) : (c));
}

__device__ double Min(double a, double b, double c)
{
    //printf("min  \n");
    return a < b ? ((a) < (c) ? (a) : (c)) : ((b) < (c) ? (b) : (c));
}

/*
* 计算all leaf nodes的包围盒
*/
__global__ void CalAllLeafBoundingBox(Node *LeafNodes,
                                      Object *triangles,
                                      vec3f *vertices,
                                      unsigned int numObjects)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x; // 当前线程的坐标
    int skip = blockDim.x * gridDim.x;               // 一个grid内的所有线程数量，作为递增量
                                    
    for (unsigned int i = tid; i < numObjects; i += skip)
    {
        //printf("node %d \n", i);
        unsigned int vertice1 = triangles[i].idx[0];
        //printf("hhhhhhhhhhhh \n");
        unsigned int vertice2 = triangles[i].idx[1];
        //printf("ning \n");
        unsigned int vertice3 = triangles[i].idx[2];
        //printf("ooooooooooooooooo \n");

        LeafNodes[i].bounding.LeftDown.x = Min(vertices[vertice1].x, vertices[vertice2].x, vertices[vertice3].x);
        LeafNodes[i].bounding.LeftDown.y = Min(vertices[vertice1].y, vertices[vertice2].y, vertices[vertice3].y);
        LeafNodes[i].bounding.LeftDown.z = Min(vertices[vertice1].z, vertices[vertice2].z, vertices[vertice3].z);

        LeafNodes[i].bounding.RightUp.x = Max(vertices[vertice1].x, vertices[vertice2].x, vertices[vertice3].x);
        LeafNodes[i].bounding.RightUp.y = Max(vertices[vertice1].y, vertices[vertice2].y, vertices[vertice3].y);
        LeafNodes[i].bounding.RightUp.z = Max(vertices[vertice1].z, vertices[vertice2].z, vertices[vertice3].z);

        LeafNodes[i].boxCaled = true;
    }
}

/*r
* double3类型 max
*/
__device__ double3 Max_double3(double3 a, double3 b)
{
    double3 temp;

    temp.x = max(a.x, b.x);
    temp.y = max(a.y, b.y);
    temp.z = max(a.z, b.z);

    return temp;
}

/*
* double3类型 min
*/
__device__ double3 Min_double3(double3 a, double3 b)
{
    double3 temp;

    temp.x = min(a.x, b.x);
    temp.y = min(a.y, b.y);
    temp.z = min(a.z, b.z);

    return temp;
}

/*
* 计算单个internal node的包围盒
*/
__device__ void CalInternalBoundingBox(Node* internal)
{
    // if(internal.ObjectIndex == 3){
    //     cout<<"here"<<endl;
    // }
    Node* leftChild = internal->LeftChild;
    Node* rightChild = internal->RightChild;

    internal->bounding.RightUp = Max_double3(leftChild->bounding.RightUp, rightChild->bounding.RightUp);
    internal->bounding.LeftDown = Min_double3(leftChild->bounding.LeftDown, rightChild->bounding.LeftDown);
}

/*
* 为每个Node计算其包围盒
* 循环从每个叶子节点向上搜索
*/
__global__ void CalculateBVHBoundingBox_Violence(Node* LeafNodes,
                                                const unsigned int numObjects)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x; // 当前线程的坐标
    int skip = blockDim.x * gridDim.x;               // 一个grid内的所有线程数量，作为递增量
                                    
    for (unsigned int i = tid; i < numObjects; i += skip)
    {
        Node* p = LeafNodes[i].parent;
        while(p)
        {
            if(!atomicAdd(&(p->boxCaled), 1))
            {       
                break;
            }
            CalInternalBoundingBox(p);
            p = p->parent;
        }
    }
}

/*
* what a fucking operation!!! cuda!!!!
*/
__global__ void CalRoot(Node* root)
{
    CalInternalBoundingBox(root);
    root->boxCaled = 2;
}

#endif