#ifndef NODE_H
#define NODE_H

#include "box.cuh"

typedef struct Node
{
    Node* LeftChild;
    Node* RightChild;
    Node* parent;
    Box bounding;
    bool isLeaf;
    int boxCaled;  //是否已经计算过该node的包围盒
    unsigned int ObjectIndex;  //obj文件中的排序， 从1开始
    unsigned int SortedObjectIndex; // sorted order 数组中的位置，从0开始
}Node;

__global__ void InitNodes(Node* Nodes, unsigned int numObjects, bool leaf)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;// 当前线程的坐标
    int skip = blockDim.x * gridDim.x; // 一个grid内的所有线程数量，作为递增量

    for(unsigned int i = tid; i < numObjects; i += skip){
        Nodes[i].LeftChild = NULL;
        Nodes[i].RightChild = NULL;
        Nodes[i].parent = NULL;
        Nodes[i].isLeaf = leaf;
        Nodes[i].boxCaled = 0;
        Nodes[i].ObjectIndex = 0;
        Nodes[i].SortedObjectIndex = 0;
    }
}

#endif  //NODE_H