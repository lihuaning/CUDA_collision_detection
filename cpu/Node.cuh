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


void CreateNodes(Node* &Nodes, const unsigned int N)
{
    Nodes = new Node[N];
    if(!Nodes)
        printf("Nodes malloc failed!!!!");
}

void InitNodes(Node* Nodes, const unsigned int N, const bool leaf)
{
    for(unsigned int i = 0; i < N; i++)
    {
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