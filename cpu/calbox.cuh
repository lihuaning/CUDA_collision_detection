#ifndef CALBOX
#define CALBOX

#include "Node.cuh"
#include "tri_contact.h"

#define Max(a, b, c) (a) > (b)? ((a) > (c)? (a) : (c)) : ((b) > (c)? (b):(c))
#define Min(a, b, c) (a) < (b)? ((a) < (c)? (a) : (c)) : ((b) < (c)? (b):(c))

/*
* 计算all leaf nodes的包围盒
*/
void CalAllLeafBoundingBox(Node* &LeafNodes,
                            const unsigned int &numObjects,
                            const std::vector<Object>& triangles,
                            const std::vector<vec3f>& vertices)
{
    for(unsigned int i = 0; i<numObjects; i++)
    { 
        unsigned int vertice1 = triangles[i].idx[0];
        unsigned int vertice2 = triangles[i].idx[1];
        unsigned int vertice3 = triangles[i].idx[2];

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
inline double3 Max_double3(double3 a, double3 b)
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
inline double3 Min_double3(double3 a, double3 b)
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
void CalInternalBoundingBox(Node &internal)
{   
    // if(internal.ObjectIndex == 3){
    //     cout<<"here"<<endl;
    // }
    Node* leftChild = internal.LeftChild;
    Node* rightChild = internal.RightChild;

    internal.bounding.RightUp = Max_double3(leftChild->bounding.RightUp, rightChild->bounding.RightUp);
    internal.bounding.LeftDown = Min_double3(leftChild->bounding.LeftDown, rightChild->bounding.LeftDown);
}

/*
* 为每个Node计算其包围盒
* 循环从每个叶子节点向上搜索
*/
void CalculateBVHBoundingBox_Violence(Node *LeafNodes, 
                                        const unsigned int &numObjects)
{
    //遍历叶子节点
    for(unsigned int i = 0; i < numObjects; i++)
    {
        Node* p = LeafNodes[i].parent;
        while(p)
        {
            // 还没有计算包围盒
            if(!p->boxCaled)
            {
                p->boxCaled += 1;
                break;
            }

            CalInternalBoundingBox(*p);
            p = p->parent;
        }
    }
}

/*
* 为每个Node计算其包围盒
* 后序遍历
* 不适合GPU并行操作
*/
void CalculateBVHBoundingBox_Post(Node *LeafNodes, Node *InternalNodes,
                                    const std::vector<Object>& triangles)
{
    //Node *stack[CAL_BOX_STACK]; //建立一个元素是Node地址的Stack
    std::stack<Node*> node_stack;
    Node* p = NULL;
    Node* pre = NULL;

    node_stack.push(InternalNodes); // push root into stack
    while(!node_stack.empty())
    {
        //一直找到最左叶子节点
        p = node_stack.top();
        while(p){
            node_stack.push(p->LeftChild);
            p = p->LeftChild;
        }

        node_stack.pop(); //将NULL退栈
        
        if(!node_stack.empty()){
            p = node_stack.top();

            //如果是叶子节点，或者其右孩子已经被访问过，那么说明
            //这个节点可以被当做中间节点
            //对internal 和 leaf 进行判断，计算包围盒
            if(p->RightChild == NULL || p->RightChild == pre)
            {
                //Leaf Node
                if(p->isLeaf){
                    ;
                }
                else{
                    //calculation internal box
                    //CalInternalBoundingBox();
                }
                node_stack.pop();
                pre = p;
                node_stack.push(NULL);
            }
            else{
                node_stack.push(p->RightChild);
            }
        }
    }
}



#endif