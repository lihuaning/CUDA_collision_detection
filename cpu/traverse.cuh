#ifndef TRAVERSE_H
#define TRAVERSE_H

#include "definition.cuh"
#include "box.cuh"
#include "Node.cuh"
#include "tri_contact.h"
#include <algorithm>


/*
* 检查两个box是否存在覆盖情况
* 将长方体映射到三个维度，只有三个维度都相交才说明长方体是相交的
*/
bool checkOverlap(const Box& query,
                    const Box& checked
                    )
{
    if(query.LeftDown.x > checked.RightUp.x || checked.LeftDown.x > query.RightUp.x){
        return false;
    }
    if(query.LeftDown.y > checked.RightUp.y || checked.LeftDown.y > query.RightUp.y){
        return false;
    }
    if(query.LeftDown.z > checked.RightUp.z || checked.LeftDown.z > query.RightUp.z){
        return false;
    }        
    return true;
}



void traverseIteractive(CollisionList& list,
                        Node* BVHTree,
                        Node* LeafNodes,
                        const Box& query,
                        const unsigned int& queryObjectIdx)
{
    // Allocate traversal stack from thread-local memory,
    // and push NULL to indicate that there are no postponed nodes.
    Node* stack[64];
    Node** stackPtr = stack;
    *stackPtr++ = NULL; // push

    // Traverse nodes starting from the root.
    Node* node = BVHTree;
    //int over = 0;
    do
    {
        // Check each child node for overlap.
        Node* childL = node->LeftChild;
        Node* childR = node->RightChild;

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
            list.push_back(potential);
        }

        if (overlapR && LeafNodes[queryObjectIdx].ObjectIndex < childR->ObjectIndex)
        {
            int2 potential;
            potential.x = queryObjectIdx;
            potential.y = childR->SortedObjectIndex;
            list.push_back(potential);
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
    }while (node != NULL);
    
    //cout<<"overlap num: "<<over<<endl;
}

/*
* 遍历所有的objects，find collision
* list是存在碰撞的包围盒，不是三角形
*/
void FindPotentialCollisions(Node* BVHTree,
                            Node* LeafNodes, 
                            CollisionList &list,
                            const unsigned int &numObjects)
{
    
    for(unsigned int i = 0; i < numObjects; i++)
    {
        traverseIteractive(list, BVHTree, LeafNodes, LeafNodes[i].bounding, i);
    }
    
}



bool share_vertice(const Object &triangle_1,
                    const Object &triangle_2)
{
    return (triangle_1.idx[0] == triangle_2.idx[0] || 
            triangle_1.idx[0] == triangle_2.idx[1] ||
            triangle_1.idx[0] == triangle_2.idx[2] ||
            triangle_1.idx[1] == triangle_2.idx[0] ||
            triangle_1.idx[1] == triangle_2.idx[1] ||
            triangle_1.idx[1] == triangle_2.idx[2] ||
            triangle_1.idx[2] == triangle_2.idx[0] ||
            triangle_1.idx[2] == triangle_2.idx[1] ||
            triangle_1.idx[2] == triangle_2.idx[2] )  ? true : false;
}



/*
* 确认box中的object是否overlap
*/
void VerifyCollisions(std::vector<unsigned int> &collision_list,
                        CollisionList &collision_pairs,
                        const CollisionList &pairs,
                        std::vector<vec3f> &vertices, 
                        std::vector<Object> &triangles
                        )
{
    
    for(unsigned int i = 0; i < (unsigned int)pairs.size(); i++)
    {
        // object index start from 1, but vector subscript start from 0
        unsigned int triangle_1 = pairs[i].x;
        unsigned int triangle_2 = pairs[i].y;
        
        vec3f P1 = vertices[triangles[triangle_1].idx[0]];
        vec3f P2 = vertices[triangles[triangle_1].idx[1]];
        vec3f P3 = vertices[triangles[triangle_1].idx[2]];

        vec3f Q1 = vertices[triangles[triangle_2].idx[0]];
        vec3f Q2 = vertices[triangles[triangle_2].idx[1]];
        vec3f Q3 = vertices[triangles[triangle_2].idx[2]];

        if(!share_vertice(triangles[triangle_1], triangles[triangle_2]) && tri_contact(P1, P2, P3, Q1, Q2, Q3))
        {
            collision_list.push_back(triangles[pairs[i].x].Index);
            collision_list.push_back(triangles[pairs[i].y].Index);
            int2 nono;
            nono.x = triangles[pairs[i].x].Index;
            nono.y = triangles[pairs[i].y].Index;
            collision_pairs.push_back(nono);
        }
    }

    sort(collision_list.begin(), collision_list.end());
    collision_list.erase(unique(collision_list.begin(), collision_list.end()), collision_list.end());
}


#endif