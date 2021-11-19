/*
* 建树  
* 建立包围盒  
*遍历BVH树
*/

#ifndef BVH_H
#define BVH_H

#include "definition.cuh"
#include "Node.cuh"
#include "object.cuh"
#include <vector>
#include  <algorithm>
#include <stdlib.h>
#include <stdio.h>
using namespace std;

inline int myclz(mortonCode code)
{
#ifdef MORTON_CODE_64
    if (!code)
        return 64;
    int e = 63;
    //1111 1111 1111 1111 0000 0000 0000 0000
    if (code & 0xFFFFFFFF00000000)
    {
        e -= 32;
        code >>= 32;
    }
    //0000 0000 0000 0000 1111 1111 0000 0000
    if (code & 0x00000000FFFF0000)
    {
        e -= 16;
        code >>= 16;
    }
    //0000 0000 0000 0000 0000 0000 1111 0000
    if (code & 0x000000000000FF00)
    {
        e -= 8;
        code >>= 8;
    }
    //0000 0000 0000 0000 0000 0000 1111 0000
    if (code & 0x00000000000000F0)
    {
        e -= 4;
        code >>= 4;
    }
    //0000 0000 0000 0000 0000 0000 0000 1100
    if (code & 0x000000000000000C)
    {
        e -= 2;
        code >>= 2;
    }
    //0000 0000 0000 0000 0000 0000 0000 0010
    if (code & 0x0000000000000002)
    {
        e -= 1;
    }

    return e;
#else
    if (!code)
        return 32;
    int e = 31;
    //1111 1111 1111 1111 0000 0000 0000 0000
    if (code & 0xFFFF0000)
    {
        e -= 16;
        code >>= 16;
    }
    //0000 0000 0000 0000 1111 1111 0000 0000
    if (code & 0x0000FF00)
    {
        e -= 8;
        code >>= 8;
    }
    //0000 0000 0000 0000 0000 0000 1111 0000
    if (code & 0x000000F0)
    {
        e -= 4;
        code >>= 4;
    }
    //0000 0000 0000 0000 0000 0000 0000 1100
    if (code & 0x0000000C)
    {
        e -= 2;
        code >>= 2;
    }
    //0000 0000 0000 0000 0000 0000 0000 0010
    if (code & 0x00000002)
    {
        e -= 1;
    }

    return e;
#endif
    
}

/*find the split index*/
int findSplit(const std::vector<mortonCode> &MortonCodes,
              const unsigned int &first,
              const unsigned int &last)
{
    mortonCode first_montor = MortonCodes[first];
    mortonCode last_montor = MortonCodes[last];

    // if same, return the middle position of them
    if (first_montor == last_montor)
    {
        return (first + last) >> 1;
    }

    int commonPrefix = myclz(first_montor ^ last_montor);

    // Use binary search to find where the next bit differs.
    // Specifically, we are looking for the highest object that
    // shares more than commonPrefix bits with the first one.

    unsigned int split = first; // initial guess
    unsigned int step = last - first;

    // ********************************
    do
    {
        step = (step + 1) >> 1;      // exponential decrease
        unsigned int newSplit = split + step; // proposed new position

        if (newSplit < last)
        {
            mortonCode split_montor = MortonCodes[newSplit];
            int splitPrefix = myclz(first_montor ^ split_montor);
            if (splitPrefix > commonPrefix)
                split = newSplit; // accept proposal
        }
    } while (step > 1);

    return split;
}

int deltaCpu(const unsigned int i, 
            const unsigned int j, 
            const std::vector<mortonCode>& MortonCodes, 
            const unsigned int n) 
{
    return ((j >= 0 && j < n) ? myclz(MortonCodes[i] ^ MortonCodes[j]) : -1);
}


/*find the range of objects in the internal node*/
int2 determineRange(const std::vector<mortonCode> &MortonCodes,
                    const unsigned int &numObjects,
                    const unsigned int &i)
{
    int d = (deltaCpu(i, i + 1, MortonCodes, numObjects) - deltaCpu(i, i - 1, MortonCodes, numObjects)) >= 0 ? 1 : -1;
    int delta_min = deltaCpu(i, i - d, MortonCodes, numObjects);
    int mlen = 2;
    //return Range(100, mlen);
    while (deltaCpu(i, i + mlen * d, MortonCodes, numObjects) > delta_min) {
        mlen <<= 1;
    }

    int l = 0;
    for (int t = mlen >> 1; t >= 1; t >>= 1) {
        if (deltaCpu(i, i + (l + t) * d, MortonCodes, numObjects) > delta_min) {
            l += t;
        }
    }
    unsigned int j = i + l * d;

    int2 nono;
    nono.x = min(i, j);
    nono.y = max(i, j);

    // 返回的范围中，前者一定是较小者，后者一定是较大者
    return nono;
}

/*generate a BVH tree by iteration*/
void GenerateBVH(std::vector<Object> &triangles,
                     const std::vector<mortonCode> &MortonCodes,
                     Node *LeafNodes, Node *InternalNodes,
                     const unsigned int &numObjects)
{

    //Construct Leaf Nodes
    for (unsigned int i = 0; i < numObjects; i++)
    {
        LeafNodes[i].ObjectIndex = triangles[i].Index;
        LeafNodes[i].SortedObjectIndex = i;
    }
    //Construct interal Nodes
    for (unsigned int i = 0; i < numObjects - 1; i++)
    {
        // Find out which range of objects the node corresponds to.
        // (This is where the magic happens!)

        int2 range = determineRange(MortonCodes, numObjects, i);
        int first = range.x;                        
        int last = range.y;

        // Determine where to split the range.
        int split = findSplit(MortonCodes, first, last);

        // Select leftChild.

        Node *left;
        if (split == first)
            left = &LeafNodes[split];
        else
            left = &InternalNodes[split];

        // Select rightChild.

        Node *right;
        if (split + 1 == last)
            right = &LeafNodes[split + 1];
        else
            right = &InternalNodes[split + 1];

        // Record parent-child relationships.
        InternalNodes[i].LeftChild = left;
        InternalNodes[i].RightChild = right;
        left->parent = &InternalNodes[i];
        right->parent = &InternalNodes[i];
    }
}

#endif //BVH_H