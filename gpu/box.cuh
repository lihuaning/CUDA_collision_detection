#ifndef BOX_H
#define BOX_H

#include <stack>
#include <vector>

#include "definition.cuh"
#include "Node.cuh"
#include "tri_contact.cuh"
#include "object.cuh"

typedef struct Box
{
    double3 RightUp;
    double3 LeftDown;
}Box;

#endif //BOX_H

