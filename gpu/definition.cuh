/*
* cuda类型的定义
* 宏预定义
*/


#ifndef DEFINITION_H
#define DEFINITION_H

#include <vector>

//#define TEST
//#define TEST_NODE
//#define TEST_NODE_BOX
//#define TEST_POTENTIAL_PAIR
//#define TEST_CUDA_DEVICE_MEMORY
//#define TEST_TRIANGLE_MONTORCODE_CAL
//#define TEST_TRIANGLE_MONTORCODE_SORTED
//#define TEST_MONTORCODE_SORTED
//#define TEST_COLLISION

//#define CHAI_HAN


#define MORTON_CODE_64
//#define MORTON_CODE_U32


//typedef unsigned int mortonCode;

typedef struct __device_builtin__ int2 int2;
typedef struct __device_builtin__ double3 double3;

typedef unsigned long long mortonCode;
typedef std::vector<int2> CollisionList;

#define CAL_BOX_STACK 64

#endif //DEFINITION_H