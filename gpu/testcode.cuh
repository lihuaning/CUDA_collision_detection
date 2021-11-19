#ifndef TEST_H
#define TEST_H

#include "definition.cuh"
#include "morton.cuh"
#include "object.cuh"
#include "Node.cuh"
#include "bvh.cuh"
#include "calbox.cuh"
#include "traverse.cuh"


__global__ void kernel(Object* triangles)
{
	for(int i = 0; i < 10; i++){
		//printf("%d \n", i);
		printf("%lld \n", triangles[i].morton);
        printf("%lf %lf %lf \n", triangles[i].cx, triangles[i].cy, triangles[i].cz);
	}
}

__global__ void kernel_morton(mortonCode* MortonCodes)
{
	for(int i = 0; i < 10; i++){
		//printf("%d \n", i);
		printf("%lld \n", MortonCodes[i]);
	}
}

__global__ void kernel_node(Node* InternalNodes)
{
	Node* p = InternalNodes;
	if(p->LeftChild)
		printf("left");
	
	if(p->RightChild)
		printf("right");

	// int i = 0;
	// while(p->LeftChild){
	// 	i++;
	// 	p = p->LeftChild;
	// }
	// printf("left depth: %d \n", i);

	// i = 0;

	// p = InternalNodes;
	// while(p->RightChild){
	// 	i++;
	// 	p = p->RightChild;
	// }
	// printf("right depth: %d \n", i);
}

__global__ void kernel_nodebox(Node* Nodes, unsigned int N)
{
	printf("%d \n", Nodes->boxCaled);
	int n = 0;
	for(int i = 0; i < N -1; i++){
		if(Nodes[i].boxCaled < 2){
			n++;
		}
	}
	printf("no cal num %d \n", n);


	printf("kernel_nodebox \n");
	for(int i = 0; i < 10; i++){
		printf("Node %d \n", i);
		printf("rightup %lf %lf %lf \n", Nodes[i].bounding.RightUp.x, Nodes[i].bounding.RightUp.y, Nodes[i].bounding.RightUp.z);
		printf("leftdown %lf %lf %lf \n", Nodes[i].bounding.LeftDown.x, Nodes[i].bounding.LeftDown.y, Nodes[i].bounding.LeftDown.z);
	}
}

__global__ void kernel_collision(int2* list, unsigned int N)
{
	for(int i = 0; i < N; i++){
		printf("pair %d \n", i);
		printf("%d %d \n", list[i].x, list[i].y);	
	}	
}
__global__ void kernel_count(int* count)
{
	printf("count: %d \n", *count);	
}

#endif