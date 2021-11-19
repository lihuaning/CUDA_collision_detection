#include <thrust/sort.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/copy.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include <vector>
#include <algorithm>
#include <string>
#include <iostream>
#include  <algorithm>
#include "definition.cuh"
#include "objread.cuh"
#include "morton.cuh"
#include "object.cuh"
#include "../common/book.h"
#include "Node.cuh"
#include "bvh.cuh"
#include "calbox.cuh"
#include "traverse.cuh"
#include "tri_contact.cuh"
#include "testcode.cuh"


using namespace std;

#define THREAD 1024
#define COLLISION_PAIR 50

int main()
{
#ifdef TEST
	double a,b,c;
	a = 4;
	b = 2;
	c = 5;
	minmaxtest<<<1,1>>>();
	return 0;
#endif

	const std::string objfile_path = "../collision detection/flag-2000-changed.obj";

 	std::vector<vec3f> vertices;   //all vertices
 	std::vector<Object> triangles;  //all objects
 	CentroidValue alltri;

 	/*****************/
     //read obj file
     /*****************/

	clock_t start_objread = clock();

 	objread(objfile_path, vertices, triangles, alltri);

	clock_t end_objread = clock();

	cout<<"objread time: "<<(double)(end_objread - start_objread)<<"ms"<<endl<<endl;

 #ifdef TEST_CUDA_DEVICE_MEMORY
 	cout << "vertice:" << endl;
 	cout << vertices[0].x << " " << vertices[0].y << " " << vertices[0].z << endl;

 	cout << "triangle:" << endl;
 	cout << triangles[0].cx << " " << triangles[0].cy << " " << triangles[0].cz << endl;

 	cout << "centroidvalue:" << endl;
 	cout << alltri.max_x << " " << triangles[250].cy << " " << triangles[250].cz << endl;

 #endif

 	/*define device name*/
 	vec3f* d_vertices;
 	Object* d_triangles;
 	CentroidValue* d_alltri;

 	//apply for device memory
 	HANDLE_ERROR( cudaMalloc( (void**)&d_vertices, vertices.size() * sizeof(vec3f) ) );
 	HANDLE_ERROR( cudaMalloc( (void**)&d_triangles, triangles.size() * sizeof(Object) ) );
 	HANDLE_ERROR( cudaMalloc( (void**)&d_alltri, sizeof(CentroidValue) ) );

 	//assign from host to device
 	HANDLE_ERROR( cudaMemcpy( d_vertices, &vertices[0], vertices.size() * sizeof(vec3f), cudaMemcpyHostToDevice) );
 	HANDLE_ERROR( cudaMemcpy( d_triangles, &triangles[0], triangles.size() * sizeof(Object), cudaMemcpyHostToDevice) );
 	HANDLE_ERROR( cudaMemcpy( d_alltri, &alltri, sizeof(CentroidValue), cudaMemcpyHostToDevice) );
	
 	/*****************/
    //calculate montor code, and sort the objects with morton code 
    /*****************/

 	unsigned int numObjects = (unsigned int)triangles.size();
	
	clock_t start_mortoncode = clock();

 	MortonCode<<<(numObjects + THREAD - 1)/THREAD, THREAD, THREAD>>>(d_triangles, *d_alltri, numObjects);

	cudaThreadSynchronize();

	clock_t end_mortoncode = clock();

	cout<<"mortoncode calculate time: "<<(double)(end_mortoncode - start_mortoncode)<<"ms"<<endl<<endl;


 #ifdef TEST_TRIANGLE_MONTORCODE_CAL
 	printf("morton \n");

	kernel<<<1, 1>>>(d_triangles);
	cudaThreadSynchronize();

	printf("morton \n");
 #endif

	/*sort the triangles with thrust::sort(a host function)*/
	clock_t start_sortObjects = clock();

	thrust::device_ptr<Object> dev_ptr = thrust::device_pointer_cast(d_triangles);

    thrust::sort(dev_ptr, dev_ptr + numObjects);

#ifdef TEST_TRIANGLE_MONTORCODE_SORTED
	 printf("kernel \n");

	 kernel<<<1, 1>>>(d_triangles);
	 cudaThreadSynchronize();

	 printf("kernel end \n");
#endif

	/*assign vector d_MortonCodes with d_triangles*/

	mortonCode* d_MortonCodes;
	HANDLE_ERROR( cudaMalloc( (void**)&d_MortonCodes, triangles.size() * sizeof(mortonCode) ) );
	AssignMortonWithObjects<<<(numObjects + THREAD - 1)/THREAD, THREAD>>>(d_triangles, d_MortonCodes, numObjects);
	
	cudaThreadSynchronize();

	clock_t end_sortObjects = clock();

	cout<<"sort objects time: "<<(double)(end_sortObjects - start_sortObjects)<<"ms"<<endl<<endl;

#ifdef TEST_MONTORCODE_SORTED
	 printf("morton code \n");

	 kernel_morton<<<1, 1>>>(d_MortonCodes);
	 cudaThreadSynchronize();

	 printf("morton code end \n");
#endif

	/*****************/
    //generate BVH tree
    /*****************/

	Node* d_LeafNodes = NULL;
	Node* d_InternalNodes = NULL;

	HANDLE_ERROR( cudaMalloc( (void**)&d_LeafNodes, numObjects * sizeof(Node) ) );
	HANDLE_ERROR( cudaMalloc( (void**)&d_InternalNodes, (numObjects - 1)  * sizeof(Node) ) );

	clock_t start_generatebvh = clock();

	InitNodes<<<(numObjects + THREAD - 1)/THREAD, THREAD>>>(d_LeafNodes, numObjects, true);
	InitNodes<<<(numObjects + THREAD - 1)/THREAD, THREAD>>>(d_InternalNodes, numObjects - 1, false);
	
	cudaThreadSynchronize();

	AssignIndexOfLeafNodes<<<(numObjects + THREAD - 1)/THREAD, THREAD>>>(d_LeafNodes, d_triangles, numObjects);

	cudaThreadSynchronize();

	GenerateBVHGPU<<<(numObjects + THREAD - 1)/THREAD, THREAD>>>(d_triangles, d_MortonCodes, d_LeafNodes, d_InternalNodes, numObjects);
	
	cudaThreadSynchronize();

	clock_t end_generatebvh = clock();

	cout<<"generate bvh time: "<<(double)(end_generatebvh - start_generatebvh)<<"ms"<<endl<<endl;


#ifdef TEST_NODE
	printf("Node \n");

	kernel_node<<<1, 1>>>(d_InternalNodes);
	cudaThreadSynchronize();

	printf("Node end \n");
#endif

	/*****************/
	//Calculate BVH box
	/*****************/
	clock_t start_calboundingbox = clock();

	CalAllLeafBoundingBox<<<(numObjects + THREAD - 1)/THREAD, THREAD>>>(d_LeafNodes, d_triangles, d_vertices, numObjects);
	
	cudaThreadSynchronize();

	CalculateBVHBoundingBox_Violence<<<(numObjects + THREAD - 1)/THREAD, THREAD>>>(d_LeafNodes, numObjects);
	
	cudaThreadSynchronize();

	CalRoot<<<1,1>>>(d_InternalNodes);
	cudaThreadSynchronize();

	clock_t end_calboundingbox = clock();

	cout<<"calculate bounding box time: "<<(double)(end_calboundingbox - start_calboundingbox)<<"ms"<<endl<<endl;

#ifdef TEST_NODE_BOX
	printf("Node box\n");

	//kernel_nodebox<<<1, 1>>>(d_LeafNodes);
	kernel_nodebox<<<1, 1>>>(d_InternalNodes, numObjects);
	cudaThreadSynchronize();

	printf("Node box end \n");
#endif

	/*****************/
	//Traverse BVH
	/*****************/
	
	//CollisionList potential_list;
	int2* collision_list = new int2[COLLISION_PAIR];
	int2* d_collision_list;

	int* d_count;
	int* count = new int;

	*count = 0;

	

	HANDLE_ERROR(cudaMalloc( (void**)&d_collision_list, COLLISION_PAIR * sizeof(int2)) );
	HANDLE_ERROR(cudaMalloc( (void**)&d_count, sizeof(int)) );

	HANDLE_ERROR(cudaMemcpy(d_count, count, sizeof(int), cudaMemcpyHostToDevice));

	clock_t start_findpotential = clock();

	FindPotentialCollisions<<<256, 256>>>(d_InternalNodes, 
											d_LeafNodes, 
											d_collision_list,
											d_vertices,
											d_triangles,  
											numObjects,
											d_count);
	
	cudaThreadSynchronize();

	clock_t end_findpotential = clock();

	cout<<"find potential collision time: "<<(double)(end_findpotential- start_findpotential)<<"ms"<<endl<<endl;

	cout<<"collision process of GPU: "<<(double)(end_findpotential- start_mortoncode)<<"ms"<<endl<<endl;
#ifdef TEST_COLLISION

	printf("collision\n");

	//kernel_nodebox<<<1, 1>>>(d_LeafNodes);
	kernel_count<<<1, 1>>>(d_count);
	cudaThreadSynchronize();
	kernel_collision<<<1, 1>>>(d_collision_list, 20);
	cudaThreadSynchronize();

	printf("collision end \n");
#endif

	/*****************/
	//verify and printf
	/*****************/

	cout<<"************************************"<<endl;
	cout<<"****************GPU*****************"<<endl;
	cout<<"************************************"<<endl;

	HANDLE_ERROR(cudaMemcpy(collision_list,d_collision_list, COLLISION_PAIR * sizeof(int2),cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(count, d_count, sizeof(int), 	cudaMemcpyDeviceToHost));
	
	for(int i = 0; i < count[0]; i++)
	{
		cout<<"#self contact found at ("<<collision_list[i].x - 1<<", "<<collision_list[i].y - 1<<")"<<endl;
	}	
	
	std::vector<int> oneshot;
	for(int i = 0; i < count[0]; i++)
	{
		oneshot.push_back(collision_list[i].x);
		oneshot.push_back(collision_list[i].y);
	}
	sort(oneshot.begin(), oneshot.end());
    oneshot.erase(unique(oneshot.begin(), oneshot.end()), oneshot.end());
	
	cout<<"totally "<<count[0]<< "colliding pairs ..."<<endl;
	for(int x : oneshot)
        cout << x - 1 << endl;
	/*****************/
    //Destroy the device memory
    /*****************/
	cudaFree(d_vertices);
	cudaFree(d_triangles);
	cudaFree(d_alltri);
	cudaFree(d_MortonCodes);
	cudaFree(d_LeafNodes);
	cudaFree(d_InternalNodes);

	return 0;	
}