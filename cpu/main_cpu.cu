#include <cuda_runtime.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>
#include  <algorithm>

#include "objread.h"
#include "morton.h"
#include "object.cuh"
//#include "../common/book.h"
#include "tri_contact.h"
#include "Node.cuh"
#include "bvh.cuh"
#include "calbox.cuh"
#include "traverse.cuh"
#include <ctime>

using namespace std;

int main()
{
#ifdef TEST
	unsigned long long i = 1;
	i <<= 62;
	cout<<myclz(i)<<endl;
	return 0;
#endif
	const std::string objfile_path = "../collision detection/flag-2000-changed.obj";

	std::vector<vec3f> vertices;   //all vertices
	std::vector<Object> triangles;  //all objects
    std::vector<mortonCode> MortonCodes; //specificly store montor codes
    
	CentroidValue alltri;
	
	/*****************/
   //read obj file
   /*****************/
    clock_t start_objread = clock();

	objread(objfile_path, vertices, triangles, alltri);

	clock_t end_objread = clock();

	cout<<"objread time: "<<(double)(end_objread - start_objread)<<"ms"<<endl<<endl;

	/*****************/
    //计算montor code
    /*****************/
   	clock_t start_mortoncode = clock();

	MortonCode(triangles, alltri);

	clock_t end_mortoncode = clock();

	cout<<"mortoncode calculate time: "<<(double)(end_mortoncode - start_mortoncode)<<"ms"<<endl<<endl;

	/*****************/
	//sort objects with the morton code
	/*****************/
	clock_t start_sortObjects = clock();

	SortObjectsWithMorton(triangles, MortonCodes);

	clock_t end_sortObjects = clock();

	cout<<"sort objects time: "<<(double)(end_sortObjects - start_sortObjects)<<"ms"<<endl<<endl;

	/*****************/
	//generate BVH
	/*****************/

	unsigned int numObjects = (unsigned int)triangles.size();

	Node* LeafNodes = NULL;
	Node* InternalNodes = NULL;
	
	CreateNodes(LeafNodes, numObjects);
	InitNodes(LeafNodes, numObjects, true);
	CreateNodes(InternalNodes, numObjects - 1);
	InitNodes(InternalNodes, numObjects - 1, false);

	clock_t start_generatebvh = clock();
	
	GenerateBVH(triangles, MortonCodes, LeafNodes, InternalNodes, numObjects);

	clock_t end_generatebvh = clock();

	cout<<"generate bvh time: "<<(double)(end_generatebvh - start_generatebvh)<<"ms"<<endl<<endl;

#ifdef TEST_NODE
	Node* p = InternalNodes;
	int i = 0;
	while(p->LeftChild){
		i++;
		p = p->LeftChild;
	}
	cout<<"left depth: "<<i<<endl;

	i = 0;

	p = InternalNodes;
	while(p->RightChild){
		i++;
		p = p->RightChild;
	}
	cout<<"right depth: "<<i<<endl;

	return 0;

	
#endif
	/*****************/
	//Calculate BVH box
	/*****************/
	clock_t start_calboundingbox = clock();
	//calculate all objects' bounding box
	CalAllLeafBoundingBox(LeafNodes, numObjects, triangles, vertices);

	//calculate all internal nodes' bounding box
	CalculateBVHBoundingBox_Violence(LeafNodes, numObjects);

	clock_t end_calboundingbox = clock();

	cout<<"calculate bounding box time: "<<(double)(end_calboundingbox - start_calboundingbox)<<"ms"<<endl<<endl;


#ifdef TEST_NODE_BOX
	// int all = 0;
	// for (unsigned int i = 0; i < numObjects-1; i++) {
	// 	if (!InternalNodes[i].boxCaled) {
	// 		cout << "no calculated internal node " << i << endl;
	// 		all++;
	// 	}
	// }
	// cout << "all " << all << endl;
	cout<<"root box:"<<endl;
	cout<<"rightup: "<<InternalNodes->bounding.RightUp.x<<" "<<InternalNodes->bounding.RightUp.y<<" "<<InternalNodes->bounding.RightUp.z<<endl;
	cout<<"leftdown: "<<InternalNodes->bounding.LeftDown.x<<" "<<InternalNodes->bounding.LeftDown.y<<" "<<InternalNodes->bounding.LeftDown.z<<endl;
	return 0;
#endif

	/*****************/
	//Traverse BVH
	/*****************/
	
	CollisionList potential_list;
	std::vector<unsigned int> collision_list;

	clock_t start_findpotential = clock();

	//find potential collisioned box pair
	FindPotentialCollisions(InternalNodes, LeafNodes, potential_list, numObjects);

	clock_t end_findpotential = clock();
#ifdef TEST_POTENTIAL_PAIR
	cout<<"potential pair"<<endl;
	cout << potential_list.size() << endl;
	return 0;
	for(int i = 0; i < potential_list.size(); i++)
	{
		cout<<potential_list[i].x<<" "<<potential_list[i].y<<endl;
	}
#endif

	/*****************/
	//verify and printf
	/*****************/

	

	CollisionList potential_pair;
	//Verify Object Collision
	VerifyCollisions(collision_list,potential_pair, potential_list, vertices, triangles);
	
	cout<<"find potential collision time: "<<(double)(end_findpotential- start_findpotential)<<"ms"<<endl<<endl;


	cout<<"collision process of CPU: "<<(double)(end_findpotential- start_mortoncode)<<"ms"<<endl<<endl;

	cout<<"************************************"<<endl;
	cout<<"****************CPU*****************"<<endl;
	cout<<"************************************"<<endl;


	for(int i = 0; i < potential_pair.size(); i++)
	{
		cout<<"#self contact found at ("<<potential_pair[i].x - 1<<", "<<potential_pair[i].y - 1<<")"<<endl;
	}
	cout<<"totally "<<collision_list.size()<< "colliding pairs ..."<<endl;


	for(int x : collision_list)
        cout << x - 1 << endl;

	/*****************/
	//Free memory
	/*****************/
	delete []LeafNodes;
	delete []InternalNodes;



	return 0;
}