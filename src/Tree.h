//
// Created by Andrew on 2025-09-12.
//

#ifndef JUNIPEREXE_TREENODE_H
#define JUNIPEREXE_TREENODE_H

#include <span>

class SimData; // forward declare to avoid circular imports with SimData.h

struct TreeNode {
    float x;
    float y;
    float z;
    float size;
    float hmax;
    int mappingStart;
    int mappingSize;
    int leftChild = -1;
    int rightChild = -1;
    int parent;
};

struct TreeData {
    TreeNode* contents = nullptr;
    int* mapping = nullptr;
    int nodeCount = 0;
    int partCount = 0;
};

struct NodeRange {
    int* data;
    int size;
};

class Tree {

    int registerNodeFromIndices(SimData& data, int index, int size);

public:
    TreeData* data;

    explicit Tree(const int count) {
        data = new TreeData();
        data->contents = new TreeNode[2 * count];
        data->mapping = new int[count];

        for (int i = 0; i < count; i++) {
            data->mapping[i] = i;
        }
        data->partCount = count;
    }

    ~Tree() {
        delete[] data->contents;
        delete[] data->mapping;
    }

    TreeNode* getNode(int nodeIndex);
    int getParent(int nodeIndex);
    int getLeftChild(int nodeIndex);
    int getRightChild(int nodeIndex);
    bool isLeaf(int nodeIndex);
    void splitLeaf(SimData& data, int nodeIndex);
    float distBetweenNodes(int nodeIndex1, int nodeIndex2);
    int getParticleCount();
    int getNodeCount();

    NodeRange getPartsFromNode(int nodeIndex);
    void build(SimData& data);
    void densityIterate(SimData& data);
};

#endif //JUNIPEREXE_TREENODE_H
