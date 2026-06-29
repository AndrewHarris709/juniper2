//
// Created by Andrew on 2025-09-12.
//

#ifndef JUNIPEREXE_TREENODE_H
#define JUNIPEREXE_TREENODE_H

#include <span>

#include "SimData.h"

#define TREE_DATA_MAP(td) \
map(to: (td).contents[0:(td).nodeCount]) \
map(to: (td).mapping[0:(td).partCount])

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

    TreeNode* getNode(int index);
    int getParent(int index);
    int getLeftChild(int index);
    int getRightChild(int index);
    bool isLeaf(int index);
    void splitLeaf(SimData& data, int index);
    int getParticleCount();
    int getNodeCount();
    std::span<int> getNodeIndices(int index);
    void build(SimData& data);
};

#endif //JUNIPEREXE_TREENODE_H
