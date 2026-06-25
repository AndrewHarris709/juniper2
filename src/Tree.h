//
// Created by Andrew on 2025-09-12.
//

#ifndef JUNIPEREXE_TREENODE_H
#define JUNIPEREXE_TREENODE_H

#include <span>

#include "SimData.h"

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

class Tree {
    TreeNode* contents;
    int* mapping;
    int nodeCount = 0;
    int partCount = 0;

    int registerNodeFromIndices(SimData& data, int index, int size);
public:
    explicit Tree(int count) {
        contents = new TreeNode[2 * count];
        mapping = new int[count];

        for (int i = 0; i < count; i++) {
            mapping[i] = i;
        }
        this->partCount = count;
    }

    ~Tree() {
        delete[] contents;
        delete[] mapping;
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
