//
// Created by Andrew on 2025-09-12.
//

#include <vector>
#include <iostream>
#include "Tree.h"

#include <stack>

#include "commons.h"
#include "ParticleOperators.cpp"

int Tree::getLeftChild(int index) {
    return this->getNode(index)->leftChild;
}

int Tree::getRightChild(int index) {
    return this->getNode(index)->rightChild;
}

bool Tree::isLeaf(int index) {
    return this->getNode(index)->leftChild == -1;
}

int Tree::getParent(int index) {
    return this->getNode(index)->parent;
}

int Tree::getParticleCount() {
    return this->data->partCount;
}

int Tree::getNodeCount() {
    return this->data->nodeCount;
}

TreeNode* Tree::getNode(int index) {
    return &this->data->contents[index];
}

std::span<int> Tree::getNodeIndices(int index) {
    TreeNode *node = getNode(index);
    int mappingStart = node->mappingStart, mappingSize = node->mappingSize;

    std::span<int> indices(this->data->mapping + mappingStart, this->data->mapping + mappingStart + mappingSize);
    return indices;
}

void Tree::build(SimData& data) {
    std::stack<int> stack;
    stack.push(0);
    registerNodeFromIndices(data, 0, this->data->partCount);

    while (!stack.empty()) {
        int newIndex = stack.top();
        TreeNode* node = this->getNode(newIndex);
        stack.pop();
        if (node->mappingSize >= 10) {
            this->splitLeaf(data, newIndex);
            stack.push(node->leftChild);
            stack.push(node->rightChild);
            std::cout << "Node " << newIndex << " separated into subnodes " << node->leftChild << " and " << node->rightChild << std::endl;
            std::cout << node->mappingSize << " particles separated into groups of " << getNode(node->leftChild)->mappingSize << " and " << getNode(node->rightChild)->mappingSize << std::endl;
        }
    }
}

void Tree::splitLeaf(SimData& data, int index) {
    if (!this->isLeaf(index)) {
        std::cout << "Warning: attempted splitLeaf() on a non-leaf node!" << std::endl;
        return;
    }

    std::span<int> indices = this->getNodeIndices(index);
    Box3f box = PartOps::getBoundingBox(data, *this, index);
    float xdiff = box.x2 - box.x1;
    float ydiff = box.y2 - box.y1;
    float zdiff = box.z2 - box.z1;

    TreeNode* thisNode = getNode(index);
    float centre = thisNode->x;
    int decomposeAxis = 0;
    if (ydiff >= xdiff && ydiff >= zdiff) {
        decomposeAxis = 1;
        centre = thisNode->y;
    } else if (zdiff >= xdiff && zdiff >= ydiff) {
        decomposeAxis = 2;
        centre = thisNode->z;
    }

    int frontInsert = 0, backInsert = indices.size() - 1;
    std::vector<int> oldIndices(indices.begin(), indices.end());
    for (int i : oldIndices) {
        if (data.xyzh[4 * i + decomposeAxis] > centre) {
            this->data->mapping[thisNode->mappingStart + frontInsert] = i;
            frontInsert++;
        } else {
            this->data->mapping[thisNode->mappingStart + backInsert] = i;
            backInsert--;
        }
    }

    int leftIndex = registerNodeFromIndices(data, thisNode->mappingStart, frontInsert);
    this->data->contents[leftIndex].parent = index;
    int rightIndex = registerNodeFromIndices(data, thisNode->mappingStart + frontInsert, indices.size() - frontInsert);
    this->data->contents[rightIndex].parent = index;

    thisNode->leftChild = leftIndex;
    thisNode->rightChild = rightIndex;
}

int Tree::registerNodeFromIndices(SimData& data, int index, int size) {
    TreeNode* node = getNode(this->getNodeCount());

    node->mappingStart = index;
    node->mappingSize = size;

    Point3f centreOfMass = PartOps::getCentreOfMass(data, *this, this->getNodeCount());
    node->x = centreOfMass.x;
    node->y = centreOfMass.y;
    node->z = centreOfMass.z;

    node->size = PartOps::getBoundingRadius(data, *this, this->getNodeCount());
    node->hmax = PartOps::getMaxSmoothingLength(data, *this, this->getNodeCount());

    this->data->nodeCount++;
    return this->data->nodeCount - 1;
}
