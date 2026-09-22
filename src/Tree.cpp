//
// Created by Andrew on 2025-09-12.
//

#include <vector>
#include <iostream>
#include <stack>

#include "commons.h"
#include "ParticleOperators.cpp"
#include "Tree.h"

int Tree::getLeftChild(int nodeIndex) {
    return this->getNode(nodeIndex)->leftChild;
}

int Tree::getRightChild(int nodeIndex) {
    return this->getNode(nodeIndex)->rightChild;
}

bool Tree::isLeaf(int index) {
    return this->getNode(index)->leftChild == -1;
}

int Tree::getParent(int nodeIndex) {
    return this->getNode(nodeIndex)->parent;
}

int Tree::getParticleCount() {
    return this->data->partCount;
}

int Tree::getNodeCount() {
    return this->data->nodeCount;
}

TreeNode* Tree::getNode(int nodeIndex) {
    return &this->data->contents[nodeIndex];
}

float Tree::distBetweenNodes(SimConfig& config, int nodeIndex1, int nodeIndex2) {
    TreeNode* node1 = this->getNode(nodeIndex1);
    TreeNode* node2 = this->getNode(nodeIndex2);

    return junipermath::distBetween(config.getBoundingBox(), node1->x, node2->x, node1->y, node2->y, node1->z, node2->z);
}

NodeRange Tree::getPartsFromNode(int nodeIndex) {
    TreeNode *node = getNode(nodeIndex);

    NodeRange range;
    range.data = this->data->mapping + node->mappingStart;
    range.size = node->mappingSize;

    return range;
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
        }
    }
}

void Tree::splitLeaf(SimData& data, int nodeIndex) {
    if (!this->isLeaf(nodeIndex)) {
        std::cout << "Warning: attempted splitLeaf() on a non-leaf node!" << std::endl;
        return;
    }

    NodeRange indices = this->getPartsFromNode(nodeIndex);
    junipermath::Box3f box = PartOps::getBoundingBox(data, *this, nodeIndex);
    float xdiff = box.x2 - box.x1;
    float ydiff = box.y2 - box.y1;
    float zdiff = box.z2 - box.z1;

    TreeNode* thisNode = getNode(nodeIndex);
    float centre = thisNode->x;
    int decomposeAxis = 0;
    if (ydiff >= xdiff && ydiff >= zdiff) {
        decomposeAxis = 1;
        centre = thisNode->y;
    } else if (zdiff >= xdiff && zdiff >= ydiff) {
        decomposeAxis = 2;
        centre = thisNode->z;
    }

    int frontInsert = 0, backInsert = indices.size - 1;
    std::vector<int> oldIndices(indices.data, indices.data + indices.size);
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
    this->data->contents[leftIndex].parent = nodeIndex;
    int rightIndex = registerNodeFromIndices(data, thisNode->mappingStart + frontInsert, indices.size - frontInsert);
    this->data->contents[rightIndex].parent = nodeIndex;

    thisNode->leftChild = leftIndex;
    thisNode->rightChild = rightIndex;
}

int Tree::registerNodeFromIndices(SimData& data, int index, int size) {
    TreeNode* node = getNode(this->getNodeCount());

    node->mappingStart = index;
    node->mappingSize = size;

    junipermath::Point3f centreOfMass = PartOps::getCentreOfMass(data, *this, this->getNodeCount());
    node->x = centreOfMass.x;
    node->y = centreOfMass.y;
    node->z = centreOfMass.z;

    node->size = PartOps::getBoundingRadius(data, *this, this->getNodeCount());
    node->hmax = PartOps::getMaxSmoothingLength(data, *this, this->getNodeCount());

    this->data->nodeCount++;
    return this->data->nodeCount - 1;
}
