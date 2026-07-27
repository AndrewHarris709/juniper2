#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "SimData.h"
#include "GPUOperators.h"

void assertAtDepth(Tree* tree, int nodeIndex, int depth, int partCountAtDepth) {
    if (depth > 0) {
        assertAtDepth(tree, tree->getLeftChild(nodeIndex), depth - 1, partCountAtDepth);
        assertAtDepth(tree, tree->getRightChild(nodeIndex), depth - 1, partCountAtDepth);
    } else {
        ASSERT_EQ(tree->getPartsFromNode(nodeIndex).size, partCountAtDepth);
        ASSERT_TRUE(tree->isLeaf(nodeIndex));
    }
}

TEST(NeighbourTest, TreeBuildTest) {
    SimData sim = SimData(".//files//kd_test.csv", ".//files//kd_test.toml");
    Tree* tree = new Tree(sim.getParticleCount());
    tree->build(sim);

    assertAtDepth(tree,0, 6, 5);
}

TEST(NeighbourTest, NeighbourFindTest) {
    SimData data = SimData(".//files//kd_test.csv", ".//files//kd_test.toml");
    Tree* tree = new Tree(data.getParticleCount());
    tree->build(data);

    int nodeIndex = 0;
    bool found = false;
    while (!found) {
        nodeIndex++;
        if (!tree->isLeaf(nodeIndex)) {
            continue;
        }

        NodeRange range = tree->getPartsFromNode(nodeIndex);
        for (int i = 0; i < range.size; i++) {
            if (range.data[i] == 208) {
                found = true;
            }
        }
    }

    GPU::NeighbourList neighbours;
    GPU::getNeighbours(tree->data, nodeIndex, data.getOffloadConfig(), data.h(208), &neighbours);

    std::vector<int> neighbourVector(neighbours.indices, neighbours.indices + neighbours.count);

    ASSERT_THAT(neighbourVector, testing::ElementsAre(209, 208, 207, 206, 205));
}
