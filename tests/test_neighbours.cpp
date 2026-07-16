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

    bool* neighbourScratch = new bool[tree->getParticleCount()];
    int *stackScratch = new int[tree->getNodeCount()];

    GPU::getNeighbours(tree->data->contents, tree->data->mapping, tree->getNodeCount(), tree->getParticleCount(), nodeIndex, 208, data.getOffloadConfig(), data.h(208), neighbourScratch, stackScratch);

    std::vector<int> neighbourList = {};
    for (int i = 0; i < tree->getParticleCount(); i++) {
        if (neighbourScratch[i]) {
            neighbourList.push_back(i);
        }
    }

    ASSERT_THAT(neighbourList, testing::ElementsAre(205, 206, 207, 208, 209));
}
