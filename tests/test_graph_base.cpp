#include <gtest/gtest.h>
#include "graph_base.h"

using namespace polygon_coverage_planning;

class SimpleGraph : public GraphBase<int, int> {
 public:
  bool addEdges() override { return true; }
  bool create() override { return true; }
  bool addEdgePublic(const EdgeId& id, const int& prop, double cost) {
    return addEdge(id, prop, cost);
  }
};

TEST(GraphBaseTest, GetEdgeProperty) {
  SimpleGraph g;
  EXPECT_TRUE(g.addNode(1));
  EXPECT_TRUE(g.addNode(2));
  EXPECT_TRUE(g.addEdgePublic({0,1}, 42, 1.0));
  const int* prop = g.getEdgeProperty({0,1});
  ASSERT_NE(prop, nullptr);
  EXPECT_EQ(*prop, 42);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
