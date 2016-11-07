#include <gtest/gtest.h>

#include <OpenVolumeMesh/System/MemoryInclude.hh>

TEST(MakeUniqueTest, MakeUniqueTest) {
  std::unique_ptr<int> foo;
  auto bar = ptr::make_unique<int>(5);
  foo = std::move(bar);

  EXPECT_EQ(*foo, 5);
  EXPECT_EQ(bar.get(), nullptr);
}

