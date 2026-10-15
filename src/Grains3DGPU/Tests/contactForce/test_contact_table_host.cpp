#include <gtest/gtest.h>

#include "ContactTable.hh"
#include "GrainsUtils.hh"

// Host-side tests for ContactHashTable.
// find / findOrInsert / remove live on ContactMemoryView<T> (the lightweight proxy
// returned by getView()), so each test that needs lookup operations uses a view.
// The raw pointers in the view remain valid for the lifetime of the table; a single
// view obtained in SetUp() stays usable after clear() and resetIndexCounter().

class ContactTableHostTest : public ::testing::Test
{
protected:
    ContactHashTable<float, MemType::HOST> table;
    ContactMemoryView<float>               view;

    void SetUp() override
    {
        table.allocate(100, 100);
        view = table.getView();
    }

    void TearDown() override
    {
        table.deallocate();
    }
};

TEST_F(ContactTableHostTest, AllocationAndDeallocation)
{
    ContactHashTable<float, MemType::HOST> t;
    EXPECT_EQ(t.getCapacity(), 0);
    EXPECT_EQ(t.getNextIndex(), 0);
    t.allocate(100, 100);
    EXPECT_EQ(t.getCapacity(), 100);
    EXPECT_EQ(t.getNextIndex(), 0);
    t.deallocate();
    EXPECT_EQ(t.getCapacity(), 0);
    EXPECT_EQ(t.getNextIndex(), 0);
}

TEST_F(ContactTableHostTest, ConstructorWithCapacity)
{
    ContactHashTable<float, MemType::HOST> t(50, 50);
    EXPECT_EQ(t.getCapacity(), 50);
    EXPECT_EQ(t.getNextIndex(), 0);
}

TEST_F(ContactTableHostTest, BasicInsertion)
{
    uint2 pair1 = make_uint2(5, 10);
    uint  index1;
    bool  success = view.findOrInsert(pair1, index1);
    EXPECT_TRUE(success);
    EXPECT_EQ(index1, 0u);
    EXPECT_EQ(table.getNextIndex(), 1u);
}

TEST_F(ContactTableHostTest, FindExisting)
{
    uint2 pair1 = make_uint2(5, 10);
    uint  index1, index2;
    view.findOrInsert(pair1, index1);
    bool found = view.find(pair1, index2);
    EXPECT_TRUE(found);
    EXPECT_EQ(index1, index2);
}

TEST_F(ContactTableHostTest, FindNonExistent)
{
    uint2 pair1 = make_uint2(5, 10);
    uint2 pair2 = make_uint2(7, 12);
    uint  index;
    view.findOrInsert(pair1, index);
    bool found = view.find(pair2, index);
    EXPECT_FALSE(found);
}

TEST_F(ContactTableHostTest, MultipleInsertions)
{
    const int numPairs = 10;
    uint      indices[numPairs];
    for(int i = 0; i < numPairs; ++i)
    {
        uint2 pair    = make_uint2(i, i + 10);
        bool  success = view.findOrInsert(pair, indices[i]);
        EXPECT_TRUE(success);
        EXPECT_EQ(indices[i], (uint)i);
    }
    EXPECT_EQ(table.getNextIndex(), (uint)numPairs);
}

TEST_F(ContactTableHostTest, DuplicateInsertion)
{
    uint2 pair = make_uint2(5, 10);
    uint  index1, index2;
    view.findOrInsert(pair, index1);
    view.findOrInsert(pair, index2);
    EXPECT_EQ(index1, index2);
    EXPECT_EQ(table.getNextIndex(), 1u);
}

TEST_F(ContactTableHostTest, Removal)
{
    uint2 pair = make_uint2(5, 10);
    uint  index;
    view.findOrInsert(pair, index);
    bool removed = view.remove(pair);
    EXPECT_TRUE(removed);
    bool found = view.find(pair, index);
    EXPECT_FALSE(found);
}

TEST_F(ContactTableHostTest, RemoveNonExistent)
{
    uint2 pair    = make_uint2(5, 10);
    bool  removed = view.remove(pair);
    EXPECT_FALSE(removed);
}

TEST_F(ContactTableHostTest, ClearOperation)
{
    for(int i = 0; i < 5; ++i)
    {
        uint2 pair = make_uint2(i, i + 10);
        uint  index;
        view.findOrInsert(pair, index);
    }
    EXPECT_EQ(table.getNextIndex(), 5u);
    table.clear();
    EXPECT_EQ(table.getNextIndex(), 0u);
    uint index;
    for(int i = 0; i < 5; ++i)
    {
        uint2 pair  = make_uint2(i, i + 10);
        bool  found = view.find(pair, index);
        EXPECT_FALSE(found);
    }
}

TEST_F(ContactTableHostTest, ResetIndexCounter)
{
    uint2 pair = make_uint2(5, 10);
    uint  index;
    view.findOrInsert(pair, index);
    EXPECT_EQ(table.getNextIndex(), 1u);
    table.resetIndexCounter();
    EXPECT_EQ(table.getNextIndex(), 0u);
}

TEST_F(ContactTableHostTest, FullTable)
{
    const int                              capacity = 10;
    ContactHashTable<float, MemType::HOST> t(capacity, capacity);
    ContactMemoryView<float>               localView = t.getView();
    for(int i = 0; i < capacity; ++i)
    {
        uint2 pair = make_uint2(i, i + 100);
        uint  index;
        bool  success = localView.findOrInsert(pair, index);
        EXPECT_TRUE(success);
    }
}

TEST_F(ContactTableHostTest, HashCollisions)
{
    uint indices[5];
    for(int i = 0; i < 5; ++i)
    {
        uint2 pair    = make_uint2(i * 1000, i * 1000 + 1);
        bool  success = view.findOrInsert(pair, indices[i]);
        EXPECT_TRUE(success);
    }
    for(int i = 0; i < 5; ++i)
    {
        for(int j = i + 1; j < 5; ++j)
            EXPECT_NE(indices[i], indices[j]);
    }
}

TEST_F(ContactTableHostTest, PairOrdering)
{
    ContactHashTable<float, MemType::HOST> t(100, 100);
    ContactMemoryView<float>               localView = t.getView();
    uint2                                  pair1     = make_uint2(5, 10);
    uint2                                  pair2     = make_uint2(10, 5);
    uint                                   index1, index2;
    localView.findOrInsert(pair1, index1);
    localView.findOrInsert(pair2, index2);
    EXPECT_NE(index1, index2);
}
