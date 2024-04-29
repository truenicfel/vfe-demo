#pragma once

#include <Eigen/Eigen>

#include <deque>
#include <array>
#include <iostream>

// A double buffered line that stores two instances of a line where one is for
// reading while the other one handles write requests. Expanding the lines is
// done on both instances.
template<int Dimensions>
class DoubleBufferedLine
{
public:
    using Vertex = Eigen::Vector<double, Dimensions>;
    using Line = std::deque<Vertex>;
    using IndexRange = std::pair<size_t, size_t>;

    // create a new line with initial seed given
    explicit DoubleBufferedLine(const Vertex& seed)
    : mWriteIndex(0)
    , mLines()
    {
        // initialize with seed
        mLines[mWriteIndex].push_back(seed);
        mLines[1 - mWriteIndex].push_back(seed);
    }

    // get size
    size_t GetSize() const {
        return mLines[1 - mWriteIndex].size();
    }

    // swaps the reading and writing lines
    // swapping is not required after prepend/append
    void Swap() {
        mWriteIndex = 1 - mWriteIndex;
    }

    // add vertex in the front (both lines)
    void Prepend(const Vertex& value) {
        GetReadLine().push_front(value);
        GetWriteLine().push_front(value);
    }

    // add vertex in the back (both lines)
    void Append(const Vertex& value) {
        GetReadLine().push_back(value);
        GetWriteLine().push_back(value);
    }

    // returns true if the line is closed
    bool IsClosed(const double& epsilon) const {
        if (GetSize() < 3)
        {
            return false;
        }
        return (Read(0) - Read(GetSize() - 1)).norm() < epsilon;
    }

    // remove the last vertex (both lines)
    void RemoveLast() {
        GetReadLine().pop_back();
        GetWriteLine().pop_back();
    }

    // remove the first vertex (both lines)
    void RemoveFirst() {
        GetReadLine().pop_front();
        GetWriteLine().pop_front();
    }

    // write value at index
    void Write(const size_t& index, const Vertex& vertex) {
        GetWriteLine().at(index) = vertex;
    }

    // read value at index
    const Vertex& Read(const size_t& index) const {
        return GetReadLine().at(index);
    }

    // get last n indices (when append, prepend, remove* are called this index range becomes invalid!)
    // first index is inclusive, second is exclusive
    IndexRange GetLastN(const size_t& n) const {
        IndexRange result;
        if (n <= GetSize()) {
            result =  std::make_pair(GetSize() - n, GetSize());
        }
        else {
            result =  std::make_pair(0, GetSize());
        }
        return result;
    }

    // get first n indices (when append, prepend, remove* are called this index range becomes invalid!)
    // first index is inclusive, second is exclusive
    IndexRange GetFirstN(const size_t& n) const {
        return std::make_pair(0, std::min(n, GetSize()));
    }

//    void Print() const
//    {
//        std::cout << "Size: " << GetSize() << std::endl;
//        std::cout << "Active: " << std::endl << "[";
//        for (int index = 0; index < GetSize(); ++index)
//        {
//            auto value = mLines[mActive]->GetVertices()->GetValue(mLeft + index);
//            std::cout << "(" << value.x() << ", " << value.y() << ", " << value.z() << "), ";
//        }
//        std::cout << "]" << std::endl;
//
//        std::cout << "Inactive: " << std::endl << "[";
//        for (int index = 0; index < GetSize(); ++index)
//        {
//            auto value = mLines[1- mActive]->GetVertices()->GetValue(mLeft + index);
//            std::cout << "(" << value.x() << ", " << value.y() << ", " << value.z() << "), ";
//        }
//        std::cout << "]" << std::endl;
//    }

    Line GetCopy() const {
        return GetReadLine();
    }

private:

    // index of the line currently in use for writing
    size_t mWriteIndex;

    // stores the two lines
    std::array<Line, 2> mLines;

    Line& GetWriteLine() {
        return mLines[mWriteIndex];
    }

    Line& GetReadLine() {
        return mLines[1 - mWriteIndex];
    }

    const Line& GetReadLine() const {
        return mLines[1 - mWriteIndex];
    }

};

using DoubleBufferedLine3 = DoubleBufferedLine<3>;
using DoubleBufferedLine2 = DoubleBufferedLine<2>;

