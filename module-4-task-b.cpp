#include <array>
#include <cassert>
#include <cstring>
#include <iostream>
#include <stack>
#include <vector>

/* We have field 8x8
 * So we have 8^2 * 8^2 * 2 = 8192 vertices in
 * out game graph
 */
enum FieldSize {
    X = 8,
    Y = 8
};

enum FieldObjects {
    WALL = 1,
    FUGITIVE = 2,
    TERMINATOR = 3
};

using comp_t  =  int8_t; //component type
using vert_t  = size_t;  //encoded vertex type

using coord_t = std::pair<comp_t , comp_t>;
using field_t = std::array<std::array<comp_t, FieldSize::X>, FieldSize::Y>;

constexpr vert_t GetBase() {
    static_assert(FieldSize::X == FieldSize::Y);
    return FieldSize::X;
}

constexpr vert_t GetGraphVSize() {
    return 2 * FieldSize::X * FieldSize::X * FieldSize::Y * FieldSize::Y;
}

using graph_vertices_t = std::vector<bool>;

struct Vertex {
    coord_t fugitiveCoord;
    coord_t terminatorCoord;
    bool fugitiveTurn;
};


Vertex ComplexCoordToVertex(uint32_t complexCoord) {
    bool fugitiveTurn = true;
    if(complexCoord % 2 != 0) {
        fugitiveTurn = false;
        --complexCoord;
    }

    uint32_t base = GetBase();
    uint32_t baseSq = base * base;
    uint32_t baseCube = baseSq * base;

    comp_t fugitiveX   = complexCoord % base;
    comp_t fugitiveY   = ((complexCoord - fugitiveX) / base) % baseSq;
    comp_t terminatorX = ((complexCoord - fugitiveX - fugitiveY * base) / baseSq) % base;
    comp_t terminatorY = (complexCoord - fugitiveX - fugitiveY * base - terminatorX * baseSq) / baseCube;

    return Vertex {
            {fugitiveX, fugitiveY},
            {terminatorX, terminatorY},
            fugitiveTurn
    };
}

uint32_t VertexToComplexCoord(const Vertex& vertex) {
    uint32_t base = GetBase();
    uint32_t result = 0;

    result += vertex.fugitiveCoord.first;
    result += vertex.fugitiveCoord.second * base;
    result += vertex.terminatorCoord.first * base * base;
    result += vertex.terminatorCoord.second * base * base * base;

    if(!vertex.fugitiveTurn) {
        ++result;
    }

    return result;
}

class Solution {
private:
    field_t field_{};

    /* @param vertices_winning: is the vertex winning position
    *  @param vertices_known: do we know about winning property of the vertex
    *  @param top_sort_used: is the certain vertex used
    */
    graph_vertices_t vertices_winning_;
    graph_vertices_t vertices_known_;

    graph_vertices_t top_sort_used_;
    std::stack<vert_t> top_sort_dfs_;

    void handleField();
    void runTopologicalSort();

    void topologicalDFS(vert_t currentVertex);
    std::array<coord_t, 8> possibleMoves_;

    std::vector<vert_t> getNeighbours(vert_t complexCoord) const;

    bool isUnderBlast(const Vertex& vertex);
    bool isInitiallyWinning(const Vertex& vertex);
    bool isInitiallyLoosing(const Vertex& vertex);

    bool isThereWallX(comp_t xComp, const Vertex &vertex) const;
    bool isThereWallY(comp_t yComp, const Vertex &vertex) const;
    bool isThereWallDiagonal(const Vertex& vertex) const;

    bool canGo(const Vertex &currentVertex, coord_t delta) const;
    static Vertex goDelta(const Vertex& currentVertex, coord_t delta) ;
public:
    explicit Solution(const field_t& field);
    bool FugitiveWins();
};

Solution::Solution(const field_t &field) : field_(field) {
    vertices_winning_.resize(GetGraphVSize(), false);
    vertices_known_.resize(GetGraphVSize(), false);
    top_sort_used_.resize(GetGraphVSize(), false);

    assert(vertices_known_.size() == GetGraphVSize());

    possibleMoves_ = {{
      {1, 0},
      {-1, 0},
      {0, 1},
      {0, -1},
      {1, 1},
      {1, -1},
      {-1, 1},
      {-1, -1},
    }};
}

void Solution::handleField() {
    /* Now we are setting game graph for
     * further topological sort
     */
    for(vert_t complexCoord = 0; complexCoord < GetGraphVSize(); ++complexCoord) {
        bool winning = isInitiallyWinning(ComplexCoordToVertex(complexCoord));
        if(winning) {
            vertices_known_[complexCoord] = true;
            vertices_winning_[complexCoord] = true;
            continue;
        }

        bool loosing = isInitiallyLoosing(ComplexCoordToVertex(complexCoord));
        if(loosing) {
            vertices_known_[complexCoord] = true;
            vertices_winning_[complexCoord] = false;
        }
    }
}

void Solution::runTopologicalSort() {
    for(vert_t complexCoord = 0; complexCoord < GetGraphVSize(); ++complexCoord) {
        if(!top_sort_used_[complexCoord]) {
            topologicalDFS(complexCoord);
        }
    }
}

void Solution::topologicalDFS(vert_t currentVertex) {
    top_sort_used_[currentVertex] = true;
    auto neigbours = getNeighbours(currentVertex);
    for(vert_t neighbour : neigbours) {
        if(!top_sort_used_[neighbour]) {
            topologicalDFS(neighbour);
        }
    }
    top_sort_dfs_.push(currentVertex);
}

bool Solution::FugitiveWins() {
    handleField();

    runTopologicalSort();
    return false;
}

bool Solution::isInitiallyWinning(const Vertex &vertex) {
    /*
     * Fugitive is winning IF:
     * 1) This is his turn
     * 2) He is on the last line
     * 3) There is no terminator blast on him
     */
    return  vertex.fugitiveTurn &&
            vertex.fugitiveCoord.second == (FieldSize::Y - 1) &&
            !isUnderBlast(vertex);
}

bool Solution::isInitiallyLoosing(const Vertex &vertex) {
    /*
     * Fugitive is definitely loosing IF:
     * 1) If he is under terminator's blast
     * 2) It is terminators turn and terminator has a way
     * to make a turn and blast the fugitive
     */

    if(isUnderBlast(vertex)) {
        return true;
    }

    return false;
}

bool Solution::isUnderBlast(const Vertex &vertex) {
    if (vertex.fugitiveCoord.first == vertex.terminatorCoord.first &&
        !isThereWallX(vertex.terminatorCoord.first, vertex))
    {
        return true;
    }

    if(vertex.fugitiveCoord.second == vertex.terminatorCoord.second &&
        !isThereWallY(vertex.terminatorCoord.second, vertex))
    {
        return true;
    }

    comp_t deltaX = std::abs(vertex.fugitiveCoord.first - vertex.terminatorCoord.first);
    comp_t deltaY = std::abs(vertex.fugitiveCoord.second - vertex.terminatorCoord.second);
    return deltaX == deltaY && !isThereWallDiagonal(vertex);
}

bool Solution::isThereWallX(comp_t xComp, const Vertex &vertex) const {
    comp_t fugitiveY = vertex.fugitiveCoord.second;
    comp_t terminatorY = vertex.terminatorCoord.second;

    auto betweenThem = [&](comp_t wallY) -> bool {
        return  ((wallY - fugitiveY > 0) && (terminatorY - wallY > 0)) ||
                ((wallY - fugitiveY < 0) && (terminatorY - wallY < 0));
    };

    for(comp_t y = 0; y < FieldSize::Y; ++y) {
        if(field_[y][xComp] == FieldObjects::WALL && betweenThem(y)) {
            return true;
        }
    }
    return false;
}

bool Solution::isThereWallY(comp_t yComp, const Vertex &vertex) const {
    comp_t fugitiveX = vertex.fugitiveCoord.first;
    comp_t terminatorX = vertex.terminatorCoord.first;

    auto betweenThem = [&](comp_t wallX) -> bool {
        return ((wallX - fugitiveX > 0) && (terminatorX - wallX > 0)) ||
               ((wallX - fugitiveX < 0) && (terminatorX - wallX < 0));
    };

    for(uint8_t x = 0; x < FieldSize::X; ++x) {
        if(field_[yComp][x] == FieldObjects::WALL && betweenThem(x)) {
            return true;
        }
    }
    return false;
}

bool Solution::isThereWallDiagonal(const Vertex &vertex) const {
    auto first  = vertex.fugitiveCoord;
    auto second = vertex.terminatorCoord;

    if(first.second > second.second) {
        std::swap(first, second);
    }

    if(second.first > first.first) {
        for(int8_t i = 0; i < second.first - first.first; ++i) {
            if(field_[first.second + i][first.first + i] == FieldObjects::WALL) {
                return true;
            }
        }
    } else {
        for(int8_t i = 0; i < second.first - first.first; ++i) {
            if(field_[first.second + i][first.first - i] == FieldObjects::WALL) {
                return true;
            }
        }
    }

    return false;
}

bool Solution::canGo(const Vertex &currentVertex, coord_t delta) const {
    comp_t currentX = 0;
    comp_t currentY = 0;

    if(currentVertex.fugitiveTurn) {
        currentX = currentVertex.fugitiveCoord.first;
        currentY = currentVertex.fugitiveCoord.second;
    } else {
        currentX = currentVertex.terminatorCoord.first;
        currentY = currentVertex.terminatorCoord.second;
    }

    comp_t deltaX = delta.first;
    comp_t deltaY = delta.second;

    if (currentX + deltaX < 0 || currentX + deltaX > FieldSize::X ||
        currentY + deltaY < 0 || currentY + deltaY > FieldSize::Y)
    {
        return false;
    }

    return field_[currentY + deltaY][currentX + deltaX] != FieldObjects::WALL;
}


std::vector<vert_t> Solution::getNeighbours(vert_t complexCoord) const {
    Vertex currentVertex = ComplexCoordToVertex(complexCoord);
    std::vector<vert_t> answer;

    Vertex newVertex = currentVertex;
    newVertex.fugitiveTurn = !newVertex.fugitiveTurn;

    for(auto move : possibleMoves_) {
        if(canGo(currentVertex, move)) {
            answer.push_back(VertexToComplexCoord(goDelta(currentVertex, move)));
        }
    }
    return answer;
}

Vertex Solution::goDelta(const Vertex &currentVertex, coord_t delta) {
    Vertex newVertex = currentVertex;
    if(newVertex.fugitiveTurn) {
        newVertex.fugitiveCoord.first += delta.first;
        newVertex.fugitiveCoord.second += delta.second;
    } else {
        newVertex.terminatorCoord.first += delta.first;
        newVertex.terminatorCoord.second += delta.second;
    }

    return newVertex;
}

int main() {
    for(size_t i = 0; i < GetGraphVSize(); ++i) {
        assert(VertexToComplexCoord(ComplexCoordToVertex(i)) == i);
    }

    field_t field;

    std::string buffer;
    for(int8_t y = 0; y < FieldSize::Y; ++y) {
        std::cin >> buffer;
        for(int8_t x = 0; x < FieldSize::X; ++x) {
            field[y][x] = buffer[x] - '0';
        }
    }

    Solution handler(field);
    bool result = handler.FugitiveWins();

    std::cout << (result ? 1 : -1);

    return 0;
}