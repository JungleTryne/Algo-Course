#include <algorithm>
#include <array>
#include <cassert>
#include <cstring>
#include <iostream>
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
    FUGITIVE [[maybe_unused]] = 2,
    TERMINATOR = 3
};

using comp_t  = int8_t; //component type
using vert_t  = size_t; //encoded vertex type

using coord_t = std::pair<comp_t , comp_t>;
using field_t = std::array<std::array<comp_t, FieldSize::X>, FieldSize::Y>;

constexpr vert_t GetBase() {
    static_assert(FieldSize::X == FieldSize::Y);
    return FieldSize::X;
}

constexpr vert_t GetGraphVSize() {
    return 2 * FieldSize::X * FieldSize::X * FieldSize::Y * FieldSize::Y;
}

constexpr vert_t GetMaxVertexCode() {
    vert_t maxBaseNum = GetBase() - 1;
    vert_t base = GetBase();
    return  maxBaseNum +
            maxBaseNum * base +
            maxBaseNum * base * base +
            maxBaseNum * base * base * base;
}

constexpr vert_t GetTurnCode() {
    return GetBase() * GetBase() * GetBase() * GetBase();
}

//char is bool actually, but I don't like std::vector<bool>
using graph_vertices_t = std::vector<char>;

struct Vertex {
    coord_t fugitiveCoord;
    coord_t terminatorCoord;
    bool fugitiveTurn;
};


constexpr Vertex ComplexCoordToVertex(uint32_t complexCoord) {
    bool fugitiveTurn = true;
    if(complexCoord > GetMaxVertexCode()) {
        fugitiveTurn = false;
        complexCoord -= GetTurnCode();
    }

    uint32_t base = GetBase();

    comp_t fugitiveX   = complexCoord % base;
    comp_t fugitiveY   = (complexCoord / base) % base;
    comp_t terminatorX = ((complexCoord / base) / base) % base;
    comp_t terminatorY = ((complexCoord / base) / base) / base;

    return Vertex {
            {fugitiveX, fugitiveY},
            {terminatorX, terminatorY},
            fugitiveTurn
    };
}

constexpr uint32_t VertexToComplexCoord(const Vertex& vertex) {
    uint32_t base = GetBase();
    uint32_t result = 0;

    result += vertex.fugitiveCoord.first;
    result += vertex.fugitiveCoord.second * base;
    result += vertex.terminatorCoord.first * base * base;
    result += vertex.terminatorCoord.second * base * base * base;

    if(!vertex.fugitiveTurn) {
        result += GetTurnCode();
    }

    return result;
}

class Solution {
private:
    field_t field_{};

    /* @param vertices_winning: is the vertex winning position
    *  @param vertices_loosing: is the vertex loosing position
    *  @param banned: is the vertex banned (because of a wall)
    */
    graph_vertices_t vertices_winning_;
    graph_vertices_t vertices_loosing_;
    graph_vertices_t banned_;

    graph_vertices_t top_sort_used_;
    std::vector<size_t> loosing_degree_;

    void handleField();
    void getDegrees();
    void banImpossible();
    void traverseTheGraph();

    void dfs(vert_t currentVert);

    std::array<coord_t, 8> possibleMoves_;

    std::vector<vert_t> getNeighbours(vert_t complexCoord) const;

    bool isTheVertexActualGame(const Vertex& vertex);

    bool isUnderBlast(const Vertex& vertex);
    bool isInitiallyWinning(const Vertex& vertex);
    bool isInitiallyLoosing(const Vertex& vertex);

    bool isThereWallX(comp_t xComp, const Vertex &vertex) const;
    bool isThereWallY(comp_t yComp, const Vertex &vertex) const;
    bool isThereWallDiagonal(const Vertex& vertex) const;

    bool canGo(const Vertex &currentVertex, coord_t delta) const;

    static Vertex goDelta(const Vertex& currentVertex, coord_t delta);
public:
    explicit Solution(const field_t& field);
    bool FugitiveWins();
};

Solution::Solution(const field_t &field) : field_(field) {
    vertices_winning_.resize(GetGraphVSize(), false);
    vertices_loosing_.resize(GetGraphVSize(), false);
    top_sort_used_.resize(GetGraphVSize(), false);
    banned_.resize(GetGraphVSize(), false);
    loosing_degree_.resize(GetGraphVSize(), 0);

    possibleMoves_ = {{
      {1, 0},
      {-1, 0},
      {0, -1},
      {0, 1},
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
        if(banned_[complexCoord]) {
            continue;
        }

        bool winning = isInitiallyWinning(ComplexCoordToVertex(complexCoord));
        if(winning) {
            vertices_winning_[complexCoord] = true;
            continue;
        }

        bool loosing = isInitiallyLoosing(ComplexCoordToVertex(complexCoord));
        if(loosing) {
            vertices_loosing_[complexCoord] = true;
        }
    }
}

bool Solution::FugitiveWins() {
    banImpossible();
    getDegrees();
    handleField();
    traverseTheGraph();

    for(size_t i = 0; i < GetGraphVSize(); ++i) {
        if(isTheVertexActualGame(ComplexCoordToVertex(i))) {
            if(vertices_winning_[i]) {
                return true;
            }
            if(vertices_loosing_[i]) {
                return false;
            }
            return false; //not reachable -> not winning
        }
    }

    assert(false);
}

void Solution::dfs(vert_t currentVertex) {
    top_sort_used_[currentVertex] = true;
    auto neighbours = getNeighbours(currentVertex);
    for(vert_t neighbour : neighbours) {
        if(!top_sort_used_[neighbour]) {
            if(vertices_loosing_[currentVertex]) {
                vertices_winning_[neighbour] = true;
                dfs(neighbour);
                continue;
            }

            --loosing_degree_[neighbour];
            if(loosing_degree_[neighbour] == 0) {
                vertices_loosing_[neighbour] = true;
                dfs(neighbour);
            }
        }
    }
}

void Solution::traverseTheGraph() {
    for(size_t i = 0; i < GetGraphVSize(); ++i) {
        if((vertices_winning_[i] || vertices_loosing_[i]) && !top_sort_used_[i]) {
            dfs(i);
        }
    }
}

bool Solution::isInitiallyWinning(const Vertex &vertex) {
    /*
     * Fugitive is winning IF:
     * 1) This is his turn
     * 2) He is on the last line
     * 3) There is no terminator blast on him
     */
    bool fugitiveTurn = vertex.fugitiveTurn;
    bool onEdge = vertex.fugitiveCoord.second == (FieldSize::Y - 1);
    bool underBlast = isUnderBlast(vertex);

    if(fugitiveTurn) {
        return onEdge && !underBlast;
    }

    /*
     * Terminator is winning IF:
     * 1) His coordinate is the same is fugitive's one
     * 2) He can blast the fugitive with extra move
     */
    if(vertex.fugitiveCoord == vertex.terminatorCoord) {
        return true;
    }

    for(auto move : possibleMoves_) {
        if(canGo(vertex, move)) {
            Vertex movedVertex = goDelta(vertex, move);
            if(isUnderBlast(movedVertex)) {
                return true;
            }
        }
    }

    return underBlast;
}

bool Solution::isInitiallyLoosing(const Vertex &vertex) {
    /*
     * Fugitive is definitely loosing IF:
     * 1) If he is under terminator's blast
     */

    if(vertex.fugitiveTurn) {
        if(vertex.fugitiveCoord == vertex.terminatorCoord) {
            return true;
        }
        return isUnderBlast(vertex);
    }

    /*
     * Terminator is definitely loosing if
     * 1) He goes to any direction and can't catch
     * the fugitive who is standing on the last line
     */

    for(auto move : possibleMoves_) {
        if(canGo(vertex, move)) {
            Vertex movedVertex = goDelta(vertex, move);
            if(isUnderBlast(movedVertex)) {
                return false;
            }
        }
    }

    return vertex.fugitiveCoord.second == (FieldSize::Y - 1);
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

    if(fugitiveY == terminatorY) {
        return false;
    }

    auto betweenThem = [&](comp_t wallY) -> bool {
        return  ((wallY - fugitiveY >= 0) && (terminatorY - wallY >= 0)) ||
                ((wallY - fugitiveY <= 0) && (terminatorY - wallY <= 0));
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
        return ((wallX - fugitiveX >= 0) && (terminatorX - wallX >= 0)) ||
               ((wallX - fugitiveX <= 0) && (terminatorX - wallX <= 0));
    };

    for(vert_t x = 0; x < FieldSize::X; ++x) {
        assert(yComp < FieldSize::Y);
        assert(x < FieldSize::X);
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
            assert(first.second + i < FieldSize::Y);
            assert(first.first + i < FieldSize::X);
            if(field_[first.second + i][first.first + i] == FieldObjects::WALL) {
                return true;
            }
        }
    } else {
        for(int8_t i = 0; i < first.first - second.first; ++i) {
            assert(first.second + i < FieldSize::Y);
            assert(first.first - i < FieldSize::X);
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

    if (currentX + deltaX < 0 || currentX + deltaX >= FieldSize::X ||
        currentY + deltaY < 0 || currentY + deltaY >= FieldSize::Y)
    {
        return false;
    }

    return field_[currentY + deltaY][currentX + deltaX] != FieldObjects::WALL;
}

std::vector<vert_t> Solution::getNeighbours(vert_t complexCoord) const {
    Vertex currentVertex = ComplexCoordToVertex(complexCoord);
    std::vector<vert_t> answer;

    Vertex newVertex = currentVertex;
    if(newVertex.fugitiveTurn) {
        newVertex.fugitiveTurn = 0;
    } else {
        newVertex.fugitiveTurn = 1;
    }

    for(auto move : possibleMoves_) {
        if(canGo(newVertex, move)) {
            Vertex vertexToPush = goDelta(newVertex, move);
            answer.push_back(VertexToComplexCoord(vertexToPush));
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

bool Solution::isTheVertexActualGame(const Vertex& vertex) {
    comp_t FugitiveX = vertex.fugitiveCoord.first;
    comp_t FugitiveY = vertex.fugitiveCoord.second;
    comp_t TerminatorX = vertex.terminatorCoord.first;
    comp_t TerminatorY = vertex.terminatorCoord.second;

    return field_[FugitiveY][FugitiveX] == 2 &&
           field_[TerminatorY][TerminatorX] == 3 &&
           vertex.fugitiveTurn;
}

void Solution::getDegrees() {
    for(size_t i = 0; i < GetGraphVSize(); ++i) {
        if(banned_[i]) {
            continue;
        }
        auto neighbours = getNeighbours(i);
        for(vert_t neighbour : neighbours) {
            loosing_degree_[neighbour]++;
        }
    }
}

void Solution::banImpossible() {
    for(size_t i = 0; i < GetGraphVSize(); ++i) {
        Vertex vertex = ComplexCoordToVertex(i);
        comp_t FugitiveX = vertex.fugitiveCoord.first;
        comp_t FugitiveY = vertex.fugitiveCoord.second;
        comp_t TerminatorX = vertex.terminatorCoord.first;
        comp_t TerminatorY = vertex.terminatorCoord.second;

        if(field_[FugitiveY][FugitiveX] == 1 || field_[TerminatorY][TerminatorX] == 1) {
            banned_[i] = true;
        }
    }
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