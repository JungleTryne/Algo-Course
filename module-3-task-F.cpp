#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <vector>
#include <utility>
#include <queue>
#include <set>

namespace VoronoiGlobals {
    inline static constexpr long double eps() {
        return 1e-12;
    };

    inline static constexpr long double subeps() {
        return 1e-11;
    }

    inline static constexpr long double infty() {
        return std::numeric_limits<long double>::infinity();
    }
}


template<typename PointType>
struct Point {
    PointType x;
    PointType y;

    explicit Point(PointType x= 0, PointType y= 0) : x(x), y(y) {};
    Point(const Point<PointType>& other) = default;
    ~Point() = default;

    bool operator==(const Point<PointType>& other) const;
    bool operator!=(const Point<PointType>& other) const;

    long double GetLength() const;
};

template<typename PointType>
bool operator<(const Point<PointType>& first, const Point<PointType>& second) {
    return std::pair<PointType, PointType>(first.x, first.y) <
           std::pair<PointType, PointType>(second.x, second.y);
}

template<typename PointType>
bool Point<PointType>::operator==(const Point<PointType>& other) const {
    if constexpr (std::is_same<PointType, long double>::value) {
        return std::abs(this->x - other.x) <= VoronoiGlobals::eps() &&
               std::abs(this->y - other.y) <= VoronoiGlobals::eps();
    } else {
        return (this->x == other.x) && (this->y == other.y);
    }

}

template<typename PointType>
bool Point<PointType>::operator!=(const Point<PointType>& other) const {
    return !(*this==other);
}

template<typename PointType>
long double Point<PointType>::GetLength() const {
    return sqrt(static_cast<long double>(x) * static_cast<long double>(x) +
                static_cast<long double>(y) * static_cast<long double>(y));
}

template<typename PointType>
Point<PointType> operator-(const Point<PointType>& one, const Point<PointType>& two) {
    Point newPoint(one.x - two.x, one.y - two.y);
    return newPoint;
}

template<typename PointType>
Point<PointType> operator+(const Point<PointType>& one, const Point<PointType>& two) {
    Point newPoint(one.x + two.x, one.y + two.y);
    return newPoint;
}

template<typename PointType>
Point<PointType> operator*(const Point<PointType>& one, PointType coefficient) {
    Point newPoint(one.x * coefficient, one.y*coefficient);
    return newPoint;
}

template<typename PointType>
Point<PointType> operator*(PointType coefficient, Point<PointType>& one) {
    return one*coefficient;
}

std::vector<Point<long double>> GetArchIntersection(const Point<long double>& first,
                                                    const Point<long double>& second, long double dir) {
    std::vector<Point<long double>> result;

    /* Case if foci are on the same line y(x) = C
    */
    if(std::abs(first.y - second.y) < VoronoiGlobals::eps()) {
        auto minToX = first.x < second.x - VoronoiGlobals::eps() ? first : second;
        auto maxToX = first.x > second.x - VoronoiGlobals::eps() ? first : second;
        long double deltaX = (maxToX.x - minToX.x) / 2;
        auto answer = minToX;
        answer.x += deltaX;
        answer.y += sqrt(dir * dir - deltaX * deltaX);
        return {answer};
    }

    /* Else just using some simple math
    */

    /* if we have focus (a,b) and directix c, then parabola formula is following:
     * (x-a)^2 + b^2 - c^2 = 2(b - c)y
     * Then we just find intersection of two parabolas
     */

    long double denominator = first.y - second.y;
    long double fixed = second.x * first.y - first.x * second.y - second.x * dir + first.x * dir;
    long double delta = sqrt((first.x * first.x - 2 * second.x * first.x + second.x * second.x + first.y * first.y +
                              second.y * second.y - 2 * first.y * second.y) * (-first.y * dir - second.y * dir + first.y * second.y +
                                                                               dir * dir));

    long double x1 = (fixed - delta) / denominator;
    long double x2 = (fixed + delta) / denominator;

    long double y1 = (pow((x1 - first.x), 2) + first.y * first.y - dir * dir) / (2 * (first.y - dir));
    long double y2 = (pow((x2 - second.x), 2) + second.y * second.y - dir * dir) / (2 * (second.y - dir));

    if (x1 > x2) {
        std::swap(x1, x2);
        std::swap(y1, y2);
    }

    result.push_back(Point<long double>(x1, y1));
    result.push_back(Point<long double>(x2, y2));
    return result;
}

Point<long double> GetCenterOfCircle(const Point<long double>& first, const Point<long double>& second, const Point<long double>& third) {
    //Getting circle's center, simplifying variables' names
    long double x1 = first.x;
    long double x2 = second.x;
    long double x3 = third.x;

    long double y1 = first.y;
    long double y2 = second.y;
    long double y3 = third.y;

    /*
     * This formula can be gotten if one expand matrix determinant formula
     * Actually circleCenterX = -B/2A, on the numerator it is -B and 2A for denominator
     */
    long double circleCenterX =
            ( (x1 * x1 + y1 * y1) * (y2 - y3) + (x2 * x2 + y2 * y2) * (y3 - y1) + (x3 * x3 + y3 * y3) * (y1 - y2) ) /
            //---------------------------------------------------------------------------------------------------
                            ( 2 * (x1 * (y2 - y3) - y1 * (x2 - x3) + x2 * y3 - x3 * y2) );

    /*
     * This formula can be gotten if one expand matrix determinant formula
     * Actually circleCenterY = -C/2A, on the numerator it is -C and 2A for denominator
     */
    long double circleCenterY =
            ( (x1 * x1 + y1 * y1) * (x3 - x2) + (x2 * x2 + y2 * y2) * (x1 - x3) + (x3 * x3 + y3 * y3) * (x2 - x1) ) /
            //---------------------------------------------------------------------------------------------------
                            ( 2 * (x1 * (y2 - y3) - y1 * (x2 - x3) + x2 * y3 - x3 * y2) );

    return Point(circleCenterX, circleCenterY);
}

class VoronoiDiagramBuilder;

class BeachLineNode;

enum class EventType {SITE, CIRCLE};
struct Event {
    long double sweep_line_coord;
    EventType type;
    bool valid;

    std::optional<Point<long double>> point;
    BeachLineNode* archToDelete = nullptr;

    explicit Event(long double sweep_line, EventType type, bool valid);
    bool operator<(const Event& other) const;
};

struct SharedPtrComparator {
    bool operator()(const std::shared_ptr<Event>& first, const std::shared_ptr<Event>& second) const {
        return first->sweep_line_coord < second->sweep_line_coord;
    }
};

class BeachLine;

class BeachLineNode {
private:
    BeachLineNode* parent_;
    std::unique_ptr<BeachLineNode> left_ = nullptr;
    std::unique_ptr<BeachLineNode> right_ = nullptr;

    int64_t height_ = 0;
    long double* sweepLineCoord_;

    std::pair<Point<long double>, Point<long double>> archs_;

    std::vector<std::shared_ptr<Event>> circleEvent_ = {};

    BeachLineNode* next_ = nullptr;
    BeachLineNode* prev_ = nullptr;

    bool IsLeaf() const { return archs_.first == archs_.second; }

    bool IsRoot() const { return !parent_; }
    int64_t GetHeight() const { return height_; }
    int64_t GetBalance() const {
        if(this->IsLeaf()) {
            return 0;
        }
        return right_->height_ - left_->height_;
    }

    long double GetValue() const;
    friend BeachLine;

    void line_print() const;

public:
    BeachLineNode* GetPrev() const { return this->prev_; }
    BeachLineNode* GetNext() const { return this->next_; }
    Point<long double> GetLeaf() const { return archs_.first; }
    BeachLineNode(BeachLineNode* parent, long double* sweepLineCoord,
                  std::pair<Point<long double>, Point<long double>> archs);

    void SetCircleEvent(std::shared_ptr<Event> circleEvent) {
        circleEvent_.push_back(circleEvent);
    }
};

class BeachLine {
private:
    long double* sweepLineCoord_ = nullptr;
    std::unique_ptr<BeachLineNode> root = nullptr;

    BeachLineNode* FindIntersectingArch(const Point<long double>& newFocus) const;
    BeachLineNode* InsertToTheLeft(BeachLineNode* node, const Point<long double>& newFocus);
    BeachLineNode* InsertToTheRight(BeachLineNode* node, const Point<long double>& newFocus);

    auto DeleteRightRight(const Point<long double>& focusToDelete,
                          BeachLineNode* superParent,
                          const Point<long double> deletedPrev,
                          const Point<long double> deletedNext) -> std::pair<BeachLineNode*, BeachLineNode*>;

    auto DeleteRightLeft(const Point<long double>& focusToDelete,
                         BeachLineNode* superParent,
                         const Point<long double> deletedPrev,
                         const Point<long double> deletedNext) -> std::pair<BeachLineNode*, BeachLineNode*>;

    auto DeleteLeftRight(const Point<long double>& focusToDelete,
                         BeachLineNode* superParent,
                         const Point<long double> deletedPrev,
                         const Point<long double> deletedNext) -> std::pair<BeachLineNode*, BeachLineNode*>;

    auto DeleteLeftLeft(const Point<long double>& focusToDelete,
                        BeachLineNode* superParent,
                        const Point<long double> deletedPrev,
                        const Point<long double> deletedNext) -> std::pair<BeachLineNode*, BeachLineNode*>;

    void Balance(BeachLineNode* currentNode);

    void RotateLeftChild(std::unique_ptr<BeachLineNode>& rotationBase);
    void RotateRightChild(std::unique_ptr<BeachLineNode>& rotationBase);

    void BalanceTreeRightRight(BeachLineNode* currentNode);
    void BalanceTreeRightLeft(BeachLineNode* currentNode);
    void BalanceTreeLeftRight(BeachLineNode* currentNode);
    void BalanceTreeLeftLeft(BeachLineNode* currentNode);

    void FindAndFixRemain(const Point<long double>& first, const Point<long double>&
            second, const Point<long double>& deletedPrev, const Point<long double>& deletedNext);

    std::unique_ptr<BeachLineNode>& ToUnique(BeachLineNode* node);

public:
    BeachLineNode* InsertNewArch(const Point<long double>& newFocus);
    std::pair<BeachLineNode*, BeachLineNode*> DeleteArch(BeachLineNode* node);
    BeachLine(long double* sweepLine) : sweepLineCoord_(sweepLine) {}

    BeachLineNode* GetFirstArch() const;
    friend VoronoiDiagramBuilder;
};

class VoronoiDiagramBuilder {
private:
    std::map<Point<long double>, std::vector<Point<long double>>> vertices;

    std::priority_queue<std::shared_ptr<Event>, std::vector<std::shared_ptr<Event>>, SharedPtrComparator> queue_;
    std::vector<Point<long double>> points_;
    long double sweepLineCoord_ = VoronoiGlobals::infty();

    void HandleSiteEvent(const Event& event);
    void HandleCircleEvent(const Event& event);

    std::shared_ptr<Event> CheckCircleEvent(BeachLineNode* first, BeachLineNode* second, BeachLineNode* third) const;

    BeachLine line{&sweepLineCoord_};

public:
    [[maybe_unused]] explicit VoronoiDiagramBuilder(const std::vector<Point<long double>>& points) : points_(points) {}
    explicit VoronoiDiagramBuilder(std::vector<Point<long double>>&& points) : points_(std::move(points)) {}
    void Build();
    long double GetAverage() const;
};

BeachLineNode::BeachLineNode(BeachLineNode *parent,
                             long double *sweepLineCoord,
                             std::pair<Point<long double>, Point<long double>> archs) :
        parent_(parent), sweepLineCoord_(sweepLineCoord),
        archs_(archs) {}

long double BeachLineNode::GetValue() const {
    if(this->IsLeaf()) {
        return this->archs_.first.x;
    } else {
        auto first_point = this->archs_.first;
        auto second_point = this->archs_.second;
        auto intersections = GetArchIntersection(first_point, second_point, *this->sweepLineCoord_);
        if(intersections.size() == 1) {
            return intersections[0].x;
        }

        if(first_point.y < second_point.y) {
            return intersections[1].x;
        }
        return intersections[0].x;

    }
}

void BeachLineNode::line_print() const {
    std::cout << "(" << this->GetLeaf().x << ", " << this->GetLeaf().y << ") ";
    if(this->next_) {
        this->next_->line_print();
    }
}


BeachLineNode*
BeachLine::FindIntersectingArch(const Point<long double> &newFocus) const {
    BeachLineNode* currentNode = root.get();
    while(!currentNode->IsLeaf()) {
        long double value = currentNode->GetValue();
        if(value <= newFocus.x - VoronoiGlobals::eps()) {
            currentNode = currentNode->right_.get();
        } else {
            currentNode = currentNode->left_.get();
        }
    }
    return currentNode;
}

BeachLineNode* BeachLine::InsertNewArch(const Point<long double> &newFocus) {
    if(!root) {
        root = std::make_unique<BeachLineNode> (
                nullptr,
                sweepLineCoord_,
                std::make_pair(newFocus, newFocus)
        );
        root->height_ = 0;
        return root.get();
    }

    BeachLineNode* ArchToBreak = FindIntersectingArch(newFocus);

    for(auto event : ArchToBreak->circleEvent_) {
        event->valid = false;
    }

    long double value = ArchToBreak->GetValue();
    if(value <= newFocus.x - VoronoiGlobals::eps()) {
        return InsertToTheRight(ArchToBreak, newFocus);
    }
    return InsertToTheLeft(ArchToBreak, newFocus);
}

BeachLineNode*
BeachLine::InsertToTheLeft(BeachLineNode *node, const Point<long double> &newFocus) {
    node->right_ = std::make_unique<BeachLineNode> (
            node,
            sweepLineCoord_,
            node->archs_
    );

    node->archs_ = std::make_pair(newFocus, node->archs_.second);

    node->left_ = std::make_unique<BeachLineNode>(
            node,
            sweepLineCoord_,
            std::make_pair(node->archs_.second, node->archs_.first)
    );

    auto leftNode = node->left_.get();

    leftNode->left_ = std::make_unique<BeachLineNode>(
            leftNode,
            sweepLineCoord_,
            std::make_pair(leftNode->archs_.first, leftNode->archs_.first)
    );

    leftNode->right_ = std::make_unique<BeachLineNode>(
            leftNode,
            sweepLineCoord_,
            std::make_pair(leftNode->archs_.second, leftNode->archs_.second)
    );

    node->right_->prev_ = leftNode->right_.get();
    leftNode->right_->prev_ = leftNode->left_.get();
    leftNode->left_->prev_ = node->prev_;

    leftNode->left_->next_ = leftNode->right_.get();
    leftNode->right_->next_ = node->right_.get();
    node->right_->next_ = node->next_;

    if(node->next_) {
        node->next_->prev_ = node->right_.get();
    }
    if(node->prev_) {
        node->prev_->next_ = leftNode->left_.get();
    }

    node->next_ = nullptr;
    node->prev_ = nullptr;

    leftNode->height_ = 1;
    auto nodeToReturn = leftNode->right_.get();
    assert(!leftNode->IsLeaf());
    auto currentNode = node;
    while(!currentNode->IsRoot()) {
        currentNode->height_ = std::max(currentNode->left_->height_, currentNode->right_->height_) + 1;
        Balance(currentNode);
        currentNode = currentNode->parent_;
    }
    currentNode->height_ = std::max(currentNode->left_->height_, currentNode->right_->height_) + 1;

    return nodeToReturn;
}

BeachLineNode*
BeachLine::InsertToTheRight(BeachLineNode *node, const Point<long double> &newFocus) {
    node->left_ = std::make_unique<BeachLineNode> (
            node,
            sweepLineCoord_,
            node->archs_
    );

    node->archs_ = std::make_pair(node->archs_.second, newFocus);

    node->right_ = std::make_unique<BeachLineNode>(
            node,
            sweepLineCoord_,
            std::make_pair(node->archs_.second, node->archs_.first)
    );

    auto rightNode = node->right_.get();

    rightNode->right_ = std::make_unique<BeachLineNode>(
            rightNode,
            sweepLineCoord_,
            std::make_pair(rightNode->archs_.second, rightNode->archs_.second)
    );

    rightNode->left_ = std::make_unique<BeachLineNode>(
            rightNode,
            sweepLineCoord_,
            std::make_pair(rightNode->archs_.first, rightNode->archs_.first)
    );

    node->left_->next_ = rightNode->left_.get();
    rightNode->left_->next_ = rightNode->right_.get();
    rightNode->right_->next_ = node->next_;

    rightNode->right_->prev_ = rightNode->left_.get();
    rightNode->left_->prev_ = node->left_.get();
    node->left_->prev_ = node->prev_;

    if(node->next_) {
        node->next_->prev_ = rightNode->right_.get();
    }
    if(node->prev_) {
        node->prev_->next_ = node->left_.get();
    }

    node->next_ = nullptr;
    node->prev_ = nullptr;

    auto nodeToReturn = rightNode->left_.get();

    rightNode->height_ = 1;
    assert(!rightNode->IsLeaf());
    auto currentNode = node;

    while(!currentNode->IsRoot()) {
        Balance(currentNode);
        currentNode->height_ = std::max(currentNode->left_->height_, currentNode->right_->height_) + 1;
        currentNode = currentNode->parent_;
    }
    currentNode->height_ = std::max(currentNode->left_->height_, currentNode->right_->height_) + 1;

    return nodeToReturn;
}

void BeachLine::BalanceTreeRightRight(BeachLineNode* currentNode) {
    RotateRightChild(this->ToUnique(currentNode));
}

void BeachLine::BalanceTreeRightLeft(BeachLineNode* currentNode) {
    RotateLeftChild(currentNode->right_);
    auto& uniqueCurrent = ToUnique(currentNode);
    RotateRightChild(uniqueCurrent);
}

void BeachLine::BalanceTreeLeftRight(BeachLineNode* currentNode) {
    RotateRightChild(currentNode->left_);
    auto& uniqueCurrent = ToUnique(currentNode);
    RotateLeftChild(uniqueCurrent);
}

void BeachLine::BalanceTreeLeftLeft(BeachLineNode* currentNode) {
    RotateLeftChild(this->ToUnique(currentNode));
}


std::shared_ptr<Event> VoronoiDiagramBuilder::CheckCircleEvent(
        BeachLineNode *first,
        BeachLineNode *second,
        BeachLineNode *third) const
{
    if(!(first && second && third)) {
        return nullptr;
    }

    auto first_point = first->GetLeaf();
    auto second_point = second->GetLeaf();
    auto third_point = third->GetLeaf();

    //checking the right turn
    if((second_point.x - first_point.x) * (third_point.y - first_point.y) -
       (third_point.x - first_point.x) * (second_point.y - first_point.y) > 0) {
        return nullptr;
    }

    auto circleCenter = GetCenterOfCircle(first_point, second_point, third_point);
    if(circleCenter.y - (circleCenter - first_point).GetLength() > sweepLineCoord_) {
        return nullptr;
    }

    long double circleRadius = (circleCenter - first_point).GetLength();
    std::shared_ptr<Event> circleEvent =
            std::make_shared<Event>(circleCenter.y - circleRadius, EventType::CIRCLE, true);
    circleEvent->point = circleCenter;
    circleEvent->archToDelete = second;

    second->SetCircleEvent(circleEvent);
    return circleEvent;
}

auto BeachLine::DeleteArch(BeachLineNode *node)
-> std::pair<BeachLineNode*, BeachLineNode*>
{
    for(auto event : node->circleEvent_) {
        event->valid = false;
    }

    BeachLineNode* parent = node->parent_;
    assert(parent->parent_);
    auto superParent = parent->parent_;

    auto focusToDelete = node->GetLeaf();
    auto deletedPrev = node->prev_->GetLeaf();
    auto deletedNext = node->next_->GetLeaf();

    if(parent->right_.get() == node) {
        if(superParent->right_.get() == parent) {
            return DeleteRightRight(focusToDelete, superParent, deletedPrev, deletedNext);
        } else {
            return DeleteLeftRight(focusToDelete, superParent, deletedPrev, deletedNext);
        }

    } else {
        if(superParent->right_.get() == parent) {
            return DeleteRightLeft(focusToDelete, superParent, deletedPrev, deletedNext);
        } else {
            return DeleteLeftLeft(focusToDelete, superParent, deletedPrev, deletedNext);
        }
    }
}

BeachLineNode *BeachLine::GetFirstArch() const {
    auto currentNode = root.get();
    while(currentNode->left_) {
        currentNode = currentNode->left_.get();
    }
    return currentNode;
}

auto BeachLine::DeleteRightRight(const Point<long double>& focusToDelete,
                                 BeachLineNode* superParent,
                                 const Point<long double> deletedPrev,
                                 const Point<long double> deletedNext) -> std::pair<BeachLineNode *, BeachLineNode *> {
    auto next = superParent->right_->right_->next_;

    superParent->right_->left_->parent_ = superParent;

    std::unique_ptr<BeachLineNode> tmp = nullptr;
    tmp.swap(superParent->right_->left_);
    tmp.swap(superParent->right_);

    auto mostRight = superParent->right_.get();
    while(mostRight->right_) {
        mostRight = mostRight->right_.get();
    }

    mostRight->next_ = next;
    if(next) {
        next->prev_ = mostRight;
    }

    FindAndFixRemain(focusToDelete, deletedNext, deletedPrev, deletedNext);

    auto currentNode = superParent;
    while(!currentNode->IsRoot()) {
        Balance(currentNode);
        currentNode->height_ = std::max(currentNode->left_->height_, currentNode->right_->height_) + 1;
        currentNode = currentNode->parent_;
    }
    currentNode->height_ = std::max(currentNode->left_->height_, currentNode->right_->height_) + 1;

    return std::make_pair(mostRight, mostRight->next_);
}

auto BeachLine::DeleteRightLeft(const Point<long double> &focusToDelete, BeachLineNode *superParent,
                                const Point<long double> deletedPrev,
                                const Point<long double> deletedNext) -> std::pair<BeachLineNode *, BeachLineNode *> {
    auto prev = superParent->right_->left_->prev_;

    superParent->right_->right_->parent_ = superParent;

    std::unique_ptr<BeachLineNode> tmp = nullptr;
    tmp.swap(superParent->right_->right_);
    tmp.swap(superParent->right_);

    auto mostLeft = superParent->right_.get();
    while(mostLeft->left_) {
        mostLeft = mostLeft->left_.get();
    }

    mostLeft->prev_ = prev;
    if(prev) {
        prev->next_ =mostLeft;
    }

    FindAndFixRemain(deletedPrev, focusToDelete, deletedPrev, deletedNext);

    auto currentNode = superParent;
    while(!currentNode->IsRoot()) {
        Balance(currentNode);
        currentNode->height_ = std::max(currentNode->left_->height_, currentNode->right_->height_) + 1;
        currentNode = currentNode->parent_;
    }
    currentNode->height_ = std::max(currentNode->left_->height_, currentNode->right_->height_) + 1;

    return std::make_pair(mostLeft->prev_, mostLeft);
}

auto BeachLine::DeleteLeftRight(const Point<long double> &focusToDelete, BeachLineNode *superParent,
                                const Point<long double> deletedPrev,
                                const Point<long double> deletedNext) -> std::pair<BeachLineNode *, BeachLineNode *> {
    auto next = superParent->left_->right_->next_;

    superParent->left_->left_->parent_ = superParent;

    std::unique_ptr<BeachLineNode> tmp = nullptr;
    tmp.swap(superParent->left_->left_);
    tmp.swap(superParent->left_);

    auto mostRight = superParent->left_.get();
    while(mostRight->right_) {
        mostRight = mostRight->right_.get();
    }

    mostRight->next_ = next;
    if(next) {
        next->prev_ = mostRight;
    }

    FindAndFixRemain(focusToDelete, deletedNext, deletedPrev, deletedNext);

    auto currentNode = superParent;
    while(!currentNode->IsRoot()) {
        Balance(currentNode);
        currentNode->height_ = std::max(currentNode->left_->height_, currentNode->right_->height_) + 1;
        currentNode = currentNode->parent_;
    }
    currentNode->height_ = std::max(currentNode->left_->height_, currentNode->right_->height_) + 1;

    return std::make_pair(mostRight, mostRight->next_);
}

auto BeachLine::DeleteLeftLeft(const Point<long double> &focusToDelete, BeachLineNode *superParent,
                               const Point<long double> deletedPrev,
                               const Point<long double> deletedNext) -> std::pair<BeachLineNode *, BeachLineNode *> {
    auto prev = superParent->left_->left_->prev_;

    superParent->left_->right_->parent_ = superParent;

    std::unique_ptr<BeachLineNode> tmp = nullptr;
    tmp.swap(superParent->left_->right_);
    tmp.swap(superParent->left_);

    auto mostLeft = superParent->left_.get();
    while(mostLeft->left_) {
        mostLeft = mostLeft->left_.get();
    }

    mostLeft->prev_ = prev;
    if(prev) {
        prev->next_ = mostLeft;
    }

    FindAndFixRemain(deletedPrev, focusToDelete, deletedPrev, deletedNext);

    auto currentNode = superParent;
    while(!currentNode->IsRoot()) {
        Balance(currentNode);
        currentNode->height_ = std::max(currentNode->left_->height_, currentNode->right_->height_) + 1;
        currentNode = currentNode->parent_;
    }
    currentNode->height_ = std::max(currentNode->left_->height_, currentNode->right_->height_) + 1;

    return std::make_pair(mostLeft->prev_, mostLeft);
}

std::unique_ptr<BeachLineNode> &BeachLine::ToUnique(BeachLineNode *node) {
    if(node->IsRoot()) {
        return root;
    }
    if(node->parent_->left_.get() == node) {
        return node->parent_->left_;
    } else if (node->parent_->right_.get() == node) {
        return node->parent_->right_;
    }
    assert(false);
}

void BeachLine::RotateLeftChild(std::unique_ptr<BeachLineNode>& rotationBase) {
    auto parent = rotationBase->parent_;

    auto leftChild = std::move(rotationBase->left_);
    rotationBase->left_ = std::move(leftChild->right_);
    leftChild->right_ = std::move(rotationBase);
    rotationBase = std::move(leftChild);

    rotationBase->right_->parent_ = rotationBase.get();
    rotationBase->parent_ = parent;

    rotationBase->right_->left_->parent_ = rotationBase->right_.get();
    rotationBase->right_->right_->parent_ = rotationBase->right_.get();

    rotationBase->right_->height_ = std::max(rotationBase->right_->left_->height_, rotationBase->right_->right_->height_) + 1;
    rotationBase->height_ = std::max(rotationBase->left_->height_, rotationBase->right_->height_) + 1;
}

void BeachLine::RotateRightChild(std::unique_ptr<BeachLineNode> &rotationBase) {
    auto parent = rotationBase->parent_;

    auto rightChild = std::move(rotationBase->right_);
    rotationBase->right_ = std::move(rightChild->left_);
    rightChild->left_ = std::move(rotationBase);
    rotationBase = std::move(rightChild);

    rotationBase->left_->parent_ = rotationBase.get();
    rotationBase->parent_ = parent;

    rotationBase->left_->left_->parent_ = rotationBase->left_.get();
    rotationBase->left_->right_->parent_ = rotationBase->left_.get();

    rotationBase->left_->height_ = std::max(rotationBase->left_->left_->height_, rotationBase->left_->right_->height_) + 1;
    rotationBase->height_ = std::max(rotationBase->left_->height_, rotationBase->right_->height_) + 1;
}

void BeachLine::Balance(BeachLineNode* currentNode) {
    if(currentNode->GetBalance() > 1) {
        if(currentNode->right_->GetBalance() < 0) {
            BalanceTreeRightLeft(currentNode);
        } else {
            BalanceTreeRightRight(currentNode);
        }
    } else if(currentNode->GetBalance() < -1) {
        if(currentNode->left_->GetBalance() < 0) {
            BalanceTreeLeftLeft(currentNode);
        } else {
            BalanceTreeLeftRight(currentNode);
        }
    }
}

void BeachLine::FindAndFixRemain(const Point<long double> &first, const Point<long double> &second,
                                 const Point<long double> &deletedPrev, const Point<long double> &deletedNext) {
    auto currentNode = root.get();
    long double wantedValueFirst = GetArchIntersection(first, second, *sweepLineCoord_)[0].x;
    long double wantedValueSecond = GetArchIntersection(first, second, *sweepLineCoord_)[1].x;

    while(!currentNode->IsLeaf() && !(currentNode->archs_.first == first && currentNode->archs_.second == second)) {
        long double value = currentNode->GetValue();
        if(value < wantedValueFirst) {
            currentNode = currentNode->right_.get();
        } else {
            currentNode = currentNode->left_.get();
        }
    }

    if(currentNode->IsLeaf()) {
        currentNode = root.get();
        while(!currentNode->IsLeaf() && !(currentNode->archs_.first == first && currentNode->archs_.second == second)) {
            long double value = currentNode->GetValue();
            if(value < wantedValueSecond) {
                currentNode = currentNode->right_.get();
            } else {
                currentNode = currentNode->left_.get();
            }
        }
    }

    if(currentNode) {
        currentNode->archs_ = std::make_pair(deletedPrev, deletedNext);
    }
}

void VoronoiDiagramBuilder::Build() {
    for(const auto& point : points_) {
        auto new_event = std::make_shared<Event>(point.y, EventType::SITE, true);
        new_event->point = point;
        queue_.push(std::move(new_event));
    }

    while(!queue_.empty()) {

        auto event = *queue_.top();
        queue_.pop();

        sweepLineCoord_ = event.sweep_line_coord;

        if(!event.valid) {
            continue;
        }
        if(event.type == EventType::SITE) {
            HandleSiteEvent(event);
        }
        else {
            HandleCircleEvent(event);
        }

    }
}

void VoronoiDiagramBuilder::HandleSiteEvent(const Event &event) {
    assert(event.point.has_value());

    BeachLineNode* newLeaf = this->line.InsertNewArch(*event.point);

    BeachLineNode* nextNewLeaf = newLeaf->GetNext();
    BeachLineNode* nextNextNewLeaf = nullptr;
    if(nextNewLeaf) {
        nextNextNewLeaf = nextNewLeaf->GetNext();
    }

    std::shared_ptr<Event> firstCircleEvent = nullptr;

    if(newLeaf && nextNewLeaf && nextNextNewLeaf) {
        firstCircleEvent = CheckCircleEvent(newLeaf, nextNewLeaf, nextNextNewLeaf);
    }

    BeachLineNode* prevNewLeaf = newLeaf->GetPrev();
    BeachLineNode* prevPrevNewLeaf = nullptr;
    if(prevNewLeaf) {
        prevPrevNewLeaf = prevNewLeaf->GetPrev();
    }

    std::shared_ptr<Event> secondCircleEvent = nullptr;

    if(newLeaf && prevNewLeaf && prevPrevNewLeaf) {
        secondCircleEvent = CheckCircleEvent(prevPrevNewLeaf, prevNewLeaf, newLeaf);
    }

    if(firstCircleEvent) {
        queue_.push(firstCircleEvent);
    }

    if(secondCircleEvent) {
        queue_.push(secondCircleEvent);
    }
}

void VoronoiDiagramBuilder::HandleCircleEvent(const Event &event) {
    auto circleCenter = event.point;

    BeachLineNode* nodeToDelete = event.archToDelete;
    assert(nodeToDelete);

    auto left = nodeToDelete->GetPrev();
    auto right = nodeToDelete->GetNext();

    vertices[left->GetLeaf()].push_back(*circleCenter);
    vertices[right->GetLeaf()].push_back(*circleCenter);
    vertices[nodeToDelete->GetLeaf()].push_back(*circleCenter);

    auto [prev, next] = line.DeleteArch(nodeToDelete);

    sweepLineCoord_ -= VoronoiGlobals::subeps();

    std::shared_ptr<Event> firstCircleEvent = nullptr;
    if(prev && prev->GetPrev() && prev->GetNext()) {
        firstCircleEvent = CheckCircleEvent(prev->GetPrev(), prev, prev->GetNext());
    }

    std::shared_ptr<Event> secondCircleEvent = nullptr;
    if(next && next->GetPrev() && next->GetNext()) {
        secondCircleEvent = CheckCircleEvent(next->GetPrev(), next, next->GetNext());
    }

    sweepLineCoord_ += VoronoiGlobals::subeps();

    if(firstCircleEvent) {
        queue_.push(firstCircleEvent);
    }

    if(secondCircleEvent) {
        queue_.push(secondCircleEvent);
    }
}

long double VoronoiDiagramBuilder::GetAverage() const {
    std::set<Point<long double>> badPoints;

    auto currentNode = line.GetFirstArch();
    while(currentNode->GetNext()) {
        badPoints.insert(currentNode->GetLeaf());
        currentNode = currentNode->GetNext();
    }

    long double answer = 0;
    uint64_t counter = 0;

    for(const auto& point : points_) {
        if(badPoints.find(point) == badPoints.end()) {
            answer += vertices.at(point).size();
            ++counter;
        }
    }

    if (counter == 0) {
        return 0;
    }
    return answer / counter;
}

bool Event::operator<(const Event &other) const {
    return sweep_line_coord < other.sweep_line_coord;
}

Event::Event(long double sweep_line, EventType type, bool valid) :
        sweep_line_coord(sweep_line), type(type), valid(valid) {}

int main() {
    long double x = 0.0;
    long double y = 0.0;

    std::vector<Point<long double>> points;

    while(std::cin >> x) {
        std::cin >> y;
        points.emplace_back(x, y);
    }

    VoronoiDiagramBuilder builder(std::move(points));
    builder.Build();
    std::cout << std::fixed << std::setprecision(7) <<  builder.GetAverage() << std::endl;
    return 0;
}