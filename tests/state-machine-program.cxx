#include <array>
#include <string>
#include <iostream>
#include <stdexcept>
#include <iomanip>
#include <cassert>
#include <cstring>

// This program tests the strategy of marking local extremes as explored left and/or right.
//
// The actual Algorithm marks LocalExtreme's as explored (in the version as of this moment
// git hash c9830f960c1cb489916ea63de79903a1ecbd1466) as follows:
//
// * In Algorithm::update_energy if check_energy_ is set and too much energy is being used
//   and last_extreme_cubic_ is set, then last_extreme_cubic_->local_extreme() is marked
//   as having been explored towards the current hdirection_.
//   After which handle_abort_hdirection is called.
//
//   This translates for our purpose into:
//
//   Case A
//   - If the current_index_ is 0 and current_direction_ is left, then mark node 0 as explored to the left,
//     reverse the current direction and jump to any even index that was already visited since the last
//     jump, including the current node (0).
//   Case B
//   - If the current_index_ is number_of_nodes - 1 and current_direction_ is right, then mark node
//     number_of_nodes - 1 as explored to the right, reverse the current direction and jump to any even
//     index that was already visited since the last jump, including the current node (number_of_nodes - 1).
//
// * In Algorithm::handle_local_extreme_jump, if the current_direction_ was not already explored,
//   and there is an adjacent local extreme in that direction, then that adjacent local extreme is
//   marked as having been explored in the direction opposite to current_direction_ and we call
//   handle_local_extreme_jump as if to continue from that neighbor.
//
//   This translates for our purpose into:
//
//   Case C
//   - in State::jump(index) if nodes_[index] is not already explored in the current_direction_,
//     and there exists a neighbor into the current_direction_ that was already visited before,
//     then that neighbor is marked as having been explored in the direction opposite to current_direction_
//     and we pretend that the jump was to this neighbor instead.
//
// * In Algorithm::handle_local_extreme, if there is an adjacent local extreme in the direction that
//   we came from, then that neighbor is marked as explored in the current_direction_.
//
//   This translates for our purpose into:
//
//   Case D
//   - in State::advance, the previous node, that we just came from, is marked as explored in the current direction.
//
// * Also in Algorithm::handle_local_extreme, if saw_minimum_ is set, then the current local extreme
//   is marked as explored in the current_direction_.
//
//   This translates for our purpose into:
//
//   Case E
//   - in State::advance, also, if saw_minimum_ is set, mark the current node as being explored
//     in the direction opposite to current_direction_.

constexpr int number_of_nodes = 9;

enum Direction { left, right };

class State;

class Node
{
 private:
  unsigned char index_ : 3;
  unsigned char visited_ : 3;
  unsigned char explored_left_ : 1;
  unsigned char explored_right_ : 1;

 public:
  Node() : index_(0x7), visited_(0), explored_left_(false), explored_right_(false) { }

  int get_index() const { return index_; }
  bool is_visited() const { return visited_; }
  bool is_visited_since_last_jump(int jump_generation) const { return visited_ == jump_generation; }
  bool is_explored(Direction dir) const { return (dir == left) ? explored_left_ : explored_right_; }

  void set_index(int index)
  {
    index_ = index;
  }

  void set_visited(int jump_generation)
  {
    visited_ = jump_generation;
  }

  void set_explored(Direction dir, State* state);
};

class State
{
 private:
  unsigned char current_index_:4;
  unsigned char current_direction_:1;
  unsigned char jump_generation_:3;
  std::array<Node, number_of_nodes> nodes_;
  bool saw_minimum_;
  bool just_jumped_;
  int str_length_;

 public:
  State(int current_index, Direction current_direction)
    : current_direction_(current_direction), jump_generation_(1), saw_minimum_(false), just_jumped_(false), str_length_(0)
  {
    for (int i = 0; i < number_of_nodes; ++i)
      nodes_[i].set_index(i);
    set_current_index(current_index);
    print("");
  }

  int get_current_index() const { return current_index_; }
  Direction get_current_direction() const { return static_cast<Direction>(current_direction_); }
  Direction get_opposite_direction() const { return static_cast<Direction>(1 - current_direction_); }
  int get_str_length() const { return str_length_; }

  bool set_current_index(int index)
  {
    current_index_ = static_cast<unsigned char>(index);
    if (current_index_ % 2 == 0)
      saw_minimum_ = true;
    if (current_index_ % 2 == 1 && nodes_[current_index_].is_visited())
      return true;
    nodes_[current_index_].set_visited(jump_generation_);
    return false;
  }

  void set_current_direction(Direction dir) { current_direction_ = dir; }

  void run();
  void jump(int index, char const* arrow);
  bool advance();
  void print(char const* arrow);
  void add_str_length(int length)
  {
    str_length_ += length;
  }
};

void Node::set_explored(Direction dir, State* state)
{
  if (dir == left)
    explored_left_ = true;
  else
    explored_right_ = true;
  std::cout << "[" << (dir == left ? "⇽" : "⇾") << get_index() << "]" << std::flush;
  state->add_str_length(4);  // 4: [ + arrow + index + ].
}

void State::print(char const* arrow)
{
  Node& current_node = nodes_[current_index_];
  std::cout << arrow;
  int len = std::strlen(arrow) == 0 ? 0 : 3;    // All our non-empty arrows are: space + arrow + space.
  bool done = false;
  if (current_node.is_explored(left))
  {
    if (current_node.is_explored(right))
    {
      std::cout << "⇿";
      done = true;
    }
    else
      std::cout << "⇽";
    ++len;
  }
  else if (current_node.is_explored(right))
  {
    std::cout << "⇾";
    ++len;
  }
  std::cout << get_current_index();
  ++len;
  if (!done)
  {
    std::cout << (current_direction_ == left ? "⇐" : "⇒") << std::flush;
    ++len;
  }
  add_str_length(len);
}

bool State::advance()
{
  Node& current_node = nodes_[current_index_];

  if (get_current_direction() == right && current_node.is_explored(right))
  {
    std::cout << std::endl;
    throw std::runtime_error("Unrecoverable state: trying to advance right while that is already explored.");
  }
  if (get_current_direction() == left && current_node.is_explored(left))
  {
    std::cout << std::endl;
    throw std::runtime_error("Unrecoverable state: trying to advance left while that is already explored.");
  }

  if (get_current_direction() == right && get_current_index() == number_of_nodes - 1)
  {
    // If we reach the end of the list, turn a call to advance into a jump to itself.
    //   Case B
    //   - If the current_index_ is number_of_nodes - 1 and current_direction_ is right, then mark node
    //     number_of_nodes - 1 as explored to the right, [reverse the current direction and jump to any even
    //     index that was already visited since the last jump], including the current node (number_of_nodes - 1).
    // Note: we only jump to ourselves here (which does reverse the direction). Jumping to other already
    // visited even nodes is handled in State::run.
    current_node.set_explored(right, this);
    jump(current_index_, " ⇸ ");
    return true;
  }
  else if (get_current_direction() == left && get_current_index() == 0)
  {
    // If we reach the end of the list, turn a call to advance into a jump to itself.
    //   Case A
    //   - If the current_index_ is 0 and current_direction_ is left, then mark node 0 as explored to the left,
    //     [reverse the current direction and jump to any even index that was already visited since the last
    //     jump], including the current node (0).
    // Note: we only jump to ourselves here (which does reverse the direction). Jumping to other already
    // visited even nodes is handled in State::run.
    current_node.set_explored(left, this);
    jump(current_index_, " ⇷ ");
    return true;
  }
  else
  {
    // Remember saw_minimum and the current index before advancing.
    bool saw_minimum = saw_minimum_;
    int prev_index = get_current_index();
    bool error;
    if (get_current_direction() == right)
      error = set_current_index(prev_index + 1);
    else
      error = set_current_index(prev_index - 1);
    print(" ⇀ ");
    if (error)
    {
      std::cout << std::endl;
      throw std::runtime_error("Visiting an odd index for the second time!");
    }
    // Case D
    //   - in State::advance, the previous node, that we just came from, is marked as explored in the current direction.
    nodes_[prev_index].set_explored(get_current_direction(), this);
    if (saw_minimum)
    {
      // Case E
      //   - in State::advance, also, if saw_minimum_ is set, mark the current node as being explored
      //     in the direction opposite to current_direction_.
      nodes_[get_current_index()].set_explored(get_opposite_direction(), this);

      // Case F
      if (current_index_ % 2 == 0)
      {
        int dir = get_current_direction() == right ? -2 : 2;
        assert(0 <= current_index_ + dir && current_index_ + dir < number_of_nodes);
        nodes_[current_index_ + dir].set_explored(get_current_direction(), this);
      }
    }
  }
  return false;
}

void State::jump(int index, char const* arrow)
{
  assert(index % 2 == 0);
  saw_minimum_ = false;
  ++jump_generation_;
  set_current_index(index);
  set_current_direction(current_direction_ == left ? right : left);
  print(arrow);

  // Are we done?
  if (nodes_[current_index_].is_explored(left) && nodes_[current_index_].is_explored(right))
    return;

  // 4) If that new node was already marked as explored in both directions the algorithm stops.
  // Otherwise we continue from there in the direction that wasn't explored yet.
  if (nodes_[current_index_].is_explored(get_current_direction()))
    current_direction_ = 1 - current_direction_;

  // Note that if at that point it is "discovered" that the next node was already visited, then we go there immediately
  // (we would anyway, but this is a special way of moving to the next node).
  bool ignored_error = false;
  if (current_direction_ == right && current_index_ < number_of_nodes - 1 && nodes_[current_index_ + 1].is_visited())
  {
    ignored_error = set_current_index(current_index_ + 1);
  }
  else if (current_direction_ == left && current_index_ > 0 && nodes_[current_index_ - 1].is_visited())
  {
    ignored_error = set_current_index(current_index_ - 1);
  }
  else
    return;     // Point 4 doesn't apply.
  assert(ignored_error);
  print(" ⇝ ");
  //   Case C
  //   - in State::jump(index) if nodes_[index] is not already explored in the current_direction_,
  //     and there exists a neighbor into the current_direction_ that was already visited before,
  //     then that neighbor is marked as having been explored in the direction opposite to current_direction_
  //     and we pretend that the jump was to this neighbor instead.
  nodes_[current_index_].set_explored(current_direction_ == left ? right : left, this);
}

void State::run()
{
  Node& current_node = nodes_[current_index_];

  if (current_node.is_explored(left) && current_node.is_explored(right))
  {
    if (current_index_ % 2 == 0)
    {
      std::cout << " DONE" << std::endl;
      return;
    }
    else
    {
      std::cout << std::endl;
      throw std::runtime_error("Unrecoverable state: both directions explored at odd index " + std::to_string(current_index_));
    }
  }

  int ways_to_advance = 1;

  // If one of these is true then Case A or B might apply; we might either set the current node
  // as explored in the current direction and then jump to another node, or we can just immediately
  // jump to another node.
  bool at_begin = (current_index_ == 0 && get_current_direction() == left &&
      !nodes_[current_index_].is_explored(get_current_direction()));
  bool at_end = (current_index_ == number_of_nodes - 1 && get_current_direction() == right &&
      !nodes_[current_index_].is_explored(get_current_direction()));

  int possibilities_per_jump = (at_begin || at_end) ? 2 : 1;

  if (current_index_ % 2 == 0)
  {
    for (int i = 0; i < number_of_nodes; i += 2)
      if (i != current_index_ && nodes_[i].is_visited_since_last_jump(jump_generation_))
        ways_to_advance += possibilities_per_jump;
  }

  for (int i = 0; i < ways_to_advance; ++i)
  {
    State new_state = *this;

    if (i == 0)
    {
      if (current_index_ % 2 == 0 && nodes_[current_index_].is_explored(get_current_direction()))
      {
        // In this case it *is* allowed to jump to ourselves, because the reason for the jump
        // is not that we just found a worse minimum, but we ran into a node that was already
        // explored into the current direction.
        if (!nodes_[current_index_].is_visited_since_last_jump(jump_generation_))
          continue;
        new_state.jump(current_index_, " ↺ ");
        just_jumped_ = true;
      }
      else
        just_jumped_ = new_state.advance();
    }
    else
    {
      int jump_target = -1;
      int k = 0;
      for (int j = 0; j < number_of_nodes; j += 2)
      {
        if (j != current_index_ && nodes_[j].is_visited_since_last_jump(jump_generation_))
        {
          if (k++ == (i - 1) / possibilities_per_jump)
          {
            jump_target = j;
            break;
          }
        }
      }
      if (i % possibilities_per_jump == 1)
      {
        if (at_begin)
        {
          std::cout << " ⇷ ";
          new_state.add_str_length(3);
          nodes_[0].set_explored(left, &new_state);
        }
        else
        {
          std::cout << " ⇸ ";
          new_state.add_str_length(3);
          nodes_[number_of_nodes - 1].set_explored(right, &new_state);
        }
      }
      new_state.jump(jump_target, " ↷ ");
      just_jumped_ = true;
    }

    try
    {
      new_state.run();
    }
    catch (std::runtime_error const& error)
    {
      std::cout << error.what() << std::endl;
    }

    if (i < ways_to_advance - 1)
      std::cout << std::setw(get_str_length()) << "";
  }
}

int main()
{
  for (int start_index = 0; start_index < number_of_nodes; ++start_index)
  {
    for (Direction start_direction : {left, right})
    {
      State initial_state(start_index, start_direction);
      initial_state.run();
      std::cout << std::endl;
    }
  }
}
