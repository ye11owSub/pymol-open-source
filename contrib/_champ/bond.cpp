#include "bond.hpp"
#include "atom.hpp"
#include "champ.hpp"
#include <pybind11/pybind11.h>

Bond::Bond(){};

void Bond::parseTag(std::string::const_iterator tag_chars_iterator,
    std::string::const_iterator it_end, bool& ok)
{
  // parse bit masks like <1> <1,2,3> <12,3,1> etc...

  unsigned int tag_mask;
  int map_index;
  char current_char;
  char next_char;
  bool not_flag = false;

  while (ok) {
    current_char = *tag_chars_iterator;
    next_char = *(std::next(tag_chars_iterator));
    if (current_char == '>') {
      tag_chars_iterator = std::next(tag_chars_iterator);
      return;
    }
    if (tag_chars_iterator == it_end) {
      ok = false;
      break;
    }
    if (current_char == ';') {
      not_flag = false;
      tag_chars_iterator = std::next(tag_chars_iterator);
    } else if (current_char == '!') {
      not_flag = true;
      tag_chars_iterator = std::next(tag_chars_iterator);
    } else if ((current_char >= '0') && (current_char <= '9')) {
      if ((next_char >= '0') && (next_char <= '9')) {
        map_index = (current_char - '0') * 10 + (next_char - '0');
        tag_chars_iterator = std::next(tag_chars_iterator, 2);
      } else {
        map_index = (current_char - '0');
        tag_chars_iterator = std::next(tag_chars_iterator);
      }
      tag_mask = 0x1;
      while (map_index) {
        tag_mask = (tag_mask << 1);
        map_index--;
      }
      if (not_flag) {
        this->not_tag |= tag_mask;
      } else {
        this->tag |= tag_mask;
      }
    } else
      tag_chars_iterator = std::next(tag_chars_iterator);
  }
}
