#ifndef _HPP_BOND
#define _HPP_BOND

#include <string>

enum class Order : unsigned int {
  Single = 0x00000001,
  Double = 0x00000002,
  Triple = 0x00000004,
  AnyOrder = 0x00000007,
  NoOrder = 0x00000000,
};
inline unsigned int operator|=(unsigned int& o1, const Order& o2)
{
  return o1 = o1 | static_cast<unsigned int>(o2);
};

/* double-bond stereochem */
enum class Direction : int {
  Specified = 0,
  Up = 1,
  Down = -1,
};

struct Bond {
public:
  int atom[2];   // connected atoms  -- directionality must reflect tree
  int pri[2];    // forward and backward lexical priorities
  int mark_tmpl; // traversal
  int mark_targ; // traversal
  int mark_read; // traversal
  int ext_index;
  unsigned int _class = 0;
  unsigned int not_class = 0;
  Direction direction =
      Direction::Specified; // 0 = specified, 1 = up, -1 = down
  unsigned int order = 0;
  unsigned int not_order = 0;
  unsigned int tag;
  unsigned int not_tag;
  unsigned int cycle;
  unsigned int not_cycle;

  Bond();
  void parseTag(std::string::const_iterator tag_chars_iterator,
      std::string::const_iterator it_end, bool& ok);
};

#endif
