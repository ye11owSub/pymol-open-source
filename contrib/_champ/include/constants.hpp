#ifndef _HPP_CONSTANTS
#define _HPP_CONSTANTS
/* cycles */
enum class Cycles : unsigned int {

  Acyclic = 0x00000001,
  Ring3 = 0x00000002,
  Ring4 = 0x00000004,
  Ring5 = 0x00000008,
  Ring6 = 0x00000010,
  Ring7 = 0x00000020,
  Ring8 = 0x00000040,
  RingN = 0x80000000, // not yet implemented

  Cyclic = 0xFFFFFFFE,
  Empty = 0,
};

inline unsigned int operator|(const unsigned int& c1, const Cycles& c2)
{
  return c1 | static_cast<unsigned int>(c2);
};

inline unsigned int operator|=(unsigned int& c1, const Cycles& c2)
{
  return c1 = c1 | static_cast<unsigned int>(c2);
};

inline unsigned int operator|(const Cycles& c1, const Cycles& c2)
{
  return static_cast<unsigned int>(c1) | static_cast<unsigned int>(c2);
};

static Cycles num_to_ring[12] = {
    Cycles::Empty,
    Cycles::Empty,
    Cycles::Empty,
    Cycles::Ring3,
    Cycles::Ring4,
    Cycles::Ring5,
    Cycles::Ring6,
    Cycles::Ring7,
    Cycles::Ring8,
    Cycles::Empty,
    Cycles::Empty,
    Cycles::Empty,
};

enum class Class : unsigned int {
  Aliphatic = 0x00000001,
  Aromatic = 0x00000002,
  Any = 0x00000003,
  Pi = 0x00000004, // non-exclusive with above
};

inline unsigned int operator|=(unsigned int& c1, const Class& c2)
{
  return c1 = c1 | static_cast<unsigned int>(c2);
};

#endif
