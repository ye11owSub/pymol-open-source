#ifndef _HPP_ATOM
#define _HPP_ATOM
#include <iterator>
#include <string>
#include <string_view>
#include "constants.hpp"

#define MAX_BOND 12
#define MAX_RING 50

enum class Atom : unsigned int {
  Any = 0xFFFFFFFF,
  H = 0x00000001,
  C = 0x00000002,
  N = 0x00000004,
  O = 0x00000008,
  Sym = 0x00000010,
  S = 0x00000020,
  P = 0x00000040,
  F = 0x00000080,
  Cl = 0x00000100,
  Br = 0x00000200,
  I = 0x00000400,
  Na = 0x00000800,
  K = 0x00001000,
  Ca = 0x00002000,
  Mg = 0x00004000,
  Zn = 0x00008000,
  Fe = 0x00010000,
  Cu = 0x00020000,
  Se = 0x00040000,
  B = 0x00080000,
  A = 0x00100000,
  E = 0x00200000,
  G = 0x00400000,
  J = 0x00800000,
  L = 0x01000000,
  M = 0x02000000,
  Q = 0x04000000,
  R = 0x08000000,
  T = 0x10000000,
  X = 0x20000000,
  Z = 0x40000000,

  NotH = 0xFFFFFFFE,
};

inline unsigned int operator|=(unsigned int& a1, const Atom& a2)
{
  return a1 = a1 | static_cast<unsigned int>(a2);
};

enum class Chirality {
  Unspecified = 0,
  Anticlock = 1,
  Clockwise = -1,
};

enum class Degree : unsigned int {
  Bond0 = 0x00000001,
  Bond1 = 0x00000002,
  Bond2 = 0x00000004,
  Bond3 = 0x00000008,
  Bond4 = 0x00000010,
  Bond5 = 0x00000020,
  Bond6 = 0x00000040,
  Bond7 = 0x00000080,
  Bond8 = 0x00000100,
};

inline unsigned int operator|=(unsigned int& a1, const Degree& a2)
{
  return a1 = a1 | static_cast<unsigned int>(a2);
};

static Degree num_to_degree[9] = {
    Degree::Bond0,
    Degree::Bond1,
    Degree::Bond2,
    Degree::Bond3,
    Degree::Bond4,
    Degree::Bond5,
    Degree::Bond6,
    Degree::Bond7,
    Degree::Bond8,
};

enum class Charge : unsigned int {
  Neutral = 0x00000001,
  Cation = 0x00000002,
  Dication = 0x00000004,
  Anion = 0x00000008,
  Dianion = 0x00000010,
  Trication = 0x00000020,
  Trianion = 0x00000040,
  Tetcation = 0x00000080,
  Tetanion = 0x00000100,
  Pentcation = 0x00000200,
  Pentanion = 0x00000400,
};

inline unsigned int operator|=(unsigned int& c1, const Charge& c2)
{
  return c1 = c1 | static_cast<unsigned int>(c2);
};

static Charge num_to_charge[6] = {
    Charge::Neutral,
    Charge::Cation,
    Charge::Dication,
    Charge::Trication,
    Charge::Tetcation,
    Charge::Pentcation,
};

/* valence */

enum class Valence : unsigned int {
  Valence0 = 0x00000001,
  Valence1 = 0x00000002,
  Valence2 = 0x00000004,
  Valence3 = 0x00000008,
  Valence4 = 0x00000010,
  Valence5 = 0x00000020,
  Valence6 = 0x00000040,
  Valence7 = 0x00000080,
  Valence8 = 0x00000100,
};

inline unsigned int operator|=(unsigned int& d1, const Valence& d2)
{
  return d1 = d1 | static_cast<unsigned int>(d2);
};

static Valence num_to_valence[9] = {
    Valence::Valence0,
    Valence::Valence1,
    Valence::Valence2,
    Valence::Valence3,
    Valence::Valence4,
    Valence::Valence5,
    Valence::Valence6,
    Valence::Valence7,
    Valence::Valence8,
};

class AtomMeta
{
public:
  bool pos_flag = false;
  bool neg_flag = false;
  bool comp_imp_hydro_flag = false; // do we need to compute implicit hydrogens?
  bool hydro_flag = false;          // are we trying to match hydrogen counts?

  Chirality stereo =
      Chirality::Unspecified; // based on lexical ordering: 0 = unspecified, 1 =
                              // anti-clockwise, -1 = clockwise

  // bitmask
  unsigned int atom_type = 0;
  unsigned int atom = 0;
  unsigned int not_atom = 0;
  unsigned int _class = 0;
  unsigned int not_class = 0;
  unsigned int degree = 0;
  unsigned int not_degree = 0;
  unsigned int cycle = 0;
  unsigned int not_cycle = 0;
  unsigned int valence = 0;
  unsigned int not_valence = 0;
  unsigned int charge = 0;
  unsigned int not_charge = 0;

  int imp_hydro;
  int tot_hydro;

  char symbol[3] = {0, 0, 0};
  char name[5] = {0, 0, 0};
  char residue[5] = {0, 0, 0};

  AtomMeta();
  void parseAliphatic(const Atom& atom_mask, const bool& imp_hyd);
  void parseAromatic(const Atom& atom_mask, const bool& imp_hyd);
  void parseBlock(
      const Atom& mask, const bool& not_atom_flag, std::string_view symbols);
  void parseBlock(
      const Atom& mask, const bool& not_atom_flag, const char& symbol);
  bool parseAtomBlock(std::string::iterator& char_iterator,
      const std::string::iterator& it_end);
  void parseString(std::string_view symbols);
  void Dump();
  void FlagDump();

  bool operator==(const AtomMeta& that)
  {
    if (

        this->atom != that.atom || this->charge != that.charge ||
        this->cycle != that.cycle || this->_class != that._class ||
        this->degree != that.degree || this->valence != that.valence ||
        this->name != that.name || this->residue != that.residue ||
        this->symbol != that.symbol || this->symbol != that.symbol ||
        // must have at least as many hydrogens as pattern...
        this->tot_hydro < that.tot_hydro)
      return false;

    return true;

    /*
    {
      if ((((!p->pos_flag) ||
               (((!p->atom) || (p->atom & a->atom)) &&
                   ((!p->charge) || (p->charge & a->charge)) &&
                   ((!p->cycle) || (p->cycle & a->cycle)) &&
                   ((!p->class) || (p->class & a->class)) &&
                   ((!p->degree) || (p->degree & a->degree)) &&
                   ((!p->valence) || (p->valence & a->valence)))) &&
              ((!p->neg_flag) ||
                  (((!p->not_atom) || (!(p->not_atom & a->atom))) &&
                      ((!p->not_charge) || (!(p->not_charge & a->charge)))
    &&
                      ((!p->not_cycle) || (!(p->not_cycle & a->cycle))) &&
                      ((!p->not_class) || (!(p->not_class & a->class))) &&
                      ((!p->not_degree) || (!(p->not_degree & a->degree)))
    &&
                      ((!p->not_valence) ||
                          (!(p->not_valence & a->valence))))))) {
        if (p->name[0])
          if (strcmp(p->name, a->name))
            return 0;
        if (p->residue[0])
          if (strcmp(p->residue, a->residue))
            return 0;
        if (p->symbol[0])
          if (strcmp(p->symbol, a->symbol))
            return 0;
        if (p->hydro_flag) {
          if (p->tot_hydro > a->tot_hydro) { // must have at least as many
                                                //hydrogens as pattern...
            return 0;
          }
        }
        return 1;
      }
      // what about implicit hydrogens?

      return 0;
  */
  }
};
#endif
