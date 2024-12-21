#ifndef _HPP_CHAMP
#define _HPP_CHAMP

#include "atom.hpp"
#include "bond.hpp"
#include <iostream>
#include <string>
#include <vector>

/* max bonds an atom can have */

enum class Sym {
  _Null,
  Atom,
  Bond,
  Openscope,
  Closescope,
  Mark,
  Openblock,
  Closeblock,
  Separator,
  Qualifier,
};

class Champ
{
public:
  std::vector<Atom> atoms;

  Champ();
  void parseTag(std::string::iterator tag_chars_iterator,
      std::string::iterator it_end, Bond bond, bool& ok);
  void memoryDump(void)
  {

    std::cout << " ChampMemory: Pat " //<< this->pats.size() << "\n"
              << " ChampMemory: Atom " << this->atoms.size() << "\n";
    //              << " ChampMemory: Bond " << this->bonds.size() << "\n"
    //              << " ChampMemory: Int " << this->Int.size()
    //              << "\n"
    //              //<< " ChampMemory: Int2 " << this->Int2.size() << "\n"
    //              //<< " ChampMemory: Int3 " << this->Int3.size() << "\n"
    //              << " ChampMemory: Tmpl " << this->tmpls.size() << "\n"
    //              << " ChampMemory: Targ " << this->targs.size() << "\n"
    //              << " ChampMemory: Scope " << this->scopes.size() << "\n"
    //              << " ChampMemory: Match " << this->matches.size() <<
    //              std::endl;
  };
};

#endif
