#include "atom.hpp"
#include "constants.hpp"
#include <cstring>
#include <iterator>

int parseNumeral(const char c)
{
  return c >= '0' && c <= '9' ? c - '0' : -1;
}

AtomMeta::AtomMeta(){};
void AtomMeta::parseAliphatic(const Atom& atom_mask, const bool& imp_hyd)
{
  this->atom |= atom_mask;
  this->pos_flag = true;
  this->comp_imp_hydro_flag = imp_hyd;
};

void AtomMeta::parseAromatic(const Atom& atom_mask, const bool& imp_hyd)
{
  this->atom |= atom_mask;
  this->atom_type |= Class::Aliphatic;
  this->pos_flag = true;
  this->comp_imp_hydro_flag = imp_hyd;
};

void AtomMeta::parseString(std::string_view symbols)
{
  this->atom |= Atom::Any;
  this->symbol[0] = symbols[0];
  this->symbol[1] = symbols[1];
  this->pos_flag = true;
}

void AtomMeta::parseBlock(
    const Atom& mask, const bool& not_atom_flag, const char& symbol)
{
  if (not_atom_flag) {
    this->not_atom |= mask;
    this->neg_flag = true;
  } else {
    this->atom |= mask;
    this->pos_flag = true;
  }
  this->hydro_flag = true;
  if (mask == Atom::Sym) {
    this->symbol[0] = symbol;
    this->symbol[1] = 0;
  }
  /* need to include code for loading symbol */
}

void AtomMeta::parseBlock(
    const Atom& mask, const bool& not_atom_flag, std::string_view symbols)
{
  if (not_atom_flag) {
    this->not_atom |= mask;
    this->neg_flag = true;
  } else {
    this->atom |= mask;
    this->pos_flag = true;
  }
  this->hydro_flag = true;
  if (mask == Atom::Sym) {
    std::strncpy(this->symbol, symbols.data(), 2);
  }
  /* need to include code for loading symbol */
}

bool AtomMeta::parseAtomBlock(
    std::string::iterator& char_iterator, const std::string::iterator& it_end)
{
  bool not_flag = false;
  bool ok = true;
  bool done = false;
  bool atom_seen = false;
  char current_char;
  char next_char;
  int num;
  int iter_step;

  this->comp_imp_hydro_flag = false;

  while (ok && !done && char_iterator != it_end) {
    iter_step = 1;

    current_char = *char_iterator;
    next_char = *(std::next(char_iterator));

    if (current_char == ']') {
      done = true;

    } else if (current_char == '!') {
      not_flag = true;
      atom_seen = false;

    } else if (current_char == ',') {
      atom_seen = false;

    } else if (current_char == ';') {
      not_flag = false;
      atom_seen = false;

    } else if (current_char == '*') { // nonstandard
      this->parseBlock(Atom::Any, not_flag, current_char);
      atom_seen = true;

    } else if (current_char == '?') { // nonstandard
      this->parseBlock(Atom::NotH, not_flag, current_char);
      atom_seen = true;

    } else if (current_char == '@') {
      num = parseNumeral(next_char);
      if (num != -1) {
        iter_step++;
      } else {
        num = 1;
        // We have to set iter_step = 0,
        // coz iterator will be updated in while loop
        iter_step = 0;
        while (char_iterator != it_end || *char_iterator != '@') {
          num++;
          char_iterator = std::next(char_iterator);
        }
      }
      this->stereo = num & 0x1 ? Chirality::Anticlock : Chirality::Clockwise;

    } else if (current_char == 'A') {
      // note there is no way to address the 'A' symbol ...
      if (std::string("cglmsu").find(next_char) != std::string::npos) {
        this->parseBlock(
            Atom::Sym, not_flag, std::string() + current_char + next_char);
        atom_seen = true;
        iter_step++;
      } else {
        if (not_flag) {
          this->neg_flag = true;
          this->not_class |= Class::Aliphatic;
        } else {
          this->pos_flag = true;
          this->_class |= Class::Aliphatic;
        }
      }

    } else if (current_char == 'a') {
      if (not_flag) {
        this->neg_flag = true;
        this->not_class |= Class::Aliphatic;
      } else {
        this->_class |= Class::Aliphatic;
        this->pos_flag = true;
      }

    } else if (current_char == 'B') {
      if (std::string("aei").find(next_char) != std::string::npos) {
        this->parseBlock(
            Atom::Sym, not_flag, std::string() + current_char + next_char);
        iter_step++;
      } else if (next_char == 'r') {
        this->parseBlock(
            Atom::Br, not_flag, std::string() + current_char + next_char);
        iter_step++;
      } else {
        this->parseBlock(Atom::B, not_flag, current_char);
      }
      atom_seen = true;

    } else if (current_char == 'C') {
      if (next_char == 'a') {
        this->parseBlock(
            Atom::Ca, not_flag, std::string() + current_char + next_char);
        iter_step++;
      } else if (next_char == 'u') {
        this->parseBlock(
            Atom::Cu, not_flag, std::string() + current_char + next_char);
        iter_step++;
      } else if (next_char == 'l') {
        this->parseBlock(
            Atom::Cl, not_flag, std::string() + current_char + next_char);
        iter_step++;
      } else if (std::string("deors").find(next_char) != std::string::npos) {
        this->parseBlock(
            Atom::Sym, not_flag, std::string() + current_char + next_char);
        iter_step++;
      } else {
        this->parseBlock(Atom::C, not_flag, current_char);
      }
      atom_seen = true;

    } else if (current_char == 'D') {
      if (next_char == 'y') {
        this->parseBlock(
            Atom::Sym, not_flag, std::string() + current_char + next_char);
        atom_seen = true;
        iter_step++;
      } else {
        num = parseNumeral(next_char);
        if (num >= 0) {
          if (not_flag) {
            this->neg_flag = true;
            this->not_degree |= num_to_degree[num];
          } else {
            this->pos_flag = true;
            this->degree |= num_to_degree[num];
          }
          iter_step++;
        } else
          ok = false;
      }

    } else if (current_char == 'E') {
      if (std::string("ru").find(next_char) != std::string::npos) {
        this->parseBlock(
            Atom::Sym, not_flag, std::string() + current_char + next_char);
        iter_step++;
      } else {
        this->parseBlock(Atom::E, not_flag, current_char);
      }
      atom_seen = true;

    } else if (current_char == 'F') {
      if (next_char == 'e') {
        this->parseBlock(
            Atom::Fe, not_flag, std::string() + current_char + next_char);
        iter_step++;
      } else if (next_char == 'r') {
        this->parseBlock(
            Atom::Sym, not_flag, std::string() + current_char + next_char);
        iter_step++;
      } else {
        this->parseBlock(Atom::F, not_flag, current_char);
      }
      atom_seen = true;

    } else if (current_char == 'G' &&
               std::string("ade").find(next_char) != std::string::npos) {
      this->parseBlock(
          Atom::Sym, not_flag, std::string() + current_char + next_char);
      atom_seen = true;
      iter_step++;

    } else if (current_char == 'H') {
      if (std::string("efgo").find(next_char) != std::string::npos) {
        this->parseBlock(
            Atom::Sym, not_flag, std::string() + current_char + next_char);
        atom_seen = true;
        iter_step++;
      } else {
        if (!atom_seen) {
          this->parseBlock(Atom::H, not_flag, current_char);
          atom_seen = true;
        } else {
          num = parseNumeral(next_char);
          if (num >= 0) {
            iter_step++;
          } else {
            num = 1;
            // We have to set iter_step = 0,
            // coz iterator will be updated in while loop
            iter_step = 0;
            while (char_iterator != it_end || *char_iterator != 'H') {
              num++;
              char_iterator = std::next(char_iterator);
            }
          }
          this->imp_hydro = num;
          this->tot_hydro = num;
          this->hydro_flag = true;
          // turn on hydrogen count matching for this atom
        }
      }

    } else if (current_char == 'I') {
      if (std::string("nr").find(next_char) != std::string::npos) {
        this->parseBlock(
            Atom::Sym, not_flag, std::string() + current_char + next_char);
        iter_step++;
      } else {
        this->parseBlock(Atom::I, not_flag, current_char);
      }
      atom_seen = true;

    } else if (current_char == 'J') {
      this->parseBlock(Atom::J, not_flag, current_char);
      atom_seen = true;

    } else if (current_char == 'K') {
      this->parseBlock(Atom::K, not_flag, current_char);
      atom_seen = true;

    } else if (current_char == 'L') {
      if (std::string("aiu").find(next_char) != std::string::npos) {
        this->parseBlock(
            Atom::Sym, not_flag, std::string() + current_char + next_char);
        iter_step++;
      } else {
        this->parseBlock(Atom::L, not_flag, current_char);
      }
      atom_seen = true;

    } else if (current_char == 'M') {
      if (std::string("on").find(next_char) != std::string::npos) {
        this->parseBlock(
            Atom::Sym, not_flag, std::string() + current_char + next_char);
        iter_step++;
      } else if (next_char == 'g') {
        this->parseBlock(
            Atom::Mg, not_flag, std::string() + current_char + next_char);
        iter_step++;
      } else {
        this->parseBlock(Atom::M, not_flag, current_char);
      }
      atom_seen = true;

    } else if (current_char == 'N') {
      if (std::string("bdi").find(next_char) != std::string::npos) {
        this->parseBlock(
            Atom::Sym, not_flag, std::string() + current_char + next_char);
        iter_step++;
      } else if (next_char == 'a') {
        this->parseBlock(
            Atom::Na, not_flag, std::string() + current_char + next_char);
        iter_step++;
      } else {
        this->parseBlock(Atom::N, not_flag, current_char);
      }
      atom_seen = true;

    } else if (current_char == 'O') {
      if (next_char == 's') {
        this->parseBlock(
            Atom::Sym, not_flag, std::string() + current_char + next_char);
        iter_step++;
      } else {
        this->parseBlock(Atom::O, not_flag, current_char);
      }
      atom_seen = true;

    } else if (current_char == 'P') {
      if (std::string("bdort").find(next_char) != std::string::npos) {
        this->parseBlock(
            Atom::Sym, not_flag, std::string() + current_char + next_char);
        iter_step++;
      } else {
        this->parseBlock(Atom::P, not_flag, current_char);
      }
      atom_seen = true;

    } else if (current_char == 'p') { // Pi system
      if (not_flag) {
        this->neg_flag = true;
        this->not_class |= Class::Pi;
      } else {
        this->pos_flag = true;
        this->_class |= Class::Pi;
      }

    } else if (current_char == 'Q') {
      this->parseBlock(Atom::Q, not_flag, current_char);
      atom_seen = true;

    } else if (current_char == 'R') {
      if (std::string("behu").find(next_char) != std::string::npos) {
        this->parseBlock(
            Atom::Sym, not_flag, std::string() + current_char + next_char);
        iter_step++;
      } else {
        this->parseBlock(Atom::R, not_flag, current_char);
      }
      atom_seen = true;

    } else if (current_char == 'r') {
      num = parseNumeral(next_char);
      if (num >= 0) {
        if (not_flag) {
          this->neg_flag = true;
          this->not_cycle |= num_to_ring[num];
        } else {
          this->pos_flag = true;
          this->cycle |= num_to_ring[num];
        }
        iter_step++;
      } else {
        if (not_flag) {
          this->neg_flag = true;
          this->not_cycle |= Cycles::Ring3 | Cycles::Ring4 | Cycles::Ring5 |
                             Cycles::Ring6 | Cycles::Ring7 | Cycles::Ring8;
        } else {
          this->pos_flag = true;
          this->cycle |= Cycles::Ring3 | Cycles::Ring4 | Cycles::Ring5 |
                         Cycles::Ring6 | Cycles::Ring7 | Cycles::Ring8;
        }
      }

    } else if (current_char == 'S') {
      if (std::string("bcimnr").find(next_char) != std::string::npos) {
        this->parseBlock(
            Atom::Sym, not_flag, std::string() + current_char + next_char);
        iter_step++;
      } else if (next_char == 'e') {
        this->parseBlock(
            Atom::Se, not_flag, std::string() + current_char + next_char);
        iter_step++;
      } else {
        this->parseBlock(Atom::S, not_flag, current_char);
      }
      atom_seen = true;

    } else if (current_char == 'T') {
      if (std::string("abeihlm").find(next_char) != std::string::npos) {
        this->parseBlock(
            Atom::Sym, not_flag, std::string() + current_char + next_char);
        iter_step++;
      } else {
        this->parseBlock(Atom::T, not_flag, current_char);
      }
      atom_seen = true;

    } else if (current_char == 'U') {
      this->parseBlock(Atom::Sym, not_flag, current_char);
      atom_seen = true;

    } else if (current_char == 'V') {
      this->parseBlock(Atom::Sym, not_flag, current_char);

    } else if (current_char == 'v') { // Pi system
      num = parseNumeral(next_char);
      if (num >= 0) {
        if (not_flag) {
          this->neg_flag = true;
          this->not_valence |= num_to_valence[num];
        } else {
          this->pos_flag = true;
          this->degree |= num_to_valence[num];
        }
        iter_step++;
      } else
        ok = false;

    } else if (current_char == 'W') {
      this->parseBlock(Atom::Sym, not_flag, current_char);
      atom_seen = true;

    } else if (current_char == 'Y') {
      if (next_char == 'b') {
        this->parseBlock(
            Atom::Sym, not_flag, std::string() + current_char + next_char);
        iter_step++;
      } else {
        this->parseBlock(Atom::Sym, not_flag, current_char);
      }
      atom_seen = true;

    } else if (current_char == 'Z') {
      if (next_char == 'r') {
        this->parseBlock(
            Atom::Sym, not_flag, std::string() + current_char + next_char);
        iter_step++;
      } else if (next_char == 'n') {
        this->parseBlock(
            Atom::Zn, not_flag, std::string() + current_char + next_char);
        iter_step++;
      } else {
        this->parseBlock(Atom::Z, not_flag, current_char);
      }
      atom_seen = true;

    } else if (current_char == '+') {
      num = parseNumeral(next_char);
      if (num >= 0) {
        iter_step++;
      } else {
        num = 1;

        // We have to set iter_step = 0,
        // coz iterator will be updated in while loop
        iter_step = 0;

        while (char_iterator != it_end || *char_iterator != '+') {
          num++;
          char_iterator = std::next(char_iterator);
        }
      }
      if (num <= 5) {
        Charge charge_by_num = num_to_charge[num];
        if (not_flag) {
          this->neg_flag = true;
          this->not_charge |= charge_by_num;
        } else {
          this->pos_flag = true;
          this->charge |= charge_by_num;
        }
      }

    } else if (current_char == '-') {
      num = parseNumeral(next_char);
      if (num >= 0) {
        iter_step++;
      } else {
        num = 1;

        // We have to set iter_step = 0,
        // coz iterator will be updated in while loop
        iter_step = 0;

        while (char_iterator != it_end || *char_iterator != '-') {
          num++;
          char_iterator = std::next(char_iterator);
        }
      }
      if (num <= 5) {
        Charge charge_by_num = num_to_charge[num];
        if (not_flag) {
          this->neg_flag = true;
          this->not_charge |= charge_by_num;
        } else {
          this->pos_flag = true;
          this->charge |= charge_by_num;
        }
      }
    } else {
      break;
    }

    char_iterator = std::next(char_iterator, iter_step);
  }

  return done;
}
