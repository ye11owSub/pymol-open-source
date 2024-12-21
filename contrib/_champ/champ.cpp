#include "champ.hpp"
#include "atom.hpp"
#include "bond.hpp"
#include "constants.hpp"
#include <pybind11/pybind11.h>

namespace py = pybind11;

Champ::Champ(){};

int smiToPat(const std::string smi)
{                         // returns root atom of list
  int mark[MAX_RING];     // ring marks 0-9
  int mark_pri[MAX_RING]; // lexical priority of mark
  int stack = 0;          // parenthetical scopes
  int base_atom = 0;
  int last_atom = 0;
  int last_bond = 0;
  int atom_list = 0;
  int bond_list = 0;
  int bond_flag = false;
  int mark_code = 0;
  int result = 0;
  bool ok = true;
  int not_bond = false;
  int lex_pri = 0;
  unsigned int bond_tags = 0;
  unsigned int bond_not_tags = 0;
  std::string::const_iterator it = smi.begin();
  std::string::const_iterator it_end = smi.begin();
  AtomMeta cur_atom = AtomMeta();
  Bond cur_bond = Bond();
  char current_char;
  char next_char;
  Sym sym;

  for (int i = 0; i < MAX_RING; i++)
    mark[i] = 0;

#define save_bond()                                                            \
  {                                                                            \
    if (last_bond) {                                                           \
      I->Bond[last_bond].link = cur_bond;                                      \
    } else {                                                                   \
      bond_list = cur_bond;                                                    \
    }                                                                          \
    last_bond = cur_bond;                                                      \
    cur_bond = ListElemNewZero(&I->Bond);                                      \
  }

#define save_atom()                                                            \
  {                                                                            \
    if (last_atom) {                                                           \
      I->Atom[last_atom].link = cur_atom;                                      \
    } else {                                                                   \
      atom_list = cur_atom;                                                    \
    }                                                                          \
    last_atom = cur_atom;                                                      \
    cur_atom = ListElemNewZero(&I->Atom);                                      \
  }

  //  PRINTFD(FB_smiles_parsing)
  //    " ChampSmiToPat: input '%s'\n",c
  //    ENDFD;

  while (it != it_end && ok) {
    lex_pri++;
    current_char = *it;
    next_char = *(std::next(it));
    //    PRINTFD(FB_smiles_parsing)
    //    " parsing: '%c' at %p\n", *c, c ENDFD;
    sym = Sym::_Null;
    /* ============ ROOT LEVEL PARSTING ============ */
    if ((current_char >= '0') && (current_char <= '9')) {
      sym = Sym::Mark;
      mark_code = current_char - '0';
      it = std::next(it);

      // standard, implicit atoms, with lowest normal valences
      // B(3), C(4), N(3,5), O(2), P(3,5), S(2,4,6), F(1), Cl(1), Br(1), I(1)
    } else if (current_char == 'C') {
      if (next_char == 'l' || next_char == 'L') {
        cur_atom.parseAliphatic(Atom::Cl, false);
        sym = Sym::Atom;
      } else {
        cur_atom.parseAliphatic(Atom::C, true);
        sym = Sym::Atom;
        // PRINTFD(FB_smiles_parsing)
        // " parsed: %p\n", c ENDFD;
      }

      // tag index/list
    } else if (current_char == '<') {
      if (bond_flag)
        cur_bond.parseTag(it, it_end, ok);
      else if (base_atom)
        cur_bond.parseTag(it, it_end, ok);
      else
        ok = false;

      sym = Sym::Qualifier;

    } else if (current_char == '*') {
      // nonstandard?
      cur_atom.parseAliphatic(Atom::Any, false);
      it = std::next(it);
      sym = Sym::Atom;

    } else if (current_char == '?') {
      // nonstandard?
      cur_atom.parseAliphatic(Atom::NotH, false);
      it = std::next(it);
      sym = Sym::Atom;

    } else if (current_char == 'H') {
      // nonstandard?
      cur_atom.parseAliphatic(Atom::H, false);
      it = std::next(it);
      sym = Sym::Atom;

    } else if (current_char == 'N' || current_char == 'n') {
      // nonstandard?
      cur_atom.parseAliphatic(Atom::N, true);
      it = std::next(it);
      sym = Sym::Atom;

    } else if (current_char == 'O' || current_char == 'o') {
      // nonstandard?
      cur_atom.parseAliphatic(Atom::O, true);
      it = std::next(it);
      sym = Sym::Atom;

    } else if (current_char == 'B') {
      if (std::string("Rr").find(next_char) != std::string::npos) {
        cur_atom.parseAliphatic(Atom::Br, true);
        it = std::next(it);
        sym = Sym::Atom;
      } else {
        cur_atom.parseAliphatic(Atom::B, true);
        it = std::next(it);
        sym = Sym::Atom;
      }

    } else if (current_char == 'P') {
      cur_atom.parseAliphatic(Atom::P, true);
      it = std::next(it);
      sym = Sym::Atom;

    } else if (current_char == 'S' || current_char == 's') {
      cur_atom.parseAliphatic(Atom::S, true);
      it = std::next(it);
      sym = Sym::Atom;

    } else if (current_char == 'F') {
      cur_atom.parseAliphatic(Atom::F, true);
      it = std::next(it);
      sym = Sym::Atom;

    } else if (current_char == 'I') {
      cur_atom.parseAliphatic(Atom::I, true);
      it = std::next(it);
      sym = Sym::Atom;

    } else if (current_char == 'c') {
      // standard implicit aromatic atoms
      cur_atom.parseAliphatic(Atom::C, true);
      it = std::next(it);
      sym = Sym::Atom;

    } else if (current_char == ';') {
      it = std::next(it);
      not_bond = false;
      sym = Sym::Qualifier;

    } else if (current_char == ',') {
      it = std::next(it);
      sym = Sym::Qualifier;

    } else if (current_char == '!') {
      it = std::next(it);
      not_bond = true;
      sym = Sym::Qualifier;

    } else if (current_char == '-') {
      it = std::next(it);
      if (not_bond)
        cur_bond.not_order |= Order::Single;
      else
        cur_bond.order |= Order::Single;

      sym = Sym::Bond;

    } else if (current_char == '/') {
      it = std::next(it);
      if (not_bond)
        cur_bond.not_order |= Order::Single;
      else
        cur_bond.order |= Order::Single;

      sym = Sym::Bond;
      cur_bond.direction = Direction::Up;

    } else if (current_char == '\\') {
      it = std::next(it);
      if (not_bond)
        cur_bond.not_order |= Order::Single;
      else
        cur_bond.order |= Order::Single;

      sym = Sym::Bond;
      cur_bond.direction = Direction::Down;

    } else if (current_char == '=') {
      it = std::next(it);
      if (not_bond)
        cur_bond.not_order |= Order::Double;
      else
        cur_bond.order |= Order::Double;

      sym = Sym::Bond;

    } else if (current_char == '#') {
      it = std::next(it);
      if (not_bond)
        cur_bond.not_order |= Order::Triple;
      else
        cur_bond.order |= Order::Triple;

      sym = Sym::Bond;

    } else if (current_char == '~') {
      it = std::next(it);
      if (not_bond) {
        cur_bond.not_order |= Order::AnyOrder;
        cur_bond.not_class |= Class::Any;
      } else {
        cur_bond.order |= Order::AnyOrder;
        cur_bond._class |= Class::Any;
      }
      sym = Sym::Bond;

    } else if (current_char == '@') {
      it = std::next(it);
      if (not_bond)
        cur_bond.not_cycle |= Cycles::Cyclic;
      else
        cur_bond.cycle |= Cycles::Cyclic;

      sym = Sym::Bond;

    } else if (current_char == ':') {
      it = std::next(it);
      if (not_bond)
        cur_bond.not_class |= Class::Aromatic;
      else
        cur_bond._class |= Class::Aromatic;

      sym = Sym::Bond;

    } else if (current_char == '.') {
      it = std::next(it);
      sym = Sym::Separator;

    } else if (current_char == '%') {
      it = std::next(it);
      sym = Sym::Separator;
    }

  case '%':
    c++;
    if (c) {
      mark_code = 10 * ((*c) - '0');
      c++;
    } /* else error */
    if (c) {
      sym = cSym_Mark;
      mark_code += (*c) - '0';
      c++;
    } /* else error */
    break;
  case '(':
    c++;
    sym = cSym_OpenScope;
    break;
  case ')':
    c++;
    sym = cSym_CloseScope;
    break;
  case '[':
    c++;
    sym = cSym_OpenBlock;
    break;
  case ']':
    c++;
    sym = cSym_CloseBlock;
    break;
  }
}
if (sym == cSym_Null) {
  PRINTFB(FB_smiles_parsing, FB_errors)
  " champ: error parsing smiles string at '%c' (char %zd) in\n champ: "
  "'%s'\n",
      *c, c - orig_c, orig_c ENDFB;
  ok = false;
}
if (ok) {
  /* =========== actions based on root level parsing ========== */
  switch (sym) {
  case cSym_OpenBlock:
    ok = ChampParseAtomBlock(I, &c, cur_atom);
  case cSym_Atom:
    /* was there a preceeding atom? if so, then form bond and save atom */
    if (base_atom) {
      PRINTFD(FB_smiles_parsing)
      " ChampSmiToPtr: saving atom %d\n", last_atom ENDFD;
      /* backward link */
      I->Bond[cur_bond].atom[0] = base_atom;
      I->Bond[cur_bond].atom[1] = cur_atom;
      I->Bond[cur_bond].pri[0] = lex_pri;
      I->Bond[cur_bond].pri[1] = lex_pri;
      if (!bond_flag) {
        if ((I->Atom[cur_atom].class & cH_Aromatic) &&
            (I->Atom[base_atom].class & cH_Aromatic))
          I->Bond[cur_bond].order =
              (cH_Single | cH_Aromatic); /* is this right? */
        else
          I->Bond[cur_bond].order = cH_Single;
      }
      I->Bond[cur_bond].tag = bond_tags;         /* save bond tags */
      I->Bond[cur_bond].not_tag = bond_not_tags; /* save bond tags */
      bond_tags = 0;
      bond_not_tags = 0;
      ok = ChampAddBondToAtom(I, cur_atom, cur_bond);
      if (ok) {
        ok = ChampAddBondToAtom(I, base_atom, cur_bond);
        save_bond();
      }
      bond_flag = false;
      not_bond = false;
    }
    base_atom = cur_atom;
    save_atom();
    break;
  case cSym_CloseBlock: /* should never be reached */
    break;
  case cSym_OpenScope: /* push base_atom onto stack */
    stack = ListElemPushInt(&I->Int, stack, base_atom);
    break;
  case cSym_CloseScope:
    if (!stack) {
      PRINTFB(FB_smiles_parsing, FB_errors)
      " champ: stack underflow for scope...\n" ENDFB;
      ok = false;
    } else {
      base_atom = I->Int[stack].value;
      stack = ListElemPop(I->Int, stack);
    }
    break;
  case cSym_Bond:
    bond_flag = true;
    break;
  case cSym_Mark:
    if (base_atom) {
      if (!mark[mark_code]) { /* opening cycle */
        mark[mark_code] = base_atom;
        mark_pri[mark_code] = lex_pri;
        bond_flag = false; /* ignore the first bond valence...we'll get it
                              from the second half of the mark*/
        not_bond = false;
      } else { /* closing cycle */
        I->Bond[cur_bond].atom[0] = base_atom;
        I->Bond[cur_bond].atom[1] = mark[mark_code];
        I->Bond[cur_bond].pri[0] = lex_pri;
        I->Bond[cur_bond].pri[1] = mark_pri[mark_code];
        if (!bond_flag) {
          I->Bond[cur_bond].order = cH_Single;
        }
        ok = ChampAddBondToAtom(I, base_atom, cur_bond);
        if (ok) {
          ok = ChampAddBondToAtom(I, mark[mark_code], cur_bond);
          save_bond();
        }
        mark[mark_code] = 0;
        bond_flag = false;
        not_bond = false;
      }
    } else {
      PRINTFB(FB_smiles_parsing, FB_errors)
      " champ:  syntax error...\n" ENDFB;
      ok = false;
    }
    break;
  case cSym_Separator:
    base_atom = 0;
    break;
  case cSym_Qualifier:
    break;
  }
}
}
if (ok && atom_list) {
  result = ListElemNewZero(&I->Pat);
  if (result) {
    I->ActivePatList = ListElemPushInt(&I->Int, I->ActivePatList, result);
    I->Pat[result].atom = atom_list;
    I->Pat[result].bond = bond_list;
  } else
    ok = false;
}
if (cur_atom)
  ChampAtomFree(I, cur_atom);
if (cur_bond)
  ChampBondFree(I, cur_bond);
if (result)
  ChampPatReindex(I, result);

PRINTFD(FB_smiles_parsing)
" ChampSmiToPtr: returning pattern %d atom_list %d bond_list %d\n", result,
    atom_list, bond_list ENDFD;

return (result);
}
PYBIND11_MODULE(contrib_with_py11, handle)
{
  handle.doc() = "example of using pybind11";
  py::class_<Champ>(handle, "Champ")
      .def(py::init())
      .def("memory_dump", &Champ::memoryDump);
};
