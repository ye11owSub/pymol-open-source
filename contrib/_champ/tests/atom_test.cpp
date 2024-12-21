#include "atom.hpp"
#include "constants.hpp"
#include <iterator>
#define CATCH_CONFIG_MAIN
#include <catch2/catch_test_macros.hpp>

using namespace std;

TEST_CASE("Atom::parseAliphatic")
{
  auto atom_meta = AtomMeta();
  atom_meta.parseAliphatic(Atom::C, true);

  REQUIRE(atom_meta.atom == static_cast<unsigned int>(Atom::C));
  REQUIRE(atom_meta.pos_flag == true);
  REQUIRE(atom_meta.comp_imp_hydro_flag == true);
}

TEST_CASE("Atom::parseAromatic")
{
  auto atom_meta = AtomMeta();
  atom_meta.parseAromatic(Atom::C, false);

  REQUIRE(atom_meta.atom == static_cast<unsigned int>(Atom::C));
  REQUIRE(
      atom_meta.atom_type == static_cast<unsigned int>(Class::Aliphatic));
  REQUIRE(atom_meta.pos_flag == true);
  REQUIRE(atom_meta.comp_imp_hydro_flag == false);
}

TEST_CASE("Atom::parseString")
{
  auto atom_meta = AtomMeta();
  atom_meta.parseString("AM");

  REQUIRE(atom_meta.atom == static_cast<unsigned int>(Atom::Any));
  REQUIRE(atom_meta.pos_flag == true);
  REQUIRE(atom_meta.symbol[0] == 'A');
  REQUIRE(atom_meta.symbol[1] == 'M');
}

TEST_CASE("Atom::parseBlock")
{
  auto atom_meta = AtomMeta();

  SECTION("not_atom == true and symbols.length == 2")
  {
    atom_meta.parseBlock(Atom::Sym, true, "AM");

    REQUIRE(atom_meta.not_atom == static_cast<unsigned int>(Atom::Sym));
    REQUIRE(atom_meta.neg_flag == true);
    REQUIRE(atom_meta.hydro_flag == true);
    REQUIRE(atom_meta.symbol[0] == 'A');
    REQUIRE(atom_meta.symbol[1] == 'M');
    REQUIRE(atom_meta.symbol[2] == 0);
  }

  SECTION("not_atom == false and symbols.length == 1")
  {
    atom_meta.parseBlock(Atom::Sym, false, 'M');

    REQUIRE(atom_meta.atom == static_cast<unsigned int>(Atom::Sym));
    REQUIRE(atom_meta.pos_flag == true);
    REQUIRE(atom_meta.hydro_flag == true);
    REQUIRE(atom_meta.symbol[0] == 'M');
    REQUIRE(atom_meta.symbol[1] == 0);
  }

  //  SECTION("not_atom == false, symbols.length == 1, type(symbols) == string")
  //  {
  //    atom_meta.parseBlock(Atom::Sym, false, "M");
  //
  //    REQUIRE(atom_meta.atom == static_cast<unsigned int>(Atom::Sym));
  //    REQUIRE(atom_meta.pos_flag == true);
  //    REQUIRE(atom_meta.hydro_flag == true);
  //    REQUIRE(atom_meta.symbol[0] == 'M');
  //    REQUIRE(atom_meta.symbol[1] == 0);
  //  }
}

TEST_CASE("Atom::parseAtomBlock")
{
  auto atom_meta = AtomMeta();

  SECTION("case ']'")
  {
    string test_string = "]";
    string::iterator test_string_iterator = test_string.begin();
    atom_meta.parseAtomBlock(test_string_iterator, test_string.end());
    REQUIRE(test_string.end() == test_string.end());
  }

  SECTION("case '*'")
  {
    const vector<pair<string, Atom>> test_cases = {
        {string("*"), Atom::Any},
        {string("?"), Atom::Any},
    };
    for (auto [test_string, expected_atom] : test_cases) {

      string::iterator test_string_iterator = test_string.begin();
      atom_meta.parseAtomBlock(test_string_iterator, test_string.end());
      REQUIRE(test_string.end() == test_string.end());

      REQUIRE(atom_meta.atom == static_cast<unsigned int>(expected_atom));
      REQUIRE(atom_meta.pos_flag == true);
      REQUIRE(atom_meta.hydro_flag == true);
      REQUIRE(atom_meta.symbol[0] == 0);
    }
  }

  const vector<pair<string, Atom>> single_atom_test_cases = {
      {string("Ac"), Atom::Sym},
      {string("B"), Atom::B},
      {string("Br"), Atom::Br},
      {string("Ba"), Atom::Sym},
      {string("C"), Atom::C},
      {string("Ca"), Atom::Ca},
      {string("Cu"), Atom::Cu},
      {string("Cl"), Atom::Cl},
      {string("Cd"), Atom::Sym},
      {string("Dy"), Atom::Sym},
      {string("E"), Atom::E},
      {string("Er"), Atom::Sym},
      {string("F"), Atom::F},
      {string("Fe"), Atom::Fe},
      {string("Fr"), Atom::Sym},
      {string("Ga"), Atom::Sym},
      {string("H"), Atom::H},
      {string("He"), Atom::Sym},
      {string("I"), Atom::I},
      {string("In"), Atom::Sym},
      {string("J"), Atom::J},
      {string("K"), Atom::K},
      {string("L"), Atom::L},
      {string("La"), Atom::Sym},
      {string("M"), Atom::M},
      {string("Mg"), Atom::Mg},
      {string("Mo"), Atom::Sym},
      {string("N"), Atom::N},
      {string("Na"), Atom::Na},
      {string("Nb"), Atom::Sym},
      {string("O"), Atom::O},
      {string("Os"), Atom::Sym},
      {string("P"), Atom::P},
      {string("Pb"), Atom::Sym},
      {string("Q"), Atom::Q},
      {string("R"), Atom::R},
      {string("Rb"), Atom::Sym},
      {string("S"), Atom::S},
      {string("Se"), Atom::Se},
      {string("Si"), Atom::Sym},
      {string("T"), Atom::T},
      {string("Ta"), Atom::Sym},
      {string("U"), Atom::Sym},
      {string("V"), Atom::Sym},
      {string("W"), Atom::Sym},
      {string("Y"), Atom::Sym},
      {string("Yb"), Atom::Sym},
      {string("Z"), Atom::Z},
      {string("Zn"), Atom::Zn},
      {string("Zr"), Atom::Sym},
  };
  for (auto [test_string, expected_atom] : single_atom_test_cases) {

    DYNAMIC_SECTION("single atom: " << test_string)
    {

      string::iterator test_string_iterator = test_string.begin();
      atom_meta.parseAtomBlock(test_string_iterator, test_string.end());
      REQUIRE(test_string.end() == test_string.end());

      REQUIRE(atom_meta.atom == static_cast<unsigned int>(expected_atom));
      REQUIRE(atom_meta.pos_flag == true);
      REQUIRE(atom_meta.hydro_flag == true);

      if (expected_atom == Atom::Sym)
        REQUIRE(atom_meta.symbol == test_string);
    }
  }
}
