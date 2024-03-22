#include "sys.h"
#include "Symbol.h"
#include <vector>
#include "debug.h"

namespace symbolic {

namespace {

struct Data {
  char const* name;
  double value;
};

std::vector<Data> symbols;

} // namespace

//static
void SymbolRegistry::register_symbol(int id, char const* name)
{
  ASSERT(0 <= id);
  if (id >= symbols.size())
    symbols.resize(id + 1);
  symbols[id].name = name;
}

//static
char const* SymbolRegistry::get_name(int id)
{
  ASSERT(0 <= id && id < symbols.size());
  return symbols[id].name;
}

//static
void SymbolRegistry::set_value(int id, double value)
{
  ASSERT(0 <= id && id < symbols.size());
  symbols[id].value = value;
}

//static
double SymbolRegistry::get_value(int id)
{
  ASSERT(0 <= id && id < symbols.size());
  return symbols[id].value;
}

} // namespace symbolic
