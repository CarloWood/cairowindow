#pragma once

#include "utils/iomanip.h"

class IOManipFullDef : public utils::iomanip::Unsticky<2>
{
 private:
  static utils::iomanip::Index s_index;

 public:
  IOManipFullDef(bool full_def) : Unsticky(s_index, full_def) { }

  static bool is_full_def(std::ostream& os) { return get_iword_from(os, s_index); }
};

extern IOManipFullDef fulldef;
