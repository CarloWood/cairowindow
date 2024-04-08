#include <iostream>
#include <vector>
#include <array>
#include <fstream>
#include <cassert>

std::string get_trailing_alpha(std::string const& text)
{
  int end = text.size() - 1;
  while (end >= 0 && std::isalpha(text[end]))
    --end;
  return text.substr(end + 1);
}

class Decoder
{
 private:
  std::string substring;

  std::array<std::string, 1000> a;
  int depth_ = 0;
  int print_no_new_lines = 0;
  bool consume_spaces = false;

 public:
  int depth() const { return depth_; }

 public:
  void add_char(char c)
  {
    if (consume_spaces && c == ' ')
      return;
    substring += c;
    consume_spaces = false;
  }

  void add_substring()
  {
    if (!substring.empty())
    {
      add(substring);
      substring.clear();
    }
  }

  void open()
  {
    add_substring();
    std::cout << '<';
    std::string word = get_trailing_alpha(a[depth_]);
    if (word == "Constant" || word == "Symbol" || word == "Sin")
      ++print_no_new_lines;
    else if (!print_no_new_lines)
    {
      std::cout << '\n';
      consume_spaces = true;
      for (int i = 0; i <= depth_; ++i)
        std::cout << "  ";
    }
    else
      ++print_no_new_lines;
    ++depth_;
    assert(depth_ < a.size());
  }

  void close()
  {
    add_substring();
    assert(depth_ > 0);
    --depth_;
    if (print_no_new_lines)
      --print_no_new_lines;
    else
    {
      std::cout << '\n';
      consume_spaces = true;
      for (int i = 0; i < depth_; ++i)
        std::cout << "  ";
    }
    std::cout << '>';
  }

  void comma()
  {
    add_substring();
    if (print_no_new_lines)
      std::cout << ", ";
    else
    {
      if (!print_no_new_lines)
      {
        std::cout << ",\n";
        consume_spaces = true;
        for (int i = 0; i < depth_; ++i)
          std::cout << "  ";
      }
      else
        std::cout << ", ";
    }
  }

  void add(std::string const& token)
  {
    std::cout << token;
    a[depth_] = token;
  }
};

int main()
{
  Decoder decoder;
  std::ifstream file("troep");

  std::string line;
  while (std::getline(file, line))
  {
    if (!line.empty() && *line.rbegin() == '\r')
      line.erase(line.end() - 1);
    if (line.rfind("cairowindow", 0) != 0)
    {
      decoder.add_substring();
      std::cout << line << '\n';
      continue;
    }
    for (char c : line)
    {
      if (c == '<')
        decoder.open();
      else if (c == '>')
        decoder.close();
      else if (c == ',')
        decoder.comma();
      else if (c != '\r')
        decoder.add_char(c);
    }
    decoder.add_char('\n');
  }

  // In case the file doesn't end with a separator.
  decoder.add_substring();
}
